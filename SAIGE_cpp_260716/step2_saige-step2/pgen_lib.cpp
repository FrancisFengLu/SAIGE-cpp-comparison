// pgen_lib.cpp -- see pgen_lib.hpp. The setup follows upstream SAIGE's
// PgenClass::loadPgen / Read (SAIGE src/PGEN.cpp, itself from pgenlibr), with
// one PgenReader per thread instead of one per file.
#include "pgen_lib.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <mutex>

#include "third_party/pgenlib/pgenlib_ffi_support.h"
#include "third_party/pgenlib/pgenlib_read.h"

namespace pgenlib_glue {

struct File {
    plink2::PgenFileInfo pgfi;
    unsigned char* pgfiAlloc = nullptr;
    uintptr_t* alleleIdxOffsets = nullptr;
    uintptr_t* nonrefFlags = nullptr;
    uint32_t maxVrecWidth = 0;
    uintptr_t pgrAllocCachelines = 0;
    std::string path;
    std::mutex initMu;   // PgrInit moves pgfi.shared_ff into the first reader
    FileFacts facts;
};

struct Reader {
    plink2::PgenReader pgr;
    unsigned char* alloc = nullptr;
    uintptr_t* genovec = nullptr;
    uintptr_t* dosagePresent = nullptr;
    uint16_t* dosageMain = nullptr;
    uint32_t rawSampleCt = 0;
    bool inited = false;
};

static std::string pglMsg(const char* buf) {
    // pgenlib's messages start with "Error: "
    std::string s(buf);
    if (s.rfind("Error: ", 0) == 0) s = s.substr(7);
    while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
    return s;
}

static void freeFileParts(File* f) {
    if (!f) return;
    plink2::PglErr reterr = plink2::kPglRetSuccess;
    plink2::CleanupPgfi(&f->pgfi, &reterr);
    if (f->pgfiAlloc) plink2::aligned_free(f->pgfiAlloc);
    if (f->alleleIdxOffsets) plink2::aligned_free(f->alleleIdxOffsets);
    if (f->nonrefFlags) plink2::aligned_free(f->nonrefFlags);
    f->pgfiAlloc = nullptr; f->alleleIdxOffsets = nullptr; f->nonrefFlags = nullptr;
}

File* openFile(const std::string& path, uint32_t rawSampleCt, uint32_t rawVariantCt,
               FileFacts& facts, std::string& err)
{
    // storage mode: byte 2 of the file
    int mode = -1;
    {
        FILE* fp = std::fopen(path.c_str(), "rb");
        if (!fp) { err = "cannot open " + path; return nullptr; }
        unsigned char hdr[3];
        if (std::fread(hdr, 1, 3, fp) != 3) { std::fclose(fp); err = "cannot read the header of " + path; return nullptr; }
        std::fclose(fp);
        if (hdr[0] != 0x6c || hdr[1] != 0x1b) { err = path + " is not a .pgen / .bed file (bad magic bytes)"; return nullptr; }
        mode = hdr[2];
    }
    File* f = new File();
    f->path = path;
    plink2::PreinitPgfi(&f->pgfi);
    char errbuf[plink2::kPglErrstrBufBlen];
    errbuf[0] = 0;
    plink2::PgenHeaderCtrl headerCtrl;
    uintptr_t pgfiAllocCachelines = 0;
    // A mode-0x01 file carries no counts; every other mode stores them and
    // pgenlib checks what we pass against them (UINT32_MAX = take the header's).
    const uint32_t vct = (mode == 0x01) ? rawVariantCt : UINT32_MAX;
    const uint32_t sct = (mode == 0x01) ? rawSampleCt : UINT32_MAX;
    if (plink2::PgfiInitPhase1(path.c_str(), nullptr, vct, sct, &headerCtrl, &f->pgfi,
                               &pgfiAllocCachelines, errbuf) != plink2::kPglRetSuccess) {
        err = pglMsg(errbuf);
        freeFileParts(f); delete f;
        return nullptr;
    }
    const uint32_t rawV = f->pgfi.raw_variant_ct;
    if (headerCtrl & 0x30) {
        if (plink2::cachealigned_malloc((rawV + 1) * sizeof(uintptr_t), &f->alleleIdxOffsets)) {
            err = "out of memory"; freeFileParts(f); delete f; return nullptr;
        }
        f->pgfi.allele_idx_offsets = f->alleleIdxOffsets;
    } else {
        f->pgfi.max_allele_ct = 2;
    }
    if ((headerCtrl & 0xc0) == 0xc0) {
        const uintptr_t words = plink2::DivUp(rawV, plink2::kBitsPerWord) + 1;
        if (plink2::cachealigned_malloc(words * sizeof(uintptr_t), &f->nonrefFlags)) {
            err = "out of memory"; freeFileParts(f); delete f; return nullptr;
        }
        f->pgfi.nonref_flags = f->nonrefFlags;
    }
    if (plink2::cachealigned_malloc(pgfiAllocCachelines * plink2::kCacheline, &f->pgfiAlloc)) {
        err = "out of memory"; freeFileParts(f); delete f; return nullptr;
    }
    if (plink2::PgfiInitPhase2(headerCtrl, 0, 0, 0, 0, rawV, &f->maxVrecWidth, &f->pgfi,
                               f->pgfiAlloc, &f->pgrAllocCachelines, errbuf) != plink2::kPglRetSuccess) {
        err = pglMsg(errbuf);
        freeFileParts(f); delete f;
        return nullptr;
    }
    const uint32_t g = f->pgfi.gflags;
    facts.mode = mode;
    facts.rawSampleCt = f->pgfi.raw_sample_ct;
    facts.rawVariantCt = rawV;
    facts.ldCompression = (g & plink2::kfPgenGlobalLdCompressionPresent) != 0;
    facts.dosage        = (g & plink2::kfPgenGlobalDosagePresent) != 0;
    facts.dosagePhase   = (g & plink2::kfPgenGlobalDosagePhasePresent) != 0;
    facts.hardcallPhase = (g & plink2::kfPgenGlobalHardcallPhasePresent) != 0;
    facts.multiallelic  = (f->pgfi.max_allele_ct > 2) ||
                          (g & plink2::kfPgenGlobalMultiallelicHardcallFound) != 0;
    f->facts = facts;
    return f;
}

void closeFile(File* f)
{
    if (!f) return;
    freeFileParts(f);
    delete f;
}

Reader* newReader(File* f, std::string& err)
{
    Reader* r = new Reader();
    plink2::PreinitPgr(&r->pgr);
    const uint32_t n = f->pgfi.raw_sample_ct;
    r->rawSampleCt = n;
    const uintptr_t mainBytes = f->pgrAllocCachelines * plink2::kCacheline;
    const uintptr_t genoBytes = plink2::NypCtToCachelineCt(n) * plink2::kCacheline;
    const uintptr_t bitBytes  = plink2::BitCtToCachelineCt(n) * plink2::kCacheline;
    const uintptr_t dosBytes  = plink2::RoundUpPow2((uintptr_t)n * sizeof(uint16_t), plink2::kCacheline);
    if (plink2::cachealigned_malloc(mainBytes + genoBytes + bitBytes + dosBytes, &r->alloc)) {
        delete r; err = "out of memory"; return nullptr;
    }
    r->genovec       = reinterpret_cast<uintptr_t*>(r->alloc + mainBytes);
    r->dosagePresent = reinterpret_cast<uintptr_t*>(r->alloc + mainBytes + genoBytes);
    r->dosageMain    = reinterpret_cast<uint16_t*>(r->alloc + mainBytes + genoBytes + bitBytes);
    plink2::PglErr rc;
    {
        std::lock_guard<std::mutex> lk(f->initMu);
        rc = plink2::PgrInit(f->path.c_str(), f->maxVrecWidth, &f->pgfi, &r->pgr, r->alloc);
    }
    if (rc != plink2::kPglRetSuccess) {
        plink2::aligned_free(r->alloc);
        delete r;
        err = "cannot open a reader on " + f->path + " (pgenlib error " + std::to_string((int)rc) + ")";
        return nullptr;
    }
    r->inited = true;
    return r;
}

void freeReader(Reader* r)
{
    if (!r) return;
    if (r->inited) {
        plink2::PglErr reterr = plink2::kPglRetSuccess;
        plink2::CleanupPgr(&r->pgr, &reterr);
    }
    if (r->alloc) plink2::aligned_free(r->alloc);
    delete r;
}

bool readAltDosage(Reader* r, uint32_t vidx, double* out, std::string& err)
{
    static const double kPairs[32] ALIGNV16 =
        PAIR_TABLE16(0.0, 1.0, 2.0, std::numeric_limits<double>::quiet_NaN());
    plink2::PgrSampleSubsetIndex pssi;
    plink2::PgrClearSampleSubsetIndex(nullptr, &pssi);
    uint32_t dosageCt = 0;
    const plink2::PglErr rc = plink2::PgrGet1D(nullptr, pssi, r->rawSampleCt, vidx, 1, &r->pgr,
                                               r->genovec, r->dosagePresent, r->dosageMain, &dosageCt);
    if (rc != plink2::kPglRetSuccess) {
        err = "pgenlib could not read variant " + std::to_string(vidx) + " (error " + std::to_string((int)rc) + ")";
        return false;
    }
    plink2::Dosage16ToDoubles(kPairs, r->genovec, r->dosagePresent, r->dosageMain,
                              r->rawSampleCt, dosageCt, out);
    return true;
}

bool readBedRow(Reader* r, uint32_t vidx, unsigned char* out, std::string& err)
{
    plink2::PgrSampleSubsetIndex pssi;
    plink2::PgrClearSampleSubsetIndex(nullptr, &pssi);
    const plink2::PglErr rc = plink2::PgrGet(nullptr, pssi, r->rawSampleCt, vidx, &r->pgr, r->genovec);
    if (rc != plink2::kPglRetSuccess) {
        err = "pgenlib could not read variant " + std::to_string(vidx) + " (error " + std::to_string((int)rc) + ")";
        return false;
    }
    // PLINK 2 code (00 hom REF, 01 het, 10 hom ALT, 11 missing) -> PLINK 1 code
    // with ALT = A1 (00 hom A1, 01 missing, 10 het, 11 hom A2):
    //   high bit' = NOT high,  low bit' = high XOR low XOR 1.
    const uint32_t n = r->rawSampleCt;
    const uint32_t nWords = (n + 31) / 32;
    const uint64_t M55 = 0x5555555555555555ULL;
    uint64_t* w = reinterpret_cast<uint64_t*>(r->genovec);
    for (uint32_t k = 0; k < nWords; k++) {
        const uint64_t x = w[k];
        const uint64_t hi = (x >> 1) & M55, lo = x & M55;
        w[k] = ((~hi & M55) << 1) | ((hi ^ lo ^ M55) & M55);
    }
    const uint32_t nb = (n + 3) / 4;
    std::memcpy(out, r->genovec, nb);
    if (n & 3u) out[nb - 1] &= (unsigned char)((1u << ((n & 3u) * 2)) - 1u);
    return true;
}

}  // namespace pgenlib_glue
