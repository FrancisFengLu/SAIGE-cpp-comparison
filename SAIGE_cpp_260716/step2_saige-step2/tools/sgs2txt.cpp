// sgs2txt -- turn the binary columnar output of `outputFormat: sgs` back into
// the exact text file the default writer would have produced.
//
//   sgs2txt <trait>.sgs [more.sgs ...]            -> writes each trait's text
//                                                    path recorded in its header
//   sgs2txt -o OUT <trait>.sgs                    -> writes OUT instead
//   sgs2txt -m markers.sgs <trait>.sgs ...        -> override the marker file
//   sgs2txt -j N ...                              -> N traits at a time
//
// The rows are emitted by the same SAIGE::outfast::format_text() the step-2
// binary itself uses, so "byte-identical" is a property of the code path, not a
// coincidence: the only question the round-trip test answers is whether the
// values survived the trip, and they do (sgs_format.hpp explains the one place
// that needed an argument, the p-value strings).
//
// A file written under sgsPrecision: fp32 is read by the same code -- the
// encoding byte carries the width -- but its values did NOT survive the trip
// intact, so its text is close to, not identical to, the original. The
// converter says so on stderr rather than leaving it to be discovered.

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "../out_fast.hpp"
#include "../sgs_format.hpp"

using namespace SAIGE;
using namespace SAIGE::sgs;

static bool slurp(const std::string& path, std::vector<char>& buf) {
    int fd = ::open(path.c_str(), O_RDONLY);
    if (fd < 0) { fprintf(stderr, "cannot open %s: %s\n", path.c_str(), strerror(errno)); return false; }
    struct stat st;
    if (fstat(fd, &st) != 0) { ::close(fd); return false; }
    buf.resize((size_t)st.st_size);
    size_t off = 0;
    while (off < buf.size()) {
        ssize_t r = ::read(fd, buf.data() + off, buf.size() - off);
        if (r <= 0) { ::close(fd); fprintf(stderr, "short read on %s\n", path.c_str()); return false; }
        off += (size_t)r;
    }
    ::close(fd);
    return true;
}

// ---- column decoders -------------------------------------------------------
template <typename T>
static void get_col_pod(Reader& r, std::vector<T>& v, size_t n) {
    v.assign(n, T());
    uint8_t e = r.u8();
    if (e == E_CONST) {
        const uint8_t* q = r.take(sizeof(T));
        if (!q) return;
        T x; memcpy(&x, q, sizeof(T));
        for (size_t i = 0; i < n; i++) v[i] = x;
    } else if (e == E_RAW) {
        const uint8_t* q = r.take(n * sizeof(T));
        if (q) memcpy(v.data(), q, n * sizeof(T));
    } else { r.bad = true; }
}

// A column the writer held as double, stored at either width.
static void get_col_f64(Reader& r, std::vector<double>& v, size_t n) {
    v.assign(n, 0.0);
    uint8_t e = r.u8();
    if (e == E_CONST) {
        const uint8_t* q = r.take(sizeof(double));
        if (!q) return;
        double x; memcpy(&x, q, sizeof(double));
        for (size_t i = 0; i < n; i++) v[i] = x;
    } else if (e == E_RAW) {
        const uint8_t* q = r.take(n * sizeof(double));
        if (q) memcpy(v.data(), q, n * sizeof(double));
    } else if (e == E_CONST32) {
        const uint8_t* q = r.take(sizeof(float));
        if (!q) return;
        float x; memcpy(&x, q, sizeof(float));
        for (size_t i = 0; i < n; i++) v[i] = (double)x;
    } else if (e == E_RAW32) {
        const uint8_t* q = r.take(n * sizeof(float));
        if (!q) return;
        const float* f = (const float*)(const void*)q;
        for (size_t i = 0; i < n; i++) v[i] = (double)f[i];
    } else { r.bad = true; }
}

static void get_col_str(Reader& r, std::vector<std::string>& v, size_t n) {
    v.assign(n, std::string());
    uint8_t e = r.u8();
    if (e == E_CONST) {
        std::string x = r.sstr();
        for (size_t i = 0; i < n; i++) v[i] = x;
    } else if (e == E_RAW) {
        for (size_t i = 0; i < n; i++) v[i] = r.sstr();
    } else { r.bad = true; }
}

static void get_col_pval(Reader& r, std::vector<std::string>& v, size_t n) {
    v.assign(n, std::string());
    uint8_t e = r.u8();
    if (e != E_PVAL && e != E_PVAL32) { r.bad = true; return; }
    const size_t w = (e == E_PVAL) ? sizeof(double) : sizeof(float);
    const uint8_t* q = r.take(n * w);
    if (!q) return;
    const double* d = (e == E_PVAL) ? (const double*)(const void*)q : nullptr;
    const float*  f = (e == E_PVAL) ? nullptr : (const float*)(const void*)q;
    char b[64];
    for (size_t i = 0; i < n; i++) {
        // The very sprintf that produced the string in score_format.hpp.
        snprintf(b, sizeof(b), "%.6E", d ? d[i] : (double)f[i]);
        v[i].assign(b);
    }
    uint32_t nexc = r.u32();
    for (uint32_t k = 0; k < nexc; k++) {
        uint32_t idx = r.u32();
        std::string s = r.sstr();
        if (idx < n) v[idx] = s;
    }
}

struct MarkerBlock {
    uint32_t nRows = 0;
    std::vector<std::string> chr, pos, mid, ref, alt;
    std::vector<double> ac, af, info;
};

static bool read_marker_block(Reader& r, MarkerBlock& B) {
    uint32_t m = r.u32();
    if (r.bad || m != BLK_MAGIC) return false;
    B.nRows = r.u32();
    const size_t n = B.nRows;
    get_col_str(r, B.chr, n); get_col_str(r, B.pos, n); get_col_str(r, B.mid, n);
    get_col_str(r, B.ref, n); get_col_str(r, B.alt, n);
    get_col_f64(r, B.ac, n);  get_col_f64(r, B.af, n);  get_col_f64(r, B.info, n);
    return !r.bad;
}

int main(int argc, char** argv) {
    std::string outOverride, markerOverride;
    int jobs = 1;
    std::vector<std::string> inputs;
    for (int i = 1; i < argc; i++) {
        std::string a = argv[i];
        if (a == "-o" && i + 1 < argc)      outOverride = argv[++i];
        else if (a == "-m" && i + 1 < argc) markerOverride = argv[++i];
        else if (a == "-j" && i + 1 < argc) jobs = atoi(argv[++i]);
        else if (a == "-h" || a == "--help") {
            fprintf(stderr, "usage: sgs2txt [-o OUT] [-m markers.sgs] [-j N] <trait>.sgs ...\n");
            return 0;
        } else inputs.push_back(a);
    }
    if (inputs.empty()) {
        fprintf(stderr, "usage: sgs2txt [-o OUT] [-m markers.sgs] [-j N] <trait>.sgs ...\n");
        return 2;
    }
    if (!outOverride.empty() && inputs.size() != 1) {
        fprintf(stderr, "-o takes exactly one input\n");
        return 2;
    }
    if (jobs < 1) jobs = 1;

    // `sgs2txt out/*.sgs` sweeps in the shared marker file too. Recognise it by
    // its magic and use it as the marker file instead of failing on it.
    {
        std::vector<std::string> keep;
        for (size_t i = 0; i < inputs.size(); i++) {
            char mg[8] = {0};
            FILE* f = fopen(inputs[i].c_str(), "rb");
            if (f) { if (fread(mg, 1, 8, f) != 8) mg[0] = 0; fclose(f); }
            if (memcmp(mg, MAGIC_MARKER, 8) == 0) {
                if (markerOverride.empty()) markerOverride = inputs[i];
            } else keep.push_back(inputs[i]);
        }
        inputs.swap(keep);
        if (inputs.empty()) { fprintf(stderr, "no trait .sgs files given\n"); return 2; }
    }

    // The marker file is shared: read it once.
    std::string markerPath = markerOverride;
    if (markerPath.empty()) {
        std::vector<char> head;
        if (!slurp(inputs[0], head)) return 1;
        Reader r(head.data(), head.size());
        const uint8_t* mg = r.take(8);
        if (!mg || memcmp(mg, MAGIC_TRAIT, 8) != 0) {
            fprintf(stderr, "%s is not a trait .sgs file\n", inputs[0].c_str());
            return 1;
        }
        r.u32(); r.u32(); r.str(); r.str(); r.str(); r.u8(); r.u8();
        markerPath = r.str();
    }
    std::vector<char> mbuf;
    if (!slurp(markerPath, mbuf)) return 1;
    std::vector<MarkerBlock> mblocks;
    bool isImputation = false;
    {
        Reader r(mbuf.data(), mbuf.size());
        const uint8_t* mg = r.take(8);
        if (!mg || memcmp(mg, MAGIC_MARKER, 8) != 0) {
            fprintf(stderr, "%s is not a marker .sgs file\n", markerPath.c_str());
            return 1;
        }
        uint32_t ver = r.u32();
        const uint32_t hf = r.u32();
        isImputation = (hf & H_IMPUTATION) != 0;
        if (ver != VERSION) { fprintf(stderr, "marker file version %u, expected %u\n", ver, VERSION); return 1; }
        if (hf & H_F32)
            fprintf(stderr, "note: %s was written with sgsPrecision: fp32; the text this "
                            "produces is NOT byte-identical to a text run\n", markerPath.c_str());
        while (!r.bad && (size_t)(r.end - r.p) > 4) {
            uint32_t peek; memcpy(&peek, r.p, 4);
            if (peek == END_MAGIC) break;
            MarkerBlock B;
            if (!read_marker_block(r, B)) { fprintf(stderr, "corrupt marker block\n"); return 1; }
            mblocks.push_back(B);
        }
    }

    int rc = 0;
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(jobs)
#endif
    for (long ii = 0; ii < (long)inputs.size(); ii++) {
        const std::string& in = inputs[(size_t)ii];
        std::vector<char> tbuf;
        if (!slurp(in, tbuf)) { rc = 1; continue; }
        Reader r(tbuf.data(), tbuf.size());
        const uint8_t* mg = r.take(8);
        if (!mg || memcmp(mg, MAGIC_TRAIT, 8) != 0) {
            fprintf(stderr, "%s is not a trait .sgs file\n", in.c_str()); rc = 1; continue;
        }
        uint32_t ver = r.u32();
        (void)r.u32();
        if (ver != VERSION) { fprintf(stderr, "%s version %u\n", in.c_str(), ver); rc = 1; continue; }
        const std::string header   = r.str();
        const std::string name     = r.str();
        const std::string ttype    = r.str();
        const bool isCond          = r.u8() != 0;
        const bool isMore          = r.u8() != 0;
        (void)r.str();                       // marker path
        const std::string textPath = r.str();
        uint32_t nc = r.u32();
        std::vector<uint8_t> cols(nc);
        for (uint32_t i = 0; i < nc; i++) cols[i] = r.u8();

        TraitMeta meta;
        meta.name = name; meta.traitType = ttype;
        meta.isCondition = isCond; meta.isMoreOutput = isMore;

        const std::string outPath = outOverride.empty() ? textPath : outOverride;
        FILE* fo = fopen(outPath.c_str(), "wb");
        if (!fo) { fprintf(stderr, "cannot write %s: %s\n", outPath.c_str(), strerror(errno)); rc = 1; continue; }
        setvbuf(fo, nullptr, _IOFBF, 1 << 22);
        fwrite(header.data(), 1, header.size(), fo);
        fputc('\n', fo);

        std::vector<double> Beta, seBeta, Tstat, varT, Beta_c, seBeta_c, Tstat_c, varT_c;
        std::vector<double> AF_case, AF_ctrl, Nch, Nche, Ncth, Nctt;
        std::vector<std::string> pval, pvalNA, pval_c, pvalNA_c;
        std::vector<char> isSPA;
        std::vector<uint32_t> N_case, N_ctrl, N;
        std::vector<double> ac, af, info;
        std::string buf;
        buf.reserve(1u << 23);

        size_t blk = 0;
        bool ok = true;
        while (!r.bad && (size_t)(r.end - r.p) > 4) {
            uint32_t peek; memcpy(&peek, r.p, 4);
            if (peek == END_MAGIC) break;
            if (blk >= mblocks.size()) { fprintf(stderr, "%s: more trait blocks than marker blocks\n", in.c_str()); ok = false; break; }
            const MarkerBlock& MB = mblocks[blk++];
            uint32_t m = r.u32();
            if (m != BLK_MAGIC) { fprintf(stderr, "%s: corrupt block\n", in.c_str()); ok = false; break; }
            const size_t n = r.u32();
            if (n != MB.nRows) { fprintf(stderr, "%s: block row count disagrees with the marker file\n", in.c_str()); ok = false; break; }
            const uint32_t flags = r.u32();

            std::vector<uint8_t> present;
            if (flags & F_HAS_PRESENT) {
                present.assign(n, 1);
                const uint8_t* q = r.take(n);
                if (q) memcpy(present.data(), q, n);
            }
            ac = MB.ac; af = MB.af; info = MB.info;
            if (flags & F_OVERRIDE_AC)   get_col_f64(r, ac, n);
            if (flags & F_OVERRIDE_AF)   get_col_f64(r, af, n);
            if (flags & F_OVERRIDE_MISS) get_col_f64(r, info, n);

            for (uint32_t i = 0; i < nc; i++) {
                switch (cols[i]) {
                    case C_BETA:   get_col_f64 (r, Beta, n); break;
                    case C_SE:     get_col_f64 (r, seBeta, n); break;
                    case C_TSTAT:  get_col_f64 (r, Tstat, n); break;
                    case C_VAR:    get_col_f64 (r, varT, n); break;
                    case C_PVAL:   get_col_pval(r, pval, n); break;
                    case C_PVALNA: get_col_pval(r, pvalNA, n); break;
                    case C_ISSPA:  get_col_pod (r, isSPA, n); break;
                    case C_BETA_C: get_col_f64 (r, Beta_c, n); break;
                    case C_SE_C:   get_col_f64 (r, seBeta_c, n); break;
                    case C_TSTAT_C:get_col_f64 (r, Tstat_c, n); break;
                    case C_VAR_C:  get_col_f64 (r, varT_c, n); break;
                    case C_PVAL_C: get_col_pval(r, pval_c, n); break;
                    case C_PVALNA_C: get_col_pval(r, pvalNA_c, n); break;
                    case C_AFCASE: get_col_f64 (r, AF_case, n); break;
                    case C_AFCTRL: get_col_f64 (r, AF_ctrl, n); break;
                    case C_NCASE:  get_col_pod (r, N_case, n); break;
                    case C_NCTRL:  get_col_pod (r, N_ctrl, n); break;
                    case C_NCASEHOM: get_col_f64(r, Nch, n); break;
                    case C_NCASEHET: get_col_f64(r, Nche, n); break;
                    case C_NCTRLHOM: get_col_f64(r, Ncth, n); break;
                    case C_NCTRLHET: get_col_f64(r, Nctt, n); break;
                    case C_N:      get_col_pod (r, N, n); break;
                    default: fprintf(stderr, "%s: unknown column code %u\n", in.c_str(), cols[i]); ok = false;
                }
                if (!ok) break;
            }
            if (!ok || r.bad) { ok = false; break; }

            // Rows the writer marked absent had pval "NA" in the text, and
            // format_text skips exactly those.
            if (!present.empty())
                for (size_t k = 0; k < n; k++) if (!present[k]) pval[k] = "NA";

            outfast::MarkerCols MC;
            MC.chr = &MB.chr; MC.pos = &MB.pos; MC.mid = &MB.mid;
            MC.ref = &MB.ref; MC.alt = &MB.alt;
            outfast::TraitCols TC;
            TC.altCounts = &ac; TC.altFreq = &af;
            TC.imputeInfo = &info; TC.missingRate = &info;
            TC.Beta = &Beta; TC.seBeta = &seBeta; TC.Tstat = &Tstat; TC.varT = &varT;
            TC.pval = &pval; TC.pvalNA = &pvalNA; TC.isSPAConverge = &isSPA;
            TC.Beta_c = &Beta_c; TC.seBeta_c = &seBeta_c;
            TC.Tstat_c = &Tstat_c; TC.varT_c = &varT_c;
            TC.pval_c = &pval_c; TC.pvalNA_c = &pvalNA_c;
            TC.AF_case = &AF_case; TC.AF_ctrl = &AF_ctrl;
            TC.N_case = &N_case; TC.N_ctrl = &N_ctrl;
            TC.N_case_hom = &Nch; TC.N_case_het = &Nche;
            TC.N_ctrl_hom = &Ncth; TC.N_ctrl_het = &Nctt;
            TC.N = &N;

            buf.clear();
            outfast::format_text(buf, meta, isImputation, MC, TC, n);
            fwrite(buf.data(), 1, buf.size(), fo);
        }
        fclose(fo);
        if (!ok || r.bad) { rc = 1; continue; }
#ifdef _OPENMP
#pragma omp critical
#endif
        fprintf(stderr, "%s -> %s\n", in.c_str(), outPath.c_str());
    }
    return rc;
}
