// pgen_lib.hpp -- the part of pgenlib (plink-ng, third_party/pgenlib, LGPL v3)
// the PGEN reader needs, behind an interface with no plink2 types in it, so
// the pgenlib headers (and their macros) stay out of every other translation
// unit.
//
// One File per open .pgen (shared, read-only after open). One Reader per
// thread: a Reader owns its FILE*, its record buffer and pgenlib's LD cache,
// so readers on different threads never touch each other. Records may be
// LD-compressed (stored as a difflist against an earlier variant); pgenlib
// resolves that inside the Reader, which is why reading variants in
// increasing order on one Reader is cheapest.
#ifndef SAIGE_PGEN_LIB_HPP
#define SAIGE_PGEN_LIB_HPP

#include <cstdint>
#include <string>

namespace pgenlib_glue {

struct File;
struct Reader;

struct FileFacts {
    int      mode = 0;               // storage mode byte (0x01, 0x02, 0x10, 0x11, ...)
    uint32_t rawSampleCt = 0;
    uint32_t rawVariantCt = 0;
    bool     ldCompression = false;  // some record is a difflist against an earlier variant
    bool     dosage = false;         // some variant stores dosages
    bool     dosagePhase = false;
    bool     hardcallPhase = false;  // phase only: hard calls are unchanged
    bool     multiallelic = false;   // some variant has more than 2 alleles
};

// Open the .pgen. rawSampleCt / rawVariantCt are the .psam / .pvar counts; a
// mode-0x01 file (the .bed layout) needs both, the other modes check them
// against the header when they are not UINT32_MAX. nullptr + err on failure.
File* openFile(const std::string& path, uint32_t rawSampleCt, uint32_t rawVariantCt,
               FileFacts& facts, std::string& err);
void closeFile(File* f);

// Thread-safe (a mutex covers the one step that touches the shared File).
Reader* newReader(File* f, std::string& err);
void freeReader(Reader* r);   // needs nothing from the File

// ALT (allele 1) dosage of every sample, in file order, as upstream SAIGE's
// PgenClass::Read: hard calls 0 / 1 / 2, NaN for missing, a stored dosage
// d / 16384. out has rawSampleCt entries.
bool readAltDosage(Reader* r, uint32_t vidx, double* out, std::string& err);

// The hard calls of every sample, in file order, re-encoded as a PLINK 1 .bed
// row with ALT as A1:  00 = hom ALT, 01 = missing, 10 = het, 11 = hom REF,
// sample i in bits 2*(i&3) of byte i>>2, (rawSampleCt+3)/4 bytes, bits past
// the last sample zero. Only for files without dosages and multiallelic
// variants (the caller checks FileFacts).
bool readBedRow(Reader* r, uint32_t vidx, unsigned char* out, std::string& err);

}  // namespace pgenlib_glue

#endif
