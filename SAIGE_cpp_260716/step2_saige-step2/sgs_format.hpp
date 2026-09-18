// sgs_format.hpp -- on-disk layout of the .sgs binary output, shared by the
// writer (out_fast.cpp) and the converter (tools/sgs2txt.cpp).
//
// Design in one paragraph: two kinds of file, both a header followed by a
// sequence of chunk blocks appended in write order and closed by a trailer.
// The MARKER file carries the columns that are a property of the marker and
// not of the trait (CHR, POS, MarkerID, Allele1, Allele2, AC_Allele2,
// AF_Allele2, MissingRate/imputationInfo) -- written once, not P times. Each
// TRAIT file carries only that trait's own columns. Block k of the marker file
// and block k of every trait file describe the same markers, in the same order,
// so the converter simply zips them.
//
// Every column is length-prefixed and tagged with a one-byte encoding, so a
// column that is constant over a block (MissingRate = 0, N = 50000, the
// all-rows-present mask) costs one value instead of nRows.
//
// Numbers are stored as the IEEE doubles the program computed; nothing is
// rounded on the way out. p-value STRINGS are stored as the double they parse
// back to whenever the string is in the canonical "%.6E" shape and its decimal
// exponent is inside +-290 (7 significant digits round-trip exactly through a
// normal double, so re-running the same sprintf reproduces the byte); anything
// else -- "NA", the "%.1fE%d" underflow form, subnormal exponents -- is kept
// verbatim in a per-block exception list. So the text is reproducible byte for
// byte, which tools/sgs2txt verifies with cmp.
//
// sgsPrecision: fp32 halves every floating-point column instead (E_RAW32 /
// E_CONST32 / E_PVAL32) and GIVES THAT GUARANTEE UP. A float carries ~7.2
// decimal digits, the text prints 6 ("%.6g") or 7 ("%.6E"), and the narrowing
// moves a value across the last printed digit's rounding boundary often enough
// to matter (0.49% of fields over a 10^6 x 128 run), and, worse, float cannot
// hold a p-value below 1.2e-38 at all while the text prints them down to
// ~1e-308 -- out_fast.hpp has the counts. The integer and flag
// columns (N, N_case, N_ctrl, Is.SPA) are untouched either way. The
// encoding byte says which width a column is, so one reader handles both and
// files written before the option existed still read. Default is fp64.
//
// Little-endian, and the reader checks that. No attempt at cross-endian
// portability: this is a scratch format for one machine's pipeline.

#ifndef SAIGE_SGS_FORMAT_HPP
#define SAIGE_SGS_FORMAT_HPP

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

namespace SAIGE {
namespace sgs {

static const char   MAGIC_MARKER[8] = {'S','A','I','G','E','S','G','M'};
static const char   MAGIC_TRAIT [8] = {'S','A','I','G','E','S','G','T'};
static const uint32_t VERSION   = 1;
static const uint32_t BLK_MAGIC = 0x314B4C42u;   // "BLK1"
static const uint32_t END_MAGIC = 0x21444E45u;   // "END!"

// Per-trait column codes, in the order the text writer emits them.
enum ColCode : uint8_t {
    C_BETA = 1, C_SE, C_TSTAT, C_VAR, C_PVAL,
    C_PVALNA, C_ISSPA,
    C_BETA_C, C_SE_C, C_TSTAT_C, C_VAR_C, C_PVAL_C, C_PVALNA_C,
    C_AFCASE, C_AFCTRL, C_NCASE, C_NCTRL,
    C_NCASEHOM, C_NCASEHET, C_NCTRLHOM, C_NCTRLHET,
    C_N
};

// Column encodings. The 32-bit forms appear only under sgsPrecision: fp32, and
// only for columns the program computed as double; uint32 and flag columns keep
// their own width.
enum Enc : uint8_t {
    E_RAW     = 0,   // nRows values back to back
    E_CONST   = 1,   // one value, repeated
    E_PVAL    = 2,   // f64[nRows] + exception list (see above)
    E_RAW32   = 3,   // nRows floats
    E_CONST32 = 4,   // one float, repeated
    E_PVAL32  = 5    // f32[nRows] + exception list
};

// Header flag bits (the u32 after VERSION in both file kinds).
enum HdrFlag : uint32_t {
    H_IMPUTATION = 1u << 0,   // the info column is imputationInfo, not MissingRate
    H_F32        = 1u << 1    // floating-point columns stored as float
};

// Block flag bits.
enum BlkFlag : uint32_t {
    F_HAS_PRESENT   = 1u << 0,   // a u8 mask follows; otherwise every row is present
    F_OVERRIDE_AC   = 1u << 1,   // this trait's AC_Allele2 differs from the marker block
    F_OVERRIDE_AF   = 1u << 2,
    F_OVERRIDE_MISS = 1u << 3
};

// ---------------------------------------------------------------- put ----
inline void put_bytes(std::string& b, const void* p, std::size_t n) {
    b.append(static_cast<const char*>(p), n);
}
inline void put_u8 (std::string& b, uint8_t  v) { b.push_back((char)v); }
inline void put_u32(std::string& b, uint32_t v) { put_bytes(b, &v, 4); }
inline void put_u64(std::string& b, uint64_t v) { put_bytes(b, &v, 8); }
inline void put_f64(std::string& b, double   v) { put_bytes(b, &v, 8); }
inline void put_f32(std::string& b, float    v) { put_bytes(b, &v, 4); }
inline void put_str(std::string& b, const std::string& s) {
    put_u32(b, (uint32_t)s.size());
    b.append(s);
}
// Short string inside a column: one length byte, or 0xFF + u32 for the rare
// marker id longer than 254 bytes.
inline void put_sstr(std::string& b, const std::string& s) {
    if (s.size() < 255) { put_u8(b, (uint8_t)s.size()); }
    else                { put_u8(b, 255); put_u32(b, (uint32_t)s.size()); }
    b.append(s);
}

// ---------------------------------------------------------------- get ----
struct Reader {
    const uint8_t* p = nullptr;
    const uint8_t* end = nullptr;
    bool bad = false;
    Reader() {}
    Reader(const void* d, std::size_t n)
        : p((const uint8_t*)d), end((const uint8_t*)d + n) {}
    bool need(std::size_t n) {
        if (bad || (std::size_t)(end - p) < n) { bad = true; return false; }
        return true;
    }
    const uint8_t* take(std::size_t n) {
        if (!need(n)) return nullptr;
        const uint8_t* q = p; p += n; return q;
    }
    uint8_t  u8 () { const uint8_t* q = take(1); return q ? *q : 0; }
    uint32_t u32() { const uint8_t* q = take(4); uint32_t v = 0; if (q) std::memcpy(&v, q, 4); return v; }
    uint64_t u64() { const uint8_t* q = take(8); uint64_t v = 0; if (q) std::memcpy(&v, q, 8); return v; }
    double   f64() { const uint8_t* q = take(8); double   v = 0; if (q) std::memcpy(&v, q, 8); return v; }
    float    f32() { const uint8_t* q = take(4); float    v = 0; if (q) std::memcpy(&v, q, 4); return v; }
    std::string str() {
        uint32_t n = u32(); const uint8_t* q = take(n);
        return q ? std::string((const char*)q, n) : std::string();
    }
    std::string sstr() {
        uint32_t n = u8(); if (n == 255) n = u32();
        const uint8_t* q = take(n);
        return q ? std::string((const char*)q, n) : std::string();
    }
};

}  // namespace sgs
}  // namespace SAIGE

#endif  // SAIGE_SGS_FORMAT_HPP
