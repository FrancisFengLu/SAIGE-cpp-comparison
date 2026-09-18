// out_fast.cpp -- see out_fast.hpp for why this exists and what was measured.

#include "out_fast.hpp"
#include "sgs_format.hpp"

#include <cerrno>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <unistd.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace SAIGE {
namespace outfast {

// ---------------------------------------------------------------- text ----
//
// libstdc++'s num_put<char>::do_put(double) builds the format "%.*g" from the
// stream's flags and precision and calls vsnprintf with the C locale, then
// widens the result. With the default flags (defaultfloat, precision 6) that is
// exactly snprintf("%.6g"). Everything else here is a verbatim copy of the
// character sequence writeOutfile_single used to push through operator<<.
static inline char* put_g(char* p, double v) {
    return p + std::snprintf(p, 40, "%.6g", v);
}
static inline char* put_u(char* p, uint32_t v) {
    return p + std::snprintf(p, 16, "%u", v);
}
static inline char* put_s(char* p, const std::string& s) {
    std::memcpy(p, s.data(), s.size());
    return p + s.size();
}

std::string header_line(const TraitMeta& meta, bool isImputation) {
    std::string h = "CHR\tPOS\tMarkerID\tAllele1\tAllele2\tAC_Allele2\tAF_Allele2\t";
    h += isImputation ? "imputationInfo\t" : "MissingRate\t";
    h += "BETA\tSE\tTstat\tvar\tp.value\t";
    const bool spa = (meta.traitType == "binary" || meta.traitType == "survival");
    if (spa) h += "p.value.NA\tIs.SPA\t";
    if (meta.isCondition) {
        h += "BETA_c\tSE_c\tTstat_c\tvar_c\tp.value_c\t";
        if (spa) h += "p.value.NA_c\t";
    }
    if (meta.traitType == "binary") {
        h += "AF_case\tAF_ctrl\tN_case\tN_ctrl";
        if (meta.isMoreOutput) h += "\tN_case_hom\tN_case_het\tN_ctrl_hom\tN_ctrl_het";
    } else if (meta.traitType == "quantitative") {
        h += "N";
    } else if (meta.traitType == "survival") {
        h += "AF_event\tAF_censor\tN_event\tN_censor";
        if (meta.isMoreOutput) h += "\tN_event_hom\tN_event_het\tN_censor_hom\tN_censor_het";
    }
    return h;
}

int format_text(std::string& buf, const TraitMeta& meta, bool isImputation,
                const MarkerCols& M, const TraitCols& T, std::size_t nRows)
{
    const bool isBinOrSurv = (meta.traitType == "binary" || meta.traitType == "survival");
    const bool isBin       = (meta.traitType == "binary");
    const bool isQuant     = (meta.traitType == "quantitative");
    const bool isSurv      = (meta.traitType == "survival");
    const bool cond        = meta.isCondition;
    const bool more        = meta.isMoreOutput;

    // Longest possible row: 22 numeric fields at <= 24 chars plus five marker
    // strings. Grown, not guessed, for the marker id.
    std::vector<char> row(1024);
    int numtest = 0;
    for (std::size_t k = 0; k < nRows; k++) {
        const std::string& pv = (*T.pval)[k];
        if (pv == "NA") continue;
        numtest++;

        // Room for every field this row can hold: 21 numeric fields at <= 24
        // characters plus their tabs fits in 1024, and the variable-length
        // ones -- the five marker strings and up to four p-value strings --
        // are added explicitly rather than assumed to be short.
        std::size_t want = 1024 + (*M.chr)[k].size() + (*M.pos)[k].size()
                         + (*M.mid)[k].size() + (*M.ref)[k].size()
                         + (*M.alt)[k].size() + pv.size();
        if (isBinOrSurv) want += (*T.pvalNA)[k].size();
        if (cond) {
            want += (*T.pval_c)[k].size();
            if (isBinOrSurv) want += (*T.pvalNA_c)[k].size();
        }
        if (row.size() < want) row.resize(want);
        char* p = row.data();

        p = put_s(p, (*M.chr)[k]); *p++ = '\t';
        p = put_s(p, (*M.pos)[k]); *p++ = '\t';
        p = put_s(p, (*M.mid)[k]); *p++ = '\t';
        p = put_s(p, (*M.ref)[k]); *p++ = '\t';
        p = put_s(p, (*M.alt)[k]); *p++ = '\t';
        p = put_g(p, (*T.altCounts)[k]); *p++ = '\t';
        p = put_g(p, (*T.altFreq)[k]);   *p++ = '\t';
        p = put_g(p, isImputation ? (*T.imputeInfo)[k] : (*T.missingRate)[k]); *p++ = '\t';
        p = put_g(p, (*T.Beta)[k]);   *p++ = '\t';
        p = put_g(p, (*T.seBeta)[k]); *p++ = '\t';
        p = put_g(p, (*T.Tstat)[k]);  *p++ = '\t';
        p = put_g(p, (*T.varT)[k]);   *p++ = '\t';
        p = put_s(p, pv);             *p++ = '\t';

        if (isBinOrSurv) {
            p = put_s(p, (*T.pvalNA)[k]); *p++ = '\t';
            const bool b = (*T.isSPAConverge)[k] != 0;
            p = put_s(p, b ? std::string("true") : std::string("false")); *p++ = '\t';
        }
        if (cond) {
            p = put_g(p, (*T.Beta_c)[k]);   *p++ = '\t';
            p = put_g(p, (*T.seBeta_c)[k]); *p++ = '\t';
            p = put_g(p, (*T.Tstat_c)[k]);  *p++ = '\t';
            p = put_g(p, (*T.varT_c)[k]);   *p++ = '\t';
            p = put_s(p, (*T.pval_c)[k]);   *p++ = '\t';
            if (isBinOrSurv) { p = put_s(p, (*T.pvalNA_c)[k]); *p++ = '\t'; }
        }
        if (isBinOrSurv) {
            p = put_g(p, (*T.AF_case)[k]); *p++ = '\t';
            p = put_g(p, (*T.AF_ctrl)[k]); *p++ = '\t';
            p = put_u(p, (*T.N_case)[k]);  *p++ = '\t';
            p = put_u(p, (*T.N_ctrl)[k]);
            if (more) {
                *p++ = '\t'; p = put_g(p, (*T.N_case_hom)[k]);
                *p++ = '\t'; p = put_g(p, (*T.N_case_het)[k]);
                *p++ = '\t'; p = put_g(p, (*T.N_ctrl_hom)[k]);
                *p++ = '\t'; p = put_g(p, (*T.N_ctrl_het)[k]);
            }
            *p++ = '\n';
        } else if (isQuant) {
            p = put_u(p, (*T.N)[k]); *p++ = '\n';
        }
        // A survival trait takes the isBinOrSurv branch above; the old writer
        // emitted no newline for any other trait type, and neither does this.
        (void)isBin; (void)isSurv;
        buf.append(row.data(), (std::size_t)(p - row.data()));
    }
    return numtest;
}

// ---------------------------------------------------------------- .sgs ----
using namespace SAIGE::sgs;

static bool write_all(int fd, const void* data, std::size_t n, std::string& err) {
    const char* p = (const char*)data;
    while (n) {
        ssize_t w = ::write(fd, p, n);
        if (w < 0) {
            if (errno == EINTR) continue;
            err = std::string("write failed: ") + std::strerror(errno);
            return false;
        }
        p += w; n -= (std::size_t)w;
    }
    return true;
}

template <typename T>
static void put_col_pod(std::string& b, const T* v, std::size_t n) {
    bool cst = (n > 0);
    for (std::size_t i = 1; i < n; i++)
        if (std::memcmp(&v[i], &v[0], sizeof(T)) != 0) { cst = false; break; }
    if (cst) { put_u8(b, E_CONST); put_bytes(b, &v[0], sizeof(T)); }
    else     { put_u8(b, E_RAW);   put_bytes(b, v, n * sizeof(T)); }
}

static void put_col_str(std::string& b, const std::vector<std::string>& v, std::size_t n) {
    bool cst = (n > 0);
    for (std::size_t i = 1; i < n; i++) if (v[i] != v[0]) { cst = false; break; }
    if (cst) { put_u8(b, E_CONST); put_sstr(b, v[0]); }
    else {
        put_u8(b, E_RAW);
        for (std::size_t i = 0; i < n; i++) put_sstr(b, v[i]);
    }
}

// "d.dddddd" 'E' ('+'|'-') two-or-more digits, decimal exponent within +-290.
// Inside that shape the string carries 7 significant digits, which round-trip
// exactly through a normal double, so storing strtod(s) and re-running
// sprintf("%.6E") gives back this very string. Everything outside it -- "NA",
// the "%.1fE%d" underflow form, subnormal exponents -- is kept verbatim.
static bool canonical_6E(const std::string& s, double& out) {
    const std::size_t n = s.size();
    if (n < 12 || n > 13) return false;
    if (s[0] < '0' || s[0] > '9') return false;
    if (s[1] != '.') return false;
    for (int i = 2; i < 8; i++) if (s[i] < '0' || s[i] > '9') return false;
    if (s[8] != 'E') return false;
    if (s[9] != '+' && s[9] != '-') return false;
    int e = 0;
    for (std::size_t i = 10; i < n; i++) {
        if (s[i] < '0' || s[i] > '9') return false;
        e = e * 10 + (s[i] - '0');
    }
    if (e > 290) return false;
    out = std::strtod(s.c_str(), nullptr);
    return true;
}

static void put_col_pval(std::string& b, const std::vector<std::string>& v, std::size_t n) {
    put_u8(b, E_PVAL);
    std::vector<double> d(n);
    std::string exc;
    uint32_t nexc = 0;
    for (std::size_t i = 0; i < n; i++) {
        double x;
        if (canonical_6E(v[i], x)) d[i] = x;
        else {
            d[i] = std::nan("");
            put_u32(exc, (uint32_t)i);
            put_sstr(exc, v[i]);
            nexc++;
        }
    }
    put_bytes(b, d.data(), n * sizeof(double));
    put_u32(b, nexc);
    b.append(exc);
}

// The per-trait column list, in text order.
static std::vector<uint8_t> col_list(const TraitMeta& meta) {
    const bool spa = (meta.traitType == "binary" || meta.traitType == "survival");
    std::vector<uint8_t> c{C_BETA, C_SE, C_TSTAT, C_VAR, C_PVAL};
    if (spa) { c.push_back(C_PVALNA); c.push_back(C_ISSPA); }
    if (meta.isCondition) {
        c.push_back(C_BETA_C); c.push_back(C_SE_C); c.push_back(C_TSTAT_C);
        c.push_back(C_VAR_C);  c.push_back(C_PVAL_C);
        if (spa) c.push_back(C_PVALNA_C);
    }
    if (spa) {
        c.push_back(C_AFCASE); c.push_back(C_AFCTRL);
        c.push_back(C_NCASE);  c.push_back(C_NCTRL);
        if (meta.isMoreOutput) {
            c.push_back(C_NCASEHOM); c.push_back(C_NCASEHET);
            c.push_back(C_NCTRLHOM); c.push_back(C_NCTRLHET);
        }
    } else if (meta.traitType == "quantitative") {
        c.push_back(C_N);
    }
    return c;
}

struct SgsSink::Impl {
    int              markerFd = -1;
    std::vector<int> traitFd;
    std::vector<TraitMeta> metas;
    std::vector<std::vector<uint8_t> > cols;
    bool isImputation = false;
    uint64_t nMarkers = 0;
    std::vector<uint64_t> nEmitted;
    std::vector<std::string> bufs;       // one scratch buffer per trait
    std::string markerBuf;
};

SgsSink::~SgsSink() { delete m_impl; }

bool SgsSink::open(const std::vector<TraitMeta>& metas, bool isImputation,
                   std::string& err)
{
    if (metas.empty()) { err = "no traits"; return false; }
    m_impl = new Impl();
    m_impl->metas = metas;
    m_impl->isImputation = isImputation;
    m_impl->nEmitted.assign(metas.size(), 0);
    m_impl->bufs.resize(metas.size());
    m_markerPath = metas[0].outFile + ".markers.sgs";

    m_impl->markerFd = ::open(m_markerPath.c_str(),
                              O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (m_impl->markerFd < 0) {
        err = "cannot open " + m_markerPath + ": " + std::strerror(errno);
        return false;
    }
    {
        std::string h;
        put_bytes(h, MAGIC_MARKER, 8);
        put_u32(h, VERSION);
        put_u32(h, isImputation ? 1u : 0u);
        if (!write_all(m_impl->markerFd, h.data(), h.size(), err)) return false;
        m_bytes += h.size();
    }

    m_impl->traitFd.assign(metas.size(), -1);
    m_impl->cols.resize(metas.size());
    for (std::size_t t = 0; t < metas.size(); t++) {
        const std::string path = metas[t].outFile + ".sgs";
        int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) { err = "cannot open " + path + ": " + std::strerror(errno); return false; }
        m_impl->traitFd[t] = fd;
        m_impl->cols[t] = col_list(metas[t]);

        std::string h;
        put_bytes(h, MAGIC_TRAIT, 8);
        put_u32(h, VERSION);
        put_u32(h, isImputation ? 1u : 0u);
        put_str(h, header_line(metas[t], isImputation));
        put_str(h, metas[t].name);
        put_str(h, metas[t].traitType);
        put_u8(h, metas[t].isCondition  ? 1 : 0);
        put_u8(h, metas[t].isMoreOutput ? 1 : 0);
        put_str(h, m_markerPath);
        put_str(h, metas[t].outFile);     // the text path this converts back to
        put_u32(h, (uint32_t)m_impl->cols[t].size());
        for (std::size_t i = 0; i < m_impl->cols[t].size(); i++)
            put_u8(h, m_impl->cols[t][i]);
        if (!write_all(fd, h.data(), h.size(), err)) return false;
        m_bytes += h.size();
    }
    m_open = true;
    return true;
}

bool SgsSink::writeChunk(const MarkerCols& M, const std::vector<TraitCols>& cols,
                         std::size_t nRows, int nThreads,
                         std::vector<int>& numtest, std::string& err)
{
    if (!m_open) { err = "sgs sink is not open"; return false; }
    Impl& I = *m_impl;
    const std::size_t P = I.metas.size();
    if (cols.size() != P) { err = "trait count mismatch"; return false; }
    numtest.assign(P, 0);
    if (nRows == 0) return true;

    // ---- marker block: written once, not P times ----
    const std::vector<double>& mInfo =
        I.isImputation ? *cols[0].imputeInfo : *cols[0].missingRate;
    {
        std::string& b = I.markerBuf;
        b.clear();
        put_u32(b, BLK_MAGIC);
        put_u32(b, (uint32_t)nRows);
        put_col_str(b, *M.chr, nRows);
        put_col_str(b, *M.pos, nRows);
        put_col_str(b, *M.mid, nRows);
        put_col_str(b, *M.ref, nRows);
        put_col_str(b, *M.alt, nRows);
        put_col_pod(b, cols[0].altCounts->data(), nRows);
        put_col_pod(b, cols[0].altFreq->data(), nRows);
        put_col_pod(b, mInfo.data(), nRows);
        if (!write_all(I.markerFd, b.data(), b.size(), err)) return false;
        m_bytes += b.size();
        I.nMarkers += nRows;
    }

    // ---- one block per trait, formatted and written in parallel ----
    std::vector<std::string> perr(P);
    const int nt = (nThreads > 0 ? nThreads : 1);
#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic) num_threads(nt)
#endif
    for (long tt = 0; tt < (long)P; tt++) {
        const std::size_t t = (std::size_t)tt;
        const TraitCols& C = cols[t];
        std::string& b = I.bufs[t];
        b.clear();

        // present mask: the text writer skips rows whose pval is "NA".
        std::vector<uint8_t> present(nRows, 1);
        int emitted = 0;
        bool allPresent = true;
        for (std::size_t k = 0; k < nRows; k++) {
            if ((*C.pval)[k] == "NA") { present[k] = 0; allPresent = false; }
            else emitted++;
        }
        numtest[t] = emitted;

        const std::vector<double>& tInfo =
            I.isImputation ? *C.imputeInfo : *C.missingRate;
        uint32_t flags = 0;
        if (!allPresent) flags |= F_HAS_PRESENT;
        if (t != 0) {
            if (std::memcmp(C.altCounts->data(), cols[0].altCounts->data(),
                            nRows * sizeof(double)) != 0) flags |= F_OVERRIDE_AC;
            if (std::memcmp(C.altFreq->data(), cols[0].altFreq->data(),
                            nRows * sizeof(double)) != 0) flags |= F_OVERRIDE_AF;
            if (std::memcmp(tInfo.data(), mInfo.data(),
                            nRows * sizeof(double)) != 0) flags |= F_OVERRIDE_MISS;
        }

        put_u32(b, BLK_MAGIC);
        put_u32(b, (uint32_t)nRows);
        put_u32(b, flags);
        if (flags & F_HAS_PRESENT) put_bytes(b, present.data(), nRows);
        if (flags & F_OVERRIDE_AC)   put_col_pod(b, C.altCounts->data(), nRows);
        if (flags & F_OVERRIDE_AF)   put_col_pod(b, C.altFreq->data(), nRows);
        if (flags & F_OVERRIDE_MISS) put_col_pod(b, tInfo.data(), nRows);

        for (std::size_t i = 0; i < I.cols[t].size(); i++) {
            switch (I.cols[t][i]) {
                case C_BETA:   put_col_pod (b, C.Beta->data(), nRows); break;
                case C_SE:     put_col_pod (b, C.seBeta->data(), nRows); break;
                case C_TSTAT:  put_col_pod (b, C.Tstat->data(), nRows); break;
                case C_VAR:    put_col_pod (b, C.varT->data(), nRows); break;
                case C_PVAL:   put_col_pval(b, *C.pval, nRows); break;
                case C_PVALNA: put_col_pval(b, *C.pvalNA, nRows); break;
                case C_ISSPA:  put_col_pod (b, C.isSPAConverge->data(), nRows); break;
                case C_BETA_C: put_col_pod (b, C.Beta_c->data(), nRows); break;
                case C_SE_C:   put_col_pod (b, C.seBeta_c->data(), nRows); break;
                case C_TSTAT_C:put_col_pod (b, C.Tstat_c->data(), nRows); break;
                case C_VAR_C:  put_col_pod (b, C.varT_c->data(), nRows); break;
                case C_PVAL_C: put_col_pval(b, *C.pval_c, nRows); break;
                case C_PVALNA_C: put_col_pval(b, *C.pvalNA_c, nRows); break;
                case C_AFCASE: put_col_pod (b, C.AF_case->data(), nRows); break;
                case C_AFCTRL: put_col_pod (b, C.AF_ctrl->data(), nRows); break;
                case C_NCASE:  put_col_pod (b, C.N_case->data(), nRows); break;
                case C_NCTRL:  put_col_pod (b, C.N_ctrl->data(), nRows); break;
                case C_NCASEHOM: put_col_pod(b, C.N_case_hom->data(), nRows); break;
                case C_NCASEHET: put_col_pod(b, C.N_case_het->data(), nRows); break;
                case C_NCTRLHOM: put_col_pod(b, C.N_ctrl_hom->data(), nRows); break;
                case C_NCTRLHET: put_col_pod(b, C.N_ctrl_het->data(), nRows); break;
                case C_N:      put_col_pod (b, C.N->data(), nRows); break;
                default: perr[t] = "unknown column code"; break;
            }
        }
        std::string e2;
        if (!write_all(I.traitFd[t], b.data(), b.size(), e2)) perr[t] = e2;
        I.nEmitted[t] += (uint64_t)emitted;
#ifdef _OPENMP
#pragma omp atomic
#endif
        m_bytes += b.size();
    }
    for (std::size_t t = 0; t < P; t++)
        if (!perr[t].empty()) { err = perr[t]; return false; }
    return true;
}

bool SgsSink::close(std::string& err)
{
    if (!m_open) return true;
    Impl& I = *m_impl;
    {
        std::string b;
        put_u32(b, END_MAGIC);
        put_u64(b, I.nMarkers);
        if (!write_all(I.markerFd, b.data(), b.size(), err)) return false;
        m_bytes += b.size();
        ::close(I.markerFd); I.markerFd = -1;
    }
    for (std::size_t t = 0; t < I.traitFd.size(); t++) {
        std::string b;
        put_u32(b, END_MAGIC);
        put_u64(b, I.nMarkers);
        put_u64(b, I.nEmitted[t]);
        if (!write_all(I.traitFd[t], b.data(), b.size(), err)) return false;
        m_bytes += b.size();
        ::close(I.traitFd[t]); I.traitFd[t] = -1;
    }
    m_open = false;
    return true;
}

}  // namespace outfast
}  // namespace SAIGE
