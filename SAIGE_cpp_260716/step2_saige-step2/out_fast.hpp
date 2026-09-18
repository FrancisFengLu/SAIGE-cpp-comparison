// out_fast.hpp -- step-2 single-assoc output: a fast byte-identical text
// writer, and an optional binary columnar format (.sgs).
//
// Why this file exists (measured 2026-09-18, 10^6 markers x 50,000 samples,
// quantitative, GPU path):
//
//   P            1      8     32    128
//   output write 4.7   36.5  159.3  614.2   seconds
//
// At P=128 the writer was 84% of the whole run.  A microbenchmark of the old
// writer on one real 1M-row output file (logs/gpu_step2/writer/wbench.cpp):
//
//   ostream -> file            3.905 s
//   ostream -> /dev/null       3.772 s      <- so disk was 0.13 s; 96% was CPU
//   snprintf into a buffer     2.901 s + 0.446 s write()
//   memcpy of pre-built rows   0.101 s + 0.437 s write()
//   one snprintf("%.6g")         332 ns     x 7 doubles per row = 2.3 s
//   one ostream << double        386 ns
//
// So the cost was number-to-string formatting, not IO.  Two answers, both here:
//
//   1. format_text() -- the same rows, formatted with snprintf into one big
//      buffer and handed to the stream as a single write().  Byte-identical by
//      construction: libstdc++'s num_put for double is itself vsnprintf with
//      "%.*g" in the C locale, so "%.6g" reproduces `stream << double` exactly.
//      Combined with running the P traits in parallel this takes the P=128
//      writer from 614 s to the disk's floor.
//
//   2. SgsSink -- .sgs, a binary columnar format that never formats a number.
//      CHR/POS/MarkerID/Allele1/Allele2/AC/AF/MissingRate are per marker, not
//      per trait, so they are stored ONCE instead of P times (at P=128 the text
//      repeats them 128 times: 4.8 GB of the 12 GB).  Losslessly convertible
//      back to the text output by tools/sgs2txt, which re-uses format_text()
//      below, so the round trip is byte-identical by construction.
//
// Rejected after measuring: zstd on the f64 columns (1.17x, 1.29x with a byte
// shuffle -- not worth the CPU), and a hand-rolled %g (the text path is disk
// bound once the traits are written in parallel, so a faster formatter buys
// nothing).

#ifndef SAIGE_OUT_FAST_HPP
#define SAIGE_OUT_FAST_HPP

#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>

#include "saige_mt.hpp"

namespace SAIGE {
namespace outfast {

// Marker-level columns: one set per chunk, shared by every trait.
struct MarkerCols {
    const std::vector<std::string>* chr = nullptr;
    const std::vector<std::string>* pos = nullptr;
    const std::vector<std::string>* mid = nullptr;
    const std::vector<std::string>* ref = nullptr;
    const std::vector<std::string>* alt = nullptr;
};

// One trait's columns for one chunk. Pointers into the caller's vectors; the
// writer never owns them. altCounts/altFreq/imputeInfo/missingRate live here
// and not in MarkerCols because with different sample sets they are per trait
// (design 4.7); the .sgs writer stores them once and only falls back to a
// per-trait override when a trait's values actually differ.
struct TraitCols {
    const std::vector<double>*      altCounts   = nullptr;
    const std::vector<double>*      altFreq     = nullptr;
    const std::vector<double>*      imputeInfo  = nullptr;
    const std::vector<double>*      missingRate = nullptr;
    const std::vector<double>*      Beta        = nullptr;
    const std::vector<double>*      seBeta      = nullptr;
    const std::vector<double>*      Tstat       = nullptr;
    const std::vector<double>*      varT        = nullptr;
    const std::vector<std::string>* pval        = nullptr;
    const std::vector<std::string>* pvalNA      = nullptr;
    const std::vector<char>*        isSPAConverge = nullptr;   // 0/1
    const std::vector<double>*      Beta_c      = nullptr;
    const std::vector<double>*      seBeta_c    = nullptr;
    const std::vector<double>*      Tstat_c     = nullptr;
    const std::vector<double>*      varT_c      = nullptr;
    const std::vector<std::string>* pval_c      = nullptr;
    const std::vector<std::string>* pvalNA_c    = nullptr;
    const std::vector<double>*      AF_case     = nullptr;
    const std::vector<double>*      AF_ctrl     = nullptr;
    const std::vector<uint32_t>*    N_case      = nullptr;
    const std::vector<uint32_t>*    N_ctrl      = nullptr;
    const std::vector<double>*      N_case_hom  = nullptr;
    const std::vector<double>*      N_ctrl_het  = nullptr;
    const std::vector<double>*      N_case_het  = nullptr;
    const std::vector<double>*      N_ctrl_hom  = nullptr;
    const std::vector<uint32_t>*    N           = nullptr;
};

// The header line openOutfile_single writes, without the trailing newline.
std::string header_line(const TraitMeta& meta, bool isImputation);

// Append rows [0, nRows) to buf exactly as the ostream loop in
// writeOutfile_single would have written them. Rows whose pval is "NA" are
// skipped, as there. Returns the number of rows emitted.
int format_text(std::string& buf, const TraitMeta& meta, bool isImputation,
                const MarkerCols& M, const TraitCols& T, std::size_t nRows);

// ---------------------------------------------------------------- .sgs ----
// One marker file plus one file per trait. Chunk blocks are appended in the
// same order to all of them, so the converter zips block k of the marker file
// with block k of every trait file.
class SgsSink {
public:
    ~SgsSink();
    // metas[t].outFile + ".sgs" per trait; metas[0].outFile + ".markers.sgs"
    // for the shared marker block. Returns false and fills err on any failure.
    bool open(const std::vector<TraitMeta>& metas, bool isImputation,
              std::string& err);
    bool isOpen() const { return m_open; }
    const std::string& markerPath() const { return m_markerPath; }

    // One chunk. cols.size() must equal the trait count. numtest[t] receives
    // the rows emitted for trait t (pval != "NA"), matching the text writer.
    // nThreads > 1 formats the per-trait blocks in parallel.
    bool writeChunk(const MarkerCols& M, const std::vector<TraitCols>& cols,
                    std::size_t nRows, int nThreads, std::vector<int>& numtest,
                    std::string& err);
    bool close(std::string& err);

    // Bytes written so far, for the run's breakdown line.
    unsigned long long bytesWritten() const { return m_bytes; }

private:
    struct Impl;
    Impl* m_impl = nullptr;
    bool  m_open = false;
    std::string m_markerPath;
    unsigned long long m_bytes = 0;
};

}  // namespace outfast
}  // namespace SAIGE

#endif  // SAIGE_OUT_FAST_HPP
