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
// Measured after the change, same sweep as before (10^6 markers x 50,000
// samples, V100, fp64, page cache dropped before every run, machine otherwise
// idle; logs/gpu_step2/writer/bench3.out):
//
//   output write (s)     P=1    P=8    P=32   P=128
//     before             4.7   36.5   159.3   614.2
//     text, this file    3.6    6.4    26.5   106.4
//     sgs                0.4    0.6     2.8     6.9
//
//   end to end (s)       P=1    P=8    P=32   P=128
//     before            66.5  100.9   238.0   737.0
//     text              66.6   71.2   109.1   265.6
//     sgs               66.4   68.1    84.1   165.1
//     TorchGWAS2 GPU    96.6   99.1   107.6   139.9
//
//   output bytes       text -> sgs:  94 -> 72 MB, 752 -> 352 MB,
//                      3.01 -> 1.31 GB, 12.03 -> 5.15 GB (2.33x at P=128)
//
// ------------------------------------------------------------- fp32 ------
// sgsPrecision: fp32 (default fp64) narrows every floating-point column to
// float. Measured on the same sweep (tools/bench_sgs_precision.sh, run
// 2026-09-18, logs/gpu_step2/writer/f32/bench4.out):
//
//                        P=1     P=8     P=32    P=128
//   .sgs bytes  fp64   72.2 MB  352 MB  1.31 GB  5.15 GB
//               fp32   44.2 MB  184 MB  0.66 GB  2.58 GB   (1.63x .. 1.99x)
//   output write (s)   0.42/0.36  1.32/0.65  2.86/1.90  7.17/6.54   fp64/fp32
//   end to end   (s)  66.48/66.50 68.48/68.47 84.20/83.68 169.4/171.7
//   sgs2txt      (s)   4.9/4.4   9.6/8.0   34.0/30.6  127.3/122.9
//
// The file halves and the wall clock does not move. At P=128 writing is 7.2 s
// of a 169 s run, and halving the bytes gives back 0.6 s of it, because the
// fp32 writer has to convert a column before it can write it where the fp64
// writer memcpy's it. The end-to-end column is inside the run-to-run spread
// (three fp64 P=128 runs came out 165.1, 169.4, and the fp32 one 171.7).
//
// What it costs, counted exactly over the P=128 run (tools/sgs_fidelity.cpp,
// cross-checked field for field against a real fp32 run at P=8):
//
//   column            fields         differ     max |dx/x|
//   AC_Allele2     128,000,000            0     -           AC is an integer
//   AF_Allele2     128,000,000            0     -           AC/1e5: 6 digits
//   MissingRate    128,000,000            0     -           constant 0 here
//   BETA           128,000,000    1,104,895     1e-5
//   SE             128,000,000    1,178,101     1e-5
//   Tstat          128,000,000    1,079,021     1e-5
//   var            128,000,000    1,012,279     1e-5
//   all %.6g       896,000,000    4,374,296     0.4882%
//   p.value        128,000,000          470     max |d(-log10 p)| 4.45e-8
//
// 3.4% of ROWS are not byte-identical. No p-value crosses 5e-8 -- but this data
// is a null simulation whose smallest p is 9.5e-9, and that is the whole
// problem with the number. float's smallest normal is 1.2e-38 and its smallest
// subnormal 1.4e-45, while the score test prints "%.6E" all the way down to
// ~1e-308, so a real hit at p = 1e-50 does not round, it VANISHES. Forced with
// tests/make_extreme_model.py (tools/check_sgs_fp32_extreme.sh): residuals
// scaled by 6 put 147 of 5,000 markers below 1.2e-38 and the fp32 file printed
// 0.000000E+00 for 81 of them; scaled by 15, 1,665 of 5,000 markers (33%) came
// back as 0.000000E+00 where the text run had p down to 5.9E-716. Ironically
// the p-values SAIGE already prints in the "%.1fE%d" underflow form survive:
// those are kept verbatim in the exception list.
//
// So fp32 buys 2x on disk, ~0 on the clock, and costs the one property the
// format was built to have. Default stays fp64.
//
// Round trip verified by tools/sgs2txt + cmp at every P: 128 traits x 1,000,000
// markers = 128,000,000 rows byte-identical at P=128, and on a file with 2%
// missing calls plus 50 markers over maxMissRate so the QC-dropped rows are
// exercised too.
//
// What is left at P=128 with sgs: 33.9 s before the main loop (28.7 s of it
// reading the 128 null models), 50.9 s read+QC+stage (the 12.5 GB .bed at the
// disk's ~207 MB/s), 40.5 s tail+finalize (the per-pair p-value formatting and
// chi-square tail), 20.0 s on the device, 6.9 s writing. The writer is no
// longer the thing to fix.
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
    // for the shared marker block. storeF32 narrows every floating-point column
    // to float (sgsPrecision: fp32), which halves the file and gives up the
    // byte-identical round trip -- see sgs_format.hpp for the measured cost.
    // Returns false and fills err on any failure.
    bool open(const std::vector<TraitMeta>& metas, bool isImputation,
              bool storeF32, std::string& err);
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
