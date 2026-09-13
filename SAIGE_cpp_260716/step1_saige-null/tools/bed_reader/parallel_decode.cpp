#include "parallel_decode.hpp"

#include <algorithm>
#include <atomic>
#include <cstring>
#include <new>
#include <stdexcept>
#include <thread>

#include <sys/mman.h>

namespace saige {

namespace {

// Packed rows are staged in fixed-size blocks rather than one growing vector.
// Two reasons, both about peak RSS on a buffer that is the size of the whole
// packed genotype matrix:
//   · a growing vector reallocs, so old + new are briefly both resident;
//   · it can only be freed as a whole, so during pass 2 the entire staging
//     copy stayed live until the last row had been written to the final
//     PackedFlat — ~2× the matrix at peak.
// With blocks, pass 2 releases each block right after copying it, so
// (staging still held) + (final store written so far) stays ≈ one matrix.
// 64 rows ≈ 0.8 MB on mid, ≈ 7 MB at UKB shape.
constexpr std::size_t kRowsPerBlock = 64;

// Blocks are mmap'd, not new[]'d. glibc returns a freed large chunk to its own
// arena rather than to the kernel — and its mmap threshold is *dynamic*, it
// rises as soon as a big mmap'd chunk is freed — so a malloc-backed block
// stayed resident after release() and peak RSS was still ~2× the packed
// matrix. munmap hands the pages back immediately, which is the whole point of
// draining block by block.
class MappedBlock {
 public:
  explicit MappedBlock(std::size_t n) : n_(n) {
    void* p = ::mmap(nullptr, n_, PROT_READ | PROT_WRITE,
                     MAP_PRIVATE | MAP_ANONYMOUS, -1, 0);
    if (p == MAP_FAILED) throw std::bad_alloc();
    p_ = static_cast<unsigned char*>(p);
  }
  ~MappedBlock() { release(); }
  MappedBlock(MappedBlock&& o) noexcept : p_(o.p_), n_(o.n_) { o.p_ = nullptr; o.n_ = 0; }
  MappedBlock& operator=(MappedBlock&& o) noexcept {
    if (this != &o) { release(); p_ = o.p_; n_ = o.n_; o.p_ = nullptr; o.n_ = 0; }
    return *this;
  }
  MappedBlock(const MappedBlock&) = delete;
  MappedBlock& operator=(const MappedBlock&) = delete;

  void release() {
    if (p_) { ::munmap(p_, n_); p_ = nullptr; n_ = 0; }
  }
  unsigned char* data() const { return p_; }

 private:
  unsigned char* p_ = nullptr;
  std::size_t    n_ = 0;
};

struct PerThreadState {
  // Pass 1 outputs, indexed by (marker - chunk_lo):
  std::vector<MarkerStats> stats;
  std::vector<char>        passQC_flags;     // 0/1 — char because std::vector<bool> can't be written concurrently
  std::vector<char>        passVR_flags;     // 0/1 — variance-ratio pool
  // Packed bytes of this thread's pass-QC markers, kRowsPerBlock rows/block:
  std::vector<MappedBlock>   packed_blocks;
  std::vector<std::size_t>   pass_orig_idx;  // length = n_pass_in_chunk (marker indices)
  // Same, for the (much smaller) variance-ratio pool:
  std::vector<MappedBlock>   vr_blocks;
  std::vector<std::size_t>   vr_orig_idx;
  std::exception_ptr         err;            // first exception (if any)
};

// Append one packed row to a block list. `count` is the number of rows already
// in it; the caller keeps the running count.
inline void push_row(std::vector<MappedBlock>& blocks,
                     std::size_t count,
                     const unsigned char* row,
                     std::size_t nbyte_out) {
  const std::size_t off = count % kRowsPerBlock;
  if (off == 0) blocks.emplace_back(kRowsPerBlock * nbyte_out);
  std::memcpy(blocks.back().data() + off * nbyte_out, row, nbyte_out);
}

// Copy a block list into `dest` starting at row `base`, releasing each block
// as soon as it has been copied.
inline void drain_blocks(std::vector<MappedBlock>& blocks,
                         std::size_t n_rows,
                         PackedFlat& dest,
                         std::size_t base,
                         std::size_t nbyte_out) {
  for (std::size_t b = 0; b < blocks.size(); ++b) {
    const std::size_t first = b * kRowsPerBlock;
    const std::size_t rows  = std::min(kRowsPerBlock, n_rows - first);
    for (std::size_t r = 0; r < rows; ++r)
      dest.write(base + first + r, blocks[b].data() + r * nbyte_out);
    blocks[b].release();  // munmap now, not at join
  }
}

} // namespace

ParallelDecodeResult parallel_decode_bed(
    BedReaderPool& reader,
    const int*     ptrsub,
    std::size_t    Nnomissing,
    std::size_t    M,
    float          min_maf,
    float          max_miss,
    int            nthreads,
    const VarRatioRule&  vr,
    const unsigned char* vr_drawn) {
  if (nthreads < 1) nthreads = 1;
  if (nthreads > reader.n_threads())
    throw std::runtime_error(
        "parallel_decode_bed: nthreads=" + std::to_string(nthreads) +
        " exceeds reader pool size " + std::to_string(reader.n_threads()));
  if (ptrsub == nullptr || Nnomissing == 0)
    throw std::runtime_error("parallel_decode_bed: empty ptrsub");

  const std::size_t nbyte_in  = reader.n_bytes_per_marker();
  const std::size_t nbyte_out = (Nnomissing + 3) / 4;

  // Static partitioning. Contiguous ranges minimize per-fd readahead churn.
  std::vector<std::size_t> chunk_lo(nthreads + 1);
  for (int t = 0; t <= nthreads; ++t) {
    chunk_lo[t] = (M * static_cast<std::size_t>(t)) / static_cast<std::size_t>(nthreads);
  }

  // Prebuild the BED byte → geno lookup once — shared read-only across threads.
  BedLut lut;
  build_bed_lookup(lut);

  std::vector<PerThreadState> state(nthreads);

  // -------- Pass 1 — parallel decode --------
  {
    std::vector<std::thread> threads;
    threads.reserve(nthreads);
    for (int t = 0; t < nthreads; ++t) {
      threads.emplace_back([&, t] {
        try {
          const std::size_t lo = chunk_lo[t];
          const std::size_t hi = chunk_lo[t + 1];
          const std::size_t chunk_size = hi - lo;
          auto& st = state[t];
          st.stats.resize(chunk_size);
          st.passQC_flags.assign(chunk_size, 0);
          st.passVR_flags.assign(chunk_size, 0);
          // Guess capacity ~ 25% pass-QC; realloc on demand is fine (this one
          // is 8 bytes per marker, not ⌈N/4⌉).
          st.pass_orig_idx.reserve(chunk_size / 4 + 1);

          std::vector<unsigned char> raw(nbyte_in);
          std::vector<unsigned char> packed(nbyte_out);

          for (std::size_t i = lo; i < hi; ++i) {
            reader.read_marker(t, i, raw.data());
            MarkerStats& s = st.stats[i - lo];
            const bool drawn =
                (vr.enabled && vr_drawn != nullptr) ? (vr_drawn[i] != 0) : false;
            bool passVR = false;
            decode_marker(raw.data(), reader.n_samples(),
                          ptrsub, Nnomissing,
                          min_maf, max_miss,
                          lut, vr, drawn, s, passVR, packed.data());
            if (s.passQC) {
              st.passQC_flags[i - lo] = 1;
              push_row(st.packed_blocks, st.pass_orig_idx.size(),
                       packed.data(), nbyte_out);
              st.pass_orig_idx.push_back(i);
            } else if (passVR) {
              st.passVR_flags[i - lo] = 1;
              push_row(st.vr_blocks, st.vr_orig_idx.size(),
                       packed.data(), nbyte_out);
              st.vr_orig_idx.push_back(i);
            }
          }
        } catch (...) {
          state[t].err = std::current_exception();
        }
      });
    }
    for (auto& th : threads) th.join();
    for (int t = 0; t < nthreads; ++t)
      if (state[t].err) std::rethrow_exception(state[t].err);
  }

  // -------- Prefix sum — where each worker's pass-QC block lands --------
  std::vector<std::size_t> t_offset(nthreads + 1, 0);
  std::vector<std::size_t> t_vr_offset(nthreads + 1, 0);
  for (int t = 0; t < nthreads; ++t) {
    t_offset[t + 1]    = t_offset[t]    + state[t].pass_orig_idx.size();
    t_vr_offset[t + 1] = t_vr_offset[t] + state[t].vr_orig_idx.size();
  }
  const std::size_t total_pass = t_offset[nthreads];
  const std::size_t total_vr   = t_vr_offset[nthreads];

  // -------- Allocate final result shape --------
  ParallelDecodeResult result;
  result.store.init(total_pass, nbyte_out);
  result.store.set_n_stored(total_pass);
  result.stats.resize(M);
  result.passQC.assign(M, false);
  result.passVR.assign(M, false);
  result.orig_plink_idx.resize(total_pass);
  if (total_vr > 0) {
    result.vr_store.init(total_vr, nbyte_out);
    result.vr_store.set_n_stored(total_vr);
  }
  result.vr_orig_idx.resize(total_vr);

  // Merge per-marker stats / passQC  (sequential, cheap — M * small stats)
  for (int t = 0; t < nthreads; ++t) {
    const std::size_t lo = chunk_lo[t];
    const std::size_t hi = chunk_lo[t + 1];
    const auto& st = state[t];
    for (std::size_t i = lo; i < hi; ++i) {
      result.stats[i]  = st.stats[i - lo];
      result.passQC[i] = (st.passQC_flags[i - lo] != 0);
      result.passVR[i] = (st.passVR_flags[i - lo] != 0);
    }
  }

  // -------- Pass 2 — parallel memcpy to PackedFlat --------
  {
    std::vector<std::thread> threads;
    threads.reserve(nthreads);
    for (int t = 0; t < nthreads; ++t) {
      threads.emplace_back([&, t] {
        try {
          auto& st = state[t];
          const std::size_t base = t_offset[t];
          for (std::size_t k = 0; k < st.pass_orig_idx.size(); ++k)
            result.orig_plink_idx[base + k] = st.pass_orig_idx[k];
          drain_blocks(st.packed_blocks, st.pass_orig_idx.size(),
                       result.store, base, nbyte_out);

          const std::size_t vr_base = t_vr_offset[t];
          for (std::size_t k = 0; k < st.vr_orig_idx.size(); ++k)
            result.vr_orig_idx[vr_base + k] = st.vr_orig_idx[k];
          drain_blocks(st.vr_blocks, st.vr_orig_idx.size(),
                       result.vr_store, vr_base, nbyte_out);
        } catch (...) {
          state[t].err = std::current_exception();
        }
      });
    }
    for (auto& th : threads) th.join();
    for (int t = 0; t < nthreads; ++t)
      if (state[t].err) std::rethrow_exception(state[t].err);
  }

  return result;
}

} // namespace saige
