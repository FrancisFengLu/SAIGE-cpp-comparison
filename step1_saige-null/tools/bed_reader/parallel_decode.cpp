#include "parallel_decode.hpp"

#include <algorithm>
#include <atomic>
#include <stdexcept>
#include <thread>

namespace saige {

namespace {

struct PerThreadState {
  // Pass 1 outputs, indexed by (marker - chunk_lo):
  std::vector<MarkerStats> stats;
  std::vector<char>        passQC_flags;     // 0/1 — char because std::vector<bool> can't be written concurrently
  // Concatenated packed bytes of this thread's pass-QC markers:
  std::vector<unsigned char> packed_local;   // size = n_pass_in_chunk * nbyte_out
  std::vector<std::size_t>   pass_orig_idx;  // length = n_pass_in_chunk (marker indices)
  std::exception_ptr         err;            // first exception (if any)
};

} // namespace

ParallelDecodeResult parallel_decode_bed(
    BedReaderPool& reader,
    const int*     ptrsub,
    std::size_t    Nnomissing,
    std::size_t    M,
    float          min_maf,
    float          max_miss,
    int            nthreads) {
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
          // Guess capacity ~ 25% pass-QC; realloc on demand is fine.
          st.packed_local.reserve((chunk_size / 4 + 1) * nbyte_out);
          st.pass_orig_idx.reserve(chunk_size / 4 + 1);

          std::vector<unsigned char> raw(nbyte_in);
          std::vector<unsigned char> packed(nbyte_out);

          for (std::size_t i = lo; i < hi; ++i) {
            reader.read_marker(t, i, raw.data());
            MarkerStats& s = st.stats[i - lo];
            decode_marker(raw.data(), reader.n_samples(),
                          ptrsub, Nnomissing,
                          min_maf, max_miss,
                          lut, s, packed.data());
            if (s.passQC) {
              st.passQC_flags[i - lo] = 1;
              st.packed_local.insert(st.packed_local.end(),
                                     packed.begin(), packed.end());
              st.pass_orig_idx.push_back(i);
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
  for (int t = 0; t < nthreads; ++t) {
    t_offset[t + 1] = t_offset[t] + state[t].pass_orig_idx.size();
  }
  const std::size_t total_pass = t_offset[nthreads];

  // -------- Allocate final result shape --------
  ParallelDecodeResult result;
  result.store.init(total_pass, nbyte_out);
  result.store.set_n_stored(total_pass);
  result.stats.resize(M);
  result.passQC.assign(M, false);
  result.orig_plink_idx.resize(total_pass);

  // Merge per-marker stats / passQC  (sequential, cheap — M * small stats)
  for (int t = 0; t < nthreads; ++t) {
    const std::size_t lo = chunk_lo[t];
    const std::size_t hi = chunk_lo[t + 1];
    const auto& st = state[t];
    for (std::size_t i = lo; i < hi; ++i) {
      result.stats[i]  = st.stats[i - lo];
      result.passQC[i] = (st.passQC_flags[i - lo] != 0);
    }
  }

  // -------- Pass 2 — parallel memcpy to PackedFlat --------
  {
    std::vector<std::thread> threads;
    threads.reserve(nthreads);
    for (int t = 0; t < nthreads; ++t) {
      threads.emplace_back([&, t] {
        try {
          const auto& st = state[t];
          const std::size_t base = t_offset[t];
          for (std::size_t k = 0; k < st.pass_orig_idx.size(); ++k) {
            result.store.write(base + k,
                               st.packed_local.data() + k * nbyte_out);
            result.orig_plink_idx[base + k] = st.pass_orig_idx[k];
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

  return result;
}

} // namespace saige
