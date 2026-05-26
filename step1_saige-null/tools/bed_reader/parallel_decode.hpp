// parallel_decode.hpp — PR-3 of PARALLEL_BED_PLAN_2026-04-22.md.
//
// Orchestrates the BED marker loop across `nthreads` std::thread workers.
// Each worker owns its own file descriptor (via BedReaderPool) and a private
// per-thread scratch space — no shared writes inside the hot loop.
//
// Design (two-pass):
//   Pass 1 (parallel decode)
//     Partition [0, M) into `nthreads` contiguous chunks. Each worker, for
//     its chunk, calls BedReaderPool::read_marker(tid, i, raw) and then
//     decode_marker(...). It accumulates:
//       - per-marker stats  (length = chunk size)
//       - per-marker passQC (length = chunk size)
//       - packed bytes of pass-QC markers in the chunk, in marker-index order
//       - original marker indices of those pass-QC markers
//
//   Prefix sum (sequential, O(nthreads))
//     Compute each worker's offset into the final PackedFlat.
//
//   Pass 2 (parallel memcpy)
//     Each worker copies its local packed bytes into PackedFlat at its
//     assigned offset; each slot writes to a disjoint region, no locks.
//
// Determinism: the marker-order iteration inside each worker is sequential,
// the chunks are contiguous disjoint ranges, and per-worker output is
// appended to PackedFlat in worker order. So PackedFlat's byte layout is
// identical to the serial pipeline's, regardless of `nthreads`. That is the
// property the regression test below checks.
#pragma once

#include "bed_reader.hpp"
#include "marker_decoder.hpp"
#include "packed_store.hpp"

#include <cstddef>
#include <vector>

namespace saige {

struct ParallelDecodeResult {
  PackedFlat                 store;          // concatenated pass-QC packed bytes
  std::vector<MarkerStats>   stats;          // length M — per-original-marker stats
  std::vector<bool>          passQC;         // length M — per-original-marker flag
  std::vector<std::size_t>   orig_plink_idx; // length n_pass — original marker i per compact slot
};

// Decode the first `M` markers of the BED in parallel using `nthreads`
// std::thread workers. `reader` must have been constructed with at least
// `nthreads` file descriptors. `ptrsub` is length Nnomissing, 1-based FAM
// indices in ascending-FAM order (matches the production `setGenoObj`).
ParallelDecodeResult parallel_decode_bed(
    BedReaderPool& reader,
    const int*     ptrsub,
    std::size_t    Nnomissing,
    std::size_t    M,
    float          min_maf,
    float          max_miss,
    int            nthreads);

} // namespace saige
