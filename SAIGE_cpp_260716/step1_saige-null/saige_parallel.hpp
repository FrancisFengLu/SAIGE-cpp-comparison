#pragma once
// OpenMP replacement for the subset of RcppParallel step 1 used
// (Worker / Split / parallelFor / parallelReduce). Worker structs keep the
// RcppParallel shape: operator()(begin, end), a splitting constructor
// Body(const Body&, Split), and join(const Body&).
//
// Thread count: RCPP_PARALLEL_NUM_THREADS if set to a positive integer, else
// omp_get_max_threads(). main.cpp sets RCPP_PARALLEL_NUM_THREADS from
// cfg.nthreads (a value already in the environment wins), so the env var keeps
// the meaning it had under RcppParallel. RCPP_PARALLEL_GRAIN_SIZE is no longer
// read.
//
// parallelReduce is deterministic: [begin, end) is cut into T contiguous
// chunks (T = thread count above, capped by the range length). Chunk 0 runs on
// the caller's body, chunk c > 0 on its own split copy, and the copies are
// joined into the caller's body in chunk order 1, 2, ..., T-1. The result
// depends only on T, not on scheduling, so repeated runs at the same thread
// count are byte-identical. With T = 1 the caller's body sees the whole range
// in one call, which is what TBB did on a one-thread arena (no stealing, so
// every sub-range landed on the same body in index order) for every worker
// whose operator() accumulates element by element.
//
// parallelFor cuts the range into blocks of at least `grain` elements and
// hands them out dynamically; workers write disjoint outputs, so the result
// does not depend on the schedule.

#include <cstddef>
#include <cstdlib>
#include <memory>
#include <vector>
#include <algorithm>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace saige {
namespace par {

struct Split {};

struct Worker {
  virtual ~Worker() {}
  virtual void operator()(std::size_t begin, std::size_t end) = 0;
};

inline int num_threads() {
  if (const char* s = std::getenv("RCPP_PARALLEL_NUM_THREADS")) {
    const int n = std::atoi(s);
    if (n > 0) return n;
  }
#ifdef _OPENMP
  return std::max(1, omp_get_max_threads());
#else
  return 1;
#endif
}

template <typename W>
inline void parallelFor(std::size_t begin, std::size_t end, W& worker,
                        std::size_t grain = 1) {
  if (end <= begin) return;
  const std::size_t n = end - begin;
  const int T = num_threads();
  if (T <= 1 || n <= std::max<std::size_t>(grain, 1)) {
    worker(begin, end);
    return;
  }
  // ~8 blocks per thread for load balance, never smaller than grain.
  std::size_t bs = (n + (std::size_t)T * 8 - 1) / ((std::size_t)T * 8);
  bs = std::max<std::size_t>(bs, std::max<std::size_t>(grain, 1));
  const long long nb = (long long)((n + bs - 1) / bs);
#pragma omp parallel for schedule(dynamic, 1) num_threads(T)
  for (long long b = 0; b < nb; ++b) {
    const std::size_t lo = begin + (std::size_t)b * bs;
    const std::size_t hi = std::min(end, lo + bs);
    worker(lo, hi);
  }
}

template <typename R>
inline void parallelReduce(std::size_t begin, std::size_t end, R& reducer,
                           std::size_t grain = 1) {
  if (end <= begin) return;
  const std::size_t n = end - begin;
  std::size_t T = (std::size_t)num_threads();
  T = std::min(T, std::max<std::size_t>(1, n / std::max<std::size_t>(grain, 1)));
  if (T <= 1) {
    reducer(begin, end);
    return;
  }
  // chunk c covers [begin + c*n/T, begin + (c+1)*n/T)
  std::vector<std::unique_ptr<R>> parts(T);
  for (std::size_t c = 1; c < T; ++c) parts[c].reset(new R(reducer, Split()));
#pragma omp parallel for schedule(static, 1) num_threads((int)T)
  for (long long c = 0; c < (long long)T; ++c) {
    const std::size_t lo = begin + (std::size_t)c * n / T;
    const std::size_t hi = begin + ((std::size_t)c + 1) * n / T;
    R& body = (c == 0) ? reducer : *parts[(std::size_t)c];
    body(lo, hi);
  }
  for (std::size_t c = 1; c < T; ++c) reducer.join(*parts[c]);
}

}  // namespace par
}  // namespace saige
