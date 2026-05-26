// bed_reader.hpp — standalone PLINK BED reader (pread-backed)
// ---------------------------------------------------------------------------
// PR-1 of PARALLEL_BED_PLAN_2026-04-22.md.
//
// Design:
//   * One open() per worker thread; each worker does positional reads
//     (::pread) on its own file descriptor. NFS-safe (no mmap page faults),
//     thread-safe by construction (no shared cursor), and each fd can carry
//     its own readahead window via posix_fadvise(SEQUENTIAL | WILLNEED).
//   * Serial use is supported by constructing with n_threads == 1.
//   * The reader only knows about the BED file; N and M are provided by the
//     caller after FAM/BIM scans. That keeps it testable without a full
//     setGenoObj wired up.
//
// This unit is intentionally isolated — it does NOT touch the existing
// `genoClass` / `genoVecofPointers` machinery. Its job is just
//     (m, out) -> out holds the raw ⌈N/4⌉ bytes of BED marker m.
// QC, standardization, flat packed storage, and the parallel fill loop
// all live in later PRs (§5.2–§5.5 of the plan).
// ---------------------------------------------------------------------------
#pragma once

#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace saige {

class BedReaderPool {
public:
  // Open the BED file, validate the 3-byte magic, and allocate `n_threads`
  // parallel file descriptors. `n_samples` is from the companion FAM scan.
  // Throws std::runtime_error on any I/O or format error.
  BedReaderPool(const std::string& bed_path,
                std::size_t n_samples,
                int n_threads = 1);

  ~BedReaderPool();

  // Copy / move disabled — owning raw fds.
  BedReaderPool(const BedReaderPool&) = delete;
  BedReaderPool& operator=(const BedReaderPool&) = delete;
  BedReaderPool(BedReaderPool&&) = delete;
  BedReaderPool& operator=(BedReaderPool&&) = delete;

  // Read the raw packed bytes of marker `m` (0-based) into `out`.
  // `out` must be preallocated with at least nbyte_per_marker() bytes.
  // `tid` is the caller's TBB/OpenMP worker index in [0, n_threads).
  // In a single-threaded call just pass tid=0.
  void read_marker(int tid, std::size_t m, unsigned char* out) const;

  // Same as above, but returns bytes in a freshly-allocated vector. Handy
  // for tests; inside the hot loop prefer the overload that writes to a
  // caller-provided buffer.
  std::vector<unsigned char> read_marker(int tid, std::size_t m) const;

  // Geometry / accessors.
  std::size_t n_samples()            const { return n_samples_; }
  std::size_t n_bytes_per_marker()   const { return nbyte_;     }
  std::size_t file_size()            const { return file_size_; }
  int         n_threads()            const { return static_cast<int>(fds_.size()); }
  // The BED file has (file_size - 3) / nbyte markers. Callers get this
  // implicitly via the BIM scan; exposed here as a sanity check.
  std::size_t n_markers_from_file_size() const {
    return (file_size_ >= 3) ? (file_size_ - 3) / nbyte_ : 0;
  }

private:
  std::string path_;
  std::size_t n_samples_;
  std::size_t nbyte_;       // ⌈n_samples / 4⌉
  std::size_t file_size_;   // size of the .bed on disk
  std::vector<int> fds_;    // one per worker thread
};

} // namespace saige
