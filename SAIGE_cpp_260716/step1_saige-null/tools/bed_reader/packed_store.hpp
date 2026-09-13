// packed_store.hpp — PR-2, §5.5 of PARALLEL_BED_PLAN_2026-04-22.md.
//
// Replace the existing `genoVecofPointers` (vector<vector<unsigned char>*>,
// allocated one-per-marker) with a single flat buffer plus a PackedSpan shim
// that keeps legacy `.at(j)` call sites compiling during migration.
//
// Shape:
//   packed_flat = concat( packed bytes of pass-QC markers in marker-index order )
//   packed_flat.size()  == n_stored * nbyte_per_marker
//   span(i).data/len    — view over the i-th stored marker's bytes
//
// This file is header-only; the sandbox test binary uses it directly.
// When we integrate with setGenoObj, genoClass can own one of these instead
// of its vector-of-pointers.
#pragma once

#include <cstddef>
#include <cstring>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace saige {

// Tiny view that walks the shim line for the current
//   genoVecofPointers[i]->at(j)
// pattern. We keep .at(j) bounds-checked just like std::vector to catch bugs
// at the moment of migration; inner loops can switch to direct .data[] later.
struct PackedSpan {
  const unsigned char* data = nullptr;
  std::size_t          len  = 0;

  unsigned char at(std::size_t j) const {
    if (j >= len)
      throw std::out_of_range("PackedSpan::at: j=" + std::to_string(j) +
                              " len=" + std::to_string(len));
    return data[j];
  }
  unsigned char operator[](std::size_t j) const { return data[j]; }
  std::size_t size() const { return len; }
};

class PackedFlat {
 public:
  // Pre-allocate for up to `capacity_markers` entries of `nbyte_per_marker`
  // bytes each. Call once before any append().
  //
  // The buffer is left UNINITIALIZED on purpose. `new unsigned char[n]`
  // (no parens) default-initializes, i.e. does not write the pages, so RSS
  // grows only as rows are actually stored. The previous assign(n, 0) touched
  // the whole buffer up front; during the parallel loader's copy-out phase
  // that stacked on top of the still-live per-thread staging buffers and made
  // peak RSS ≈ 2× the packed genotype matrix. Every slot below n_stored() is
  // fully overwritten by append()/write(), and span() refuses to hand out a
  // slot at or beyond n_stored(), so no caller can observe the garbage.
  void init(std::size_t capacity_markers, std::size_t nbyte_per_marker) {
    nbyte_    = nbyte_per_marker;
    capacity_ = capacity_markers;
    n_stored_ = 0;
    bytes_    = capacity_markers * nbyte_per_marker;
    data_.reset(bytes_ ? new unsigned char[bytes_] : nullptr);
  }

  // Append one marker's packed bytes. Returns the compact index assigned to
  // it (== `n_stored()` before the call). Throws if capacity exceeded.
  std::size_t append(const unsigned char* bytes) {
    if (n_stored_ >= capacity_)
      throw std::runtime_error("PackedFlat::append beyond capacity " +
                               std::to_string(capacity_));
    std::memcpy(data_.get() + n_stored_ * nbyte_, bytes, nbyte_);
    return n_stored_++;
  }

  // Directly fill the i-th slot (for parallel post-scan writes).
  void write(std::size_t compact_idx, const unsigned char* bytes) {
    if (compact_idx >= capacity_)
      throw std::out_of_range(
          "PackedFlat::write idx=" + std::to_string(compact_idx) +
          " capacity=" + std::to_string(capacity_));
    std::memcpy(data_.get() + compact_idx * nbyte_, bytes, nbyte_);
  }

  // Set the final stored count (use after parallel writes done via write()).
  void set_n_stored(std::size_t n) {
    if (n > capacity_)
      throw std::out_of_range("PackedFlat::set_n_stored > capacity");
    n_stored_ = n;
  }

  PackedSpan span(std::size_t compact_idx) const {
    if (compact_idx >= n_stored_)
      throw std::out_of_range(
          "PackedFlat::span idx=" + std::to_string(compact_idx) +
          " n_stored=" + std::to_string(n_stored_));
    return { data_.get() + compact_idx * nbyte_, nbyte_ };
  }

  // Accessors.
  std::size_t n_stored() const { return n_stored_; }
  std::size_t nbyte()    const { return nbyte_; }
  std::size_t capacity() const { return capacity_; }
  std::size_t bytes()    const { return bytes_; }
  const unsigned char* raw() const { return data_.get(); }

 private:
  // unique_ptr, not vector: vector has no way to allocate without value-
  // initializing. Move-only is fine — nothing copies a PackedFlat, it is
  // always std::move()d into its final owner.
  std::unique_ptr<unsigned char[]> data_;
  std::size_t nbyte_    = 0;
  std::size_t n_stored_ = 0;
  std::size_t capacity_ = 0;
  std::size_t bytes_    = 0;
};

} // namespace saige
