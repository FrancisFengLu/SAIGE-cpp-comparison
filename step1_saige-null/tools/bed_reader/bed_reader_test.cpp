// bed_reader_test.cpp — standalone smoke test for BedReaderPool.
//
// Purpose: verify that ::pread over the UKB BED produces bytes byte-identical
// to the ifstream-based sequential reader that `setGenoObj` uses today. We
// don't link against SAIGE_step1_fast.cpp so the production binary is
// untouched during PR-1.
//
// Usage:
//   ./bed_reader_test <bedfile> <famfile> [marker_idx ...]
//
// What it checks:
//   1. Header magic + mode byte parsing.
//   2. A round-trip read at markers 0 and (M-1): first/last bytes match the
//      values produced by a plain ifstream at the same offsets.
//   3. If you pass extra marker indices on the CLI, it prints their first
//      5 bytes so you can eyeball any marker.
//
// Build (example):
//   g++ -std=c++17 -O2 -o bed_reader_test bed_reader.cpp bed_reader_test.cpp
//
// The file is NOT added to the Makefile automatically. Build manually with
// the one-liner above during PR-1 validation; we'll wire it into a
// `make test_bed` target in a later PR.

#include "bed_reader.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <sys/stat.h>
#include <vector>

namespace {

std::size_t count_lines(const std::string& path) {
  std::ifstream f(path);
  if (!f) throw std::runtime_error("failed to open " + path);
  std::size_t n = 0;
  std::string junk;
  while (std::getline(f, junk)) ++n;
  return n;
}

std::vector<unsigned char>
ifstream_read_marker(const std::string& bed_path, std::size_t m,
                     std::size_t nbyte) {
  std::ifstream f(bed_path, std::ios::binary);
  if (!f) throw std::runtime_error("failed to open " + bed_path);
  f.seekg(static_cast<std::streamoff>(3 + m * nbyte));
  std::vector<unsigned char> buf(nbyte);
  f.read(reinterpret_cast<char*>(buf.data()),
         static_cast<std::streamsize>(nbyte));
  if (static_cast<std::size_t>(f.gcount()) != nbyte)
    throw std::runtime_error("short read at marker " + std::to_string(m));
  return buf;
}

void print_first5(const char* label, std::size_t m,
                  const std::vector<unsigned char>& v) {
  std::printf("[%s] marker %zu first5: %02x %02x %02x %02x %02x\n",
              label, m, v[0], v[1], v[2], v[3], v[4]);
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    std::fprintf(stderr,
                 "usage: %s <bedfile> <famfile> [marker_idx ...]\n",
                 argv[0]);
    return 2;
  }
  const std::string bed_path = argv[1];
  const std::string fam_path = argv[2];

  // FAM → N; BED size → M.
  const std::size_t N = count_lines(fam_path);
  struct stat st{};
  if (::stat(bed_path.c_str(), &st) != 0) {
    std::perror(("stat " + bed_path).c_str());
    return 1;
  }
  const std::size_t nbyte = (N + 3) / 4;
  const std::size_t M = (static_cast<std::size_t>(st.st_size) - 3) / nbyte;

  std::printf("BED  : %s  (%lld bytes)\n", bed_path.c_str(),
              static_cast<long long>(st.st_size));
  std::printf("FAM  : %s  N=%zu\n", fam_path.c_str(), N);
  std::printf("nbyte per marker = %zu\n", nbyte);
  std::printf("M (markers) = %zu\n", M);
  std::printf("\n");

  try {
    saige::BedReaderPool reader(bed_path, N, /*n_threads=*/1);
    if (reader.n_markers_from_file_size() != M)
      throw std::runtime_error(
          "M mismatch: file says " +
          std::to_string(reader.n_markers_from_file_size()) +
          ", we computed " + std::to_string(M));

    // --- check marker 0 ---
    {
      auto ref = ifstream_read_marker(bed_path, 0, nbyte);
      auto got = reader.read_marker(0, 0);
      print_first5("ifstream", 0, ref);
      print_first5("pread   ", 0, got);
      if (ref != got) {
        std::fprintf(stderr, "MISMATCH at marker 0\n");
        return 1;
      }
      std::printf("marker 0: OK (%zu bytes match)\n\n", nbyte);
    }

    // --- check marker M-1 ---
    {
      auto ref = ifstream_read_marker(bed_path, M - 1, nbyte);
      auto got = reader.read_marker(0, M - 1);
      print_first5("ifstream", M - 1, ref);
      print_first5("pread   ", M - 1, got);
      if (ref != got) {
        std::fprintf(stderr, "MISMATCH at marker M-1\n");
        return 1;
      }
      std::printf("marker M-1 (%zu): OK (%zu bytes match)\n\n", M - 1, nbyte);
    }

    // --- user-specified markers for visual inspection ---
    for (int k = 3; k < argc; ++k) {
      const std::size_t m = std::strtoull(argv[k], nullptr, 10);
      if (m >= M) {
        std::fprintf(stderr, "skip m=%zu (M=%zu)\n", m, M);
        continue;
      }
      auto ref = ifstream_read_marker(bed_path, m, nbyte);
      auto got = reader.read_marker(0, m);
      print_first5("ifstream", m, ref);
      print_first5("pread   ", m, got);
      if (ref != got) {
        std::fprintf(stderr, "MISMATCH at marker %zu\n", m);
        return 1;
      }
      std::printf("marker %zu: OK\n", m);
    }

    std::printf("\nbed_reader_test: all checks passed\n");
    return 0;
  } catch (const std::exception& e) {
    std::fprintf(stderr, "bed_reader_test: %s\n", e.what());
    return 1;
  }
}
