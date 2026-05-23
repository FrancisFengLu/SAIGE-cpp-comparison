// bed_reader.cpp — see bed_reader.hpp for design notes.
#include "bed_reader.hpp"

#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

#include <cerrno>
#include <cstring>
#include <sstream>

namespace saige {

namespace {
constexpr unsigned char kBedMagic0 = 0x6C;
constexpr unsigned char kBedMagic1 = 0x1B;
constexpr unsigned char kBedModeSnpMajor = 0x01;
constexpr std::size_t  kBedHeaderBytes = 3;

std::string errno_str(const std::string& what) {
  std::ostringstream os;
  os << what << ": " << std::strerror(errno) << " (errno=" << errno << ")";
  return os.str();
}
} // namespace

BedReaderPool::BedReaderPool(const std::string& bed_path,
                             std::size_t n_samples,
                             int n_threads)
    : path_(bed_path),
      n_samples_(n_samples),
      nbyte_((n_samples + 3) / 4),   // ⌈N/4⌉
      file_size_(0) {
  if (n_threads < 1) n_threads = 1;
  if (n_samples_ == 0)
    throw std::runtime_error("BedReaderPool: n_samples must be > 0");

  // Probe file size.
  {
    struct stat st{};
    if (::stat(path_.c_str(), &st) != 0)
      throw std::runtime_error(errno_str("stat(" + path_ + ")"));
    if (st.st_size < static_cast<off_t>(kBedHeaderBytes + nbyte_))
      throw std::runtime_error("BED file too small for given N: " + path_);
    file_size_ = static_cast<std::size_t>(st.st_size);
  }

  // Open and validate magic once on fd 0; then open the remaining fds.
  fds_.reserve(static_cast<std::size_t>(n_threads));
  for (int t = 0; t < n_threads; ++t) {
    int fd = ::open(path_.c_str(), O_RDONLY);
    if (fd < 0)
      throw std::runtime_error(errno_str("open(" + path_ + ")"));

    if (t == 0) {
      unsigned char hdr[kBedHeaderBytes] = {0, 0, 0};
      ssize_t got = ::pread(fd, hdr, kBedHeaderBytes, 0);
      if (got != static_cast<ssize_t>(kBedHeaderBytes)) {
        int e = errno;
        ::close(fd);
        throw std::runtime_error("BED header read failed: " +
                                 std::string(std::strerror(e)));
      }
      if (hdr[0] != kBedMagic0 || hdr[1] != kBedMagic1) {
        ::close(fd);
        throw std::runtime_error("BED magic mismatch in " + path_ +
                                 " (got " + std::to_string(hdr[0]) + "," +
                                 std::to_string(hdr[1]) + ")");
      }
      if (hdr[2] != kBedModeSnpMajor) {
        ::close(fd);
        throw std::runtime_error(
            "BED is not SNP-major (mode byte " + std::to_string(hdr[2]) +
            ", only 0x01 supported)");
      }
    }

    // Hint the kernel: we'll read linearly and we want the data warm.
    ::posix_fadvise(fd, 0, 0, POSIX_FADV_SEQUENTIAL);
    ::posix_fadvise(fd, 0, 0, POSIX_FADV_WILLNEED);

    fds_.push_back(fd);
  }
}

BedReaderPool::~BedReaderPool() {
  for (int fd : fds_)
    if (fd >= 0) ::close(fd);
}

void BedReaderPool::read_marker(int tid, std::size_t m,
                                unsigned char* out) const {
  if (tid < 0 || tid >= static_cast<int>(fds_.size()))
    throw std::out_of_range("BedReaderPool::read_marker: bad tid " +
                            std::to_string(tid));
  if (out == nullptr)
    throw std::invalid_argument("BedReaderPool::read_marker: null out buffer");

  const off_t off = static_cast<off_t>(kBedHeaderBytes) +
                    static_cast<off_t>(m) * static_cast<off_t>(nbyte_);
  if (off < 0 || static_cast<std::size_t>(off) + nbyte_ > file_size_)
    throw std::out_of_range(
        "BedReaderPool::read_marker: marker " + std::to_string(m) +
        " out of BED file bounds (" + std::to_string(file_size_) + " bytes)");

  std::size_t done = 0;
  while (done < nbyte_) {
    ssize_t got = ::pread(fds_[tid], out + done, nbyte_ - done, off + done);
    if (got < 0) {
      if (errno == EINTR) continue;
      throw std::runtime_error(errno_str("pread(" + path_ + ")"));
    }
    if (got == 0)
      throw std::runtime_error("pread short at marker " + std::to_string(m));
    done += static_cast<std::size_t>(got);
  }
}

std::vector<unsigned char>
BedReaderPool::read_marker(int tid, std::size_t m) const {
  std::vector<unsigned char> buf(nbyte_);
  read_marker(tid, m, buf.data());
  return buf;
}

} // namespace saige
