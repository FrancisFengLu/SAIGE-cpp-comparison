// gpu_matvec_cpu_stub.cpp — linked when SAIGE is built without CUDA.
// Implements the facade from gpu_matvec.hpp as a no-op so saige-null always
// links, regardless of whether nvcc was present at build time.
#include "gpu_matvec.hpp"

namespace saige::gpu {

bool available() { return false; }

Handle* create(const saige::PackedFlat&,
               const std::vector<float>&,
               const std::vector<float>&,
               int,
               int) {
  return nullptr;
}

Handle* create_rows(const unsigned char* const*,
                    std::size_t,
                    std::size_t,
                    const std::vector<float>&,
                    const std::vector<float>&,
                    int,
                    int) {
  return nullptr;
}

bool matvec(Handle*, const float*, float*) { return false; }

bool matvec_mat_available(const Handle*) { return false; }

bool matvec_mat(Handle*, const float*, int, float*) { return false; }

void destroy(Handle*) {}

int tier(const Handle*) { return 0; }

}  // namespace saige::gpu
