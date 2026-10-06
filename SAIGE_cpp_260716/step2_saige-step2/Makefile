UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  CXX ?= clang++
else
  CXX ?= g++
endif

# W1-2 (same fix as step1 Phase 1/2): `?=` was silently defeated by conda
# activate scripts that export CXXFLAGS (-march=nocona -O2, no -fopenmp!), so
# objects were built WITHOUT -O3/-fopenpm — OpenMP pragmas compiled to no-ops.
# `:=` makes the Makefile assignment win over the environment (make only lets
# the environment override `?=`, or any assignment under `make -e`).
# -march=native / -funroll-loops / -DARMA_NO_DEBUG mirror the step1 flags.
CXXFLAGS := -std=c++17 -O3 -march=native -funroll-loops -DARMA_NO_DEBUG -DARMA_USE_SUPERLU=1 -DARMA_64BIT_WORD=1 -MMD -MP -fopenmp -pthread

# REQUIRED_DEFS are always appended, protecting the mandatory defines even if
# CXXFLAGS is force-overridden (make -e / make CXXFLAGS=...).
REQUIRED_DEFS := -DARMA_USE_SUPERLU=1 -DARMA_64BIT_WORD=1
override CXXFLAGS += $(REQUIRED_DEFS)

PKG_CONFIG   := $(shell command -v pkg-config 2>/dev/null)

ifeq ($(strip $(PKG_CONFIG)),)
  ARMA_CFLAGS :=
  ARMA_LIBS   :=
  YAML_CFLAGS :=
  YAML_LIBS   :=
  HTS_CFLAGS  :=
  HTS_LIBS    :=
  ZSTD_CFLAGS :=
  ZSTD_LIBS   :=
else
  ARMA_CFLAGS := $(shell "$(PKG_CONFIG)" --cflags armadillo 2>/dev/null)
  ARMA_LIBS   := $(shell "$(PKG_CONFIG)" --libs armadillo 2>/dev/null)
  YAML_CFLAGS := $(shell "$(PKG_CONFIG)" --cflags yaml-cpp 2>/dev/null)
  YAML_LIBS   := $(shell "$(PKG_CONFIG)" --libs yaml-cpp 2>/dev/null)
  HTS_CFLAGS  := $(shell "$(PKG_CONFIG)" --cflags htslib 2>/dev/null)
  HTS_LIBS    := $(shell "$(PKG_CONFIG)" --libs htslib 2>/dev/null)
  ZSTD_CFLAGS := $(shell "$(PKG_CONFIG)" --cflags libzstd 2>/dev/null)
  ZSTD_LIBS   := $(shell "$(PKG_CONFIG)" --libs libzstd 2>/dev/null)
endif

ifeq ($(UNAME_S),Darwin)
  INCLUDES := -I/usr/local/include \
              -I/opt/homebrew/opt/boost/include \
              -I/opt/homebrew/include \
              -I/opt/homebrew/Cellar/zstd/1.5.7/include
else
  INCLUDES := -I/usr/local/include
  # W1-2: headers (armadillo, boost, yaml-cpp, ...) live in the conda env on
  # the build machines; previously they were only found via the -isystem that
  # the conda-exported CXXFLAGS carried. Now that CXXFLAGS is `:=`, add the
  # env include path explicitly (same pattern as step1's Makefile).
  ifneq ($(strip $(CONDA_PREFIX)),)
    INCLUDES += -isystem $(CONDA_PREFIX)/include
  endif
endif

ifneq ($(strip $(YAML_CFLAGS)),)
  INCLUDES += $(YAML_CFLAGS)
endif
ifneq ($(strip $(ARMA_CFLAGS)),)
  INCLUDES += $(ARMA_CFLAGS)
endif
ifneq ($(strip $(HTS_CFLAGS)),)
  INCLUDES += $(HTS_CFLAGS)
endif
ifneq ($(strip $(ZSTD_CFLAGS)),)
  INCLUDES += $(ZSTD_CFLAGS)
endif

ifeq ($(UNAME_S),Darwin)
  LIBS := -lyaml-cpp \
          -larmadillo \
          -L/opt/homebrew/opt/openblas/lib -lopenblas \
          -llapack \
          -L/opt/homebrew/opt/superlu/lib -lsuperlu \
          -L/opt/homebrew/lib -lhts \
          -L/opt/homebrew/Cellar/zstd/1.5.7/lib -lzstd \
          -lz \
          -fopenmp -pthread
else
  LIBS := -lyaml-cpp \
          -larmadillo \
          -lopenblas \
          -llapack \
          -lsuperlu \
          -L/usr/local/lib -lhts \
          -lzstd \
          -lz \
          -lsqlite3 \
          -fopenmp -pthread
  # W1-2: link against the conda env libs and bake the rpath (step1 pattern).
  ifneq ($(strip $(CONDA_PREFIX)),)
    LIBS += -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
  endif
endif

ifneq ($(strip $(YAML_LIBS)),)
  LIBS += $(YAML_LIBS)
endif
ifneq ($(strip $(ARMA_LIBS)),)
  LIBS += $(ARMA_LIBS)
endif
ifneq ($(strip $(HTS_LIBS)),)
  LIBS += $(HTS_LIBS)
endif
ifneq ($(strip $(ZSTD_LIBS)),)
  LIBS += $(ZSTD_LIBS)
endif

# ---- optional: glibc libmvec's 4-wide double erfc (score_vec.hpp) ---------
# Used by the mtVecQuantStats block path for the chi-square(1) upper tail.
# Probe-gated so the build stays portable: on anything that is not glibc 2.35+
# on x86-64 with AVX2 the probe fails, SAIGE_HAVE_MVEC_ERFC stays undefined and
# score_vec.hpp uses the scalar libm erfc instead (18.4 ns vs 5.1 ns, i.e. 13 ns
# out of a ~1500 ns per-pair tail). Set USE_MVEC=0 to skip the probe.
USE_MVEC ?= 1
ifeq ($(USE_MVEC),1)
  MVEC_OK := $(shell $(CXX) -O0 -std=c++17 -march=native tools/mvec_probe.cpp \
                     -o /dev/null -lmvec > /dev/null 2>&1 && echo yes)
  ifeq ($(MVEC_OK),yes)
    override CXXFLAGS += -DSAIGE_HAVE_MVEC_ERFC=1
    LIBS += -lmvec
  endif
endif

# ---- optional stage timers for the multi-trait loop (phase_timing.hpp) ----
# DEFAULT OFF. Without PHASE_TIMING=1 every PT_* macro expands to ((void)0),
# phase_timing.cpp compiles to an empty translation unit, and the binary is
# byte-for-byte the production one. Used to measure how the multi-trait
# single-variant loop splits into the dense fast-path scan and the scalar
# slow-path corrections (SPA / ER / Firth).
#   make clean && make PHASE_TIMING=1 TARGET=saige-step2.phase
# The object files are shared between the two configurations, so a `make clean`
# between builds is required.
PHASE_TIMING ?= 0
ifeq ($(PHASE_TIMING),1)
  override CXXFLAGS += -DMT_PHASE_TIMING
endif

# ---- optional GPU path for step 2 (gpu/gpu_step2.hpp) ----------------------
# DEFAULT OFF. Without USE_CUDA=1 the stub translation unit is linked, so
# saige::gpu2::available() always says no, main()'s GPU gate can never be
# satisfied, and the binary behaves exactly as it did before the GPU work
# existed. Nothing in main.cpp is #ifdef'd on CUDA.
#   make USE_CUDA=1            build with the CUDA kernels (needs nvcc)
#   make USE_CUDA=1 SM=70      target a specific compute capability (V100 = 70)
USE_CUDA ?= 0
INCLUDES += -Igpu

ifeq ($(USE_CUDA),1)
  NVCC      ?= nvcc
  NVCC_PATH := $(shell command -v $(NVCC) 2>/dev/null)
  ifeq ($(strip $(NVCC_PATH)),)
    $(error USE_CUDA=1 but nvcc is not on PATH; set NVCC=/path/to/nvcc)
  endif
  CUDA_HOME ?= $(patsubst %/bin/,%,$(dir $(NVCC_PATH)))
  SM        ?= 70
  NVCCFLAGS ?= -O3 -std=c++17 -arch=sm_$(SM) -lineinfo -Xcompiler -fPIC
  GPU_OBJ   := gpu/gpu_step2.cu.o gpu/gpu_spa.cu.o gpu/spa_gpu/spa_gpu.cu.o gpu/gpu_firth.cu.o gpu/gpu_sparse.cu.o gpu/gpu_er.cu.o
  LIBS      += -L$(CUDA_HOME)/lib64 -lcudart -lcublas -Wl,-rpath,$(CUDA_HOME)/lib64
else
  GPU_OBJ   := gpu/gpu_step2_stub.o gpu/gpu_spa_stub.o gpu/spa_gpu/spa_gpu_stub.o gpu/gpu_firth_stub.o gpu/gpu_sparse_stub.o gpu/gpu_er_stub.o
endif

MAKEFLAGS += --no-builtin-rules

SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)
DEPS    := $(OBJECTS:.o=.d)

TARGET  ?= saige-step2

all: $(TARGET) tools/sgs2txt

# sgs2txt converts `outputFormat: sgs` output back to the exact text file.
# It links the very out_fast.o the step-2 binary uses, so the rows it emits
# come out of the same format_text().
tools/sgs2txt: tools/sgs2txt.cpp out_fast.o
	$(CXX) $(CXXFLAGS) $(INCLUDES) tools/sgs2txt.cpp out_fast.o $(LIBS) -o $@

# The two analysis tools behind the sgsPrecision measurement in out_fast.hpp.
# Not in `all`: they answer a question, they are not part of a run.
#   sgs_fidelity  counts, per column, the fields an fp32 .sgs would print
#                 differently, read off an fp64 .sgs
#   out_fieldcmp  the same count between two actual text outputs -- the fp32
#                 round-trip check, since `cmp` can only say "differs"
.PHONY: sgs-tools
sgs-tools: tools/sgs_fidelity tools/out_fieldcmp

tools/sgs_fidelity: tools/sgs_fidelity.cpp sgs_format.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) tools/sgs_fidelity.cpp -o $@

tools/out_fieldcmp: tools/out_fieldcmp.cpp
	$(CXX) $(CXXFLAGS) tools/out_fieldcmp.cpp -o $@

# The device saddlepoint kernel against the production CPU SPA, pair by pair
# (USE_CUDA=1 only; links every step-2 object but main.o plus the CUDA objects).
#   make USE_CUDA=1 SM=70 tools/spa_gpu_check
tools/spa_gpu_check: tools/spa_gpu_check.cpp $(filter-out main.o,$(OBJECTS)) $(GPU_OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) tools/spa_gpu_check.cpp $(filter-out main.o,$(OBJECTS)) $(GPU_OBJ) $(LIBS) -o $@

# Accuracy + cost of score_vec.hpp's chi-square(1) tail against boost's.
tools/mtvec_tail_accuracy: tools/mtvec_tail_accuracy.cpp score_vec.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) tools/mtvec_tail_accuracy.cpp $(filter -lmvec,$(LIBS)) -o $@

$(TARGET): $(OBJECTS) $(GPU_OBJ)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(GPU_OBJ) $(LIBS) -fopenmp -pthread -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

gpu/gpu_step2.cu.o: gpu/gpu_step2.cu gpu/gpu_step2.hpp
	$(NVCC) $(NVCCFLAGS) -Igpu -c $< -o $@

# The saddlepoint kernel: no FMA contraction, so a*b+c rounds the way the CPU
# build's does (gcc -std=c++17 implies -ffp-contract=off); the kernel is exp /
# log bound and the difference in speed is not measurable.
gpu/gpu_spa.cu.o: gpu/gpu_spa.cu gpu/gpu_spa.hpp gpu/gpu_step2.hpp
	$(NVCC) $(NVCCFLAGS) --fmad=false -Igpu -c $< -o $@

# Firth on the device (gpu/gpu_firth.hpp, config key gpuFirth): the same
# one-block-per-pair shape as the SPA kernels, --fmad=false for the same reason.
gpu/gpu_firth.cu.o: gpu/gpu_firth.cu gpu/gpu_firth.hpp gpu/gpu_step2.hpp
	$(NVCC) $(NVCCFLAGS) --fmad=false -Igpu -c $< -o $@

# The within-block cross terms of g' Sigma^-1 g for sparse-GRM binary traits
# (gpu/gpu_sparse.hpp, config key gpuSparse).
gpu/gpu_sparse.cu.o: gpu/gpu_sparse.cu gpu/gpu_sparse.hpp gpu/gpu_step2.hpp
	$(NVCC) $(NVCCFLAGS) -Igpu -c $< -o $@

# The exact test for binary traits on the device (gpu/gpu_er.hpp, config key
# gpuER): --fmad=false like the SPA / Firth kernels; the only fused operations are
# the explicit fma() of the glibc exp / log ports (gpu_er_glibc.cuh).
gpu/gpu_er.cu.o: gpu/gpu_er.cu gpu/gpu_er.hpp gpu/gpu_er_glibc.cuh gpu/gpu_er_glibc_tables.cuh
	$(NVCC) $(NVCCFLAGS) --fmad=false -Igpu -c $< -o $@

# Device ER against er_binary.cpp pair by pair, and the device exp / log against
# the host's bit by bit (USE_CUDA=1 only).
#   make USE_CUDA=1 SM=70 tools/er_gpu_check
tools/er_gpu_check: tools/er_gpu_check.cpp er_binary.o $(GPU_OBJ)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -I. tools/er_gpu_check.cpp er_binary.o $(GPU_OBJ) $(LIBS) -o $@

# The SPA GPU library (gpu/spa_gpu/, from branch spa-gpu-lib, SPA_GPU_LIB.md):
# same contract as gpu_spa.hpp, carriers-only scans, one pass per Newton step,
# Boost's 53-bit erfc with the underflow deferred. Config key gpuSpaImpl picks
# it (lib, the default) or gpu_spa.cu (own).
gpu/spa_gpu/spa_gpu.cu.o: gpu/spa_gpu/spa_gpu.cu gpu/spa_gpu/spa_gpu.hpp gpu/spa_gpu/erfc_boost53.cuh gpu/spa_gpu/erfc_boost53_consts.hpp
	$(NVCC) $(NVCCFLAGS) --fmad=false -Igpu -c $< -o $@

clean:
	rm -f $(OBJECTS) $(DEPS) gpu/*.o gpu/*.d gpu/spa_gpu/*.o gpu/spa_gpu/*.d $(TARGET) tools/sgs2txt \
	      tools/sgs_fidelity tools/out_fieldcmp tools/mtvec_tail_accuracy tools/er_gpu_check

test_region: $(TARGET)
	./$(TARGET) config_compare_region.yaml

test_ldmat: $(TARGET)
	./$(TARGET) config_ldmat.yaml

# LOCO (leave-one-chromosome-out) regression tests.
# Optionally set BASELINE=/path/to/pre-change/saige-step2 to also assert that a
# LOCO=false run is byte-identical to the pre-LOCO build.
# tests/saige-step2.pre-loco is a checked-in copy of the build immediately
# before the LOCO change, used as that baseline by default.
BASELINE ?= $(wildcard tests/saige-step2.pre-loco)

test: $(TARGET)
	bash tests/run_loco_tests.sh ./$(TARGET) $(BASELINE)

.PHONY: all clean test test_region test_ldmat

-include $(DEPS)
