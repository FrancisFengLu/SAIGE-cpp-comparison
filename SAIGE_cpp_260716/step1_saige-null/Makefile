# Platform detection
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  CXX ?= clang++
else
  CXX ?= g++
endif

CXXFLAGS ?= -std=c++17 -O2 -DARMA_64BIT_WORD=1 -DRCPP_PARALLEL_USE_TBB=1 -MMD -MP -pthread -fopenmp

# REQUIRED_DEFS are always appended to CXXFLAGS, even when the environment
# (e.g. a conda/mamba activate script) pre-exports CXXFLAGS — in which case the
# `?=` above is a no-op and -DARMA_64BIT_WORD would otherwise be dropped.
# ARMA_64BIT_WORD is mandatory: sparse n×n GRM matrices have n_rows*n_cols
# (e.g. 158598² ≈ 2.5e10) exceeding the 32-bit arma::uword max (~4.29e9),
# which triggers "SpMat::init(): requested size is too large" in the
# sparse-GRM-fit / cate-VR path.
REQUIRED_DEFS := -DARMA_64BIT_WORD=1 -DRCPP_PARALLEL_USE_TBB=1
CXXFLAGS += $(REQUIRED_DEFS)

RSCRIPT      := $(shell command -v Rscript 2>/dev/null)
PKG_CONFIG   := $(shell command -v pkg-config 2>/dev/null)

# -------------------------
# R / Rcpp includes (auto-detect via Rscript if available)
# -------------------------
ifeq ($(strip $(RSCRIPT)),)
  R_HOME       :=
  RCPP_INC     :=
  RCPPARM_INC  :=
  RCPPARAL_INC :=
  RCPPEIG_INC  :=
else
  R_HOME       := $(shell "$(RSCRIPT)" -e 'cat(R.home())')
  RCPP_INC     := $(shell "$(RSCRIPT)" -e 'cat(system.file("include", package="Rcpp"))')
  RCPPARM_INC  := $(shell "$(RSCRIPT)" -e 'cat(system.file("include", package="RcppArmadillo"))')
  RCPPARAL_INC := $(shell "$(RSCRIPT)" -e 'cat(system.file("include", package="RcppParallel"))')
  RCPPEIG_INC  := $(shell "$(RSCRIPT)" -e 'cat(system.file("include", package="RcppEigen"))')
endif

# -------------------------
# pkg-config (Armadillo / yaml-cpp)
# -------------------------
ifeq ($(strip $(PKG_CONFIG)),)
  ARMA_CFLAGS :=
  ARMA_LIBS   :=
  YAML_CFLAGS :=
  YAML_LIBS   :=
else
  ARMA_CFLAGS := $(shell "$(PKG_CONFIG)" --cflags armadillo 2>/dev/null)
  ARMA_LIBS   := $(shell "$(PKG_CONFIG)" --libs armadillo 2>/dev/null)
  YAML_CFLAGS := $(shell "$(PKG_CONFIG)" --cflags yaml-cpp 2>/dev/null)
  YAML_LIBS   := $(shell "$(PKG_CONFIG)" --libs yaml-cpp 2>/dev/null)
endif

# -------------------------
# Include paths
# -------------------------
ifeq ($(UNAME_S),Darwin)
  INCLUDES := -I. -I/usr/local/include \
              -I/opt/homebrew/opt/boost/include \
              -I/opt/homebrew/include \
              -I/opt/homebrew/opt/eigen/include/eigen3
else
  INCLUDES := -I. -I/usr/local/include \
              -I/usr/include/eigen3
  # Conda/mamba env (if active, add its include paths)
  ifneq ($(strip $(CONDA_PREFIX)),)
    INCLUDES += -isystem $(CONDA_PREFIX)/include \
                -I$(CONDA_PREFIX)/include \
                -I$(CONDA_PREFIX)/include/eigen3
  endif
endif

# Optional R headers
ifneq ($(strip $(R_HOME)),)
  INCLUDES += -I"$(R_HOME)/include"
endif
ifneq ($(strip $(RCPP_INC)),)
  INCLUDES += -I"$(RCPP_INC)"
endif
ifneq ($(strip $(RCPPARM_INC)),)
  INCLUDES += -I"$(RCPPARM_INC)"
endif
ifneq ($(strip $(RCPPARAL_INC)),)
  INCLUDES += -I"$(RCPPARAL_INC)"
endif
ifneq ($(strip $(RCPPEIG_INC)),)
  INCLUDES += -I"$(RCPPEIG_INC)"
endif

# pkg-config cflags
ifneq ($(strip $(YAML_CFLAGS)),)
  INCLUDES += $(YAML_CFLAGS)
endif
ifneq ($(strip $(ARMA_CFLAGS)),)
  INCLUDES += $(ARMA_CFLAGS)
endif

# -------------------------
# Libraries / linker flags
# -------------------------
ifeq ($(UNAME_S),Darwin)
  LIBS := -lyaml-cpp \
          -larmadillo \
          -L/opt/homebrew/opt/openblas/lib -lopenblas \
          -llapack \
          -L/opt/homebrew/opt/superlu/lib -lsuperlu \
          -L/opt/homebrew/opt/tbb/lib -ltbb
else
  LIBS := -pthread -fopenmp \
          -lyaml-cpp \
          -larmadillo \
          -lopenblas \
          -llapack \
          -lsuperlu \
          -ltbb \
          -ldl
  # Conda/mamba env (if active, add its library paths)
  ifneq ($(strip $(CONDA_PREFIX)),)
    LIBS := -pthread -fopenmp \
            -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib \
            -lyaml-cpp -larmadillo -lopenblas -llapack -lsuperlu -ltbb -ldl
  endif
endif

# pkg-config libs
ifneq ($(strip $(YAML_LIBS)),)
  LIBS += $(YAML_LIBS)
endif
ifneq ($(strip $(ARMA_LIBS)),)
  LIBS += $(ARMA_LIBS)
endif

# Optional R library
ifneq ($(strip $(R_HOME)),)
  LIBS += -L"$(R_HOME)/lib" -lR
  LIBS += -Wl,-rpath,"$(R_HOME)/lib"
  # libR has versioned pcre2 refs (PCRE2_10.47); link pcre2 explicitly so the
  # env's libpcre2-8 resolves them (fixes "undefined reference pcre2_*@PCRE2_10.47").
  LIBS += -lpcre2-8
endif

# Optional RcppParallel runtime
RCPPARAL_LIBDIR := $(shell "$(RSCRIPT)" -e 'cat(system.file("libs", package="RcppParallel"))' 2>/dev/null)
RCPPARAL_SO     := $(shell ls "$(RCPPARAL_LIBDIR)/RcppParallel.so" 2>/dev/null)
ifneq ($(strip $(RCPPARAL_LIBDIR)),)
  LIBS += -Wl,-rpath,"$(RCPPARAL_LIBDIR)"
endif
ifneq ($(strip $(RCPPARAL_SO)),)
  LIBS += "$(RCPPARAL_SO)"
endif

# -------------------------
# Sources / targets
# -------------------------
SOURCES := $(wildcard *.cpp)

# Pull in the standalone BED reader library units (bed_reader.cpp,
# marker_decoder.cpp, parallel_decode.cpp) — but NOT their test binaries
# which each have their own main().
BED_READER_DIR  := tools/bed_reader
BED_READER_TESTS := $(BED_READER_DIR)/bed_reader_test.cpp \
                    $(BED_READER_DIR)/bed_pipeline_test.cpp \
                    $(BED_READER_DIR)/bed_parallel_test.cpp
BED_READER_SRC   := $(filter-out $(BED_READER_TESTS), $(wildcard $(BED_READER_DIR)/*.cpp))
SOURCES         += $(BED_READER_SRC)
INCLUDES        += -I$(BED_READER_DIR)

# ---- GPU backend (optional, tier-1 cuBLAS sgemv) --------------------------
# Enable if nvcc is present AND tools/gpu_matvec/gpu_matvec.cu exists.
# When enabled: -DSAIGE_WITH_CUDA in CXX flags, nvcc compiles the .cu file,
# link is done via nvcc with the SYSTEM toolchain. When disabled the stub
# TU is linked so saige-null still exports saige::gpu::{available,…} symbols
# as no-ops.
GPU_MATVEC_DIR  := tools/gpu_matvec
GPU_MATVEC_TEST := $(GPU_MATVEC_DIR)/gpu_matvec_test.cpp
NVCC           ?= /usr/bin/nvcc
GPU_SM         ?= sm_86
HAS_CUDA       := $(shell test -x $(NVCC) && test -f $(GPU_MATVEC_DIR)/gpu_matvec.cu && echo yes)
INCLUDES       += -I$(GPU_MATVEC_DIR)

ifeq ($(HAS_CUDA),yes)
  CXXFLAGS       += -DSAIGE_WITH_CUDA
  GPU_EXTRA_SRC  := $(filter-out $(GPU_MATVEC_TEST) $(GPU_MATVEC_DIR)/gpu_matvec_cpu_stub.cpp $(GPU_MATVEC_DIR)/gpu_matvec.cu, \
                                  $(wildcard $(GPU_MATVEC_DIR)/*.cpp))
  GPU_CU_OBJ     := $(GPU_MATVEC_DIR)/gpu_matvec.cu.o
  GPU_LIBS       := -L/lib/x86_64-linux-gnu -Xlinker -rpath=/lib/x86_64-linux-gnu -lcublas -lcudart
else
  # stub path: include the noop implementation so saige::gpu::available() links
  GPU_EXTRA_SRC  := $(GPU_MATVEC_DIR)/gpu_matvec_cpu_stub.cpp
  GPU_CU_OBJ     :=
  GPU_LIBS       :=
endif
SOURCES         += $(GPU_EXTRA_SRC)

OBJECTS := $(SOURCES:.cpp=.o)
DEPS    := $(OBJECTS:.o=.d)

TARGET  ?= saige-null

all: $(TARGET)

$(TARGET): $(OBJECTS) $(GPU_CU_OBJ)
ifeq ($(HAS_CUDA),yes)
	# Final link via nvcc on the SYSTEM toolchain so CUDA runtime symbols
	# resolve correctly (libcublas/libcudart from /lib/x86_64-linux-gnu).
	# Transform host-compiler flags nvcc doesn't understand:
	#   -pthread, -fopenmp    → -Xcompiler=<flag>
	#   -Wl,-rpath,FOO        → -Xlinker -rpath=FOO    (GNU ld accepts =path form)
	#   -Wl,<other>           → -Xlinker <flag>
	PATH=/usr/bin:/bin $(NVCC) -std=c++17 -arch=$(GPU_SM) -ccbin /usr/bin/g++ \
	    $(OBJECTS) $(GPU_CU_OBJ) \
	    $(shell echo "$(LIBS)" | sed \
	        -e 's/-pthread/-Xcompiler=-pthread/g' \
	        -e 's/-fopenmp/-Xcompiler=-fopenmp/g' \
	        -e 's/-Wl,-rpath,\([^ ]*\)/-Xlinker -rpath=\1/g') \
	    $(GPU_LIBS) -o $@
else
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LIBS) -o $@
endif

$(GPU_MATVEC_DIR)/gpu_matvec.cu.o: $(GPU_MATVEC_DIR)/gpu_matvec.cu $(GPU_MATVEC_DIR)/gpu_matvec.hpp
	$(NVCC) -std=c++17 -O2 -arch=$(GPU_SM) -I$(BED_READER_DIR) -I$(GPU_MATVEC_DIR) -c $< -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Profile build: debug symbols + light optimization for meaningful callgrind output
profile: CXXFLAGS := $(filter-out -O2,$(CXXFLAGS)) -O1 -g
profile: clean $(TARGET)

# Run under callgrind (requires CONFIG= to be set, e.g. make callgrind CONFIG=config_test.yaml)
CONFIG ?= config_test.yaml
callgrind: profile
	valgrind --tool=callgrind \
	         --callgrind-out-file=callgrind.out.%p \
	         --collect-jumps=yes \
	         --cache-sim=yes \
	         --branch-sim=yes \
	         ./$(TARGET) -c $(CONFIG)

clean:
	rm -f $(OBJECTS) $(DEPS) $(TARGET)

.PHONY: all clean profile callgrind

-include $(DEPS)
