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
  GPU_OBJ   := gpu/gpu_step2.cu.o
  LIBS      += -L$(CUDA_HOME)/lib64 -lcudart -lcublas -Wl,-rpath,$(CUDA_HOME)/lib64
else
  GPU_OBJ   := gpu/gpu_step2_stub.o
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

$(TARGET): $(OBJECTS) $(GPU_OBJ)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(GPU_OBJ) $(LIBS) -fopenmp -pthread -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

gpu/gpu_step2.cu.o: gpu/gpu_step2.cu gpu/gpu_step2.hpp
	$(NVCC) $(NVCCFLAGS) -Igpu -c $< -o $@

clean:
	rm -f $(OBJECTS) $(DEPS) gpu/*.o gpu/*.d $(TARGET) tools/sgs2txt

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
