# Standalone build for the GPU matvec unit (G1–G5 of GPU_INTEGRATION_PLAN).
# Not linked into saige-null from here — that wiring lives in the main
# Makefile. This file just builds a test binary we can iterate on without
# touching the production build.
#
# Targets:
#   make                  → build gpu_matvec_test (if CUDA), else "no-CUDA" stub exit
#   make test             → run the synthetic Au_cpu ≡ Au_gpu check
#   make clean
#
# CUDA knobs:
#   NVCC   : path to nvcc   (probed under /usr/local/cuda*/bin then /usr/bin)
#   SM     : target arch    (default sm_70 = V100; sm_86/sm_89 for Ampere/Ada)
#   CUDA_LIBS : extra -L paths for libcublas/libcudart if needed

CXX      ?= g++
NVCC     ?= $(firstword $(wildcard /usr/local/cuda/bin/nvcc /usr/local/cuda-*/bin/nvcc /usr/bin/nvcc) nvcc-not-found)
SM       ?= sm_70

# The CUDA path is taken only if (a) nvcc is present AND (b) the kernel
# source exists. G1 ships just the scaffold; G2 adds gpu_matvec.cu and
# flips us into the nvcc build automatically.
HAS_NVCC := $(shell test -x $(NVCC) && test -f gpu_matvec.cu && echo yes)

# Append flags instead of `?=` so conda-env-provided CXXFLAGS don't
# overwrite the minimum we need.
CXXFLAGS  += -std=c++17 -O2 -Wall -Wextra -pthread -I../bed_reader
NVCCFLAGS += -std=c++17 -O2 -arch=$(SM) --compiler-options "-Wall -Wextra" -I../bed_reader

# cuBLAS and cudart are in /lib/x86_64-linux-gnu on this Ubuntu box (Debian
# packages libcublas12, libcudart12). The conda-env linker has its own
# sysroot and does NOT search /lib/x86_64-linux-gnu, so point -L there.
CUDA_HOME   := $(patsubst %/bin/nvcc,%,$(NVCC))
CUDA_LIBDIR ?= $(firstword $(dir $(wildcard $(CUDA_HOME)/lib64/libcudart.so /lib/x86_64-linux-gnu/libcudart.so)) /lib/x86_64-linux-gnu)
CUDA_LIBS   ?= -L$(CUDA_LIBDIR) -Wl,-rpath,$(CUDA_LIBDIR) -lcublas -lcudart

ifeq ($(HAS_NVCC),yes)
  CXXFLAGS  += -DSAIGE_WITH_CUDA
endif

all: gpu_matvec_test

# --- Build rules ---
# All .cu units in this dir (gpu_matvec.cu dispatches into the others).
CU_OBJS := $(patsubst %.cu,%.cu.o,$(wildcard *.cu))

gpu_matvec_test: gpu_matvec_test.o \
                 $(if $(filter yes,$(HAS_NVCC)),$(CU_OBJS),gpu_matvec_cpu_stub.o) \
                 ../bed_reader/bed_reader.o \
                 ../bed_reader/marker_decoder.o
ifeq ($(HAS_NVCC),yes)
	# Force the final link to use the SYSTEM toolchain (ld/libc/libpthread
	# from /usr), not the conda-env ones. Object files compile fine with
	# either, but link resolution of CUDA runtime + libc startup objects
	# needs the system linker that nvcc expects.
	PATH=/usr/bin:/bin $(NVCC) -std=c++17 -arch=$(SM) \
	    -ccbin /usr/bin/g++ \
	    $^ -o $@ -lcublas -lcudart -lpthread
else
	$(CXX) $(CXXFLAGS) $^ -o $@ -lpthread
endif

%.cu.o: %.cu gpu_matvec.hpp ../bed_reader/packed_store.hpp
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

gpu_matvec_cpu_stub.o: gpu_matvec_cpu_stub.cpp gpu_matvec.hpp ../bed_reader/packed_store.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

gpu_matvec_test.o: gpu_matvec_test.cpp gpu_matvec.hpp ../bed_reader/packed_store.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Side-builds in case someone wants just the tools/bed_reader deps available.
../bed_reader/bed_reader.o:
	$(MAKE) -C ../bed_reader bed_reader.o
../bed_reader/marker_decoder.o:
	$(MAKE) -C ../bed_reader marker_decoder.o

test: gpu_matvec_test
	./gpu_matvec_test

clean:
	rm -f *.o *.cu.o gpu_matvec_test

# Tell make that .cu files are terminal sources; kills the spurious
# "Circular gpu_matvec.cu <- gpu_matvec.cu.o" warning from implicit rules.
.SUFFIXES: .cu

.PHONY: all test clean
