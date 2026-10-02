# gpu/spa_gpu -- the device SPA library and its test driver.
#
#   make              spa_gpu.cu.o (the library object) + spa_gpu_test
#   make SM=89        another compute capability (V100 = 70, L4 = 89)
#   make stub         spa_gpu_stub.o for a build without CUDA
#   make consts       regenerate erfc_boost53_consts.hpp (needs the host's Boost headers' values; see gen_erfc_consts.cpp)
#
# The test driver links the step-2 tree's own spa.cpp / spa_binary.cpp /
# UTIL.cpp (compiled here, into build/, with the flags the step-2 Makefile
# uses), so the CPU side of every comparison is the shipped SPA.
#
# Run from this directory with the saige-build conda env active
# (CONDA_PREFIX set), as the step-2 build is.

MAKEFLAGS += --no-builtin-rules
NVCC      ?= nvcc
SM        ?= 70
NVCCFLAGS ?= -O3 -std=c++17 -arch=sm_$(SM) --fmad=false -lineinfo -Xcompiler -fPIC
CUDA_HOME ?= $(patsubst %/bin/,%,$(dir $(shell command -v $(NVCC))))

CXX      ?= g++
CXXFLAGS := -std=c++17 -O3 -march=native -funroll-loops -DARMA_NO_DEBUG -DARMA_USE_SUPERLU=1 -DARMA_64BIT_WORD=1 -fopenmp -pthread
INCLUDES := -I. -I.. -I../.. -I/usr/local/include
LIBS     := -larmadillo -lopenblas -llapack -lsuperlu -fopenmp -pthread
ifneq ($(strip $(CONDA_PREFIX)),)
  INCLUDES += -isystem $(CONDA_PREFIX)/include
  LIBS     += -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib
endif

STEP2 := ../..
CPU_SRCS := $(STEP2)/spa.cpp $(STEP2)/spa_binary.cpp $(STEP2)/UTIL.cpp
CPU_OBJS := $(patsubst $(STEP2)/%.cpp,build/%.o,$(CPU_SRCS))

all: spa_gpu.cu.o spa_gpu_test

spa_gpu.cu.o: spa_gpu.cu spa_gpu.hpp erfc_boost53_consts.hpp
	$(NVCC) $(NVCCFLAGS) -c $< -o $@

stub: spa_gpu_stub.o
spa_gpu_stub.o: spa_gpu_stub.cpp spa_gpu.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

build/%.o: $(STEP2)/%.cpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

build/spa_gpu_test.o: spa_gpu_test.cpp spa_gpu.hpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

spa_gpu_test: build/spa_gpu_test.o spa_gpu.cu.o $(CPU_OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -L$(CUDA_HOME)/lib64 -lcudart -Wl,-rpath,$(CUDA_HOME)/lib64 -o $@

consts: gen_erfc_consts.cpp
	$(CXX) -std=c++17 -O0 $< -o build/gen_erfc_consts && build/gen_erfc_consts > erfc_boost53_consts.hpp

# ---- the same gate against the integrator's kernel -------------------------
# INTEG_DIR holds its gpu_spa.cu + gpu_spa.hpp (a snapshot, or a `git show`
# of the committed files). integ_shim/ supplies the reducer accessors it
# links against, so no part of the reducer is needed.
INTEG_DIR ?= /opt/saige/logs/spa-gpu-lib/integ_snapshot

build/integ_gpu_spa.cu.o: $(INTEG_DIR)/gpu_spa.cu $(INTEG_DIR)/gpu_spa.hpp integ_shim/gpu_step2.hpp
	@mkdir -p build
	$(NVCC) $(NVCCFLAGS) -Iinteg_shim -I$(INTEG_DIR) -c $< -o $@

build/spa_gpu_test_integ.o: spa_gpu_test.cpp spa_gpu.hpp $(INTEG_DIR)/gpu_spa.hpp integ_shim/gpu_step2.hpp
	@mkdir -p build
	$(CXX) $(CXXFLAGS) -DWITH_INTEGRATOR -Iinteg_shim -I$(INTEG_DIR) $(INCLUDES) -c $< -o $@

spa_gpu_test_integ: build/spa_gpu_test_integ.o spa_gpu.cu.o build/integ_gpu_spa.cu.o $(CPU_OBJS)
	$(CXX) $(CXXFLAGS) $^ $(LIBS) -L$(CUDA_HOME)/lib64 -lcudart -Wl,-rpath,$(CUDA_HOME)/lib64 -o $@

clean:
	rm -rf build spa_gpu.cu.o spa_gpu_stub.o spa_gpu_test spa_gpu_test_integ

.PHONY: all stub consts clean
