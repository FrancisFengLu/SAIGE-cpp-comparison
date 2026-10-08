# Top-level build: the command-line front end plus the two engines, in bin/.
#
#   make                          CPU build, -march=native
#   make USE_CUDA=1 SM=70         GPU build for compute capability 7.0 (V100)
#   make USE_CUDA=1 SM="70 80 90" one GPU build for V100, A100 and H100
#   make ARCH=x86-64-v3           portable CPU code (any AVX2 machine) instead of native
#   make PROGRAM=other-name       name of the front end (docs use saige-gpu-cpp)
#
# Result:
#   bin/$(PROGRAM)    the front end:  $(PROGRAM) step1|step2|sgs2txt ...
#   bin/saige-null    step-1 engine   (saige-null -c step1.yaml)
#   bin/saige-step2   step-2 engine   (saige-step2 step2.yaml)
#   bin/sgs2txt       .sgs -> text converter
#
# Needs the build environment of docs/index.md active (conda activate ...).
# Changing ARCH, USE_CUDA or SM triggers a clean rebuild of both engines (stale
# objects from another configuration break the link or the ABI).

PROGRAM  ?= saige-gpu-cpp
ARCH     ?= native
USE_CUDA ?= 0
SM       ?= 70
CUDA_HOME ?= /usr/local/cuda
NVCC     ?= $(CUDA_HOME)/bin/nvcc
JOBS     ?= 8
BINDIR   ?= bin

S1DIR := step1_saige-null
S2DIR := step2_saige-step2

VERSION := $(shell git describe --always --dirty --abbrev=8 2>/dev/null || echo unknown)
comma := ,
BUILD_DESC := ARCH=$(ARCH) $(if $(filter 1,$(USE_CUDA)),USE_CUDA=1 SM=$(SM),USE_CUDA=0)

ifeq ($(USE_CUDA),1)
  S1ARGS := NVCC=$(NVCC) SM="$(SM)"
  S2ARGS := USE_CUDA=1 NVCC=$(NVCC) SM="$(SM)" CUDA_HOME=$(CUDA_HOME)
else
  S1ARGS := NVCC=none
  S2ARGS := USE_CUDA=0
endif

# front end: yaml-cpp only
CLI_CXX      ?= $(CXX)
CLI_CXXFLAGS := -std=c++17 -O2 -Wall
CLI_INC      := $(if $(CONDA_PREFIX),-isystem $(CONDA_PREFIX)/include)
CLI_LIBS     := $(if $(CONDA_PREFIX),-L$(CONDA_PREFIX)/lib -Wl$(comma)-rpath$(comma)$(CONDA_PREFIX)/lib) -lyaml-cpp

STAMP := .build-config
CONFIG_STR := ARCH=$(ARCH) USE_CUDA=$(USE_CUDA) SM=$(SM) NVCC=$(NVCC)

.PHONY: all engines step1 step2 clean FORCE

all: engines $(BINDIR)/$(PROGRAM)
	@echo "built $(BINDIR)/: $(PROGRAM) saige-null saige-step2 sgs2txt ($(BUILD_DESC))"

# Rewrites the stamp only when the configuration changed; a change cleans both engines.
$(STAMP): FORCE
	@if [ "$$(cat $(STAMP) 2>/dev/null)" != "$(CONFIG_STR)" ]; then \
	  echo "build configuration changed ($(CONFIG_STR)): cleaning the engines"; \
	  $(MAKE) -C $(S1DIR) clean NVCC=none > /dev/null; \
	  $(MAKE) -C $(S2DIR) clean > /dev/null; \
	  echo "$(CONFIG_STR)" > $(STAMP); \
	fi

engines: step1 step2

step1: $(STAMP)
	$(MAKE) -C $(S1DIR) -j$(JOBS) ARCH=$(ARCH) $(S1ARGS)
	@mkdir -p $(BINDIR)
	cp -f $(S1DIR)/saige-null $(BINDIR)/saige-null

step2: $(STAMP)
	$(MAKE) -C $(S2DIR) -j$(JOBS) ARCH=$(ARCH) $(S2ARGS)
	@mkdir -p $(BINDIR)
	cp -f $(S2DIR)/saige-step2 $(BINDIR)/saige-step2
	cp -f $(S2DIR)/tools/sgs2txt $(BINDIR)/sgs2txt

$(BINDIR)/$(PROGRAM): cli/saige_cli.cpp Makefile $(STAMP)
	@mkdir -p $(BINDIR)
	$(CLI_CXX) $(CLI_CXXFLAGS) $(CLI_INC) \
	    -DSAIGE_CLI_NAME='"$(PROGRAM)"' -DSAIGE_CLI_VERSION='"$(VERSION)"' \
	    -DSAIGE_CLI_BUILD='"$(BUILD_DESC)"' \
	    cli/saige_cli.cpp $(CLI_LIBS) -o $@

clean:
	$(MAKE) -C $(S1DIR) clean NVCC=none
	$(MAKE) -C $(S2DIR) clean
	rm -rf $(BINDIR) $(STAMP)

FORCE:
