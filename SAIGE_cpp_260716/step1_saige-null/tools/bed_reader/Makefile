# Standalone build for the new BED reader unit (PR-1 of PARALLEL_BED_PLAN).
# Kept separate from the main saige-null Makefile so the production build
# is never affected until we're ready to integrate.
#
#   make                 → build bed_reader_test
#   make test            → build + run against the UKB BED + mapped FAM
#   make clean           → remove build artefacts

CXX      ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pthread

# Use the conda env compiler + libs when active, matching saige-null's build.
ifneq ($(strip $(CONDA_PREFIX)),)
  CXX      := $(CONDA_PREFIX)/bin/x86_64-conda-linux-gnu-c++
  CXXFLAGS += -isystem $(CONDA_PREFIX)/include
endif

SRC_SMOKE    := bed_reader.cpp bed_reader_test.cpp
SRC_PIPELINE := bed_reader.cpp marker_decoder.cpp bed_pipeline_test.cpp
SRC_PARALLEL := bed_reader.cpp marker_decoder.cpp parallel_decode.cpp bed_parallel_test.cpp
SRC_MASK     := bed_reader.cpp marker_decoder.cpp parallel_decode.cpp mask_stats.cpp mask_stats_test.cpp

OBJ_SMOKE    := $(SRC_SMOKE:.cpp=.o)
OBJ_PIPELINE := $(SRC_PIPELINE:.cpp=.o)
OBJ_PARALLEL := $(SRC_PARALLEL:.cpp=.o)
OBJ_MASK     := $(SRC_MASK:.cpp=.o)

all: bed_reader_test bed_pipeline_test bed_parallel_test mask_stats_test

bed_reader_test: $(OBJ_SMOKE)
	$(CXX) $(CXXFLAGS) $^ -o $@

bed_pipeline_test: $(OBJ_PIPELINE)
	$(CXX) $(CXXFLAGS) $^ -o $@

bed_parallel_test: $(OBJ_PARALLEL)
	$(CXX) $(CXXFLAGS) $^ -o $@ -lpthread

# Scheme C stage A — gates A1/A2 of optimization/missing_mt/SCHEME_C_DESIGN.md.
mask_stats_test: $(OBJ_MASK)
	$(CXX) $(CXXFLAGS) $^ -o $@ -lpthread

%.o: %.cpp bed_reader.hpp marker_decoder.hpp packed_store.hpp parallel_decode.hpp mask_stats.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Default smoke test against the UKB data used this session.
UKB_BED   ?= /media/leelabsg-storage0/DATA/UKBB/cal/pruned/output/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.bed
UKB_FAM   ?= /media/leelabsg-storage0/seokho/ukb_allchr_v2_newID_passedQC_white.British_geno0.05_poly_500_50_0.2.pruned.mapped.fam
# Use the cov_input.csv that saige-null already wrote during covariate_offset.
# This file contains the exact 961 IIDs that the production run used (post
# NA-filter + FAM overlap), so our pipeline-test builds ptrsubSampleInGeno
# from the same sample set the reference log was computed against.
UKB_IIDS  ?= /media/leelabsg-storage0/seokho/seokho92-SAIGE-cpp-comparison/output/ukb_ldl_dense_output_cov_input.csv
UKB_REF   ?= /media/leelabsg-storage0/seokho/seokho92-SAIGE-cpp-comparison/code_copy/cpp_standalone/benchmark_results/ukb_ldl/cpp_covT/stdout.log

test: bed_reader_test
	./bed_reader_test $(UKB_BED) $(UKB_FAM) 0 1 2 100 340446

pipeline-test: bed_pipeline_test
	./bed_pipeline_test $(UKB_BED) $(UKB_FAM) $(UKB_IIDS) iid $(UKB_REF) 20

# PR-3 parallel correctness + scaling harness.
# n_markers=0 means the full BED (340 447 markers, full ~5 min serial pass).
# Default: just 1000 markers for a fast correctness check + indicative scaling.
PAR_N_MARKERS ?= 1000
PAR_THREADS   ?= 1,4,16
parallel-test: bed_parallel_test
	./bed_parallel_test $(UKB_BED) $(UKB_FAM) $(UKB_IIDS) iid $(UKB_REF) \
	    $(PAR_N_MARKERS) $(PAR_THREADS)

clean:
	rm -f $(OBJ_SMOKE) $(OBJ_PIPELINE) $(OBJ_PARALLEL) $(OBJ_MASK) \
	      bed_reader_test bed_pipeline_test bed_parallel_test mask_stats_test

# Scheme C stage A gate. MASK_CASE picks the tolerance-gate case whose sample
# sets the run uses; qc and fill are the two built to make per-trait QC lists
# and per-trait fill values diverge.
MASK_SCALE ?= small
MASK_CASE  ?= qc
MASK_DATA  ?= /opt/saige/logs/mt_gate/data
MASK_GENO  ?= /opt/saige/data/$(MASK_SCALE)
MASK_ARGS  ?=
mask-test: mask_stats_test
	./mask_stats_test $(MASK_GENO).bed $(MASK_GENO).fam \
	    $(MASK_DATA)/$(MASK_SCALE)/$(MASK_CASE).pheno.txt $(MASK_ARGS)

# Full stage-A gate: every case, VR off and on, plus a synthetic P=8 run where
# each trait drops its own 1-8% (the biobank-like regime none of the gate cases
# covers). Non-zero exit on the first failure.
MASK_CASES ?= qc fill tiny block nomiss indep vrcate
mask-gate: mask_stats_test
	@set -e; for c in $(MASK_CASES); do \
	  for v in "" "--vr 20,-1" "--vr 20.5,200.5"; do \
	    echo "=== $(MASK_SCALE)/$$c  VR[$$v] ==="; \
	    ./mask_stats_test $(MASK_GENO).bed $(MASK_GENO).fam \
	        $(MASK_DATA)/$(MASK_SCALE)/$$c.pheno.txt --threads 8 --reps 1 $$v 2>/dev/null \
	      | grep -E "A1 mismatches|A2 mismatches|keep-rule|pack_if_keep|missing-cell|PASS|FAIL"; \
	  done; \
	done; \
	echo "=== $(MASK_SCALE)/synth 8 traits x 1-8% ==="; \
	./mask_stats_test $(MASK_GENO).bed $(MASK_GENO).fam \
	    $(MASK_DATA)/$(MASK_SCALE)/$(firstword $(MASK_CASES)).pheno.txt \
	    --threads 8 --reps 1 --vr 20,-1 --synth 8,8,7 2>/dev/null \
	  | grep -E "A1 mismatches|A2 mismatches|keep-rule|pack_if_keep|missing-cell|PASS|FAIL"

# What the scheme-C decode outputs cost. Interleaved reps + min, so the numbers
# survive a busy machine; run it on an idle one anyway.
MASK_REPS ?= 5
mask-price: mask_stats_test
	@for c in $(MASK_CASES); do \
	  echo "=== $(MASK_SCALE)/$$c ==="; \
	  ./mask_stats_test $(MASK_GENO).bed $(MASK_GENO).fam \
	      $(MASK_DATA)/$(MASK_SCALE)/$$c.pheno.txt --threads 8 \
	      --reps $(MASK_REPS) --price-only --vr 20,-1 2>/dev/null; \
	done

.PHONY: all test pipeline-test parallel-test mask-test mask-gate mask-price clean
