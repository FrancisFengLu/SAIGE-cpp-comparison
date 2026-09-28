# step1bench -- standalone. No R, no Rcpp, no linkage against the step1/step2
# trees. The only file shared with them is block_sigma.{hpp,cpp}, copied in.
#
# Build env:  source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh
# CHOLMOD 5.3.1 + CXSparse 4.4.1 are installed into that conda env
# (saige-build), pinned to the same versions TorchGWAS2 links.

CXX      ?= g++

# `:=` on purpose: conda activate scripts export CXXFLAGS, which would
# otherwise win over `?=` and silently drop everything below.
CXXFLAGS := -std=c++17 -O3 -march=native -funroll-loops \
            -DARMA_NO_DEBUG -DARMA_64BIT_WORD=1 \
            -MMD -MP -pthread -fopenmp

PREFIX   := $(CONDA_PREFIX)
ifeq ($(strip $(PREFIX)),)
  $(error CONDA_PREFIX is empty -- source /opt/saige/logs/tg2_step2/scripts/env_cpp.sh first)
endif

INCLUDES := -I. -isystem $(PREFIX)/include
LDFLAGS  := -L$(PREFIX)/lib -Wl,-rpath,$(PREFIX)/lib
LDLIBS   := -larmadillo -lsuperlu -lcholmod -lcxsparse -lopenblas -lgomp -lpthread -lm

SRCS := main.cpp aireml.cpp bench_io.cpp sigma_solver.cpp \
        solver_block.cpp solver_cholmod.cpp solver_superlu.cpp block_sigma.cpp
OBJS := $(SRCS:.cpp=.o)
DEPS := $(OBJS:.o=.d)
BIN  := step1bench

all: $(BIN)

$(BIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@ $(LDFLAGS) $(LDLIBS)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJS) $(DEPS) $(BIN)

.PHONY: all clean
-include $(DEPS)
