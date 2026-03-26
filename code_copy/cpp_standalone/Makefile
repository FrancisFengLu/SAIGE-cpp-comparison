# Platform detection
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  CXX ?= clang++
else
  CXX ?= g++
endif

CXXFLAGS ?= -std=c++17 -O2 -DARMA_64BIT_WORD=1 -DRCPP_PARALLEL_USE_TBB=0 -MMD -MP -pthread

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
  LIBS := -pthread \
          -lyaml-cpp \
          -larmadillo \
          -lopenblas \
          -llapack \
          -lsuperlu \
          -ltbb
  # Conda/mamba env (if active, add its library paths)
  ifneq ($(strip $(CONDA_PREFIX)),)
    LIBS := -pthread \
            -L$(CONDA_PREFIX)/lib -Wl,-rpath,$(CONDA_PREFIX)/lib \
            -lyaml-cpp -larmadillo -lopenblas -llapack -lsuperlu -ltbb
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
OBJECTS := $(SOURCES:.cpp=.o)
DEPS    := $(OBJECTS:.o=.d)

TARGET  ?= saige-null

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LIBS) -o $@

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
