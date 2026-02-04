CXX      ?= clang++
CXXFLAGS ?= -std=c++17 -O3 -DARMA_64BIT_WORD=1 -DRCPP_PARALLEL_USE_TBB=0 -MMD -MP

RSCRIPT      := $(shell command -v Rscript 2>/dev/null)
PKG_CONFIG   := $(shell command -v pkg-config 2>/dev/null)

ifeq ($(strip $(RSCRIPT)),)
  R_HOME       :=
  RCPP_INC     :=
  RCPPARM_INC  :=
  RCPPARAL_INC :=
else
  R_HOME       := $(shell "$(RSCRIPT)" -e 'cat(R.home())')
  RCPP_INC     := $(shell "$(RSCRIPT)" -e 'cat(system.file("include", package="Rcpp"))')
  RCPPARM_INC  := $(shell "$(RSCRIPT)" -e 'cat(system.file("include", package="RcppArmadillo"))')
  RCPPARAL_INC := $(shell "$(RSCRIPT)" -e 'cat(system.file("include", package="RcppParallel"))')
  RCPPEIG_INC  := $(shell "$(RSCRIPT)" -e 'cat(system.file("include", package="RcppEigen"))')
endif

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

INCLUDES := -I/usr/local/include -I/usr/include/eigen3 -I/opt/homebrew/opt/boost/include -I/opt/homebrew/include -I/opt/homebrew/var/homebrew/tmp/.cellar/eigen/5.0.1/include/eigen3

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

ifneq ($(strip $(YAML_CFLAGS)),)
  INCLUDES += $(YAML_CFLAGS)
endif
ifneq ($(strip $(ARMA_CFLAGS)),)
  INCLUDES += $(ARMA_CFLAGS)
endif

LIBS := -lyaml-cpp -larmadillo -L/opt/homebrew/opt/openblas/lib -lopenblas -llapack -L/opt/homebrew/opt/superlu/lib -lsuperlu -L/opt/homebrew/opt/tbb/lib -ltbb

ifneq ($(strip $(YAML_LIBS)),)
  LIBS += $(YAML_LIBS)
endif
ifneq ($(strip $(ARMA_LIBS)),)
  LIBS += $(ARMA_LIBS)
endif

ifneq ($(strip $(R_HOME)),)
  LIBS += -L"$(R_HOME)/lib" -lR
endif

SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)
DEPS    := $(OBJECTS:.o=.d)

TARGET  ?= saige-null

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LIBS) -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(DEPS) $(TARGET)

.PHONY: all clean

-include $(DEPS)
