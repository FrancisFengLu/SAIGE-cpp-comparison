UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S),Darwin)
  CXX ?= clang++
else
  CXX ?= g++
endif

CXXFLAGS ?= -std=c++17 -O3 -DARMA_USE_SUPERLU=1 -DARMA_64BIT_WORD=1 -MMD -MP -fopenmp -pthread

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
          -lz
else
  LIBS := -lyaml-cpp \
          -larmadillo \
          -lopenblas \
          -llapack \
          -lsuperlu \
          -L/usr/local/lib -lhts \
          -lzstd \
          -lz \
          -lsqlite3
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

SOURCES := $(wildcard *.cpp)
OBJECTS := $(SOURCES:.cpp=.o)
DEPS    := $(OBJECTS:.o=.d)

TARGET  ?= saige-step2

all: $(TARGET)

$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LIBS) -fopenmp -pthread -o $@

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

clean:
	rm -f $(OBJECTS) $(DEPS) $(TARGET)

test_region: $(TARGET)
	./$(TARGET) config_compare_region.yaml

test_ldmat: $(TARGET)
	./$(TARGET) config_ldmat.yaml

.PHONY: all clean test_region test_ldmat

-include $(DEPS)
