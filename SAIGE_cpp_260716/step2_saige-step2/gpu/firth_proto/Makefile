# Same optimisation flags as the production step-2 Makefile so the verbatim legacy
# fit evaluates with the same contraction/vectorisation as saige-step2.
CXX ?= g++
CXXFLAGS := -std=c++17 -O3 -march=native -funroll-loops -DARMA_NO_DEBUG -DARMA_64BIT_WORD=1
ARMA_CFLAGS := $(shell pkg-config --cflags armadillo 2>/dev/null)
ARMA_LIBS   := $(shell pkg-config --libs armadillo 2>/dev/null)
firth_ref: firth_ref.cpp firth_moment.hpp
	$(CXX) $(CXXFLAGS) $(ARMA_CFLAGS) firth_ref.cpp -o $@ $(ARMA_LIBS)
clean:
	rm -f firth_ref
