#pragma once
// Armadillo with the configuration step 1 was built with under RcppArmadillo,
// without R. Include this instead of <armadillo> / <RcppArmadillo.h>.
//
// Before: every TU picked up RcppArmadillo's bundled Armadillo 15.4.2 (its
// -I path shadowed the conda env's, which is -isystem). Its config.hpp is the
// stock one with WRAPPER / ARPACK / SUPERLU off; the conda env's Armadillo
// (same 15.4.2) has them on, so they are switched back off here.
// SAIGE_step1_fast.cpp defines ARMA_USE_SUPERLU itself and keeps SuperLU.
//
// What RcppArmadillo added on top and what became of it:
//   ARMA_USE_LAPACK / ARMA_USE_BLAS / ARMA_USE_OPENMP   kept (below)
//   ARMA_DONT_USE_WRAPPER                               kept (BLAS/LAPACK called directly)
//   ARMA_32BIT_WORD unless ARMA_64BIT_WORD              kept (Makefile sets ARMA_64BIT_WORD)
//   ARMA_RNG_ALT (R's RNG for arma::randu/randn)        dropped: step 1 never calls Armadillo's RNG
//   ARMA_COUT/CERR_STREAM = Rcpp::Rcout/Rcerr           dropped: std::cout / std::cerr
//   ARMA_EXTRA_*_PROTO/MEAT (SEXP constructors)         dropped: no R objects

#ifndef ARMA_USE_LAPACK
#define ARMA_USE_LAPACK
#endif
#ifndef ARMA_USE_BLAS
#define ARMA_USE_BLAS
#endif
#ifndef ARMA_DONT_USE_WRAPPER
#define ARMA_DONT_USE_WRAPPER
#endif
#ifndef ARMA_DONT_USE_ARPACK
#define ARMA_DONT_USE_ARPACK
#endif
#if !defined(ARMA_USE_SUPERLU) && !defined(ARMA_DONT_USE_SUPERLU)
#define ARMA_DONT_USE_SUPERLU
#endif
#if defined(_OPENMP) && !defined(ARMA_USE_OPENMP) && !defined(ARMA_DONT_USE_OPENMP)
#define ARMA_USE_OPENMP 1
#endif
#ifndef ARMA_DONT_PRINT_OPENMP_WARNING
#define ARMA_DONT_PRINT_OPENMP_WARNING 1
#endif
#if !defined(ARMA_64BIT_WORD) && !defined(ARMA_32BIT_WORD)
#define ARMA_32BIT_WORD 1
#endif

#include <armadillo>
