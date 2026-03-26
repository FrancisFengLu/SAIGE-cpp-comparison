#include <RcppArmadillo.h>
#include "saige_shims.hpp"
#include <cstdio>
#include <unistd.h>

// Forward to legacy non-const functions by copying to locals.
// This avoids const-binding errors and centralizes the friction.

namespace saige {

arma::fvec getPCG1ofSigmaAndVector(const arma::fvec& w,
                                   const arma::fvec& tau,
                                   const arma::fvec& v,
                                   int maxiterPCG, float tolPCG)
{
    { const char _m[] = "[SHIM] PCG enter\n"; write(2, _m, sizeof(_m)-1); }
    fprintf(stderr, "[SHIM] w.n=%zu tau.n=%zu v.n=%zu\n",
            (size_t)w.n_elem, (size_t)tau.n_elem, (size_t)v.n_elem); fflush(stderr);
    arma::fvec wc = w, tc = tau, vc = v;  // make non-const locals
    fprintf(stderr, "[SHIM] copies done, calling global PCG\n"); fflush(stderr);
    return ::getPCG1ofSigmaAndVector(wc, tc, vc, maxiterPCG, tolPCG);
}

arma::fvec getPCG1ofSigmaAndVector_LOCO(const arma::fvec& w,
                                        const arma::fvec& tau,
                                        const arma::fvec& v,
                                        int maxiterPCG, float tolPCG)
{
    arma::fvec wc = w, tc = tau, vc = v;  // make non-const locals
    return ::getPCG1ofSigmaAndVector_LOCO(wc, tc, vc, maxiterPCG, tolPCG);
}

} // namespace saige
