#include <RcppArmadillo.h>
#include "saige_shims.hpp"

// Forward to legacy non-const functions by copying to locals.
// This avoids const-binding errors and centralizes the friction.

namespace saige {

arma::fvec getPCG1ofSigmaAndVector(const arma::fvec& w,
                                   const arma::fvec& tau,
                                   const arma::fvec& v,
                                   int maxiterPCG, float tolPCG)
{
    arma::fvec wc = w, tc = tau, vc = v;  // make non-const locals
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
