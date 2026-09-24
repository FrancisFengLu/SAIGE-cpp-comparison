// mvec_probe.cpp -- build-time probe for glibc libmvec's 4-wide double erfc.
// The Makefile compiles and links this; if it succeeds it defines
// SAIGE_HAVE_MVEC_ERFC and adds -lmvec. See score_vec.hpp.
#if !defined(__x86_64__)
#error "not x86-64"
#endif
#include <immintrin.h>
extern "C" __m256d _ZGVdN4v_erfc(__m256d);
int main()
{
    __m256d z = _mm256_set1_pd(1.0);
    __m256d r = _ZGVdN4v_erfc(z);
    double o[4];
    _mm256_storeu_pd(o, r);
    return (o[0] > 0.0 && o[0] < 1.0) ? 0 : 1;
}
