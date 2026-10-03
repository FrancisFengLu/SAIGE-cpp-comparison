#include <boost/math/special_functions/erf.hpp>
#include <boost/math/tools/rational.hpp>
#include <cstdio>
#include <cmath>
#include "erfc_boost53.cuh"
using namespace saige::spa_gpu::erfc53;
int main(){
  std::printf("POLY_METHOD=%d RATIONAL_METHOD=%d\n", BOOST_MATH_POLY_METHOD, BOOST_MATH_RATIONAL_METHOD);
  // compare the polynomial evaluation scheme directly on P2 (6) and Q2 (7)
  static const double P2[6] = {ERF_P2[0],ERF_P2[1],ERF_P2[2],ERF_P2[3],ERF_P2[4],ERF_P2[5]};
  static const double Q2[7] = {ERF_Q2[0],ERF_Q2[1],ERF_Q2[2],ERF_Q2[3],ERF_Q2[4],ERF_Q2[5],ERF_Q2[6]};
  static const double P1[5] = {ERF_P1[0],ERF_P1[1],ERF_P1[2],ERF_P1[3],ERF_P1[4]};
  int bad6=0,bad7=0,bad5=0,n=0;
  for(double x=0.0;x<1.0;x+=1e-5){ n++;
    if(boost::math::tools::evaluate_polynomial(P2,x)!=poly6(P2,x)) bad6++;
    if(boost::math::tools::evaluate_polynomial(Q2,x)!=poly7(Q2,x)) bad7++;
    if(boost::math::tools::evaluate_polynomial(P1,x)!=poly5(P1,x)) bad5++; }
  std::printf("poly scheme mismatches over %d x: N=5 %d  N=6 %d  N=7 %d\n", n, bad5, bad6, bad7);
  // per-branch erfc exactness
  struct B{const char*nm; double lo,hi; long n=0,ex=0; double z0=0,h0=0,p0=0;} bs[]={{"<0.5",-6,0.5},{"[0.5,1.5)",0.5,1.5},{"[1.5,2.5)",1.5,2.5},{"[2.5,4.5)",2.5,4.5},{"[4.5,28)",4.5,28}};
  for(double x=-6;x<28;x+=1e-4){ double h=boost::math::erfc(x), p=erfImp53(x,true);
    for(auto&b:bs) if(x>=b.lo&&x<b.hi){ b.n++; if(h==p) b.ex++; else if(b.z0==0){b.z0=x;b.h0=h;b.p0=p;} } }
  for(auto&b:bs) std::printf("  %-10s n=%ld exact=%ld  first mismatch z=%.6f host=%.17g port=%.17g\n", b.nm,b.n,b.ex,b.z0,b.h0,b.p0);
  // the compensated exp piece alone vs boost's expression for z=3.3
  double z=3.3; int expon; double hi=std::floor(std::ldexp(std::frexp(z,&expon),26)); hi=std::ldexp(hi,expon-26); double lo=z-hi; double sq=z*z; double err=((hi*hi-sq)+2*hi*lo)+lo*lo;
  std::printf("z=3.3 hi=%.17g lo=%.17g err_sqr=%.17g exp(-sq)*exp(-err)/z=%.17g\n", hi, lo, err, std::exp(-sq)*std::exp(-err)/z);
  return 0; }
