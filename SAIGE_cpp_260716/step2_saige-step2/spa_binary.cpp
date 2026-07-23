// Standalone port of SAIGE/src/SPA_binary.cpp
// Saddlepoint approximation for binary traits

#include <armadillo>
#include <iostream>
#include <fstream>
#include <cassert>
#include <stdexcept>
#include <memory>
#include <sstream>
#include <time.h>
#include <stdint.h>
#include <cmath>
#include <limits>
#include <boost/math/distributions/normal.hpp>
#include "spa_binary.hpp"
#include "UTIL.hpp"


// SPA thread_local scratch (2026-06-18 refactor, re-applied 2026-06-23):
// replace per-call arma::vec locals with thread_local buffers sized once per
// thread. Called many times per SPA-firing marker (Newton-Raphson iterations);
// each call under the mallopt fix was triggering 1–4 mmap/munmap syscalls +
// ~750 pf each. Bit-identical results.

double Korg_Binom(double t1, arma::vec & mu, arma::vec & g)
{
	thread_local arma::vec _kb_temp;
	const arma::uword N = mu.n_elem;
	if (_kb_temp.n_elem != N) _kb_temp.set_size(N);
	_kb_temp = arma::log(1 - mu + mu % arma::exp(g * t1));
	return arma::sum(_kb_temp);
}


double K1_adj_Binom(double t1, arma::vec & mu, arma::vec & g, double q)
{
	thread_local arma::vec _k1_temp1, _k1_temp2, _k1_temp3;
	const arma::uword N = mu.n_elem;
	if (_k1_temp1.n_elem != N) { _k1_temp1.set_size(N); _k1_temp2.set_size(N); _k1_temp3.set_size(N); }
	_k1_temp1 = (1 - mu) % arma::exp(-g * t1) + mu;
	_k1_temp2 = mu % g;
	_k1_temp3 = _k1_temp2 / _k1_temp1;
	return arma::sum(_k1_temp3) - q;
}


double K2_Binom(double t1, arma::vec & mu, arma::vec & g)
{
	thread_local arma::vec _k2_temp0, _k2_temp1, _k2_temp2, _k2_temp3;
	const arma::uword N = mu.n_elem;
	if (_k2_temp0.n_elem != N) { _k2_temp0.set_size(N); _k2_temp1.set_size(N); _k2_temp2.set_size(N); _k2_temp3.set_size(N); }
	_k2_temp0 = arma::exp(-g * t1);
	_k2_temp1 = (1 - mu) % _k2_temp0 + mu;
	_k2_temp1 = arma::pow(_k2_temp1, 2);
	_k2_temp2 = arma::pow(g, 2) % _k2_temp0;
	_k2_temp2 = (1 - mu) % mu % _k2_temp2;
	_k2_temp3 = _k2_temp2 / _k2_temp1;
	return sum_arma1(_k2_temp3);
}


RootResult getroot_K1_Binom(double init, arma::vec & mu, arma::vec & g, double q, double tol, int maxiter){
	double root;
	int niter;
	bool Isconverge;
	double K1_eval, K2_eval, t, tnew, newK1;
	double prevJump;
	double gpos = arma::accu( g.elem( find(g > 0) ) );
	double gneg = arma::accu( g.elem( find(g < 0) ) );
	if(q >= gpos || q <= gneg){
		root = std::numeric_limits<double>::infinity();
		niter = 0;
		Isconverge = true;
	} else{
		t = init;
		K1_eval = K1_adj_Binom(t,mu,g,q);
		prevJump = std::numeric_limits<double>::infinity();
		int rep = 1;
		bool conv = true;
		while(rep <= maxiter){
			K2_eval = K2_Binom(t,mu,g);
			// P2 fix (2026-05-09): K2_eval ≈ 0 when MAC ≈ 0 (no carriers — g=0 vector
			// → temp2=0 everywhere → sum=0). Without this guard, K1_eval/K2_eval
			// would produce ±inf (not NaN, so the std::isnan check below misses it),
			// leading to SIGFPE downstream when inf values feed into chi-square /
			// pchisq calls with FE_INVALID enabled. Treat as non-convergence → skip
			// SPA for this marker; the score-test p-value (uncorrected) is used.
			if(!std::isfinite(K2_eval) || std::abs(K2_eval) < 1e-15){
				conv = false;
				break;
			}
			tnew = t-K1_eval/K2_eval;
			if(!std::isfinite(tnew)){
				conv = false;
				break;
			}

			if(std::abs(tnew-t)<tol){
				conv = true;
				break;
			}

			if(rep == maxiter)
                        {
                                conv = false;
                                break;
                        }

			newK1 = K1_adj_Binom(tnew,mu,g,q);
                        if(arma::sign(K1_eval) != arma::sign(newK1))
                        {
                                if(std::abs(tnew-t) > (prevJump-tol))
                                {
                                        tnew = t + (arma::sign(newK1-K1_eval))*prevJump/2;
                                        newK1 = K1_adj_Binom(tnew,mu,g,q);
                                        prevJump = prevJump/2;
                                } else {
                                        prevJump = std::abs(tnew-t);
                                }
                        }

			rep = rep + 1;
			t = tnew;
                        K1_eval = newK1;
		}
		root=t;
		niter=rep;
		Isconverge=conv;
	}
	return RootResult{root, niter, Isconverge};
}



SaddleResult Get_Saddle_Prob_Binom(double zeta, arma::vec & mu, arma::vec & g, double q, bool logp)
{
	double k1 = Korg_Binom(zeta, mu, g);
	double k2 = K2_Binom(zeta, mu, g);
	double temp1, w, v, Ztest, pval;
	double negative_infinity = - std::numeric_limits<double>::infinity();

	temp1 = zeta * q - k1;
	bool isSaddle = false;

        bool flagrun=false;
	if(std::isfinite(k1) && std::isfinite(k2) && temp1 >= 0 && k2 >= 0){
		 w = arma::sign(zeta) * std::sqrt(2 *temp1);
		 v = zeta *  std::sqrt(k2);
		 if(w != 0){
			flagrun = true;
		 }
	}

	if(flagrun)
	{
		Ztest = w + (1/w) * std::log(v/w);

		boost::math::normal norm_dist(0,1);
		double pval0;
	        if(Ztest > 0){
			if(logp){
				// R::pnorm(Ztest,0,1,false,true) = log(1 - Phi(Ztest))
				pval0 = std::log(boost::math::cdf(complement(norm_dist, Ztest)));
			}else{
				pval0 = boost::math::cdf(complement(norm_dist, Ztest));
			}
                        pval = pval0;
                }else {
			if(logp){
				// R::pnorm(Ztest,0,1,true,true) = log(Phi(Ztest))
				pval0 = std::log(boost::math::cdf(norm_dist, Ztest));
			}else{
				pval0 = boost::math::cdf(norm_dist, Ztest);
			}
                        pval = -pval0;
                }

		isSaddle = true;
	} else {
			if(logp)
			{
				pval =  negative_infinity;
			}else{
				pval= 0;
			}
	}
	return SaddleResult{pval, isSaddle};
}



SPAResult SPA_binary(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp){
	double p1, p2, pval;
	bool Isconverge = true;
	RootResult outuni1 = getroot_K1_Binom(0, mu, g, q, tol);
	RootResult outuni2 = getroot_K1_Binom(0, mu, g, qinv, tol);

	if(outuni1.Isconverge && outuni2.Isconverge)
	{
		SaddleResult getSaddle = Get_Saddle_Prob_Binom(outuni1.root, mu, g, q, logp);

		if(getSaddle.isSaddle){
			p1 = getSaddle.pval;
		}else{
			if(logp){
				p1 = pval_noadj-std::log(2);
			}else{
				p1 = pval_noadj/2;
			}
		}
		SaddleResult getSaddle2 = Get_Saddle_Prob_Binom(outuni2.root, mu, g, qinv, logp);
		if(getSaddle2.isSaddle){
                        p2 = getSaddle2.pval;
                }else{
			if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
		}

		if(logp)
		{
			pval = add_logp(p1,p2);
		} else {
			pval = std::abs(p1)+std::abs(p2);
		}
		Isconverge=true;
	}else {
			pval = pval_noadj;
			Isconverge=false;
		}
	return SPAResult{pval, Isconverge};
}


// 2026-06-18 refactor: same thread_local scratch treatment as the non-fast variants.
// Note: temps are sized to muNB.n_elem (the "below-cutoff" subset used by SPA fast path),
// not mu.n_elem.

// SPA fast variants — same thread_local scratch treatment; temps are sized to
// muNB.n_elem (the "below-cutoff" subset used by SPA fast path).

double Korg_fast_Binom(double t1, arma::vec & mu, arma::vec & g, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma)
{
	thread_local arma::vec _kbf_temp;
	const arma::uword N = muNB.n_elem;
	if (_kbf_temp.n_elem != N) _kbf_temp.set_size(N);
	_kbf_temp = arma::log(1 - muNB + muNB % (arma::exp(gNB * t1)));
	return arma::sum(_kbf_temp) + NAmu * t1 + 0.5 * NAsigma * pow(t1, 2);
}


double K1_adj_fast_Binom(double t1, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma)
{
	thread_local arma::vec _k1f_temp1, _k1f_temp2, _k1f_temp4;
	const arma::uword N = muNB.n_elem;
	if (_k1f_temp1.n_elem != N) { _k1f_temp1.set_size(N); _k1f_temp2.set_size(N); _k1f_temp4.set_size(N); }
	_k1f_temp1 = (1 - muNB) % arma::exp(-gNB * t1) + muNB;
	_k1f_temp2 = muNB % gNB;
	double temp3 = NAmu + NAsigma * t1;
	_k1f_temp4 = _k1f_temp2 / _k1f_temp1;
	return arma::sum(_k1f_temp4) + temp3 - q;
}



double K2_fast_Binom(double t1, arma::vec & mu, arma::vec & g, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma)
{
	thread_local arma::vec _k2f_temp0, _k2f_temp1, _k2f_temp2, _k2f_temp3;
	const arma::uword N = muNB.n_elem;
	if (_k2f_temp0.n_elem != N) { _k2f_temp0.set_size(N); _k2f_temp1.set_size(N); _k2f_temp2.set_size(N); _k2f_temp3.set_size(N); }
	_k2f_temp0 = arma::exp(-gNB * t1);
	_k2f_temp1 = (1 - muNB) % _k2f_temp0 + muNB;
	_k2f_temp1 = pow(_k2f_temp1, 2);
	_k2f_temp2 = arma::pow(gNB, 2) % _k2f_temp0;
	_k2f_temp2 = (1 - muNB) % muNB % _k2f_temp2;
	_k2f_temp3 = _k2f_temp2 / _k2f_temp1;
	return sum_arma1(_k2f_temp3) + NAsigma;
}


RootResult getroot_K1_fast_Binom(double init, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, double tol, int maxiter){
	double root;
	int niter;
	bool Isconverge;
	double K1_eval, K2_eval, t, tnew, newK1;
	double prevJump;
	double gpos = arma::accu( g.elem( find(g > 0) ) );
	double gneg = arma::accu( g.elem( find(g < 0) ) );
	if(q >= gpos || q <= gneg){
		root = std::numeric_limits<double>::infinity();
		niter = 0;
		Isconverge = true;
	} else{
		t = init;
		K1_eval = K1_adj_fast_Binom(t,mu,g,q,gNA,gNB,muNA,muNB,NAmu, NAsigma);
		prevJump = std::numeric_limits<double>::infinity();
		int rep = 1;
		bool conv = true;
		while(rep <= maxiter){
			K2_eval = K2_fast_Binom(t,mu,g, gNA,gNB,muNA,muNB,NAmu, NAsigma);
			tnew = t-K1_eval/K2_eval;
			if(std::isnan(tnew)){
				conv = false;
				break;
			}

			if(std::abs(tnew-t)<tol){
				conv = true;
				break;
			}

			if(rep == maxiter)
                        {
                                conv = false;
                                break;
                        }

			newK1 = K1_adj_fast_Binom(tnew,mu,g,q, gNA,gNB,muNA,muNB,NAmu, NAsigma);
                        if((K1_eval * newK1) < 0)
                        {
                                if(std::abs(tnew-t) > (prevJump-tol))
                                {
                                        tnew = t + (arma::sign(newK1-K1_eval))*prevJump/2;
                                        newK1 = K1_adj_fast_Binom(tnew,mu,g,q, gNA,gNB,muNA,muNB,NAmu, NAsigma);
                                        prevJump = prevJump/2;
                                } else {
                                        prevJump = std::abs(tnew-t);
                                }
                        }

			rep = rep + 1;
			t = tnew;
                        K1_eval = newK1;
		}
		root=t;
		niter=rep;
		Isconverge=conv;
	}
	return RootResult{root, niter, Isconverge};
}



SaddleResult Get_Saddle_Prob_fast_Binom(double zeta, arma::vec & mu, arma::vec & g, double q, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, bool logp)
{
	double k1 = Korg_fast_Binom(zeta, mu, g, gNA,gNB,muNA,muNB,NAmu, NAsigma);
	double k2 = K2_fast_Binom(zeta, mu, g, gNA,gNB,muNA,muNB,NAmu, NAsigma);
	double temp1, w, v, Ztest, pval;
	double negative_infinity = - std::numeric_limits<double>::infinity();
	temp1 = zeta * q - k1;
	bool isSaddle = false;


	bool flagrun=false;
        if(std::isfinite(k1) && std::isfinite(k2) && temp1 >= 0 && k2 >= 0){
                 w = arma::sign(zeta) * std::sqrt(2 *temp1);
                 v = zeta *  std::sqrt(k2);
                 if(w != 0){
                        flagrun = true;
                 }
        }


	if(flagrun)
	{
		Ztest = w + (1/w) * std::log(v/w);

		boost::math::normal norm_dist(0,1);
                double pval0;

		if(Ztest > 0){
			if(logp){
				pval0 = std::log(boost::math::cdf(complement(norm_dist, Ztest)));
			}else{
				pval0 = boost::math::cdf(complement(norm_dist, Ztest));
			}
			pval=pval0;
		}else {
			if(logp){
				pval0 = std::log(boost::math::cdf(norm_dist, Ztest));
			}else{
				pval0 = boost::math::cdf(norm_dist, Ztest);
			}
			pval= -pval0;
		}
		isSaddle = true;
	}else{
		if(logp)
		{
			pval =  negative_infinity;
		}else {
			pval=0;
		}
	}
	return SaddleResult{pval, isSaddle};
}



SPAResult SPA_binary_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, double tol){
	double p1, p2, pval;
	bool Isconverge = true;
	RootResult outuni1 = getroot_K1_fast_Binom(0, mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
	RootResult outuni2 = getroot_K1_fast_Binom(0, mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);

	if(outuni1.Isconverge && outuni2.Isconverge)
	{
		SaddleResult getSaddle = Get_Saddle_Prob_fast_Binom(outuni1.root, mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
		if(getSaddle.isSaddle){
			p1 = getSaddle.pval;
		}else{
		        if(logp){
				p1 = pval_noadj-std::log(2);
			}else{
				p1 = pval_noadj/2;
			}
		}

		SaddleResult getSaddle2 = Get_Saddle_Prob_fast_Binom(outuni2.root, mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
		if(getSaddle2.isSaddle){
			p2 = getSaddle2.pval;
		}else{
			if(logp){
				p2 = pval_noadj-std::log(2);
			}else{
				p2 = pval_noadj/2;
			}
		}

		// NOTE: Removed debug print statements from original:
		// std::cout << "p1  first " << p1 << "p2 " << p2 << std::endl;
		// std::cout << "HEREHERE " << std::endl;
		if(logp){
			pval = add_logp(p1,p2);
		} else {
			pval = std::abs(p1)+std::abs(p2);
		}
		Isconverge=true;
	}else {
			pval = pval_noadj;
			Isconverge=false;
		}
	return SPAResult{pval, Isconverge};
}
