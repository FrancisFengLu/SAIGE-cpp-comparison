// Standalone port of SAIGE/src/SPA.cpp
// SPA dispatcher: routes to binary (and eventually survival) SPA functions

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
#include <boost/math/distributions/normal.hpp>
#include "spa_binary.hpp"
// TODO: #include "spa_survival.hpp" when survival trait is ported
#include "UTIL.hpp"


void SPA(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, double & pval, bool & isSPAConverge){
        double p1, p2;
        bool Isconverge = true;
	RootResult outuni1;
	RootResult outuni2;
        if( traitType == "binary"){
          outuni1 = getroot_K1_Binom(0, mu, g, q, tol);
          outuni2 = getroot_K1_Binom(0, mu, g, qinv, tol);
        }else if(traitType == "survival"){
          // TODO: Port survival SPA (getroot_K1_Poi)
          throw std::runtime_error("survival trait SPA not yet implemented in standalone");
        }

        SaddleResult getSaddle;
        SaddleResult getSaddle2;
        if(outuni1.Isconverge && outuni2.Isconverge)
        {
                if( traitType == "binary"){
                  getSaddle = Get_Saddle_Prob_Binom(outuni1.root, mu, g, q, logp);
                  getSaddle2 = Get_Saddle_Prob_Binom(outuni2.root, mu, g, qinv, logp);
                }else if(traitType == "survival"){
                  // TODO: Port survival SPA (Get_Saddle_Prob_Poi)
                  throw std::runtime_error("survival trait SPA not yet implemented in standalone");
                }


                if(getSaddle.isSaddle){
                        p1 = getSaddle.pval;
                }else{
			Isconverge = false;
                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }
                if(getSaddle2.isSaddle){
                        p2 = getSaddle2.pval;
                }else{
			Isconverge = false;
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
        }else {
                        pval = pval_noadj;
                        Isconverge=false;
                }
        isSPAConverge = Isconverge;
}


void SPA_fast(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, bool logp, arma::vec & gNA, arma::vec & gNB, arma::vec & muNA, arma::vec & muNB, double NAmu, double NAsigma, double tol, std::string traitType, double & pval, bool & isSPAConverge){

        double p1, p2;
        bool Isconverge = true;
	RootResult outuni1;
        RootResult outuni2;

        if( traitType == "binary"){
          outuni1 = getroot_K1_fast_Binom(0, mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
          outuni2 = getroot_K1_fast_Binom(0, mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, tol);
        }else if(traitType == "survival"){
          // TODO: Port survival SPA (getroot_K1_fast_Poi)
          throw std::runtime_error("survival trait SPA_fast not yet implemented in standalone");
        }

        SaddleResult getSaddle;
        SaddleResult getSaddle2;
        if(outuni1.Isconverge && outuni2.Isconverge)
        {
          if( traitType == "binary"){
                getSaddle  = Get_Saddle_Prob_fast_Binom(outuni1.root, mu, g, q, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
                getSaddle2 = Get_Saddle_Prob_fast_Binom(outuni2.root, mu, g, qinv, gNA,gNB,muNA,muNB,NAmu, NAsigma, logp);
          }else if(traitType == "survival"){
                // TODO: Port survival SPA
                throw std::runtime_error("survival trait SPA_fast not yet implemented in standalone");
          }
                if(getSaddle.isSaddle){
                        p1 = getSaddle.pval;
                }else{
			Isconverge = false;
                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }

                if(getSaddle2.isSaddle){
                        p2 = getSaddle2.pval;
                }else{
			Isconverge = false;
                        if(logp){
                                p2 = pval_noadj-std::log(2);
                        }else{
                                p2 = pval_noadj/2;
                        }
                }

                if(logp){
                        pval = add_logp(p1,p2);
                }else {
                        pval = std::abs(p1)+std::abs(p2);
                }
        }else {
                        pval = pval_noadj;
                        Isconverge=false;
        }
	isSPAConverge = Isconverge;
}




double SPA_pval(arma::vec & mu, arma::vec & g, double q, double qinv, double pval_noadj, double tol, bool logp, std::string traitType, bool & isSPAConverge){
       double pval;
        double p1, p2;
        bool Isconverge = true;
        RootResult outuni1;
        RootResult outuni2;
        if( traitType == "binary"){
          outuni1 = getroot_K1_Binom(0, mu, g, q, tol);
          outuni2 = getroot_K1_Binom(0, mu, g, qinv, tol);
        }else if(traitType == "survival"){
          // TODO: Port survival SPA
          throw std::runtime_error("survival trait SPA_pval not yet implemented in standalone");
        }


        double outuni1root = outuni1.root;
        double outuni2root = outuni2.root;
        bool Isconverge1 = outuni1.Isconverge;
        bool Isconverge2 = outuni2.Isconverge;

        SaddleResult getSaddle;
        SaddleResult getSaddle2;
        if(outuni1.Isconverge && outuni2.Isconverge)
        {
                if( traitType == "binary"){
                  getSaddle = Get_Saddle_Prob_Binom(outuni1.root, mu, g, q, logp);
                  getSaddle2 = Get_Saddle_Prob_Binom(outuni2.root, mu, g, qinv, logp);
                }else if(traitType == "survival"){
                  // TODO: Port survival SPA
                  throw std::runtime_error("survival trait SPA_pval not yet implemented in standalone");
                }

                if(getSaddle.isSaddle){
                        p1 = getSaddle.pval;
                }else{
                        Isconverge = false;
                        if(logp){
                                p1 = pval_noadj-std::log(2);
                        }else{
                                p1 = pval_noadj/2;
                        }
                }
                if(getSaddle2.isSaddle){
                        p2 = getSaddle2.pval;
                }else{
                        Isconverge = false;
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
        }else {
                        pval = pval_noadj;
                        Isconverge=false;
                }
        isSPAConverge = Isconverge;

	return(pval);
}
