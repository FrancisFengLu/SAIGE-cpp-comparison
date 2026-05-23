// Standalone port of UTIL.hpp for Step 2
// Shared utility functions

#ifndef UTIL_HPP
#define UTIL_HPP

#include <armadillo>
#include <sys/time.h>
#include <unordered_map>
#include <string>
#include <vector>
#include <random>

const static std::unordered_map<std::string,int> string_to_case{
   {"best_guess",1},
   {"mean",2},
   {"minor",3}
};

double getWeights(std::string t_kernel,
                  double t_freq,
                  arma::vec t_wBeta);

void imputeGeno(arma::vec& GVec,
                double freq,
                std::vector<uint32_t> posMissingGeno);

double getInnerProd(arma::mat& x1Mat, arma::mat& x2Mat);


bool imputeGenoAndFlip(arma::vec& t_GVec,
                       double & t_altFreq,
		       double & t_altCount,
                       std::vector<uint32_t> &  t_indexForMissing,
                       std::string t_impute_method,
                       double t_dosage_zerod_cutoff,
                       double t_dosage_zerod_MAC_cutoff,
                       double & t_MAC,
		       std::vector<uint> & t_indexZero,
                       std::vector<uint> & t_indexNonZero);

arma::vec getTime();

void printTime(arma::vec t1, arma::vec t2, std::string message);

double getinvStd(double t_freq);

// Standalone replacement for Rcpp::rbinom
// Generates n random 0/1 values with p=0.5
arma::vec nb(unsigned int n);

double sum_arma1(arma::vec& X);

double add_logp(double p1, double p2);

#endif
