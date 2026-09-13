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
#include <cmath>
#include <locale>
#include <ios>

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

// ---------------------------------------------------------------------------
// R-compatible spelling of non-finite doubles in the association output.
//
// C++ iostreams write "inf" / "-inf" / "nan"; R writes "Inf" / "-Inf" / "NaN".
// This is not cosmetic: read.table()/fread() parse "Inf" as infinity but turn
// "inf" into NA (or a factor level), so an output file that spells them the C++
// way silently loses those cells for anyone reading it back in R -- which is
// every downstream SAIGE user.  Non-finite cells are rare but real: SE_Burden
// is Inf whenever Pvalue_Burden is exactly 1 (R does abs(BETA/qnorm(p/2)) with
// qnorm(0.5) == 0, SAIGE_SPATest_Region_Func.R:343,351).
//
// Installing this as a num_put facet rather than wrapping ~40 individual
// stream inserts keeps the finite path byte-for-byte identical -- finite values
// are delegated straight back to the standard facet with the stream's own
// precision/format flags -- and makes it impossible for a newly added output
// column to miss the treatment.
class RNumPut : public std::num_put<char> {
protected:
    static iter_type emit(iter_type out, const char* s) {
        while (*s) { *out = *s; ++out; ++s; }
        return out;
    }
    iter_type do_put(iter_type out, std::ios_base& str, char_type fill,
                     double v) const override {
        if (std::isnan(v)) return emit(out, "NaN");
        if (std::isinf(v)) return emit(out, v > 0 ? "Inf" : "-Inf");
        return std::num_put<char>::do_put(out, str, fill, v);
    }
    iter_type do_put(iter_type out, std::ios_base& str, char_type fill,
                     long double v) const override {
        if (std::isnan(v)) return emit(out, "NaN");
        if (std::isinf(v)) return emit(out, v > 0 ? "Inf" : "-Inf");
        return std::num_put<char>::do_put(out, str, fill, v);
    }
};

// Call once, before any std::ofstream is constructed.
void installRNumericLocale();

#endif
