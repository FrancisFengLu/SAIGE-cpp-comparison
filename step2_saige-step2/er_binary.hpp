// Standalone port of SAIGE Efficient Resampling (ER) for binary traits
// Ported from:
//   SAIGE/src/ER_binary_func.cpp
//   SAIGE/src/Binary_HyperGeo.cpp / .hpp
//   SAIGE/src/Binary_ComputeExact.cpp / .hpp
//   SAIGE/src/Binary_global.cpp / .hpp
//   SAIGE/src/Binary_resampling.cpp
//   SAIGE/src/Binary_ComputeExactSKATO.cpp
//
// All Rcpp/R dependencies removed; algorithm logic preserved exactly.

#ifndef ER_BINARY_HPP
#define ER_BINARY_HPP

#include <armadillo>
#include <vector>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <cstdio>
#include <numeric>

// ============================================================
//  Global helper functions (from Binary_global.cpp/.hpp)
// ============================================================
namespace ER {

void * SL_calloc(size_t num, size_t size);
void   SL_free(void * ptr);
double SL_runif_double();
int    SL_runif_INT(int max);
void   SL_setseed(int seed);
void   SL_out();
void   SL_Sample(int k, int n, std::vector<int> & y, std::vector<int> & x);
void   SL_GetSample(int n, int k, std::vector<int> & y, std::vector<int> & x);
void   SL_GetPermu(int n, std::vector<int> & y, std::vector<int> & x);
void   SL_Binary_Boot1(int n, int ncase, std::vector<double> & pcase,
                        std::vector<int> & buf1, std::vector<int> & buf2,
                        std::vector<int> & z_one, int *err);

// ============================================================
//  HyperGeo class (from Binary_HyperGeo.cpp/.hpp)
// ============================================================
class HyperGeo {
public:
    HyperGeo() {};
    ~HyperGeo();

    int    Run(int k, int ngroup, int ncase, int * group, double * weight);
    int    Get_lprob(double * prob);
    int    Print();
    double lCombinations(int n, int k);

protected:
    double GetLogProb(int idx, int i);
    int    SaveProb(double lprob, int ncase_used);
    int    Recursive(double lprob, int idx, int ncase_used);

    int m_ngroup;
    int m_ncase;
    std::vector<int>    m_group;
    std::vector<double> m_lweight;
    std::vector<double> m_kprob;
    int m_k;

    std::vector<double *> m_probtbl;

    double m_ref;
};

// ============================================================
//  ComputeExact class (from Binary_ComputeExact.cpp/.hpp)
// ============================================================
class ComputeExact {
public:
    ComputeExact();

    int Init(std::vector<int> resarray, int nres, int * nres_k,
             double * Z0, double *Z1, int k, int m, int total,
             int * total_k, double *prob_k, double * odds, double * p1,
             int * IsExact, double epsilon, bool IsSmallmemory=0);

    int Run(int test_type);
    int GetPvalues(double * pval, double * pval_same, double * prob_k, double * minP = NULL);
    int PrintPval();

protected:
    int SaveParam(double * Z0, double *Z1, int k, int m, int total,
                  int * total_k, double *prob_k, double * odds, double * p1,
                  int * IsExact, double epsilon, bool IsSmallmemory);
    virtual double CalTestStat(int k, int * array, bool is_save=true,
                               bool is_minIdx = false, int * minIdx = NULL);
    virtual double CalTestStat_INV(int k, int * array, bool is_save=true,
                                   bool is_minIdx = false, int * minIdx = NULL);

    int CalFisherProb(int k, std::vector<int> & array);
    int SKAT_Exact_Recurse(int k, std::vector<int> & array, int cell, int start, int end);
    int SKAT_Resampling(int k, std::vector<int> & array);
    int SKAT_Resampling_Random(int k, std::vector<int> & array);

    int CalFisherProb_INV(int k, std::vector<int> & array);
    int SKAT_Exact_Recurse_INV(int k, std::vector<int> & array, int cell, int start, int end);

protected:
    std::vector<double> m_fprob;
    std::vector<double> m_teststat;
    std::vector<double> m_Z0;
    std::vector<double> m_Z1;

    std::vector<double> m_teststat_one;
    std::vector<double> m_teststat_Z0;
    std::vector<double> m_teststat_Z1;

    int m_k;
    int m_m;
    int m_total;
    double m_pprod;
    std::vector<int>    m_total_k;
    std::vector<double> m_prob_k;
    std::vector<double> m_odds;
    std::vector<double> m_p1;
    std::vector<double> m_p1_inv;
    std::vector<double> m_Q;
    std::vector<double> m_denomi;
    std::vector<int>    m_IsExact;

    std::vector<double> m_logOdds;

    int m_idx;
    std::vector<int> m_temp_x;
    std::vector<int> m_temp_x1;

    // Output
    std::vector<double> m_pval;
    std::vector<double> m_pval_same;

    bool m_IsSmallmemory;

    // Investigate possible min p-value
    double m_LargestQ;
    double m_minP;

    // for debug
    std::vector<double> m_pr1_debug;

    // precision
    double m_epsilon;
};

// ============================================================
//  SKAT_Exact global interface function
// ============================================================
void SKAT_Exact(std::vector<int> & resarray, int nres, int * nres_k,
                double * Z0, double *Z1, int k, int m, int total,
                int * total_k, double *prob_k, double * odds, double * p1,
                int * IsExact, double * pval, double *pval_same, double *minP,
                int test_type, double epsilon);

// ============================================================
//  ER_binary_func top-level functions (from ER_binary_func.cpp)
// ============================================================
void GetProb_new(int k, int ngroup, int ncase, int* group, double* weight, double* prob);

void SKATExactBin_ComputeProb_Group(arma::uvec & idx, arma::uvec & idxCompVec,
                                    arma::vec & pi1, uint32_t n, uint32_t ncase,
                                    int type_group, std::vector<double> & prob);

int  fact(int n);
int  n_choose_r(int n, int r);
void Get_Total_K(int k, std::vector<int> & n_total_k);

void SKATExactBin_ComputProb_New(arma::uvec & idx, arma::uvec & idxCompVec,
                                 arma::vec & pi1, uint32_t n, uint32_t ncase,
                                 int NResampling, int ExactMax, int test_type,
                                 int type_group,
                                 std::vector<double> & prob,
                                 std::vector<int> & IsExactVec,
                                 std::vector<int> & n_total_k,
                                 int & n_total, bool & Is_ExactP);

void Get_Res_Arrays(arma::mat & res_out, arma::uvec & idx,
                    std::vector<int> & resarray, int & nres,
                    std::vector<int> & nres_k);

double SKATExactBin_Work(arma::mat & Z, arma::vec & res, arma::vec & pi1,
                         uint32_t ncase, arma::uvec & idx, arma::uvec & idxCompVec,
                         arma::mat & res_out,
                         int NResampling, int ExactMax, double epsilon, int test_type);

} // namespace ER

#endif // ER_BINARY_HPP
