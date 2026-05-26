// Standalone port of SAIGE Efficient Resampling (ER) for binary traits
// Ported from:
//   SAIGE/src/ER_binary_func.cpp
//   SAIGE/src/Binary_HyperGeo.cpp
//   SAIGE/src/Binary_ComputeExact.cpp
//   SAIGE/src/Binary_global.cpp
//   SAIGE/src/Binary_resampling.cpp
//
// All Rcpp/R dependencies removed; algorithm logic preserved exactly.

#include "er_binary.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>
#include <numeric>
#include <algorithm>
#include <cfloat>
#include <limits>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MIN_SIM1 0

namespace ER {

// ============================================================
//  Global helper functions (from Binary_global.cpp)
// ============================================================

void * SL_calloc(size_t num, size_t size){
    void * re = calloc(num, size);
    if(re == NULL){
        printf("memory allocation error!");
        return(re);
    }
    return(re);
}

void SL_free(void * ptr){
    if(ptr != NULL){
        free(ptr);
        ptr = NULL;
    }
}

// ============================================================
//  Binary_resampling functions (from Binary_resampling.cpp)
//  Standalone branch: uses rand()/srand()
// ============================================================

double SL_runif_double(){
    double val;
    val = ((double)rand()/(double)RAND_MAX);
    return val;
}

int SL_runif_INT(int max){
    int val;
    val = rand() % max;
    return val;
}

void SL_setseed(int seed){
    srand(seed);
}

void SL_out(){
    // no-op in standalone mode
}

void SL_Sample(int k, int n, std::vector<int> & y, std::vector<int> & x){
    int i, j;
    for (i = 0; i < n; i++){
        x[i] = i;
    }
    for (i = 0; i < k; i++) {
        j = SL_runif_INT(n);
        y[i] = x[j];
        x[j] = x[--n];
    }
}

void SL_GetSample(int n, int k, std::vector<int> & y, std::vector<int> & x){
    int i, j;
    for (i = 0; i < n; i++){
        x[i] = i;
    }
    for (i = 0; i < k; i++) {
        j = SL_runif_INT(n);
        y[i] = x[j];
        x[j] = x[--n];
    }
}

void SL_GetPermu(int n, std::vector<int> & y, std::vector<int> & x){
    int i, j;
    for (i = 0; i < n; i++){
        x[i] = y[i];
    }
    for (i = 0; i < n; i++) {
        j = SL_runif_INT(n);
        y[i] = x[j];
        x[j] = x[--n];
    }
}

void SL_Binary_Boot1(int n, int ncase, std::vector<double> & pcase,
                     std::vector<int> & buf1, std::vector<int> & buf2,
                     std::vector<int> & z_one, int *err){
    int i, i1, j, k, n1, ncase1;
    double temp;

    SL_GetSample(n, n, buf1, buf2);
    ncase1 = 0;
    n1 = n;
    for(k=0; k<500; k++){
        i1 = 0;
        for(i=0; i<n1; i++){
            j = buf1[i];
            temp = SL_runif_double();
            if(temp <= pcase[j]){
                z_one[j] = 1;
                ncase1++;
            } else {
                buf2[i1] = j;
                i1++;
            }
            if(ncase1 == ncase)
                break;
        }

        if(ncase1 == ncase){
            break;
        } else if(ncase1 > ncase) {
            *err = -1;
            return;
        } else {
            n1 = n - ncase1;
            memcpy(buf1.data(), buf2.data(), sizeof(int) * n1);
        }
    }

    if(ncase != ncase1){
        *err = -1;
        return;
    }

    *err = 1;
    return;
}

// ============================================================
//  HyperGeo class (from Binary_HyperGeo.cpp)
// ============================================================

HyperGeo::~HyperGeo(){
    double * ptemp;
    for(size_t i=0; i<m_probtbl.size(); i++){
        ptemp = m_probtbl[i];
        SL_free(ptemp);
    }
}

int HyperGeo::Run(int k, int ngroup, int ncase, int * group, double * weight){
    int i;
    m_ngroup = ngroup;
    m_ncase = ncase;
    m_k = k;

    for(i=0; i<ngroup; i++){
        m_group.push_back(group[i]);
        m_lweight.push_back(log(weight[i]));
    }

    for(i=0; i<=m_k; i++){
        m_kprob.push_back(0);
    }

    // generate prob table
    double * temp;
    double lch;
    m_ref = 0;
    for(i=0; i<ngroup; i++){
        if(i < ngroup - 1){
            temp = (double *) SL_calloc(m_group[i] + 1, sizeof(double));
            for(int j=0; j<=m_group[i]; j++){
                lch = lCombinations(m_group[i], j);
                temp[j] = lch + m_lweight[i]*j;
            }
        } else {
            temp = (double *) SL_calloc(k+1, sizeof(double));
            for(int j=0; j<=k; j++){
                lch = lCombinations(m_group[i], m_ncase - j);
                temp[j] = lch + m_lweight[i] * ((double) m_ncase - j);
                m_ref = MAX(temp[j], m_ref);
            }
        }
        m_probtbl.push_back(temp);
    }

    Recursive(0, 0, 0);

    return 1;
}

double HyperGeo::GetLogProb(int idx, int i){
    return m_probtbl[idx][i];
}

int HyperGeo::SaveProb(double lprob, int ncase_used){
    m_kprob[ncase_used] += exp(lprob - m_ref);
    return 1;
}

int HyperGeo::Recursive(double lprob, int idx, int ncase_used){
    int i;
    double lprob1;

    int num = m_group[idx];

    if(idx == m_ngroup - 1){
        lprob1 = GetLogProb(idx, ncase_used);
        SaveProb(lprob + lprob1, ncase_used);
        return 1;
    }

    for(i=0; i<=num; i++){
        if(ncase_used + i <= m_ncase){
            lprob1 = GetLogProb(idx, i);
            Recursive(lprob + lprob1, idx+1, ncase_used + i);
        }
    }

    return 0;
}

int HyperGeo::Get_lprob(double * prob){
    double sum1 = 0;
    for(int i=0; i<=m_k; i++){
        sum1 += m_kprob[i];
    }
    for(int i=0; i<=m_k; i++){
        prob[i] = m_kprob[i] / sum1;
    }
    return 1;
}

int HyperGeo::Print(){
    double sum1 = 0;
    double prob;
    for(int i=0; i<=m_k; i++){
        sum1 += m_kprob[i];
    }
    for(int i=0; i<=m_k; i++){
        prob = m_kprob[i] / sum1;
        printf("%d:[%e][%e]\n", i, prob, m_kprob[i]);
    }
    return 1;
}

double HyperGeo::lCombinations(int n, int k){
    if (k > n)
        return 0;
    double r = 0;
    // Standalone mode: compute log-combinations directly
    for (int d = 1; d <= k; ++d){
        r += log(n--);
        r -= log(d);
    }
    return r;
}

// ============================================================
//  ComputeExact class (from Binary_ComputeExact.cpp)
// ============================================================

double ComputeExact::CalTestStat(int k, int * array, bool is_save, bool is_minIdx, int * minIdx){
    int i, j, l, temp;
    double stat = 0;
    m_teststat_one = m_teststat_Z0;

    for(i=0; i<k; i++){
        l = array[i];
        temp = l * m_m;
        for(j=0; j<m_m; j++){
            m_teststat_one.at(j) += m_Z1.at(temp + j) - m_Z0.at(temp + j);
        }
    }

    for(j=0; j<m_m; j++){
        stat += m_teststat_one.at(j) * m_teststat_one.at(j);
    }

    if(is_save){
        m_teststat.at(m_idx) = stat;
    }
    return stat;
}

double ComputeExact::CalTestStat_INV(int k, int * array, bool is_save, bool is_minIdx, int * minIdx){
    int i, j, l, temp;
    double stat = 0;
    m_teststat_one = m_teststat_Z1;

    for(i=0; i<k; i++){
        l = array[i];
        temp = l * m_m;
        for(j=0; j<m_m; j++){
            m_teststat_one.at(j) += m_Z0.at(temp + j) - m_Z1.at(temp + j);
        }
    }

    for(j=0; j<m_m; j++){
        stat += m_teststat_one.at(j) * m_teststat_one.at(j);
    }

    if(is_save){
        m_teststat.at(m_idx) = stat;
    }
    return stat;
}

int ComputeExact::CalFisherProb(int k, std::vector<int> & array){
    int i, l;
    double temp = 1;
    for(i=0; i<k; i++){
        l = array.at(i);
        temp = temp * m_odds.at(l);
    }
    m_fprob.at(m_idx) = temp;
    m_denomi.at(k) = m_denomi.at(k) + temp;
    return 0;
}

int ComputeExact::CalFisherProb_INV(int k, std::vector<int> & array){
    int i, l, k1;
    k1 = m_k - k;
    double temp = m_pprod;
    for(i=0; i<k; i++){
        l = array.at(i);
        temp = temp / m_odds.at(l);
    }
    m_fprob.at(m_idx) = temp;
    m_denomi.at(k1) = m_denomi.at(k1) + temp;
    return 0;
}

int ComputeExact::SKAT_Exact_Recurse(int k, std::vector<int> & array, int cell, int start, int end){
    int i;
    if(k == cell){
        CalTestStat(k, array.data());
        CalFisherProb(k, array);
        m_idx++;
    } else {
        for(i=start; i<end; i++){
            array.at(cell) = i;
            SKAT_Exact_Recurse(k, array, cell+1, i+1, end);
        }
    }
    return 0;
}

int ComputeExact::SKAT_Exact_Recurse_INV(int k, std::vector<int> & array, int cell, int start, int end){
    int i;
    if(k == cell){
        CalTestStat_INV(k, array.data());
        CalFisherProb_INV(k, array);
        m_idx++;
    } else {
        for(i=start; i<end; i++){
            array.at(cell) = i;
            SKAT_Exact_Recurse_INV(k, array, cell+1, i+1, end);
        }
    }
    return 0;
}

int ComputeExact::SKAT_Resampling(int k, std::vector<int> & array){
    int k1 = m_k - k;
    if(k <= m_k/2 + 1){
        for(int i=0; i<m_total_k.at(k); i++){
            SL_Sample(k, m_k, m_temp_x, array);
            CalTestStat(k, m_temp_x.data());
            CalFisherProb(k, m_temp_x);
            m_idx++;
        }
    } else {
        for(int i=0; i<m_total_k.at(k); i++){
            SL_Sample(k1, m_k, m_temp_x, array);
            CalTestStat_INV(k1, m_temp_x.data());
            CalFisherProb_INV(k1, m_temp_x);
            m_idx++;
        }
    }
    return 1;
}

int ComputeExact::SKAT_Resampling_Random(int k, std::vector<int> & array){
    int err;
    int k1 = m_k - k;
    if(k <= m_k/2 + 1){
        for(int i=0; i<m_total_k.at(k); i++){
            SL_Binary_Boot1(m_k, k, m_p1, array, m_temp_x1, m_temp_x, &err);
            CalFisherProb(k, m_temp_x);
            m_fprob.at(m_idx) = 1;
            m_denomi.at(k) = m_denomi.at(k) + 1;
            m_idx++;
        }
    } else {
        for(int i=0; i<m_total_k.at(k); i++){
            SL_Binary_Boot1(m_k, k1, m_p1_inv, array, m_temp_x1, m_temp_x, &err);
            CalFisherProb_INV(k1, m_temp_x);
            m_fprob.at(m_idx) = 1;
            m_denomi.at(k) = m_denomi.at(k) + 1;
            m_idx++;
        }
    }
    return 1;
}

ComputeExact::ComputeExact(){
    m_pprod = 1;
}

int ComputeExact::Init(std::vector<int> resarray, int nres, int * nres_k,
                        double * Z0, double *Z1, int k, int m, int total,
                        int * total_k, double *prob_k, double * odds, double * p1,
                        int * IsExact, double epsilon, bool IsSmallmemory){
    int i, idx, k1;
    int * array;
    double Q;

    SaveParam(Z0, Z1, k, m, total, total_k, prob_k, odds, p1, IsExact, epsilon, IsSmallmemory);

    idx = 0;
    for(i=0; i<nres; i++){
        array = (resarray.data() + idx);
        k1 = nres_k[i];
        idx += k1;
        Q = CalTestStat(k1, array, false);
        m_Q.push_back(Q);
    }

    return 1;
}

int ComputeExact::SaveParam(double * Z0, double *Z1, int k, int m, int total,
                             int * total_k, double *prob_k, double * odds, double * p1,
                             int * IsExact, double epsilon, bool IsSmallmemory){
    int i, j;
    m_idx = 0;

    m_k = k;
    m_m = m;
    m_total = total;
    m_IsSmallmemory = IsSmallmemory;
    m_epsilon = epsilon;

    // Init _k vectors
    m_pprod = 1;
    for(i=0; i<=k; i++){
        m_total_k.push_back(total_k[i]);
        m_prob_k.push_back(prob_k[i]);
        m_denomi.push_back(0);
        m_IsExact.push_back(IsExact[i]);

        if(i < k){
            m_p1.push_back(p1[i]);
            m_odds.push_back(odds[i]);
            m_pprod = m_pprod * odds[i];
            m_p1_inv.push_back(1 - p1[i]);
        }
    }

    m_Z0.resize(m_k * m_m);
    m_Z1.resize(m_k * m_m);
    m_teststat_Z0.resize(m_m);
    m_teststat_Z1.resize(m_m);

    memcpy(m_Z0.data(), Z0, sizeof(double) * m_k * m_m);
    memcpy(m_Z1.data(), Z1, sizeof(double) * m_k * m_m);

    /* generate prob matrix */
    for(i=0; i<m_k; i++){
        int idx = i * m_m;
        for(j=0; j<m_m; j++){
            m_teststat_Z0.at(j) += m_Z0.at(idx + j);
            m_teststat_Z1.at(j) += m_Z1.at(idx + j);
        }
        // debug
        m_pr1_debug.push_back(0);
    }

    if(!m_IsSmallmemory){
        m_fprob.resize(m_total);
        m_teststat.resize(m_total);
    } else {
        m_fprob.clear();
        m_teststat.clear();
    }

    m_teststat_one.resize(m_m);
    m_temp_x.resize(m_k);
    m_temp_x1.resize(m_k);

    return 1;
}

int ComputeExact::Run(int test_type){
    int i, j, idx, l;
    std::vector<int> array(m_k);
    SL_setseed(1);

    for(i=0; i<m_k+1; i++){
        if(m_IsExact.at(i) == 1){
            if(i <= m_k/2 + 1){
                SKAT_Exact_Recurse(i, array, 0, 0, m_k);
            } else {
                SKAT_Exact_Recurse_INV(m_k - i, array, 0, 0, m_k);
            }
        } else if(m_total_k.at(i) < MIN_SIM1 && test_type == 3){
            SKAT_Resampling_Random(i, array);
        } else {
            SKAT_Resampling(i, array);
        }
    }

    // Get probability
    idx = 0;
    double total_prob_sum = 0;
    for(i=0; i<m_k+1; i++){
        for(j=idx; j<idx + m_total_k.at(i); j++){
            m_fprob.at(j) = m_fprob.at(j) / m_denomi.at(i) * m_prob_k.at(i);
            total_prob_sum += m_fprob.at(j);
        }
        idx = idx + m_total_k.at(i);
    }
    idx = 0;
    for(i=0; i<m_k+1; i++){
        m_prob_k.at(i) = 0;
        for(j=idx; j<idx + m_total_k.at(i); j++){
            m_fprob.at(j) = m_fprob.at(j) / total_prob_sum;
            m_prob_k.at(i) += m_fprob.at(j);
        }
        idx = idx + m_total_k.at(i);
    }

    double temp1 = 0;
    for(l=0; l<(int)m_Q.size(); l++){
        double n_num = 0;
        double n_same = 0;

        for(i=0; i<m_total; i++){
            temp1 = m_Q.at(l) - m_teststat.at(i);
            if(fabs(temp1) <= m_epsilon){
                temp1 = 0;
            }
            if(temp1 <= 0){
                n_num += m_fprob.at(i);
                if(temp1 == 0){
                    n_same += m_fprob.at(i);
                }
            }
        }

        m_pval.push_back(n_num);
        m_pval_same.push_back(n_same);
    }

    m_LargestQ = m_teststat.at(0);
    m_minP = 0;
    for(i=0; i<m_total; i++){
        temp1 = m_LargestQ - m_teststat.at(i);
        if(fabs(temp1) <= m_epsilon){
            temp1 = 0;
        }
        if(temp1 < 0){
            m_LargestQ = m_teststat.at(i);
            m_minP = m_fprob.at(i);
        } else if(temp1 == 0){
            m_minP += m_fprob.at(i);
        }
    }

    return 1;
}

int ComputeExact::GetPvalues(double * pval, double * pval_same, double * prob_k, double * minP){
    int i;
    for(i=0; i<(int)m_pval.size(); i++){
        pval[i] = m_pval[i];
        pval_same[i] = m_pval_same[i];
    }
    for(i=0; i<m_k+1; i++){
        prob_k[i] = m_prob_k[i];
    }
    if(minP != NULL){
        *minP = m_minP;
    }
    return 1;
}

int ComputeExact::PrintPval(){
    for(int i=0; i<(int)m_pval.size(); i++){
        printf("[%e][%e]\n", m_pval[i], m_pval_same[i]);
    }
    printf("MinP: [%e]\n", m_minP);
    return 1;
}

// ============================================================
//  SKAT_Exact global interface (from Binary_global.cpp)
// ============================================================
void SKAT_Exact(std::vector<int> & resarray, int nres, int * nres_k,
                double * Z0, double *Z1, int k, int m, int total,
                int * total_k, double *prob_k, double * odds, double * p1,
                int * IsExact, double * pval, double *pval_same, double *minP,
                int test_type, double epsilon){

    class ComputeExact exact;

    exact.Init(resarray, nres, nres_k, Z0, Z1, k, m, total, total_k, prob_k, odds, p1, IsExact, epsilon);
    exact.Run(test_type);
    exact.GetPvalues(pval, pval_same, prob_k, minP);
}

// ============================================================
//  ER_binary_func functions (from ER_binary_func.cpp)
// ============================================================

void GetProb_new(int k, int ngroup, int ncase, int* group, double* weight, double* prob){
    HyperGeo geo;
    geo.Run(k, ngroup, ncase, group, weight);
    geo.Get_lprob(prob);
}

void SKATExactBin_ComputeProb_Group(arma::uvec & idx, arma::uvec & idxCompVec,
                                    arma::vec & pi1, uint32_t n, uint32_t ncase,
                                    int type_group, std::vector<double> & prob){
    uint32_t k = idx.n_elem;
    int ngroup1 = 10;  // use the default value as ER will only be used for variants with MAC <= 10
    arma::vec p1 = pi1(idx);
    arma::vec p2 = pi1(idxCompVec);

    arma::uvec id_temp = arma::find(p1 >= 1);

    if(id_temp.n_elem > 0){
        for(unsigned int j = 0; j < id_temp.n_elem; j++){
            unsigned int id_temp_j = id_temp(j);
            p1(id_temp_j) = 0.999;
        }
    }

    std::vector<double> weight;
    std::vector<int> group;
    arma::uvec a1Vec, a2Vec, IDX;
    double a1, a2, p1temp, p2temp, oddtemp, p2oddtemp;
    arma::vec p1tempVec;

    for(int i = 0; i < ngroup1; i++){
        a1 = double(i) / ngroup1;
        a2 = double(i + 1) / ngroup1;
        if((i+1) < ngroup1){
            a1Vec = arma::find(p1 >= a1);
            a2Vec = arma::find(p1 < a2);
        } else {
            a1Vec = arma::find(p1 >= a1);
            a2Vec = arma::find(p1 <= a2);
        }
        IDX = arma::intersect(a1Vec, a2Vec);

        if(IDX.n_elem > 0){
            p1tempVec = p1(IDX);
            p1temp = arma::mean(p1tempVec);
            oddtemp = p1temp / (1 - p1temp);
            weight.push_back(oddtemp);
            group.push_back(IDX.n_elem);
        }
    }
    p2temp = arma::mean(p2);
    p2oddtemp = p2temp / (1 - p2temp);
    weight.push_back(p2oddtemp);

    for(size_t i = 0; i < weight.size(); i++){
        weight[i] = weight[i] / p2oddtemp;
    }
    group.push_back(n - k);

    int ngroup = group.size();
    int ncasei = int(ncase);
    GetProb_new(k, ngroup, ncasei, &group[0], &weight[0], &prob[0]);
}

int fact(int n) {
    if (n == 0 || n == 1){
        return 1;
    } else {
        return n * fact(n - 1);
    }
}

// P2 fix (2026-05-09): iterative C(n,r) to avoid 32-bit factorial overflow.
// Original implementation `fact(n)/(fact(r)*fact(n-r))` overflowed at n>=13
// (13! > INT_MAX), produced 0 in the denominator, and triggered SIGFPE
// (integer divide-by-zero) inside ER (efficient resampling) on rare binary
// markers with >=13 carriers. The iterative form keeps intermediates in
// long-long range while the final C(n,r) returns to int (typical values are
// small for SAIGE's ER use-case: k carriers, k ~ MAC, well below 100).
int n_choose_r(int n, int r){
    if (r < 0 || r > n) return 0;
    if (r == 0 || r == n) return 1;
    if (r > n - r) r = n - r;          // exploit symmetry C(n,r) = C(n, n-r)
    long long result = 1;
    for (int i = 0; i < r; ++i) {
        result *= static_cast<long long>(n - i);
        result /= static_cast<long long>(i + 1);
        if (result > static_cast<long long>(std::numeric_limits<int>::max())) {
            return std::numeric_limits<int>::max();   // saturate, don't overflow
        }
    }
    return static_cast<int>(result);
}

void Get_Total_K(int k, std::vector<int> & n_total_k){
    for(int i = 0; i <= k; i++){
        n_total_k[i] = n_choose_r(k, i);
    }
}

void SKATExactBin_ComputProb_New(arma::uvec & idx, arma::uvec & idxCompVec,
                                 arma::vec & pi1, uint32_t n, uint32_t ncase,
                                 int NResampling, int ExactMax, int test_type,
                                 int type_group,
                                 std::vector<double> & prob,
                                 std::vector<int> & IsExactVec,
                                 std::vector<int> & n_total_k,
                                 int & n_total, bool & Is_ExactP){

    uint32_t k = idx.n_elem;

    SKATExactBin_ComputeProb_Group(idx, idxCompVec, pi1, n, ncase, type_group, prob);

    Get_Total_K(k, n_total_k);
    // P2 fix (2026-05-09): use long-long accumulator. With n_choose_r now
    // returning correct (non-overflowed) C(k,i), the sum of (k+1) entries can
    // exceed INT_MAX for moderate k (e.g. k=20, sum=2^20=1M is fine, but
    // k=31 → 2^31 already overflows). Without this fix, an overflowed n_total
    // wraps negative, the cap below is skipped, and SaveParam tries to
    // resize() a vector to a huge unsigned value → throws std::length_error.
    {
        long long ll_total = 0;
        for (size_t i = 0; i < n_total_k.size(); ++i) ll_total += static_cast<long long>(n_total_k[i]);
        if (ll_total > static_cast<long long>(NResampling) || ll_total > std::numeric_limits<int>::max()) {
            n_total = NResampling + 1;  // sentinel forces the cap branch below
        } else {
            n_total = static_cast<int>(ll_total);
        }
    }
    Is_ExactP = true;
    if(n_total > NResampling){
        for(int i = 0; i <= (int)k; i++){
            if(n_total_k[i] > ExactMax){
                n_total_k[i] = int(ceil(NResampling * prob[i]));
                IsExactVec[i] = 0;
            }
        }
        Is_ExactP = false;
    }
    n_total = std::accumulate(n_total_k.begin(), n_total_k.end(), 0);
}

void Get_Res_Arrays(arma::mat & res_out, arma::uvec & idx,
                    std::vector<int> & resarray, int & nres,
                    std::vector<int> & nres_k){
    arma::vec res_out_colvec;
    for(int i = 0; i < nres; i++){
        res_out_colvec = res_out.col(i);
        arma::uvec res_out_i = arma::find(res_out_colvec > 0);
        arma::uvec res_out_i_s = arma::sort(res_out_i);
        int res_out_i_s_k = res_out_i_s.n_elem;
        nres_k[i] = res_out_i_s.n_elem;
        for(int k = 0; k < res_out_i_s_k; k++){
            resarray.push_back(res_out_i_s(k));
        }
    }
}

std::vector<std::vector<double>> mat_to_std_vec(arma::mat &A){
    std::vector<std::vector<double>> V(A.n_rows);
    typedef std::vector<double> stdvec;
    for (size_t i = 0; i < A.n_rows; ++i) {
        V[i] = arma::conv_to<stdvec>::from(A.row(i));
    }
    return V;
}

double SKATExactBin_Work(arma::mat & Z, arma::vec & res, arma::vec & pi1,
                         uint32_t ncase, arma::uvec & idx, arma::uvec & idxCompVec,
                         arma::mat & res_out,
                         int NResampling, int ExactMax, double epsilon, int test_type){

    uint32_t n = res.n_elem;
    arma::vec p1 = pi1(idx);
    arma::vec p2 = pi1(idxCompVec);

    arma::mat Z_1 = Z.rows(idx);
    arma::mat Z1temp = (Z_1 % (-p1)).t();
    arma::vec Z0 = arma::vectorise(Z1temp);

    arma::mat Z1temp2 = (Z_1 % (1 - p1)).t();
    arma::vec Z1 = arma::vectorise(Z1temp2);

    uint32_t m = Z_1.n_cols;
    int k = idx.n_elem;
    std::vector<int> n_total_k(k+1, 0);

    std::vector<double> prob(k+1, 0.0);
    std::vector<int> IsExactVec(k+1, 1);

    int n_total = 0;
    bool Is_ExactP = false;
    SKATExactBin_ComputProb_New(idx, idxCompVec, pi1, n, ncase, NResampling, ExactMax,
                                test_type, 2, prob, IsExactVec, n_total_k, n_total, Is_ExactP);

    double p1mean = arma::mean(p1);
    arma::vec p1_adj = p1 / p1mean;
    arma::vec odds = p1 / (1 - p1);
    int test_type_new = 1;

    if(res_out.is_empty()){
        res_out = res(idx);
    } else {
        arma::vec res_out2 = arma::join_cols(res(idx), res_out(idx));
        res_out.resize(res_out2.n_elem);
        res_out = res_out2;
    }

    std::vector<int> resarray;
    int nres = res_out.n_cols;
    std::vector<int> nres_k(nres, 0);

    Get_Res_Arrays(res_out, idx, resarray, nres, nres_k);
    std::vector<double> pval(nres, 0.0);
    std::vector<double> pval1(nres, 0.0);
    double minP = 100;

    typedef std::vector<double> stdvec;
    stdvec Z1std = arma::conv_to<stdvec>::from(Z1);
    stdvec Z0std = arma::conv_to<stdvec>::from(Z0);
    stdvec oddsstd = arma::conv_to<stdvec>::from(odds);
    stdvec p1_adjstd = arma::conv_to<stdvec>::from(p1_adj);

    SKAT_Exact(resarray, nres, &nres_k[0], &Z0std[0], &Z1std[0], k, m, n_total,
               &n_total_k[0], &prob[0], &oddsstd[0], &p1_adjstd[0], &IsExactVec[0],
               &pval[0], &pval1[0], &minP, test_type_new, epsilon);

    double pvalue = pval[0] - pval1[0] / 2;
    return(pvalue);
}

} // namespace ER
