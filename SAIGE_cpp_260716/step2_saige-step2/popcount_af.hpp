// popcount_af.hpp -- case / control allele sums from the packed 2-bit column.
//
// The CPU multi-trait loop (main.cpp, mainMarkerMT) reports AF_case / AF_ctrl
// for every (marker, binary trait) pair by summing the marker's dosage vector
// over the trait's case indices and again over its control indices: two
// index-driven gathers of O(N) each, per pair. With the PLINK fused decode the
// marker's 2-bit codes are still in hand when the block is finalised, and the
// dosage of a sample is a pure function of its code (the 4-entry table
// FusedMarkerStats::fd), so the same sum is
//
//     sum_case g  =  sum_c  n_case[c] * fd[c]
//
// with n_case[c] the number of cases carrying code c -- three popcounts per
// 64-bit word of (column AND case-mask), where the mask has 11 in every case
// sample's field and 00 elsewhere.
//
// Exactness. The gather adds doubles one sample at a time in index order, so
// the question is whether the counts formula reproduces that floating-point
// sum bit for bit, not whether it equals it as real numbers.
//   * Every non-missing dosage is 0, 1 or 2. When the missing dosage
//     fd[MISSING] is an integer too (best_guess / minor imputation, a cleaned
//     value, or no missing call among the trait's samples), every partial sum
//     of the gather is an integer below 2^53 and no addition rounds; the
//     counts formula is a sum of exact integer products and rounds nowhere
//     either. Identical.
//   * With mean imputation fd[MISSING] = 2 * AF is not an integer, and once
//     it has been added the running sum carries a fraction; from then on an
//     integer addition rounds exactly when it carries the sum past the next
//     power of two (the grid doubles and the fraction is rounded to it, ties
//     to even -- so WHICH sample crosses matters, because the tie-break reads
//     the parity of the partial sum). pcReplaySum walks the column word by
//     word and adds a whole word's integer dosage in one step whenever that
//     step provably cannot round (the sum is still an integer, or stays
//     below the next power of two); a word that holds a missing case field
//     or would cross a power of two is replayed field by field, which IS the
//     gather over those 32 samples. So it performs the same rounding
//     operations on the same operands in the same order. Identical, given
//     that the gather visits samples in ascending index order -- true for
//     arma::find's case_indices, and for a trait with its own sample list
//     only when its union positions ascend, which the caller checks once.
//
// The tests in tools/popcount_af_test.cpp compare both paths against the
// literal gather on random columns, masks and tables, including ties.
#ifndef SAIGE_POPCOUNT_AF_HPP
#define SAIGE_POPCOUNT_AF_HPP

#include <cmath>
#include <cstdint>
#include <cstring>

namespace SAIGE {

// PLINK .bed codes, as PlinkClass has them.
constexpr unsigned PC_HOM_ALT = 0x0, PC_MISSING = 0x1, PC_HET = 0x2, PC_HOM_REF = 0x3;
constexpr uint64_t PC_M55 = 0x5555555555555555ULL;

// Words needed for n samples at 2 bits each.
inline int pcWords(int t_n) { return (t_n + 31) / 32; }

// Sets the 2-bit field of every listed sample to 11. t_mask must hold
// t_nWords zeroed words.
template <class Idx>
inline void pcBuildMask(const Idx* t_pos, std::size_t t_n, uint64_t* t_mask)
{
    for (std::size_t k = 0; k < t_n; ++k) {
        const uint64_t u = static_cast<uint64_t>(t_pos[k]);
        t_mask[u >> 5] |= 3ULL << (2 * (u & 31));
    }
}

// Whether t_pos[k] is strictly increasing in k.
template <class Idx>
inline bool pcAscending(const Idx* t_pos, std::size_t t_n)
{
    for (std::size_t k = 1; k < t_n; ++k)
        if (!(t_pos[k] > t_pos[k - 1])) return false;
    return true;
}

// n[c] = number of masked samples whose code is c, c = 0..3. t_nMasked is the
// number of samples in the mask (the fourth count is the remainder).
inline void pcCountCodes(const uint64_t* t_col, const uint64_t* t_mask, int t_nWords,
                         uint64_t t_nMasked, uint64_t t_n[4])
{
    uint64_t n11 = 0, n10 = 0, n01 = 0;
    for (int w = 0; w < t_nWords; ++w) {
        const uint64_t x  = t_col[w] & t_mask[w];
        const uint64_t lo = x & PC_M55;
        const uint64_t hi = (x >> 1) & PC_M55;
        n11 += (uint64_t)__builtin_popcountll(hi & lo);
        n10 += (uint64_t)__builtin_popcountll(hi & ~lo);
        n01 += (uint64_t)__builtin_popcountll(lo & ~hi);
    }
    t_n[PC_HOM_REF] = n11;
    t_n[PC_HET]     = n10;
    t_n[PC_MISSING] = n01;
    t_n[PC_HOM_ALT] = t_nMasked - n11 - n10 - n01;
}

inline bool pcIsInteger(double t_x) { return std::isfinite(t_x) && t_x == std::floor(t_x); }

// The counts formula. Exact (see the header comment) when fd[MISSING] is an
// integer or n[MISSING] == 0; the caller decides that.
inline double pcCountsSum(const uint64_t t_n[4], const double* t_fd)
{
    return (double)t_n[PC_HOM_REF] * t_fd[PC_HOM_REF]
         + (double)t_n[PC_HET]     * t_fd[PC_HET]
         + (double)t_n[PC_MISSING] * t_fd[PC_MISSING]
         + (double)t_n[PC_HOM_ALT] * t_fd[PC_HOM_ALT];
}

// The smallest power of two strictly above t_s (t_s > 0).
inline double pcNextPow2Above(double t_s)
{
    int e;
    (void)std::frexp(t_s, &e);      // t_s = m * 2^e, m in [0.5, 1)
    return std::ldexp(1.0, e);
}

// Exact replay of
//     s = 0; for u ascending: if mask[u]: s += fd[code[u]]
// in double arithmetic. See the header comment for why the word-level
// shortcuts cannot change a rounding.
inline double pcReplaySum(const uint64_t* t_col, const uint64_t* t_mask, int t_nWords,
                          const double* t_fd)
{
    double s = 0.0;
    for (int w = 0; w < t_nWords; ++w) {
        const uint64_t mk = t_mask[w];
        if (mk == 0) continue;
        const uint64_t x   = t_col[w] & mk;
        const uint64_t lo  = x & PC_M55;
        const uint64_t hi  = (x >> 1) & PC_M55;
        const uint64_t mis = lo & ~hi;
        if (mis == 0) {
            // Every dosage in this word is 0, 1 or 2: an integer step K.
            const uint64_t n11 = (uint64_t)__builtin_popcountll(hi & lo);
            const uint64_t n10 = (uint64_t)__builtin_popcountll(hi & ~lo);
            const uint64_t nmk = (uint64_t)__builtin_popcountll(mk & PC_M55);
            const uint64_t n00 = nmk - n11 - n10;
            const double K = (double)n11 * t_fd[PC_HOM_REF]
                           + (double)n10 * t_fd[PC_HET]
                           + (double)n00 * t_fd[PC_HOM_ALT];
            if (s == std::floor(s)) {
                // integer + integer below 2^53: no step rounds
                s += K;
                continue;
            }
            // s carries a fraction: every partial sum s + P_j (P_j integer)
            // lies on s's grid and is representable as long as it stays
            // below the next power of two, so no step rounds if the whole
            // word's step does not reach it.
            if (s + K < pcNextPow2Above(s)) {
                s += K;
                continue;
            }
        }
        // A missing case field, or a power-of-two crossing somewhere in this
        // word: do what the gather does, sample by sample, in order.
        for (int j = 0; j < 32; ++j) {
            if ((mk >> (2 * j)) & 1ULL) s += t_fd[(x >> (2 * j)) & 3ULL];
        }
    }
    return s;
}

// hom / het classification of a dosage, as main.cpp's finalize does it.
inline bool pcIsHom(double t_d) { return t_d >= 1.5 && t_d <= 2.0; }
inline bool pcIsHet(double t_d) { return t_d >= 0.5 && t_d < 1.5; }

inline int pcSumPairFromCounts(const uint64_t nc[4], const uint64_t no[4],
                               const uint64_t* t_col,
                               const uint64_t* t_caseMask, const uint64_t* t_ctrlMask, int t_nWords,
                               const double* t_fd, bool t_ascending,
                               double& t_sumCase, double& t_sumCtrl,
                               bool t_moreOutput,
                               uint32_t& t_caseHom, uint32_t& t_caseHet,
                               uint32_t& t_ctrlHom, uint32_t& t_ctrlHet);

// The whole per-(marker, trait) computation. Returns
//   0  not done: a mean-imputed missing call sits in the case or control set
//      and the trait's sample positions do not ascend, so the sequential sum
//      cannot be replayed; the caller must gather
//   1  both sums by the counts formula
//   2  at least one sum by the exact replay
// t_nCase / t_nCtrl are the mask populations. The four hom / het counts are
// written only when t_moreOutput.
// t_total (config key mtPopcountCtrlFromTotal; nullptr otherwise) is the
// column's code count over ALL of the trait's samples. When the case and
// control masks partition those samples the control counts are the total
// minus the case counts -- integers, no rounding -- and the second popcount
// pass is skipped. The caller asserts the partition (N_case + N_ctrl == n and
// the trait's sample list is the column's); the replay path, which needs the
// control mask itself, is unchanged.
inline int pcSumPair(const uint64_t* t_col,
                     const uint64_t* t_caseMask, const uint64_t* t_ctrlMask, int t_nWords,
                     uint64_t t_nCase, uint64_t t_nCtrl, const double* t_fd, bool t_ascending,
                     double& t_sumCase, double& t_sumCtrl,
                     bool t_moreOutput,
                     uint32_t& t_caseHom, uint32_t& t_caseHet,
                     uint32_t& t_ctrlHom, uint32_t& t_ctrlHet,
                     const uint64_t* t_total = nullptr)
{
    uint64_t nc[4], no[4];
    pcCountCodes(t_col, t_caseMask, t_nWords, t_nCase, nc);
    if (t_total) {
        for (unsigned c = 0; c < 4; ++c) no[c] = t_total[c] - nc[c];
    } else {
        pcCountCodes(t_col, t_ctrlMask, t_nWords, t_nCtrl, no);
    }
    return pcSumPairFromCounts(nc, no, t_col, t_caseMask, t_ctrlMask, t_nWords, t_fd, t_ascending,
                               t_sumCase, t_sumCtrl, t_moreOutput,
                               t_caseHom, t_caseHet, t_ctrlHom, t_ctrlHet);
}

// The decision and the sums, given the case / control code counts from
// wherever they were computed (pcCountCodes here, count_codes on the device in
// gpu_step2.cu). Same return codes as pcSumPair.
inline int pcSumPairFromCounts(const uint64_t nc[4], const uint64_t no[4],
                               const uint64_t* t_col,
                               const uint64_t* t_caseMask, const uint64_t* t_ctrlMask, int t_nWords,
                               const double* t_fd, bool t_ascending,
                               double& t_sumCase, double& t_sumCtrl,
                               bool t_moreOutput,
                               uint32_t& t_caseHom, uint32_t& t_caseHet,
                               uint32_t& t_ctrlHom, uint32_t& t_ctrlHet)
{
    const bool missInt  = pcIsInteger(t_fd[PC_MISSING]);
    const bool caseInt  = missInt || nc[PC_MISSING] == 0;
    const bool ctrlInt  = missInt || no[PC_MISSING] == 0;
    if ((!caseInt || !ctrlInt) && !t_ascending) return 0;

    t_sumCase = caseInt ? pcCountsSum(nc, t_fd) : pcReplaySum(t_col, t_caseMask, t_nWords, t_fd);
    t_sumCtrl = ctrlInt ? pcCountsSum(no, t_fd) : pcReplaySum(t_col, t_ctrlMask, t_nWords, t_fd);

    if (t_moreOutput) {
        uint64_t ch = 0, ce = 0, oh = 0, oe = 0;
        for (unsigned c = 0; c < 4; ++c) {
            if (pcIsHom(t_fd[c]))      { ch += nc[c]; oh += no[c]; }
            else if (pcIsHet(t_fd[c])) { ce += nc[c]; oe += no[c]; }
        }
        t_caseHom = (uint32_t)ch; t_caseHet = (uint32_t)ce;
        t_ctrlHom = (uint32_t)oh; t_ctrlHet = (uint32_t)oe;
    }
    return (caseInt && ctrlInt) ? 1 : 2;
}

}  // namespace SAIGE

#endif  // SAIGE_POPCOUNT_AF_HPP
