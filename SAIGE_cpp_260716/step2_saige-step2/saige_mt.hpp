// saige_mt.hpp — multi-trait (P phenotypes over one genotype stream) step 2.
//
// See MULTITRAIT_DESIGN.md. Phase 0 scope: the data structures and the config
// front end only. Nothing here is constructed when P == 1; the single-trait
// path in main() is untouched and must stay byte-identical (design section 5).
//
// Layout (design section 1.1) is mixed SoA:
//   - small per-trait quantities (scalars, p-vectors, p x p matrices) -> AoS,
//     indexed by internal trait index (TraitMeta / std::vector<arma::mat>);
//   - large per-trait quantities (N-vectors, N x p matrices) -> SoA, stacked
//     side by side into wide matrices so one GEMM streams Gb exactly once;
//   - P untouched SAIGEClass instances for the scalar fallback path, so a
//     fallback pair runs literally the same code as today's single-trait run.

#ifndef SAIGE_MT_HPP
#define SAIGE_MT_HPP

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <armadillo>
#include <yaml-cpp/yaml.h>

#include "null_model_loader.hpp"

namespace SAIGE {

enum class TraitKind { Binary, Quantitative, Survival };

// Throws on anything other than binary / quantitative / survival. Used when a
// multi-trait run classifies its models -- an unrecognised type there would
// silently pick the wrong batch kernel, so it must be fatal.
TraitKind   traitKindFromString(const std::string& t_traitType);
// Non-throwing form. The single-trait path has never validated traitType (an
// unknown string just falls through to the quantitative branches), so the P==1
// TraitMeta is built with this and keeps that behaviour exactly.
bool        tryTraitKindFromString(const std::string& t_traitType, TraitKind& t_out);
const char* traitKindName(TraitKind t_kind);

// ------------------------------------------------------------------
// Config front end
// ------------------------------------------------------------------

// YAML keys that may be written either at the top level (applies to every
// trait) or inside one `models:` entry (applies to that trait and beats the
// top-level value). Design section 4.1. Each has an explicit "was it written"
// flag because "absent" and "false"/"0.0" are different: absent means keep
// whatever the trait's own nullmodel.json says.
struct ModelOverrides {
    bool   has_is_Firth_beta    = false;
    bool   is_Firth_beta        = false;
    bool   has_pCutoffforFirth  = false;
    double pCutoffforFirth      = 0.0;
    bool   has_isnoadjCov       = false;
    bool   isnoadjCov           = false;
    bool   has_cateVarRatioMinMACVecExclude = false;
    std::vector<double> cateVarRatioMinMACVecExclude;
    bool   has_cateVarRatioMaxMACVecInclude = false;
    std::vector<double> cateVarRatioMaxMACVecInclude;

    bool any() const {
        return has_is_Firth_beta || has_pCutoffforFirth || has_isnoadjCov ||
               has_cateVarRatioMinMACVecExclude || has_cateVarRatioMaxMACVecInclude;
    }
};

// One entry of `models:` — or the single implied entry built from the legacy
// scalar modelFile / varianceRatioFile / outputFile trio.
struct MTModelSpec {
    std::string traitName;          // label for logs/errors; defaults to the model dir name
    std::string modelFile;          // step 1 output dir (*.arma + nullmodel.json)
    std::string varianceRatioFile;
    std::string outputFile;
    ModelOverrides ov;
};

// Reads the `models:` sequence, or falls back to the scalar trio. Throws with a
// message naming the offending model/key on: both forms present, `models:` not
// a sequence, empty `models:`, a model entry missing one of the three paths, or
// two models writing to the same outputFile. Order of the returned vector is
// the order written in the config == the order of the output files.
std::vector<MTModelSpec> parseModelSpecs(const YAML::Node& t_config);

// Reads the override keys out of a node (top level or one model entry).
ModelOverrides parseModelOverrides(const YAML::Node& t_node);

// ------------------------------------------------------------------
// Per-trait constants — built once, read-only inside the marker loop
// ------------------------------------------------------------------

struct TraitMeta {
    std::string name;
    std::string modelDir, vrFile, outFile;
    // Verbatim string out of nullmodel.json. The output writers branch on this
    // and not on `kind`, so their behaviour is unchanged for any value the
    // loader accepts, recognised or not.
    std::string traitType;
    TraitKind   kind = TraitKind::Quantitative;
    int    p      = 0;      // covariate count including the intercept
    int    colOff = 0;      // first column of this trait inside Xstack / Astack
    int    binOff = -1;     // first column inside WXstack; -1 when not binary
    int    binIdx = -1;     // trait index inside MU2bin; -1 when not binary
    double tau0            = 1.0;
    double SPA_Cutoff      = 2.0;
    bool   is_Firth_beta   = false;
    double pCutoffforFirth = 0.0;
    bool   isFastTest      = false;
    double pval_cutoff_for_fastTest = 0.0;
    bool   isnoadjCov      = false;
    bool   flagSparseGRM   = false;
    bool   isCondition     = false;
    bool   isMoreOutput    = false;
    bool   locoApplied     = false;   // did this trait really read chr<N>/ ?
    bool   batchable = false;         // static gating result, design section 3.1
    int    outIdx = 0;                // position in the config's `models:` order
    // Different sample sets (design section 4.7). Column of this trait inside
    // MTContext::MASKq; -1 unless the trait is quantitative AND its sample
    // list differs from the union's.
    int    maskIdx = -1;
};

// Where one trait's samples sit inside the union sample list (design 4.7).
// Built once; read-only in the marker loop.
struct MTTraitSamples {
    // The trait's sampleIDs are the union's, element for element. Such a
    // trait's genotype vector IS the block's union column, so it takes exactly
    // the arithmetic of a same-sample-set run and none of the corrections.
    bool sameAsUnion = true;
    int  n = 0;                            // this trait's sample count
    std::vector<arma::uword> pos;          // [n]  union index of the trait's k-th sample
    std::vector<arma::uword> comp;         // union indices NOT in the trait, ascending
};

// Built once after all P null models are loaded; read-only from then on.
// NEVER constructed when P == 1.
struct MTContext {
    int N = 0, P = 0, nBin = 0;
    int sumP = 0, sumPbin = 0, sumPqnt = 0, qOff = 0;   // qOff == sumPbin

    // Design section 4.6 invariant I5: one run == one chromosome. The per-trait
    // constants below are only constant because the chromosome never changes
    // mid-run; assert against this if that ever stops being true.
    std::string locoChrom;
    bool        locoEnabled = false;

    arma::mat Xstack;    // N x sumP     internal order (binary traits first)
    arma::mat Astack;    // N x sumP     block t = XVX_inv_XV of trait t
    arma::mat WXstack;   // N x sumPbin  block t = mu2_t % X_t (binary only)
    arma::mat RES;       // N x P        column t = res_t
    arma::mat MU2bin;    // N x nBin
    // No CCM (the case/control indicator stack of design 2.3): the case and
    // control allele counts are still summed per pair in index order, because
    // a GEMM re-associates that sum and would move the printed AF columns in
    // their last digit on mean-imputed dosages. See the commit that added the
    // kernel.

    std::vector<arma::mat> XVX;    // P matrices, p_t x p_t
    std::vector<arma::vec> S_a;    // P vectors, p_t
    std::vector<TraitMeta> meta;   // internal order; TraitMeta::outIdx carries
                                   // the position in the config's models: list

    // ---- different sample sets (design section 4.7) ----
    // N above is the size of the union. Every stacked matrix is union-length:
    // a trait's rows sit at its union positions and the rows of samples the
    // trait does not have are exact zeros, so any inner product of a stack
    // column with a union-length genotype column only ever sees that trait's
    // own samples.
    bool sampleSetsDiffer = false;          // at least one trait is !sameAsUnion
    std::vector<std::string> unionIDs;      // [N], R's union_vector order
    std::vector<MTTraitSamples> samp;       // [P], internal order
    // Per-trait sums over the trait's samples (union-length columns summed in
    // union order), the constant half of the flip correction in scoreTestBatchMT.
    std::vector<arma::vec> sumA;            // p_t   column sums of A_t
    std::vector<arma::vec> sumW;            // p_t   of mu2_t % X_t (binary) / X_t (quantitative)
    std::vector<double>    sumR;            // sum res_t
    std::vector<double>    sumM;            // sum mu2_t (binary) / n_t (quantitative)
    arma::mat MASKq;    // N x nMask  0/1 sample indicator of each quantitative
                        // trait whose sample list differs from the union

    std::vector<int> batchTraits;       // batchable, binary first
    std::vector<int> batchQuantTraits;  // batchable and quantitative
    std::vector<int> scalarTraits;      // not batchable

    // ---- folded covariate projection for quantitative traits ----
    // Config key mtFoldQuantProj. For a quantitative trait V_t = (1/tau0_t)*I,
    // so XVX_inv_XV_t = V_t X (X'V_t X)^{-1} lies in the column space of X --
    // and every trait in a run shares the same covariate matrix. Both wide
    // (sum_t p_t)-column GEMMs of scoreTestBatchMT therefore collapse into one
    // p-column GEMM Z0 = Xref' G, with
    //     A_t' G = K_t' Z0      K_t the p x p matrix with A_t = Xref K_t
    //     X_t' G = Z0           (X_t is Xref bit for bit)
    // K_t is FITTED to the stored A_t (least squares against Xref) rather than
    // rebuilt from tau0 and XVX_inv: step 1 computes XVX_inv_XV in fp32 and
    // nullmodel.json rounds tau0 to 6 significant digits, so the analytic form
    // is 5-10x further from the number the wide path actually contracts.
    // foldResid records how far the fit is from the stored matrix; a trait
    // whose residual is above mtFoldResidTol keeps the wide path.
    bool foldQuant = false;             // at least one trait takes the fold
    int  foldRef   = -1;                // internal index of the trait whose X block is Xref
    int  foldRefCol0 = 0;               // Xref's first column inside Xstack
    int  foldRefP    = 0;               // Xref's column count
    std::vector<char>      foldable;    // [P] 1 = this trait takes the fold
    std::vector<arma::mat> foldK;       // [P] p x p, only for foldable traits
    std::vector<double>    foldResid;   // [P] max|A_t - Xref K_t| / max|A_t|, -1 when not tried
    std::vector<int>       foldTraits;  // internal indices that took the fold

    // ---- block-at-a-time (marker, trait) statistics ----
    // Config key mtVecQuantStats. When set, a quantitative trait's block tail
    // runs emitBlockResultsVecQuant instead of the per-pair
    // format_score_result: one branch-free pass over the block's B columns for
    // Beta / seBeta / Tstat / var / StdStat, one vectorised erfc for the
    // chi-square(1) upper tail, and one std::to_chars pass for the "%.6E"
    // string. Pairs whose p-value would fall below MT_VEC_EXACT_BELOW_P, and
    // every degenerate pair, still go through format_score_result itself, so
    // the tail is bit-identical. Binary traits are never touched: their tail
    // gates SPA / ER / Firth on each pair's own statistic. See score_vec.hpp.
    bool vecQuantStats = false;
};

// Largest relative fit residual a trait may have and still take the fold.
// The floor is step 1's fp32 rounding of XVX_inv_XV, ~6e-8 relative; anything
// that is not a scalar-V model misses the column space of X outright and
// lands near 1.0, so this threshold separates the two cases by four orders.
constexpr double MT_FOLD_RESID_TOL = 1e-5;

// The batch kernel's intermediates: one per thread, reused across blocks
// (grow-only, never reallocated per block). The block's genotype matrix itself
// is NOT here -- the caller owns it, because it is filled marker by marker
// during the read and only then handed to the kernel.
struct MTScratch {
    arma::mat Gb2;     // N x B   = Gb % Gb
    arma::mat Zall;    // sumP x B
    arma::mat GWbin;   // sumPbin x B
    arma::mat GWqnt;   // sumPqnt x B
    arma::mat GR;      // B x P
    arma::mat G2Mu2;   // B x nBin
    arma::vec Gsq;     // B
    // Different sample sets only (design 4.7).
    arma::mat GMu2;    // B x nBin   Gb' MU2bin      (flip correction, binary)
    arma::mat GMask1;  // B x nMask  Gb' MASKq       (flip correction, quantitative)
    arma::mat GMask2;  // B x nMask  Gb2' MASKq      (sum of g^2 over the trait's samples)
    arma::mat MissA, MissWbin, MissWqnt;   // (cols) x B  stack rows summed over a column's missing cells
    arma::mat MissR, MissMu2, MissMask;    // B x (P | nBin | nMask), same
    arma::mat Zc, Wc;                      // p x B   one trait's corrected Z / GW
    arma::vec Rc, Qc;                      // B
    // Folded quantitative projection only (MTContext::foldQuant).
    arma::mat Z0;      // p x B   Xref' G, shared by every folded trait
    arma::mat Zf;      // p x B   one folded trait's K_t' Z0
    // Block tail only (MTContext::vecQuantStats). Length B, grow-only.
    std::vector<double> evVar1, evStat, evZ, evP;
};

// Per-(block column, trait) description of how the trait's own genotype vector
// relates to the block's union column (design section 4.7). With g the union
// column and g_t the vector a single-trait run of trait t would build, on the
// trait's samples
//     g_t = a * g + b + d * [cell is a missing genotype]
// exactly: a = +1, b = 0 when the trait and the union agree on the flip,
// a = -1, b = 2 when they do not, and d is the difference of the two imputed
// values. The kernel only reads the traits whose sample list differs from the
// union's; everyone else takes the unadjusted arithmetic.
struct MTBlockAdj {
    arma::mat a, b, d;    // B x P
    arma::mat q;          // B x P  g_t^2 - (a*g + b)^2 on a missing cell
    arma::umat nMiss;     // B x P  missing genotypes among the trait's samples
    std::vector<std::vector<arma::uword>> miss;   // [B] union indices of the column's missing cells

    void resize(int t_B, int t_P);
};

// ------------------------------------------------------------------
// Multi-trait model set: ordering, validation, static gating
// ------------------------------------------------------------------

// Internal trait order (design section 1.2): binary first, then quantitative,
// then everything else; stable within each group so the config order decides
// ties. Returns a permutation of [0, P) holding config indices. For P == 1
// this is always {0}, so the single-trait path sees no reordering at all.
std::vector<int> mtInternalOrder(const std::vector<NullModelData>& t_nms);

// Hard-fails when the P models cannot legitimately share one genotype stream
// (design sections 4.3 / 4.7): different impute_method; a model whose
// sampleIDs list is empty, does not have n entries, or repeats an ID while the
// lists differ between models; or -- only when t_requireSameSamples -- sample
// lists that differ in content or order. Warns -- does not fail -- when the
// models disagree on whether LOCO really applied. `t_names` supplies the trait
// label used in the messages and must be the same length as t_nms.
// Returns true when the sample lists are not all identical (content and order).
bool validateMTModels(const std::vector<NullModelData>& t_nms,
                      const std::vector<std::string>& t_names,
                      bool t_requireSameSamples);

// R's union_vector (readInGLMM.R ReadModel_multiTrait): every model's
// sampleIDs, first occurrence wins, models taken in config order. When all
// lists are identical this is model 0's list, element for element.
std::vector<std::string> mtUnionSampleIDs(const std::vector<NullModelData>& t_nms);

// Static per-trait gate, evaluated once before the marker loop (design
// section 3.1). False means the trait runs the per-pair scalar path for every
// marker; it never changes per marker.
bool isBatchable(const TraitMeta& t_meta);

// Why isBatchable() said no; "-" when it said yes. For the startup table.
const char* batchableReason(const TraitMeta& t_meta);

// Design section 4.5: the table has to be printed, or a numerical problem in a
// P = 64 run cannot be localised to a trait.
void printMTGateTable(const std::vector<TraitMeta>& t_meta, bool t_locoEnabled);

// ------------------------------------------------------------------
// Batch kernel
// ------------------------------------------------------------------

// Fills t_meta's colOff / binOff / binIdx and builds the stacked matrices.
// t_meta must already be in internal order; t_order[t] is the config index of
// internal trait t, so t_nms[t_order[t]] is that trait's loaded model.
// Every trait gets a column block, batchable or not: the few extra columns cost
// p_t*N doubles and keep the offsets a single uniform rule.
// t_unionIDs is mtUnionSampleIDs(t_nms); the stacks are built at its length
// with each trait's rows placed at its union positions (design 4.7).
void buildMTContext(MTContext& t_ctx,
                    const std::vector<NullModelData>& t_nms,
                    const std::vector<int>& t_order,
                    std::vector<TraitMeta>& t_meta,
                    bool t_locoEnabled,
                    const std::string& t_locoChrom,
                    const std::vector<std::string>& t_unionIDs,
                    bool t_foldQuantProj = false);

// One marker block's normal-approximation results, B x P. Only the columns of
// the trait set passed to scoreTestBatchMT are written.
struct MTBlockResult {
    arma::mat Beta, seBeta, Tstat, var1, var2, StdStat;
    // The p-value as format_score_result returns it: linear, or natural log
    // when pvalIsLog is set for that pair. This is the number the scalar path
    // gates Firth on, so the fallback decision uses it and not a re-parse of
    // the printed string.
    arma::mat pvalRaw;
    std::vector<std::vector<std::string>> pvalStr;    // [P][B]
    std::vector<std::vector<char>>        pvalIsLog;  // [P][B]

    void resize(int t_B, int t_P);
};

// Normal-approximation score test for one marker block against several traits
// at once (design section 2). Computes exactly the quantities
// SAIGEClass::scoreTestFast computes -- S, var2, var1, and the p-value through
// the same format_score_result -- but with the per-trait cancellations of
// design section 2.2 applied, so neither B nor gtilde is ever materialised.
//
//   t_Gb        N x B, imputed / QC'd / flipped; one column per marker
//   t_j0, t_j1  the half-open column range to score. The caller packs the
//               high-MAC markers from column 0 up and the low-MAC ones from
//               column B down, so both groups are contiguous and each gets one
//               call with its own trait set (design section 3.2) without ever
//               copying or masking the block.
//   t_traitSet  internal trait indices to score; other columns are untouched
//   t_VR        B x P, the per-pair variance ratio (only t_traitSet read)
//
// NOT bit-identical to scoreTestFast: the scalar version sums over the carrier
// samples only, this one sums over all N (the non-carriers contribute exact
// zeros but change the association order). Algebraically equal; the caller
// must treat the difference as last-bit rounding, never as a licence to skip a
// fallback.
//
//   t_adj       per-(column, trait) map from the union column to the trait's
//               own genotype vector; read only for traits whose sample list
//               differs from the union's (ctx.samp[t].sameAsUnion == false),
//               may be nullptr when there are none.
void scoreTestBatchMT(const MTContext& t_ctx,
                      const std::vector<int>& t_traitSet,
                      const arma::mat& t_Gb,
                      int t_j0, int t_j1,
                      const arma::mat& t_VR,
                      const MTBlockAdj* t_adj,
                      MTScratch& t_scr,
                      MTBlockResult& t_out);

// Same kernel, but the sample-space reductions in t_scr were produced
// elsewhere -- by the GPU path in gpu/gpu_step2.hpp -- so only the O(p^2) /
// O(P) tail runs here. Restricted to quantitative, batchable traits whose
// sample list is the union's, because that is the case in which the four
// prefilled quantities are the whole of scoreTestBatchMT's sample-space work.
// See the definition in saige_mt.cpp for the exact contract on t_scr; it
// throws rather than silently reading the wrong rows if the trait set does not
// match.
void scoreTestBatchMTQuantPre(const MTContext& t_ctx,
                              const std::vector<int>& t_traitSet,
                              int t_j0, int t_j1,
                              const arma::mat& t_VR,
                              MTScratch& t_scr,
                              MTBlockResult& t_out);


#ifdef MTFOLD_PROF
// Defined in saige_mt.cpp; see the MTF_TIC/MTF_TOC block there. cpu-seconds
// summed across threads, not wall clock.
extern double g_mtfProfZall, g_mtfProfGW, g_mtfProfZ0, g_mtfProfGR;
#endif

#ifdef MTVEC_PROF
// Defined in saige_mt.cpp. cpu-seconds summed across threads.
//   g_mtvProfEmit   the whole per-(block, trait) tail, either path
//   g_mtvProfStat   block path: variance ratio, Beta/seBeta/Tstat/var/StdStat,
//                   the vectorised erfc, and the stores
//   g_mtvProfFmt    block path: the "%.6E" pass (plus the rare boost fallback)
//   g_mtvProfPairs  pairs that took the block path
//   g_mtvProfFall   pairs handed back to format_score_result
// Set MTVEC_PROF_NOFMT=1 in the environment to skip the string pass entirely;
// the output is then garbage, but the difference is the formatting cost.
extern double g_mtvProfEmit, g_mtvProfStat, g_mtvProfFmt;
extern long   g_mtvProfPairs, g_mtvProfFall;
extern bool   g_mtvProfNoFmt;
#endif

}  // namespace SAIGE

#endif  // SAIGE_MT_HPP
