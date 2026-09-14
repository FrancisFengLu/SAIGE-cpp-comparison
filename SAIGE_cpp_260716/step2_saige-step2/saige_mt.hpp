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
    int    binIdx = -1;     // trait index inside MU2bin / CCM; -1 when not binary
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
    int    nCase = 0, nCtrl = 0;
    bool   batchable = false;         // static gating result, design section 3.1
    int    outIdx = 0;                // position in the config's `models:` order
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
    arma::mat CCM;       // N x 2*nBin   [case_0, ctrl_0, case_1, ctrl_1, ...]

    std::vector<arma::mat> XVX;    // P matrices, p_t x p_t
    std::vector<arma::vec> S_a;    // P vectors, p_t
    std::vector<TraitMeta> meta;   // internal order
    std::vector<int> outOrder;     // internal index -> config order

    std::vector<int> batchTraits;       // batchable, binary first
    std::vector<int> batchQuantTraits;  // batchable and quantitative
    std::vector<int> scalarTraits;      // not batchable
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
// (design section 4.3): different sample IDs (content OR order), different n,
// or different impute_method. Warns -- does not fail -- when the models
// disagree on whether LOCO really applied. `t_names` supplies the trait label
// used in the messages and must be the same length as t_nms.
void validateMTModels(const std::vector<NullModelData>& t_nms,
                      const std::vector<std::string>& t_names);

// Static per-trait gate, evaluated once before the marker loop (design
// section 3.1). False means the trait runs the per-pair scalar path for every
// marker; it never changes per marker.
bool isBatchable(const TraitMeta& t_meta);

// Why isBatchable() said no; "-" when it said yes. For the startup table.
const char* batchableReason(const TraitMeta& t_meta);

// Design section 4.5: the table has to be printed, or a numerical problem in a
// P = 64 run cannot be localised to a trait.
void printMTGateTable(const std::vector<TraitMeta>& t_meta, bool t_locoEnabled);

// One per thread, reused across blocks (grow-only, never reallocated per block).
struct MTScratch {
    arma::mat Gb;      // N x B
    arma::mat Gb2;     // N x B   = Gb % Gb
    arma::mat Zall;    // sumP x B
    arma::mat GWbin;   // sumPbin x B
    arma::mat GWqnt;   // sumPqnt x B
    arma::mat GR;      // B x P
    arma::mat G2Mu2;   // B x nBin
    arma::vec Gsq;     // B
    arma::mat AC;      // B x 2*nBin
    arma::mat VR;      // B x P    per-pair variance ratio
    std::vector<std::pair<int,int>> fb;   // fallback queue: (column in block, internal trait)
};

}  // namespace SAIGE

#endif  // SAIGE_MT_HPP
