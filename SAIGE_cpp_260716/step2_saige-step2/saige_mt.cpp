// saige_mt.cpp — multi-trait step 2, Phase 0 (config front end only).
// See MULTITRAIT_DESIGN.md.

#include "saige_mt.hpp"

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <stdexcept>

#include "score_format.hpp"

// Opt-in stage timers for the sample-space GEMMs. Not compiled into the
// shipped binary; build a throwaway copy with -DMTFOLD_PROF to split the
// marker-loop cost between Zall / GWqnt / the folded Z0 / GR.
#ifdef MTFOLD_PROF
#include <omp.h>
namespace SAIGE {
double g_mtfProfZall = 0.0, g_mtfProfGW = 0.0, g_mtfProfZ0 = 0.0, g_mtfProfGR = 0.0;
}
#define MTF_TIC()      const double mtfT0 = omp_get_wtime()
#define MTF_TOC(v)     do { const double mtfD = omp_get_wtime() - mtfT0; \
                            _Pragma("omp atomic") v += mtfD; } while (0)
#else
#define MTF_TIC()      do {} while (0)
#define MTF_TOC(v)     do {} while (0)
#endif

namespace SAIGE {

bool tryTraitKindFromString(const std::string& t_traitType, TraitKind& t_out)
{
    if (t_traitType == "binary")       { t_out = TraitKind::Binary;       return true; }
    if (t_traitType == "quantitative") { t_out = TraitKind::Quantitative; return true; }
    if (t_traitType == "survival")     { t_out = TraitKind::Survival;     return true; }
    return false;
}

TraitKind traitKindFromString(const std::string& t_traitType)
{
    TraitKind k;
    if (tryTraitKindFromString(t_traitType, k)) return k;
    throw std::runtime_error("Unsupported traitType: '" + t_traitType +
                             "' (expected binary / quantitative / survival)");
}

const char* traitKindName(TraitKind t_kind)
{
    switch (t_kind) {
        case TraitKind::Binary:       return "binary";
        case TraitKind::Quantitative: return "quantitative";
        case TraitKind::Survival:     return "survival";
    }
    return "unknown";
}

// Default trait label: the last path component of the model directory, with a
// trailing slash tolerated. Only used for logs and error messages.
static std::string defaultTraitName(const std::string& t_modelFile)
{
    std::filesystem::path p(t_modelFile);
    if (p.has_filename() && !p.filename().empty()) return p.filename().string();
    if (p.has_parent_path()) return p.parent_path().filename().string();
    return t_modelFile;
}

ModelOverrides parseModelOverrides(const YAML::Node& t_node)
{
    ModelOverrides ov;
    if (!t_node || !t_node.IsMap()) return ov;

    if (t_node["is_Firth_beta"]) {
        ov.has_is_Firth_beta = true;
        ov.is_Firth_beta = t_node["is_Firth_beta"].as<bool>();
    }
    if (t_node["pCutoffforFirth"]) {
        ov.has_pCutoffforFirth = true;
        ov.pCutoffforFirth = t_node["pCutoffforFirth"].as<double>();
    }
    if (t_node["isnoadjCov"]) {
        ov.has_isnoadjCov = true;
        ov.isnoadjCov = t_node["isnoadjCov"].as<bool>();
    }
    if (t_node["cateVarRatioMinMACVecExclude"] &&
        t_node["cateVarRatioMinMACVecExclude"].IsSequence()) {
        ov.has_cateVarRatioMinMACVecExclude = true;
        for (const auto& v : t_node["cateVarRatioMinMACVecExclude"])
            ov.cateVarRatioMinMACVecExclude.push_back(v.as<double>());
    }
    if (t_node["cateVarRatioMaxMACVecInclude"] &&
        t_node["cateVarRatioMaxMACVecInclude"].IsSequence()) {
        ov.has_cateVarRatioMaxMACVecInclude = true;
        for (const auto& v : t_node["cateVarRatioMaxMACVecInclude"])
            ov.cateVarRatioMaxMACVecInclude.push_back(v.as<double>());
    }
    return ov;
}

std::vector<MTModelSpec> parseModelSpecs(const YAML::Node& t_config)
{
    const bool hasModels = static_cast<bool>(t_config["models"]);
    const bool hasScalarModel  = static_cast<bool>(t_config["modelFile"]);
    const bool hasScalarVR     = static_cast<bool>(t_config["varianceRatioFile"]);
    const bool hasScalarOut    = static_cast<bool>(t_config["outputFile"]);
    const bool hasAnyScalar    = hasScalarModel || hasScalarVR || hasScalarOut;

    std::vector<MTModelSpec> specs;

    if (hasModels) {
        // Never guess which form the user meant.
        if (hasAnyScalar) {
            std::string which;
            if (hasScalarModel) which += " modelFile";
            if (hasScalarVR)    which += " varianceRatioFile";
            if (hasScalarOut)   which += " outputFile";
            throw std::runtime_error(
                "Config sets both `models:` and the single-trait key(s)" + which +
                ". They are mutually exclusive: use `models:` for one or more "
                "traits, or the three scalar keys for exactly one.");
        }
        const YAML::Node& models = t_config["models"];
        if (!models.IsSequence()) {
            throw std::runtime_error("Config key `models` must be a sequence of "
                                     "maps (one per trait).");
        }
        if (models.size() == 0) {
            throw std::runtime_error("Config key `models` is empty; it needs at "
                                     "least one model.");
        }
        for (std::size_t i = 0; i < models.size(); ++i) {
            const YAML::Node& m = models[i];
            const std::string where = "models[" + std::to_string(i) + "]";
            if (!m.IsMap()) {
                throw std::runtime_error(where + " must be a map with modelFile / "
                                         "varianceRatioFile / outputFile.");
            }
            MTModelSpec s;
            if (!m["modelFile"])
                throw std::runtime_error(where + " missing required key: modelFile");
            if (!m["varianceRatioFile"])
                throw std::runtime_error(where + " missing required key: varianceRatioFile");
            if (!m["outputFile"])
                throw std::runtime_error(where + " missing required key: outputFile");
            s.modelFile         = m["modelFile"].as<std::string>();
            s.varianceRatioFile = m["varianceRatioFile"].as<std::string>();
            s.outputFile        = m["outputFile"].as<std::string>();
            s.traitName = m["traitName"] ? m["traitName"].as<std::string>()
                                         : defaultTraitName(s.modelFile);
            s.ov = parseModelOverrides(m);
            specs.push_back(std::move(s));
        }
    } else {
        // Legacy single-trait form. Keep the three original error messages
        // verbatim so existing configs fail exactly as before.
        if (!hasScalarModel)
            throw std::runtime_error("Config missing required key: modelFile");
        if (!hasScalarVR)
            throw std::runtime_error("Config missing required key: varianceRatioFile");
        if (!hasScalarOut)
            throw std::runtime_error("Config missing required key: outputFile");
        MTModelSpec s;
        s.modelFile         = t_config["modelFile"].as<std::string>();
        s.varianceRatioFile = t_config["varianceRatioFile"].as<std::string>();
        s.outputFile        = t_config["outputFile"].as<std::string>();
        s.traitName = t_config["traitName"] ? t_config["traitName"].as<std::string>()
                                            : defaultTraitName(s.modelFile);
        // No per-model overrides in the legacy form: the top-level keys are
        // applied by main() exactly where they always were, so this stays empty
        // and the P == 1 legacy path is bit-identical.
        specs.push_back(std::move(s));
    }

    // Two traits writing to one file would interleave silently.
    std::set<std::string> seenOut;
    for (std::size_t i = 0; i < specs.size(); ++i) {
        if (!seenOut.insert(specs[i].outputFile).second) {
            throw std::runtime_error("models[" + std::to_string(i) +
                                     "] (" + specs[i].traitName + ") reuses outputFile '" +
                                     specs[i].outputFile + "' already claimed by an "
                                     "earlier model; every trait needs its own file.");
        }
    }
    return specs;
}

// ------------------------------------------------------------------
// Model set: ordering, validation, static gating
// ------------------------------------------------------------------

std::vector<int> mtInternalOrder(const std::vector<NullModelData>& t_nms)
{
    const int P = static_cast<int>(t_nms.size());
    std::vector<int> order;
    order.reserve(P);
    // Pass 1 binary, pass 2 quantitative, pass 3 the rest (survival, or a
    // traitType string the loader accepted but we do not classify). Stable by
    // construction: each pass walks the config order.
    for (int i = 0; i < P; ++i) if (t_nms[i].traitType == "binary")       order.push_back(i);
    for (int i = 0; i < P; ++i) if (t_nms[i].traitType == "quantitative") order.push_back(i);
    for (int i = 0; i < P; ++i) {
        if (t_nms[i].traitType != "binary" && t_nms[i].traitType != "quantitative")
            order.push_back(i);
    }
    return order;
}

bool validateMTModels(const std::vector<NullModelData>& t_nms,
                      const std::vector<std::string>& t_names,
                      bool t_requireSameSamples)
{
    const std::size_t P = t_nms.size();
    if (P == 0) throw std::runtime_error("validateMTModels: no models");
    if (t_names.size() != P)
        throw std::runtime_error("validateMTModels: name count does not match model count");
    if (P == 1) return false;   // nothing to cross-check; never reached on the P == 1 path anyway

    const NullModelData& ref = t_nms[0];
    bool differ = false;
    for (std::size_t i = 1; i < P; ++i) {
        const NullModelData& m = t_nms[i];
        const std::string where = "model '" + t_names[i] + "' (models[" +
                                  std::to_string(i) + "])";

        if (m.impute_method != ref.impute_method) {
            throw std::runtime_error(
                where + " uses impute_method='" + m.impute_method + "' but '" +
                t_names[0] + "' uses '" + ref.impute_method +
                "'. Imputation happens once per marker and is shared by every "
                "trait, so the models must agree on it.");
        }
        // Order matters, not just membership: two lists holding the same IDs in
        // a different order are different sample lists (design 4.3 / I2).
        if (m.sampleIDs != ref.sampleIDs) {
            if (t_requireSameSamples) {
                std::string detail;
                if (m.sampleIDs.size() != ref.sampleIDs.size()) {
                    detail = " lists " + std::to_string(m.sampleIDs.size()) +
                             " sample IDs but '" + t_names[0] + "' lists " +
                             std::to_string(ref.sampleIDs.size());
                } else {
                    std::size_t k = 0;
                    while (m.sampleIDs[k] == ref.sampleIDs[k]) ++k;
                    detail = " has a different sample list from '" + t_names[0] +
                             "': position " + std::to_string(k) + " is '" +
                             m.sampleIDs[k] + "' vs '" + ref.sampleIDs[k] + "'";
                }
                throw std::runtime_error(
                    where + detail + ". mtRequireSameSamples: true requires "
                    "identical sample IDs in identical order.");
            }
            differ = true;
        } else if (m.n != ref.n) {
            throw std::runtime_error(
                where + " has n=" + std::to_string(m.n) + " but '" + t_names[0] +
                "' has n=" + std::to_string(ref.n) + " although both list the "
                "same sample IDs.");
        }
    }

    if (differ) {
        // The union / per-trait index maps are built from sampleIDs, so every
        // list has to be present, match its model's n, and name each sample
        // once -- a repeated ID would make "the trait's k-th sample" ambiguous.
        for (std::size_t i = 0; i < P; ++i) {
            const NullModelData& m = t_nms[i];
            const std::string where = "model '" + t_names[i] + "' (models[" +
                                      std::to_string(i) + "])";
            if (m.sampleIDs.empty() || (int)m.sampleIDs.size() != m.n) {
                throw std::runtime_error(
                    where + " has n=" + std::to_string(m.n) + " but lists " +
                    std::to_string(m.sampleIDs.size()) + " sample IDs. With "
                    "different sample sets across models every model must list "
                    "exactly its n sample IDs.");
            }
            std::set<std::string> seen;
            for (const auto& id : m.sampleIDs) {
                if (!seen.insert(id).second) {
                    throw std::runtime_error(
                        where + " lists sample ID '" + id + "' more than once.");
                }
            }
        }
    }

    // LOCO mix: legal (the loader silently falls back to the genome-wide fit
    // when `chrom` is not in that model's loco_chroms, matching R), but silent
    // is exactly what makes it dangerous across P models.
    std::vector<std::size_t> applied, notApplied;
    for (std::size_t i = 0; i < P; ++i)
        (t_nms[i].loco_applied ? applied : notApplied).push_back(i);
    if (!applied.empty() && !notApplied.empty()) {
        std::cout << "WARNING: " << notApplied.size() << " of " << P
                  << " models fell back to the genome-wide fit for this chromosome "
                     "while " << applied.size() << " used their chr<N>/ files."
                  << std::endl;
        std::cout << "WARNING:   genome-wide fallback:";
        for (std::size_t i : notApplied) std::cout << " " << t_names[i];
        std::cout << std::endl;
    }

    // A sparse GRM per trait means P copies of m_spSigmaMat resident at once.
    std::size_t nSparse = 0;
    for (std::size_t i = 0; i < P; ++i) if (t_nms[i].dimNum > 0) ++nSparse;
    if (nSparse > 0) {
        std::cout << "WARNING: " << nSparse << " of " << P
                  << " models carry a sparse GRM; each is held in full for the "
                     "duration of the run." << std::endl;
    }
    return differ;
}

std::vector<std::string> mtUnionSampleIDs(const std::vector<NullModelData>& t_nms)
{
    std::vector<std::string> out;
    std::unordered_set<std::string> seen;
    for (const auto& m : t_nms) {
        for (const auto& id : m.sampleIDs) {
            if (seen.insert(id).second) out.push_back(id);
        }
    }
    return out;
}

bool isBatchable(const TraitMeta& t_meta)
{
    return batchableReason(t_meta)[0] == '-';
}

const char* batchableReason(const TraitMeta& t_meta)
{
    if (t_meta.kind == TraitKind::Survival)  return "survival";
    if (t_meta.isCondition)                  return "isCondition=true";
    if (t_meta.isnoadjCov)                   return "isnoadjCov=true";
    // First-pass flagSparseGRM_cur, mirroring main()'s per-marker ctx: with
    // isFastTest the first pass is forced onto the dense path, which is the one
    // the batch kernel implements.
    if (!t_meta.isFastTest && t_meta.flagSparseGRM) return "sparseGRM first pass";
    return "-";
}

void printMTGateTable(const std::vector<TraitMeta>& t_meta, bool t_locoEnabled)
{
    const std::size_t P = t_meta.size();
    std::cout << "===== Multi-trait: " << P << " models =====" << std::endl;
    std::cout << "  idx  name                 type            p  batch";
    if (t_locoEnabled) std::cout << "  loco";
    std::cout << "  reason" << std::endl;
    int nBin = 0, nQnt = 0, nBatch = 0;
    for (std::size_t t = 0; t < P; ++t) {
        const TraitMeta& M = t_meta[t];
        std::cout << "  " << std::setw(3) << t << "  " << std::left << std::setw(20)
                  << M.name.substr(0, 20) << std::right << " " << std::setw(13)
                  << M.traitType.substr(0, 13) << " " << std::setw(3) << M.p
                  << "  " << std::setw(5) << (M.batchable ? "yes" : "no");
        if (t_locoEnabled) std::cout << "  " << std::setw(4) << (M.locoApplied ? "chr" : "gw");
        std::cout << "  " << batchableReason(M) << std::endl;
        if (M.batchable) {
            ++nBatch;
            if (M.kind == TraitKind::Binary) ++nBin; else ++nQnt;
        }
    }
    std::cout << "  batch traits: " << nBatch << " (binary " << nBin
              << " / quantitative " << nQnt << "), scalar traits: "
              << (P - nBatch) << std::endl;
}


// ------------------------------------------------------------------
// Batch kernel
// ------------------------------------------------------------------

void buildMTContext(MTContext& t_ctx,
                    const std::vector<NullModelData>& t_nms,
                    const std::vector<int>& t_order,
                    std::vector<TraitMeta>& t_meta,
                    bool t_locoEnabled,
                    const std::string& t_locoChrom,
                    const std::vector<std::string>& t_unionIDs,
                    bool t_foldQuantProj)
{
    const int P = static_cast<int>(t_meta.size());
    if (P == 0 || static_cast<int>(t_order.size()) != P)
        throw std::runtime_error("buildMTContext: trait count mismatch");

    t_ctx.P = P;
    t_ctx.locoEnabled = t_locoEnabled;
    t_ctx.locoChrom = t_locoChrom;

    // ---- sample maps (design 4.7) ----
    // A trait whose list equals the union's, element for element, keeps the
    // identity map and is flagged sameAsUnion; when every trait is, this whole
    // block reduces to the same-sample-set layout that existed before.
    t_ctx.unionIDs = t_unionIDs;
    // Models with no sampleIDs list (hand-built; step 1 and rda_to_arma.R both
    // write one) give an empty union. validateMTModels only lets that through
    // when every list is empty and every n agrees, so each trait is the
    // same-set identity layout at that n -- which is what N was before sample
    // maps existed. Taking the union's length there would make N zero.
    t_ctx.N = t_unionIDs.empty() ? t_nms[t_order[0]].n
                                 : static_cast<int>(t_unionIDs.size());
    t_ctx.samp.assign(P, MTTraitSamples());
    t_ctx.sampleSetsDiffer = false;
    {
        std::unordered_map<std::string, arma::uword> unionPos;
        unionPos.reserve(t_unionIDs.size() * 2);
        for (std::size_t u = 0; u < t_unionIDs.size(); ++u) unionPos[t_unionIDs[u]] = u;
        for (int t = 0; t < P; t++) {
            const NullModelData& nm = t_nms[t_order[t]];
            MTTraitSamples& S = t_ctx.samp[t];
            S.sameAsUnion = (nm.sampleIDs == t_unionIDs);
            if (S.sameAsUnion) {
                S.n = t_ctx.N;
                continue;
            }
            t_ctx.sampleSetsDiffer = true;
            S.n = static_cast<int>(nm.sampleIDs.size());
            S.pos.resize(nm.sampleIDs.size());
            std::vector<char> inTrait(t_unionIDs.size(), 0);
            for (std::size_t k = 0; k < nm.sampleIDs.size(); ++k) {
                auto it = unionPos.find(nm.sampleIDs[k]);
                if (it == unionPos.end())
                    throw std::runtime_error("buildMTContext: sample '" + nm.sampleIDs[k] +
                                             "' of trait '" + t_meta[t].name +
                                             "' is not in the union sample list");
                S.pos[k] = it->second;
                inTrait[it->second] = 1;
            }
            for (std::size_t u = 0; u < t_unionIDs.size(); ++u)
                if (!inTrait[u]) S.comp.push_back(u);
        }
    }

    // ---- offsets ----
    // Internal order is binary, then quantitative, then the rest, so each
    // group's stacked columns are one contiguous range and every GEMM below is
    // a single BLAS call.
    int colOff = 0, binOff = 0, binIdx = 0, maskIdx = 0;
    for (int t = 0; t < P; t++) {
        TraitMeta& M = t_meta[t];
        M.colOff = colOff;
        colOff += M.p;
        if (M.kind == TraitKind::Binary) {
            M.binOff = binOff;
            M.binIdx = binIdx;
            binOff += M.p;
            binIdx += 1;
        } else {
            M.binOff = -1;
            M.binIdx = -1;
        }
        M.maskIdx = -1;
        if (M.kind == TraitKind::Quantitative && !t_ctx.samp[t].sameAsUnion)
            M.maskIdx = maskIdx++;
    }
    t_ctx.sumP    = colOff;
    t_ctx.sumPbin = binOff;
    t_ctx.nBin    = binIdx;
    t_ctx.qOff    = binOff;
    t_ctx.sumPqnt = 0;
    for (int t = 0; t < P; t++)
        if (t_meta[t].kind == TraitKind::Quantitative) t_ctx.sumPqnt += t_meta[t].p;

    const arma::uword N = static_cast<arma::uword>(t_ctx.N);

    // ---- stacks, filled in place (never via a per-trait temporary) ----
    t_ctx.Xstack.set_size(N, t_ctx.sumP);
    t_ctx.Astack.set_size(N, t_ctx.sumP);
    t_ctx.RES.set_size(N, P);
    if (t_ctx.sumPbin > 0) t_ctx.WXstack.set_size(N, t_ctx.sumPbin);
    if (t_ctx.nBin > 0)    t_ctx.MU2bin.set_size(N, t_ctx.nBin);
    if (maskIdx > 0)       t_ctx.MASKq.zeros(N, maskIdx);
    t_ctx.XVX.assign(P, arma::mat());
    t_ctx.S_a.assign(P, arma::vec());
    t_ctx.sumA.assign(P, arma::vec());
    t_ctx.sumW.assign(P, arma::vec());
    t_ctx.sumR.assign(P, 0.0);
    t_ctx.sumM.assign(P, 0.0);

    for (int t = 0; t < P; t++) {
        const NullModelData& nm = t_nms[t_order[t]];
        const TraitMeta& M = t_meta[t];
        const MTTraitSamples& S = t_ctx.samp[t];
        const arma::uword nT = S.sameAsUnion ? N : S.pos.size();
        const arma::uword a = static_cast<arma::uword>(M.colOff);
        const arma::uword b = a + static_cast<arma::uword>(M.p) - 1;
        if (nm.X.n_rows != nT || nm.X.n_cols != static_cast<arma::uword>(M.p) ||
            nm.XVX_inv_XV.n_rows != nT || nm.XVX_inv_XV.n_cols != static_cast<arma::uword>(M.p) ||
            nm.res.n_elem != nT || (M.kind == TraitKind::Binary && nm.mu2.n_elem != nT))
            throw std::runtime_error("buildMTContext: X / XVX_inv_XV / res / mu2 shape mismatch for trait '" +
                                     M.name + "'");
        t_ctx.XVX[t] = nm.XVX;
        t_ctx.S_a[t] = nm.S_a;
        if (S.sameAsUnion) {
            t_ctx.Xstack.cols(a, b) = nm.X;
            t_ctx.Astack.cols(a, b) = nm.XVX_inv_XV;
            t_ctx.RES.col(t) = nm.res;
            if (M.kind == TraitKind::Binary) {
                const arma::uword wa = static_cast<arma::uword>(M.binOff);
                const arma::uword wb = wa + static_cast<arma::uword>(M.p) - 1;
                // mu2_t scaled into each covariate column: this is what turns
                // sum_i mu2_i g_i B_i into one GEMM plus a p-length contraction.
                t_ctx.WXstack.cols(wa, wb) = nm.X.each_col() % nm.mu2;
                t_ctx.MU2bin.col(M.binIdx) = nm.mu2;
            }
        } else {
            // Embed: the trait's k-th row goes to union row pos[k]; every other
            // row is an exact zero, so a GEMM against a union-length genotype
            // column only ever collects this trait's own samples.
            t_ctx.Xstack.cols(a, b).zeros();
            t_ctx.Astack.cols(a, b).zeros();
            t_ctx.RES.col(t).zeros();
            for (arma::uword c = 0; c < (arma::uword)M.p; ++c) {
                double* xd = t_ctx.Xstack.colptr(a + c);
                double* ad = t_ctx.Astack.colptr(a + c);
                const double* xs = nm.X.colptr(c);
                const double* as = nm.XVX_inv_XV.colptr(c);
                for (arma::uword k = 0; k < nT; ++k) { xd[S.pos[k]] = xs[k]; ad[S.pos[k]] = as[k]; }
            }
            {
                double* rd = t_ctx.RES.colptr(t);
                for (arma::uword k = 0; k < nT; ++k) rd[S.pos[k]] = nm.res[k];
            }
            if (M.kind == TraitKind::Binary) {
                const arma::uword wa = static_cast<arma::uword>(M.binOff);
                t_ctx.MU2bin.col(M.binIdx).zeros();
                double* md = t_ctx.MU2bin.colptr(M.binIdx);
                for (arma::uword k = 0; k < nT; ++k) md[S.pos[k]] = nm.mu2[k];
                for (arma::uword c = 0; c < (arma::uword)M.p; ++c) {
                    t_ctx.WXstack.col(wa + c).zeros();
                    double* wd = t_ctx.WXstack.colptr(wa + c);
                    const double* xs = nm.X.colptr(c);
                    for (arma::uword k = 0; k < nT; ++k) wd[S.pos[k]] = xs[k] * nm.mu2[k];
                }
            }
            if (M.maskIdx >= 0) {
                double* kd = t_ctx.MASKq.colptr(M.maskIdx);
                for (arma::uword k = 0; k < nT; ++k) kd[S.pos[k]] = 1.0;
            }
            // Constant half of the flip correction: sums over the trait's
            // samples of every column the kernel contracts against.
            t_ctx.sumA[t] = arma::sum(t_ctx.Astack.cols(a, b), 0).t();
            t_ctx.sumR[t] = arma::accu(t_ctx.RES.col(t));
            if (M.kind == TraitKind::Binary) {
                const arma::uword wa = static_cast<arma::uword>(M.binOff);
                t_ctx.sumW[t] = arma::sum(t_ctx.WXstack.cols(wa, wa + M.p - 1), 0).t();
                t_ctx.sumM[t] = arma::accu(t_ctx.MU2bin.col(M.binIdx));
            } else {
                t_ctx.sumW[t] = arma::sum(t_ctx.Xstack.cols(a, b), 0).t();
                t_ctx.sumM[t] = static_cast<double>(nT);
            }
        }
    }

    // ---- static gating partitions (design sections 3.1 / 3.2) ----
    t_ctx.batchTraits.clear();
    t_ctx.batchQuantTraits.clear();
    t_ctx.scalarTraits.clear();
    for (int t = 0; t < P; t++) {
        if (!t_meta[t].batchable) { t_ctx.scalarTraits.push_back(t); continue; }
        t_ctx.batchTraits.push_back(t);
        if (t_meta[t].kind != TraitKind::Binary) t_ctx.batchQuantTraits.push_back(t);
    }
    t_ctx.meta = t_meta;

    // ---- folded covariate projection (config mtFoldQuantProj) ----
    // See the block comment on MTContext::foldQuant. Everything here is a
    // one-off O(N p^2) per trait; the marker loop only ever reads foldable /
    // foldK / foldRefCol0.
    t_ctx.foldQuant   = false;
    t_ctx.foldRef     = -1;
    t_ctx.foldRefCol0 = 0;
    t_ctx.foldRefP    = 0;
    t_ctx.foldable.assign(P, 0);
    t_ctx.foldK.assign(P, arma::mat());
    t_ctx.foldResid.assign(P, -1.0);
    t_ctx.foldTraits.clear();
    if (t_foldQuantProj) {
        for (int t : t_ctx.batchQuantTraits) {
            if (!t_ctx.samp[t].sameAsUnion) continue;   // its stack block is padded with zeros
            t_ctx.foldRef = t;
            break;
        }
        if (t_ctx.foldRef >= 0) {
            const TraitMeta& R  = t_ctx.meta[t_ctx.foldRef];
            const arma::uword pr = static_cast<arma::uword>(R.p);
            const arma::uword ar = static_cast<arma::uword>(R.colOff);
            t_ctx.foldRefCol0 = R.colOff;
            t_ctx.foldRefP    = R.p;
            const arma::mat Xr  = t_ctx.Xstack.cols(ar, ar + pr - 1);   // N x p
            const arma::mat XtX = Xr.t() * Xr;
            for (int t : t_ctx.batchQuantTraits) {
                const TraitMeta& M = t_ctx.meta[t];
                if (!t_ctx.samp[t].sameAsUnion) continue;
                if (static_cast<arma::uword>(M.p) != pr) continue;
                const arma::uword a = static_cast<arma::uword>(M.colOff);
                // X_t' G is only Xref' G if the two blocks are the same bits.
                bool sameX = true;
                for (arma::uword c = 0; c < pr && sameX; ++c) {
                    const double* u = t_ctx.Xstack.colptr(a + c);
                    const double* v = t_ctx.Xstack.colptr(ar + c);
                    if (u != v && std::memcmp(u, v, N * sizeof(double)) != 0) sameX = false;
                }
                if (!sameX) continue;
                const arma::mat At = t_ctx.Astack.cols(a, a + pr - 1);
                arma::mat K;
                if (!arma::solve(K, XtX, Xr.t() * At, arma::solve_opts::no_approx)) continue;
                const double scale = arma::abs(At).max();
                const double resid = arma::abs(At - Xr * K).max();
                const double rel   = (scale > 0.0) ? resid / scale : resid;
                t_ctx.foldResid[t] = rel;
                if (!(rel <= MT_FOLD_RESID_TOL)) continue;
                t_ctx.foldable[t] = 1;
                t_ctx.foldK[t]    = K;
                t_ctx.foldTraits.push_back(t);
            }
            t_ctx.foldQuant = !t_ctx.foldTraits.empty();
        }
    }
}

void MTBlockAdj::resize(int t_B, int t_P)
{
    const arma::uword B = static_cast<arma::uword>(t_B);
    const arma::uword P = static_cast<arma::uword>(t_P);
    a.set_size(B, P); b.set_size(B, P); d.set_size(B, P); q.set_size(B, P);
    nMiss.set_size(B, P);
    miss.resize(B);
}

void MTBlockResult::resize(int t_B, int t_P)
{
    const arma::uword B = static_cast<arma::uword>(t_B);
    const arma::uword P = static_cast<arma::uword>(t_P);
    // set_size on unchanged dimensions is a no-op in armadillo, so calling this
    // once per block reuses the buffers rather than reallocating them.
    Beta.set_size(B, P);   seBeta.set_size(B, P);
    Tstat.set_size(B, P);  var1.set_size(B, P);
    var2.set_size(B, P);   StdStat.set_size(B, P);
    pvalRaw.set_size(B, P);
    pvalStr.resize(P);
    pvalIsLog.resize(P);
    for (int t = 0; t < t_P; t++) {
        pvalStr[t].resize(B);
        pvalIsLog[t].assign(B, 0);
    }
}

// Non-owning view over a contiguous column range of a stacked matrix. Column
// major storage means cols [c0, c1) are one contiguous block, so this aliases
// rather than copies -- which is the whole point of stacking.
static inline arma::mat colView(const arma::mat& t_m, int t_c0, int t_c1)
{
    return arma::mat(const_cast<double*>(t_m.colptr(static_cast<arma::uword>(t_c0))),
                     t_m.n_rows, static_cast<arma::uword>(t_c1 - t_c0),
                     /*copy_aux_mem=*/false, /*strict=*/true);
}

// Sum of the rows listed in t_idx of the stack columns [t_c0, t_c1), written
// to t_out.col(t_j) for the first form and t_out.row(t_j) for the second.
// Outer loop over stack columns so each inner loop walks one column with
// increasing row index.
static inline void sumRowsIntoCol(const arma::mat& t_stack, int t_c0, int t_c1,
                                  const std::vector<arma::uword>& t_idx,
                                  arma::mat& t_out, arma::uword t_j)
{
    for (int s = t_c0; s < t_c1; ++s) {
        const double* col = t_stack.colptr(static_cast<arma::uword>(s));
        double acc = 0.0;
        for (arma::uword k : t_idx) acc += col[k];
        t_out(static_cast<arma::uword>(s - t_c0), t_j) = acc;
    }
}
static inline void sumRowsIntoRow(const arma::mat& t_stack,
                                  const std::vector<arma::uword>& t_idx,
                                  arma::mat& t_out, arma::uword t_j)
{
    for (arma::uword s = 0; s < t_stack.n_cols; ++s) {
        const double* col = t_stack.colptr(s);
        double acc = 0.0;
        for (arma::uword k : t_idx) acc += col[k];
        t_out(t_j, s) = acc;
    }
}

// The per-trait tail both batch kernels share: S and var2 for the block's
// columns turn into the printed row through the SAME format_score_result the
// scalar path calls, so the degenerate branches (var1 <= DBL_MIN, non-finite
// stat, p underflow to the "%.1fE%d" form) cannot drift apart between the two.
// Extracted, not duplicated, when the GPU kernel was added -- a second copy of
// this loop is exactly the kind of thing that silently diverges.
static void emitBlockResults(int t_t, int t_j0,
                             const arma::vec& t_S, const arma::vec& t_var2,
                             const arma::mat& t_VR, MTBlockResult& t_out)
{
    std::vector<std::string>& pstr = t_out.pvalStr[t_t];
    std::vector<char>&        plog = t_out.pvalIsLog[t_t];
    const arma::uword B = t_S.n_elem;
    for (arma::uword j = 0; j < B; j++) {
        const arma::uword jo = static_cast<arma::uword>(t_j0) + j;
        const double v2 = t_var2[j];
        const double v1 = v2 * t_VR(jo, t_t);
        double Beta, seBeta, pval, TstatOut, var1Out, var2Out;
        bool islogp = false;
        std::string pvalStr;
        format_score_result(t_S[j], v1, v2, Beta, seBeta, pvalStr, pval,
                            islogp, TstatOut, var1Out, var2Out);
        t_out.Beta(jo, t_t)    = Beta;
        t_out.seBeta(jo, t_t)  = seBeta;
        t_out.Tstat(jo, t_t)   = TstatOut;
        t_out.var1(jo, t_t)    = var1Out;
        t_out.var2(jo, t_t)    = var2Out;
        t_out.StdStat(jo, t_t) = std::fabs(t_S[j]) / std::sqrt(v1);
        t_out.pvalRaw(jo, t_t) = pval;
        pstr[jo] = pvalStr;
        plog[jo] = islogp ? 1 : 0;
    }
}

void scoreTestBatchMT(const MTContext& t_ctx,
                      const std::vector<int>& t_traitSet,
                      const arma::mat& t_Gb,
                      int t_j0, int t_j1,
                      const arma::mat& t_VR,
                      const MTBlockAdj* t_adj,
                      MTScratch& t_scr,
                      MTBlockResult& t_out)
{
    const arma::uword B = static_cast<arma::uword>(t_j1 - t_j0);
    if (t_traitSet.empty() || t_j1 <= t_j0) return;
    const arma::mat Gv = colView(t_Gb, t_j0, t_j1);   // N x B, no copy

    // Column ranges actually needed. A trait that is not in t_traitSet but lies
    // between two that are gets computed and ignored; that is cheaper than
    // splitting the GEMM.
    const int INTMAX = std::numeric_limits<int>::max();
    const bool foldOn = t_ctx.foldQuant;
    bool anyBin = false, anyQnt = false;
    bool anyFold = false, anyWideA = false;
    bool anyAdj = false, anyAdjBin = false, anyAdjQnt = false;
    int c0 = INTMAX, c1 = 0, b0 = INTMAX, b1 = 0, q0 = INTMAX, q1 = 0;
    for (int t : t_traitSet) {
        const TraitMeta& M = t_ctx.meta[t];
        const bool adj = !t_ctx.samp.empty() && !t_ctx.samp[t].sameAsUnion;
        // A folded trait reads neither Astack nor Xstack in the marker loop --
        // its whole sample-space contribution is the shared Z0 -- so it must
        // not widen either GEMM's column range. buildMTContext only ever marks
        // a quantitative, batchable, same-sample-set trait foldable, so `adj`
        // is false here by construction.
        const bool fold = foldOn && t_ctx.foldable[t];
        if (M.kind == TraitKind::Binary) anyBin = true; else anyQnt = true;
        if (fold) { anyFold = true; continue; }
        anyWideA = true;
        c0 = std::min(c0, M.colOff); c1 = std::max(c1, M.colOff + M.p);
        if (M.kind == TraitKind::Binary) {
            b0 = std::min(b0, M.binOff); b1 = std::max(b1, M.binOff + M.p);
            anyAdjBin = anyAdjBin || adj;
        } else {
            q0 = std::min(q0, M.colOff); q1 = std::max(q1, M.colOff + M.p);
            anyAdjQnt = anyAdjQnt || adj;
        }
        anyAdj = anyAdj || adj;
    }
    if (anyAdj && t_adj == nullptr)
        throw std::runtime_error("scoreTestBatchMT: traits with their own sample list need t_adj");

    // Z = A^T G, one GEMM for every trait in the set (design section 2.3).
    if (anyWideA) {
        MTF_TIC();
        t_scr.Zall = colView(t_ctx.Astack, c0, c1).t() * Gv;        // (c1-c0) x B
        MTF_TOC(SAIGE::g_mtfProfZall);
    }

    // g^2 feeds both the binary sum_i mu2_i g_i^2 and the quantitative g'g.
    t_scr.Gb2 = Gv % Gv;                                            // N x B

    if (anyBin) {
        t_scr.GWbin  = colView(t_ctx.WXstack, b0, b1).t() * Gv;     // (b1-b0) x B
        t_scr.G2Mu2  = t_scr.Gb2.t() * t_ctx.MU2bin;                // B x nBin
    }
    if (anyQnt) {
        if (q0 != INTMAX) {
            MTF_TIC();
            t_scr.GWqnt = colView(t_ctx.Xstack, q0, q1).t() * Gv;   // (q1-q0) x B
            MTF_TOC(SAIGE::g_mtfProfGW);
        }
        t_scr.Gsq    = arma::sum(t_scr.Gb2, 0).t();                 // B
    }
    // The one shared covariate GEMM: p columns, not sum_t p_t, and it stands in
    // for the folded traits' Astack block AND their Xstack block at once.
    if (anyFold) {
        MTF_TIC();
        t_scr.Z0 = colView(t_ctx.Xstack, t_ctx.foldRefCol0,
                           t_ctx.foldRefCol0 + t_ctx.foldRefP).t() * Gv;   // p x B
        MTF_TOC(SAIGE::g_mtfProfZ0);
    }
    {
        MTF_TIC();
        t_scr.GR = Gv.t() * t_ctx.RES;                              // B x P
        MTF_TOC(SAIGE::g_mtfProfGR);
    }

    // ---- different sample sets (design 4.7) ----
    // Only traits whose sample list differs from the union's read any of this.
    bool anyFlipAdj = false;
    if (anyAdj) {
        for (int t : t_traitSet) {
            if (t_ctx.samp[t].sameAsUnion) continue;
            for (arma::uword j = 0; j < B; ++j)
                if (t_adj->b(t_j0 + j, t) != 0.0) { anyFlipAdj = true; break; }
            if (anyFlipAdj) break;
        }
        const arma::uword nMask = t_ctx.MASKq.n_cols;
        if (anyAdjBin && anyFlipAdj) t_scr.GMu2 = Gv.t() * t_ctx.MU2bin;       // B x nBin
        if (anyAdjQnt && nMask > 0) {
            t_scr.GMask2 = t_scr.Gb2.t() * t_ctx.MASKq;                          // B x nMask
            if (anyFlipAdj) t_scr.GMask1 = Gv.t() * t_ctx.MASKq;                 // B x nMask
        }
        // Stack rows summed over each column's missing cells. Rows outside a
        // trait are zero in its stack columns, so summing over every missing
        // cell of the union column collects exactly the trait's own missing
        // cells.
        t_scr.MissA.zeros(static_cast<arma::uword>(c1 - c0), B);
        t_scr.MissR.zeros(B, t_ctx.RES.n_cols);
        if (anyAdjBin) {
            t_scr.MissWbin.zeros(static_cast<arma::uword>(b1 - b0), B);
            t_scr.MissMu2.zeros(B, t_ctx.MU2bin.n_cols);
        }
        if (anyAdjQnt) {
            t_scr.MissWqnt.zeros(static_cast<arma::uword>(q1 - q0), B);
            if (nMask > 0) t_scr.MissMask.zeros(B, nMask);
        }
        for (arma::uword j = 0; j < B; ++j) {
            const std::vector<arma::uword>& mv = t_adj->miss[t_j0 + j];
            if (mv.empty()) continue;
            sumRowsIntoCol(t_ctx.Astack, c0, c1, mv, t_scr.MissA, j);
            sumRowsIntoRow(t_ctx.RES, mv, t_scr.MissR, j);
            if (anyAdjBin) {
                sumRowsIntoCol(t_ctx.WXstack, b0, b1, mv, t_scr.MissWbin, j);
                sumRowsIntoRow(t_ctx.MU2bin, mv, t_scr.MissMu2, j);
            }
            if (anyAdjQnt) {
                sumRowsIntoCol(t_ctx.Xstack, q0, q1, mv, t_scr.MissWqnt, j);
                if (nMask > 0) sumRowsIntoRow(t_ctx.MASKq, mv, t_scr.MissMask, j);
            }
        }
    }

    for (int t : t_traitSet) {
        const TraitMeta& M = t_ctx.meta[t];
        if (foldOn && t_ctx.foldable[t]) {
            // Same three contractions as the !adj branch below, with
            //   A_t' G  ->  K_t' Z0        (Astack block == Xref K_t)
            //   X_t' G  ->  Z0             (Xstack block == Xref, bit for bit)
            t_scr.Zf = t_ctx.foldK[t].t() * t_scr.Z0;                // p x B
            const arma::mat& Z_t = t_scr.Zf;
            arma::rowvec zxz = arma::sum(Z_t % (t_ctx.XVX[t] * Z_t), 0);
            arma::rowvec saz = t_ctx.S_a[t].t() * Z_t;
            arma::rowvec gwz = arma::sum(t_scr.Z0 % Z_t, 0);
            arma::vec S    = (t_scr.GR.col(t) - saz.t()) / M.tau0;
            arma::vec var2 = zxz.t() * M.tau0 + t_scr.Gsq - 2.0 * gwz.t();
            emitBlockResults(t, t_j0, S, var2, t_VR, t_out);
            continue;
        }
        const arma::uword r0 = static_cast<arma::uword>(M.colOff - c0);
        const arma::uword r1 = r0 + static_cast<arma::uword>(M.p) - 1;
        const arma::uword w0 = (M.kind == TraitKind::Binary)
                                   ? static_cast<arma::uword>(M.binOff - b0)
                                   : static_cast<arma::uword>(M.colOff - q0);
        const arma::mat& GW = (M.kind == TraitKind::Binary) ? t_scr.GWbin : t_scr.GWqnt;
        const bool adj = anyAdj && !t_ctx.samp[t].sameAsUnion;

        arma::rowvec zxz, saz, gwz;
        arma::vec S, var2;
        if (!adj) {
            const auto Z_t = t_scr.Zall.rows(r0, r1);                   // p x B

            // All three contractions are O(p*B); vectorised, never per pair.
            zxz = arma::sum(Z_t % (t_ctx.XVX[t] * Z_t), 0);
            saz = t_ctx.S_a[t].t() * Z_t;
            gwz = arma::sum(GW.rows(w0, w0 + M.p - 1) % Z_t, 0);

            S = (t_scr.GR.col(t) - saz.t()) / M.tau0;
            if (M.kind == TraitKind::Binary) {
                var2 = zxz.t() + t_scr.G2Mu2.col(M.binIdx) - 2.0 * gwz.t();
            } else {
                var2 = zxz.t() * M.tau0 + t_scr.Gsq - 2.0 * gwz.t();
            }
        } else {
            // The trait's own genotype vector is g_t = a*g + b + d*[missing]
            // on its samples (MTBlockAdj). Every kernel input is linear or
            // quadratic in g, so each is mapped exactly:
            //   L(g_t)   = a*L(g) + b*L(1_t) + d*L(e_miss)       L = A', W', res'
            //   Q_v(g_t) = Q_v(g) + 2ab*L_v(g) + b^2*sum_t v + q*sum_miss v
            // where v is mu2_t (binary) or the trait's indicator (quantitative)
            // and q = g_t^2 - (a*g + b)^2 on a missing cell. a^2 = 1 has been
            // used. Terms are added only when they are not identically zero,
            // so which other pairs share the block cannot move a result.
            const arma::uword p = static_cast<arma::uword>(M.p);
            const arma::mat& MissW = (M.kind == TraitKind::Binary) ? t_scr.MissWbin : t_scr.MissWqnt;
            const arma::vec& sA = t_ctx.sumA[t];
            const arma::vec& sW = t_ctx.sumW[t];
            t_scr.Zc.set_size(p, B);
            t_scr.Wc.set_size(p, B);
            t_scr.Rc.set_size(B);
            t_scr.Qc.set_size(B);
            for (arma::uword j = 0; j < B; ++j) {
                const arma::uword jo = static_cast<arma::uword>(t_j0) + j;
                const double a = t_adj->a(jo, t), b = t_adj->b(jo, t);
                const double d = t_adj->d(jo, t), q = t_adj->q(jo, t);
                const bool flip = (a < 0.0), shift = (b != 0.0);
                const bool miss = (t_adj->nMiss(jo, t) > 0);
                for (arma::uword r = 0; r < p; ++r) {
                    double z = t_scr.Zall(r0 + r, j);
                    double w = GW(w0 + r, j);
                    if (flip)  { z = -z; w = -w; }
                    if (shift) { z += b * sA[r]; w += b * sW[r]; }
                    if (miss)  { z += d * t_scr.MissA(r0 + r, j); w += d * MissW(w0 + r, j); }
                    t_scr.Zc(r, j) = z;
                    t_scr.Wc(r, j) = w;
                }
                double R = t_scr.GR(j, t);
                if (flip)  R = -R;
                if (shift) R += b * t_ctx.sumR[t];
                if (miss)  R += d * t_scr.MissR(j, t);
                t_scr.Rc[j] = R;
                double Q;
                if (M.kind == TraitKind::Binary) {
                    Q = t_scr.G2Mu2(j, M.binIdx);
                    if (shift) Q += 2.0 * a * b * t_scr.GMu2(j, M.binIdx) + b * b * t_ctx.sumM[t];
                    if (miss)  Q += q * t_scr.MissMu2(j, M.binIdx);
                } else {
                    Q = t_scr.GMask2(j, M.maskIdx);
                    if (shift) Q += 2.0 * a * b * t_scr.GMask1(j, M.maskIdx) + b * b * t_ctx.sumM[t];
                    if (miss)  Q += q * t_scr.MissMask(j, M.maskIdx);
                }
                t_scr.Qc[j] = Q;
            }
            const arma::mat& Z_t = t_scr.Zc;
            zxz = arma::sum(Z_t % (t_ctx.XVX[t] * Z_t), 0);
            saz = t_ctx.S_a[t].t() * Z_t;
            gwz = arma::sum(t_scr.Wc % Z_t, 0);

            S = (t_scr.Rc - saz.t()) / M.tau0;
            if (M.kind == TraitKind::Binary) {
                var2 = zxz.t() + t_scr.Qc - 2.0 * gwz.t();
            } else {
                var2 = zxz.t() * M.tau0 + t_scr.Qc - 2.0 * gwz.t();
            }
        }

        emitBlockResults(t, t_j0, S, var2, t_VR, t_out);
    }
}


// ---------------------------------------------------------------------------
// GPU-fed variant of scoreTestBatchMT: the sample-space reductions are already
// in t_scr, so this does only the O(p^2) / O(P) tail.
//
// Scope, enforced by the caller's gate (main.cpp, mainMarkerMT): every trait in
// t_traitSet is quantitative, batchable, and has the union's sample list. That
// is the case in which scoreTestBatchMT's body reduces to Zall / GWqnt / Gsq /
// GR -- no MU2bin, no MASKq, no MTBlockAdj -- which is exactly what the GPU
// computes. Anything else keeps the CPU kernel.
//
// t_scr must hold, for the block columns [t_j0, t_j1) indexed 0-based:
//   Zall   sumP x B    Astack^T G
//   GWqnt  sumP x B    Xstack^T G
//   GR     B x P       G^T RES      (only the t_traitSet columns are read)
//   Gsq    B           colsum(G % G)
// with STACK rows, i.e. row r is stack column r -- so the traits in t_traitSet
// must start at stack column 0. Throws if they do not, rather than reading the
// wrong rows.
void scoreTestBatchMTQuantPre(const MTContext& t_ctx,
                              const std::vector<int>& t_traitSet,
                              int t_j0, int t_j1,
                              const arma::mat& t_VR,
                              MTScratch& t_scr,
                              MTBlockResult& t_out)
{
    const arma::uword B = static_cast<arma::uword>(t_j1 - t_j0);
    if (t_traitSet.empty() || t_j1 <= t_j0) return;

    int c0 = std::numeric_limits<int>::max();
    for (int t : t_traitSet) {
        const TraitMeta& M = t_ctx.meta[t];
        if (M.kind != TraitKind::Quantitative)
            throw std::runtime_error("scoreTestBatchMTQuantPre: non-quantitative trait");
        if (!t_ctx.samp.empty() && !t_ctx.samp[t].sameAsUnion)
            throw std::runtime_error("scoreTestBatchMTQuantPre: trait has its own sample list");
        c0 = std::min(c0, M.colOff);
    }
    if (c0 != 0)
        throw std::runtime_error("scoreTestBatchMTQuantPre: trait set does not start at stack column 0");

    for (int t : t_traitSet) {
        const TraitMeta& M = t_ctx.meta[t];
        const arma::uword r0 = static_cast<arma::uword>(M.colOff);
        const arma::uword r1 = r0 + static_cast<arma::uword>(M.p) - 1;
        const auto Z_t = t_scr.Zall.rows(r0, r1);                    // p x B

        // Identical expressions to scoreTestBatchMT's non-adjusted branch.
        arma::rowvec zxz = arma::sum(Z_t % (t_ctx.XVX[t] * Z_t), 0);
        arma::rowvec saz = t_ctx.S_a[t].t() * Z_t;
        arma::rowvec gwz = arma::sum(t_scr.GWqnt.rows(r0, r1) % Z_t, 0);

        arma::vec S    = (t_scr.GR.col(t) - saz.t()) / M.tau0;
        arma::vec var2 = zxz.t() * M.tau0 + t_scr.Gsq - 2.0 * gwz.t();
        if (S.n_elem != B || var2.n_elem != B)
            throw std::runtime_error("scoreTestBatchMTQuantPre: prefilled scratch has the wrong width");

        emitBlockResults(t, t_j0, S, var2, t_VR, t_out);
    }
}

}  // namespace SAIGE
