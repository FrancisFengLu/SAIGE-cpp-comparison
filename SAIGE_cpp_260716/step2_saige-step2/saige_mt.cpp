// saige_mt.cpp — multi-trait step 2, Phase 0 (config front end only).
// See MULTITRAIT_DESIGN.md.

#include "saige_mt.hpp"

#include <algorithm>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

#include "score_format.hpp"

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

void validateMTModels(const std::vector<NullModelData>& t_nms,
                      const std::vector<std::string>& t_names)
{
    const std::size_t P = t_nms.size();
    if (P == 0) throw std::runtime_error("validateMTModels: no models");
    if (t_names.size() != P)
        throw std::runtime_error("validateMTModels: name count does not match model count");
    if (P == 1) return;   // nothing to cross-check; never reached on the P == 1 path anyway

    const NullModelData& ref = t_nms[0];
    for (std::size_t i = 1; i < P; ++i) {
        const NullModelData& m = t_nms[i];
        const std::string where = "model '" + t_names[i] + "' (models[" +
                                  std::to_string(i) + "])";

        if (m.n != ref.n) {
            throw std::runtime_error(
                where + " has n=" + std::to_string(m.n) + " but '" + t_names[0] +
                "' has n=" + std::to_string(ref.n) +
                ". Multi-trait testing requires every model to be fitted on the "
                "same samples.");
        }
        if (m.sampleIDs.size() != ref.sampleIDs.size()) {
            throw std::runtime_error(
                where + " lists " + std::to_string(m.sampleIDs.size()) +
                " sample IDs but '" + t_names[0] + "' lists " +
                std::to_string(ref.sampleIDs.size()) + ".");
        }
        // Order matters, not just membership: the genotype reader builds ONE
        // sample-position map from model 0, so a permuted ID list would make
        // every other trait read a permuted genotype vector.
        for (std::size_t k = 0; k < m.sampleIDs.size(); ++k) {
            if (m.sampleIDs[k] != ref.sampleIDs[k]) {
                throw std::runtime_error(
                    where + " has a different sample list from '" + t_names[0] +
                    "': position " + std::to_string(k) + " is '" + m.sampleIDs[k] +
                    "' vs '" + ref.sampleIDs[k] +
                    "'. Multi-trait testing requires identical sample IDs in "
                    "identical order (subset traits are not supported).");
            }
        }
        if (m.impute_method != ref.impute_method) {
            throw std::runtime_error(
                where + " uses impute_method='" + m.impute_method + "' but '" +
                t_names[0] + "' uses '" + ref.impute_method +
                "'. Imputation happens once per marker and is shared by every "
                "trait, so the models must agree on it.");
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
                    const std::string& t_locoChrom)
{
    const int P = static_cast<int>(t_meta.size());
    if (P == 0 || static_cast<int>(t_order.size()) != P)
        throw std::runtime_error("buildMTContext: trait count mismatch");

    t_ctx.P = P;
    t_ctx.N = t_nms[t_order[0]].n;
    t_ctx.locoEnabled = t_locoEnabled;
    t_ctx.locoChrom = t_locoChrom;

    // ---- offsets ----
    // Internal order is binary, then quantitative, then the rest, so each
    // group's stacked columns are one contiguous range and every GEMM below is
    // a single BLAS call.
    int colOff = 0, binOff = 0, binIdx = 0;
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
    t_ctx.XVX.assign(P, arma::mat());
    t_ctx.S_a.assign(P, arma::vec());

    for (int t = 0; t < P; t++) {
        const NullModelData& nm = t_nms[t_order[t]];
        const TraitMeta& M = t_meta[t];
        const arma::uword a = static_cast<arma::uword>(M.colOff);
        const arma::uword b = a + static_cast<arma::uword>(M.p) - 1;
        if (nm.X.n_rows != N || nm.X.n_cols != static_cast<arma::uword>(M.p) ||
            nm.XVX_inv_XV.n_rows != N || nm.XVX_inv_XV.n_cols != static_cast<arma::uword>(M.p))
            throw std::runtime_error("buildMTContext: X / XVX_inv_XV shape mismatch for trait '" +
                                     M.name + "'");
        t_ctx.Xstack.cols(a, b) = nm.X;
        t_ctx.Astack.cols(a, b) = nm.XVX_inv_XV;
        t_ctx.RES.col(t) = nm.res;
        t_ctx.XVX[t] = nm.XVX;
        t_ctx.S_a[t] = nm.S_a;
        if (M.kind == TraitKind::Binary) {
            const arma::uword wa = static_cast<arma::uword>(M.binOff);
            const arma::uword wb = wa + static_cast<arma::uword>(M.p) - 1;
            // mu2_t scaled into each covariate column: this is what turns
            // sum_i mu2_i g_i B_i into one GEMM plus a p-length contraction.
            t_ctx.WXstack.cols(wa, wb) = nm.X.each_col() % nm.mu2;
            t_ctx.MU2bin.col(M.binIdx) = nm.mu2;
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

void scoreTestBatchMT(const MTContext& t_ctx,
                      const std::vector<int>& t_traitSet,
                      const arma::mat& t_Gb,
                      int t_j0, int t_j1,
                      const arma::mat& t_VR,
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
    bool anyBin = false, anyQnt = false;
    int c0 = INTMAX, c1 = 0, b0 = INTMAX, b1 = 0, q0 = INTMAX, q1 = 0;
    for (int t : t_traitSet) {
        const TraitMeta& M = t_ctx.meta[t];
        c0 = std::min(c0, M.colOff); c1 = std::max(c1, M.colOff + M.p);
        if (M.kind == TraitKind::Binary) {
            anyBin = true;
            b0 = std::min(b0, M.binOff); b1 = std::max(b1, M.binOff + M.p);
        } else {
            anyQnt = true;
            q0 = std::min(q0, M.colOff); q1 = std::max(q1, M.colOff + M.p);
        }
    }

    // Z = A^T G, one GEMM for every trait in the set (design section 2.3).
    t_scr.Zall = colView(t_ctx.Astack, c0, c1).t() * Gv;            // (c1-c0) x B

    // g^2 feeds both the binary sum_i mu2_i g_i^2 and the quantitative g'g.
    t_scr.Gb2 = Gv % Gv;                                            // N x B

    if (anyBin) {
        t_scr.GWbin  = colView(t_ctx.WXstack, b0, b1).t() * Gv;     // (b1-b0) x B
        t_scr.G2Mu2  = t_scr.Gb2.t() * t_ctx.MU2bin;                // B x nBin
    }
    if (anyQnt) {
        t_scr.GWqnt  = colView(t_ctx.Xstack, q0, q1).t() * Gv;      // (q1-q0) x B
        t_scr.Gsq    = arma::sum(t_scr.Gb2, 0).t();                 // B
    }
    t_scr.GR = Gv.t() * t_ctx.RES;                                  // B x P

    for (int t : t_traitSet) {
        const TraitMeta& M = t_ctx.meta[t];
        const arma::uword r0 = static_cast<arma::uword>(M.colOff - c0);
        const arma::uword r1 = r0 + static_cast<arma::uword>(M.p) - 1;
        const auto Z_t = t_scr.Zall.rows(r0, r1);                   // p x B

        // All three contractions are O(p*B); vectorised, never per pair.
        const arma::rowvec zxz = arma::sum(Z_t % (t_ctx.XVX[t] * Z_t), 0);
        const arma::rowvec saz = t_ctx.S_a[t].t() * Z_t;

        arma::rowvec gwz;
        if (M.kind == TraitKind::Binary) {
            const arma::uword w0 = static_cast<arma::uword>(M.binOff - b0);
            gwz = arma::sum(t_scr.GWbin.rows(w0, w0 + M.p - 1) % Z_t, 0);
        } else {
            const arma::uword w0 = static_cast<arma::uword>(M.colOff - q0);
            gwz = arma::sum(t_scr.GWqnt.rows(w0, w0 + M.p - 1) % Z_t, 0);
        }

        const arma::vec S = (t_scr.GR.col(t) - saz.t()) / M.tau0;
        arma::vec var2;
        if (M.kind == TraitKind::Binary) {
            var2 = zxz.t() + t_scr.G2Mu2.col(M.binIdx) - 2.0 * gwz.t();
        } else {
            var2 = zxz.t() * M.tau0 + t_scr.Gsq - 2.0 * gwz.t();
        }

        std::vector<std::string>& pstr = t_out.pvalStr[t];
        std::vector<char>&        plog = t_out.pvalIsLog[t];
        for (arma::uword j = 0; j < B; j++) {
            const arma::uword jo = static_cast<arma::uword>(t_j0) + j;
            const double v2 = var2[j];
            const double v1 = v2 * t_VR(jo, t);
            double Beta, seBeta, pval, TstatOut, var1Out, var2Out;
            bool islogp = false;
            std::string pvalStr;
            // Same function the scalar path calls, so the degenerate branches
            // (var1 <= DBL_MIN, non-finite stat, p underflow to the "%.1fE%d"
            // form) cannot drift apart.
            format_score_result(S[j], v1, v2, Beta, seBeta, pvalStr, pval,
                                islogp, TstatOut, var1Out, var2Out);
            t_out.Beta(jo, t)    = Beta;
            t_out.seBeta(jo, t)  = seBeta;
            t_out.Tstat(jo, t)   = TstatOut;
            t_out.var1(jo, t)    = var1Out;
            t_out.var2(jo, t)    = var2Out;
            t_out.StdStat(jo, t) = std::fabs(S[j]) / std::sqrt(v1);
            t_out.pvalRaw(jo, t) = pval;
            pstr[jo] = pvalStr;
            plog[jo] = islogp ? 1 : 0;
        }
    }
}

}  // namespace SAIGE
