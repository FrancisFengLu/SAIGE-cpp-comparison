// saige_mt.cpp — multi-trait step 2, Phase 0 (config front end only).
// See MULTITRAIT_DESIGN.md.

#include "saige_mt.hpp"

#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <set>
#include <stdexcept>

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

}  // namespace SAIGE
