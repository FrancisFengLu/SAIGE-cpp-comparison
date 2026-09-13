// saige_mt.cpp — multi-trait step 2, Phase 0 (config front end only).
// See MULTITRAIT_DESIGN.md.

#include "saige_mt.hpp"

#include <filesystem>
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

}  // namespace SAIGE
