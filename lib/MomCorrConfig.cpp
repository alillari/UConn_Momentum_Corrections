#include "MomCorrConfig.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <unordered_set>

namespace {

bool HasKey(const json& object, const std::string& key) {
    return object.find(key) != object.end();
}

PhiHandling PhiHandlingFromString(const std::string& input) {
    if (input == "None") return PhiHandling::None;
    if (input == "CLAS12_CD_Standard") return PhiHandling::CLAS12_CD_Standard;
    if (input == "CLAS12_FD_El_Standard") return PhiHandling::CLAS12_FD_El_Standard;

    throw std::runtime_error("Unknown phi handling '" + input + "'");
}

std::vector<std::string> StringVectorFromJson(const json& value,
                                              const std::string& path) {
    if (!value.is_array()) {
        throw std::runtime_error(path + " must be an array of strings");
    }

    std::vector<std::string> out;
    for (const auto& item : value) {
        if (!item.is_string()) {
            throw std::runtime_error(path + " must contain only strings");
        }
        out.push_back(item.get<std::string>());
    }
    return out;
}

} // namespace

MomCorrConfig::MomCorrConfig(const std::string& jsonFile) {
    LoadJSON(jsonFile);
}

void MomCorrConfig::ValidateKnownKeys(const json& object,
                                      const std::vector<std::string>& allowed,
                                      const std::string& path) const {
    if (!object.is_object()) {
        throw std::runtime_error(path + " must be a JSON object");
    }

    for (const auto& item : object.items()) {
        if (std::find(allowed.begin(), allowed.end(), item.key()) == allowed.end()) {
            throw std::runtime_error(path + " has unknown key '" + item.key() + "'");
        }
    }
}

AxisConfig MomCorrConfig::ParseAxis(const json& j, const std::string& path) const {
    ValidateKnownKeys(j, {"bins", "low", "high"}, path);

    AxisConfig axis;
    axis.bins = j.at("bins").get<int>();
    axis.low = j.at("low").get<double>();
    axis.high = j.at("high").get<double>();
    ValidateAxis(axis, path);
    return axis;
}

void MomCorrConfig::ValidateAxis(const AxisConfig& axis, const std::string& path) const {
    if (axis.bins <= 0) {
        throw std::runtime_error(path + ".bins must be positive");
    }
    if (!std::isfinite(axis.low) || !std::isfinite(axis.high) || axis.high <= axis.low) {
        throw std::runtime_error(path + " must have finite low/high with high > low");
    }
}

void MomCorrConfig::LoadJSON(const std::string& filename) {
    std::ifstream f(filename);

    if (!f.is_open()) {
        throw std::runtime_error("MomCorrConfig: cannot open config file " + filename);
    }

    json j;
    f >> j;
    storedJSON_ = j;

    ValidateKnownKeys(
        j,
        {"schema_version", "dataset", "analysis", "detector", "particles",
         "observables", "histograms", "outputs"},
        "config"
    );

    schemaVersion_ = j.at("schema_version").get<int>();
    if (schemaVersion_ != 2) {
        throw std::runtime_error("Only schema_version 2 is supported");
    }

    const auto& dataset = j.at("dataset");
    ValidateKnownKeys(dataset, {"tree"}, "dataset");
    treeName_ = dataset.at("tree").get<std::string>();
    if (treeName_.empty()) {
        throw std::runtime_error("dataset.tree must not be empty");
    }

    const auto& a = j.at("analysis");
    ValidateKnownKeys(
        a,
        {"experiment", "specifier", "beam_energy", "bending", "channel",
         "momentum_bin"},
        "analysis"
    );
    experiment_ = a.at("experiment").get<std::string>();
    specifier_ = a.at("specifier").get<std::string>();
    beamEnergy_ = a.at("beam_energy").get<double>();
    bending_ = a.at("bending").get<std::string>();
    channel_ = a.at("channel").get<std::string>();
    defaultMomBin_ = a.value("momentum_bin", 0.05);
    if (!std::isfinite(beamEnergy_) || beamEnergy_ <= 0.0) {
        throw std::runtime_error("analysis.beam_energy must be positive");
    }

    const auto& detector = j.at("detector");
    ValidateKnownKeys(detector, {"adapter"}, "detector");
    detectorAdapter_ = detector.at("adapter").get<std::string>();
    if (detectorAdapter_ != "CLAS12") {
        throw std::runtime_error("Unsupported detector.adapter '" + detectorAdapter_ + "'");
    }

    const auto& outputs = j.at("outputs");
    ValidateKnownKeys(outputs, {"directory"}, "outputs");
    outputDirectory_ = outputs.at("directory").get<std::string>();
    if (outputDirectory_.empty()) {
        throw std::runtime_error("outputs.directory must not be empty");
    }

    const auto& particles = j.at("particles");
    if (!particles.is_array() || particles.empty()) {
        throw std::runtime_error("particles must be a non-empty array");
    }

    std::unordered_set<std::string> roles;
    std::unordered_set<std::string> names;
    for (const auto& p : particles) {
        ValidateKnownKeys(
            p,
            {"role", "name", "pid", "mass", "px", "py", "pz", "sector_branch",
             "detector", "sectors", "momentum", "phi"},
            "particles[]"
        );

        const std::string role = p.at("role").get<std::string>();
        const std::string name = p.at("name").get<std::string>();
        if (role.empty() || name.empty()) {
            throw std::runtime_error("particle role/name must not be empty");
        }
        if (!roles.insert(role).second) {
            throw std::runtime_error("Duplicate particle role '" + role + "'");
        }
        if (!names.insert(name).second) {
            throw std::runtime_error("Duplicate particle name '" + name + "'");
        }

        const auto& mom = p.at("momentum");
        ValidateKnownKeys(mom, {"low", "high", "bin_width"}, "particle.momentum");
        const double pLow = mom.at("low").get<double>();
        const double pHigh = mom.at("high").get<double>();
        const double pBin = mom.at("bin_width").get<double>();
        if (!std::isfinite(pLow) || !std::isfinite(pHigh) ||
            !std::isfinite(pBin) || pHigh <= pLow || pBin <= 0.0) {
            throw std::runtime_error("particle.momentum must have high > low and bin_width > 0");
        }

        const auto& phi = p.at("phi");
        ValidateKnownKeys(phi, {"enabled", "handling"}, "particle.phi");
        (void)PhiHandlingFromString(phi.at("handling").get<std::string>());
    }

    const auto& observables = j.at("observables");
    if (!observables.is_array() || observables.empty()) {
        throw std::runtime_error("observables must be a non-empty array");
    }

    observables_.clear();
    std::unordered_set<std::string> observableOutputs;
    for (const auto& o : observables) {
        ValidateKnownKeys(
            o,
            {"name", "type", "output", "detected_roles", "roles",
             "missing_pid", "missing_mass"},
            "observables[]"
        );

        ObservableConfig obs;
        obs.name = o.at("name").get<std::string>();
        obs.type = o.at("type").get<std::string>();
        obs.output = o.value("output", obs.name);
        if (obs.name.empty() || obs.type.empty() || obs.output.empty()) {
            throw std::runtime_error("observable name/type/output must not be empty");
        }
        if (observableOutputs.find(obs.output) != observableOutputs.end()) {
            throw std::runtime_error("Duplicate observable output '" + obs.output + "'");
        }

        if (HasKey(o, "detected_roles")) {
            obs.detectedRoles = StringVectorFromJson(
                o.at("detected_roles"), "observable.detected_roles");
            for (const auto& role : obs.detectedRoles) {
                if (roles.find(role) == roles.end()) {
                    throw std::runtime_error("observable references unknown role '" + role + "'");
                }
            }
        }

        if (HasKey(o, "roles")) {
            for (const auto& roleItem : o.at("roles").items()) {
                const std::string role = roleItem.value().get<std::string>();
                if (obs.type == "difference") {
                    if (observableOutputs.find(role) == observableOutputs.end()) {
                        throw std::runtime_error(
                            "difference observable references unknown prior output '" +
                            role + "'"
                        );
                    }
                } else if (roles.find(role) == roles.end()) {
                    throw std::runtime_error("observable references unknown role '" + role + "'");
                }
                obs.roles[roleItem.key()] = role;
            }
        }

        obs.missingPid = o.value("missing_pid", 0);
        obs.missingMass = o.value("missing_mass", 0.0);
        observables_.push_back(obs);
        observableOutputs.insert(obs.output);
    }

    const auto& histograms = j.at("histograms");
    if (!histograms.is_array() || histograms.empty()) {
        throw std::runtime_error("histograms must be a non-empty array");
    }

    histograms_.clear();
    std::unordered_set<std::string> histNames;
    for (const auto& h : histograms) {
        ValidateKnownKeys(
            h,
            {"name", "particle_role", "x", "y", "x_axis", "y_axis",
             "split_sector", "split_phi"},
            "histograms[]"
        );

        HistogramConfig hist;
        hist.name = h.at("name").get<std::string>();
        hist.particleRole = h.at("particle_role").get<std::string>();
        hist.x = h.at("x").get<std::string>();
        hist.y = h.at("y").get<std::string>();
        hist.xAxis = ParseAxis(h.at("x_axis"), "histogram.x_axis");
        hist.yAxis = ParseAxis(h.at("y_axis"), "histogram.y_axis");
        hist.splitSector = h.value("split_sector", true);
        hist.splitPhi = h.value("split_phi", true);

        if (hist.name.empty() || hist.x.empty() || hist.y.empty()) {
            throw std::runtime_error("histogram name/x/y must not be empty");
        }
        if (!histNames.insert(hist.name).second) {
            throw std::runtime_error("Duplicate histogram name '" + hist.name + "'");
        }
        if (roles.find(hist.particleRole) == roles.end()) {
            throw std::runtime_error("histogram references unknown particle_role '" +
                                     hist.particleRole + "'");
        }
        if (hist.x != "p") {
            throw std::runtime_error("Only histogram x='p' is supported in this refactor");
        }
        if (observableOutputs.find(hist.y) == observableOutputs.end()) {
            throw std::runtime_error("histogram y references unknown observable output '" +
                                     hist.y + "'");
        }

        histograms_.push_back(hist);
    }

    mmLow_ = 0.0;
    mmHigh_ = 0.0;
    mmBin_ = 0.0;
    dpLow_ = 0.0;
    dpHigh_ = 0.0;
    dpBin_ = 0.0;
}

json MomCorrConfig::GetStoredJSON() const { return storedJSON_; }

int MomCorrConfig::GetSchemaVersion() const { return schemaVersion_; }
std::string MomCorrConfig::GetTreeName() const { return treeName_; }
std::string MomCorrConfig::GetOutputDirectory() const { return outputDirectory_; }
std::string MomCorrConfig::GetDetectorAdapter() const { return detectorAdapter_; }

double MomCorrConfig::GetBeamEnergy() const { return beamEnergy_; }

double MomCorrConfig::GetMissingMassLow() const { return mmLow_; }
double MomCorrConfig::GetMissingMassHigh() const { return mmHigh_; }
double MomCorrConfig::GetMissingMassBinWidth() const { return mmBin_; }

double MomCorrConfig::GetDpLow() const { return dpLow_; }
double MomCorrConfig::GetDpHigh() const { return dpHigh_; }
double MomCorrConfig::GetDpBinWidth() const { return dpBin_; }

double MomCorrConfig::GetDefaultMomentumBin() const { return defaultMomBin_; }

std::string MomCorrConfig::GetSpecifier() const { return specifier_; }
std::string MomCorrConfig::GetBending() const { return bending_; }
std::string MomCorrConfig::GetChannel() const { return channel_; }
std::string MomCorrConfig::GetExperiment() const { return experiment_; }

std::vector<MomCorrParticle> MomCorrConfig::BuildParticles() const {
    std::vector<MomCorrParticle> createdParticles;

    const auto& particles = storedJSON_.at("particles");

    for (const auto& p : particles) {
        PhiHandling phiMode =
            PhiHandlingFromString(p.at("phi").at("handling").get<std::string>());

        MomCorrParticle particle(
            p.at("name").get<std::string>(),
            p.at("pid").get<int>(),
            p.at("mass").get<double>(),
            p.at("px").get<std::string>(),
            p.at("py").get<std::string>(),
            p.at("pz").get<std::string>(),
            p.at("sector_branch").get<std::string>(),
            p.at("detector").get<std::string>(),
            p.at("sectors").get<std::vector<int>>(),
            p.at("momentum").at("low").get<double>(),
            p.at("momentum").at("high").get<double>(),
            p.at("momentum").at("bin_width").get<double>(),
            phiMode,
            p.at("phi").at("enabled").get<bool>()
        );
        particle.SetRole(p.at("role").get<std::string>());

        createdParticles.push_back(particle);
    }

    return createdParticles;
}

const MomCorrParticle& MomCorrConfig::ParticleByRole(
    const std::vector<MomCorrParticle>& particles,
    const std::string& role) const {
    for (const auto& particle : particles) {
        if (particle.GetRole() == role) {
            return particle;
        }
    }

    throw std::runtime_error("No particle configured for role '" + role + "'");
}

const std::vector<ObservableConfig>& MomCorrConfig::GetObservables() const {
    return observables_;
}

const std::vector<HistogramConfig>& MomCorrConfig::GetHistograms() const {
    return histograms_;
}
