#include "MomCorrConfig.h"

#include <fstream>
#include <stdexcept>

namespace {

PhiHandling PhiHandlingFromString(const std::string& input) {
    if (input == "None") return PhiHandling::None;
    if (input == "CLAS12_CD_Standard") return PhiHandling::CLAS12_CD_Standard;
    if (input == "CLAS12_FD_El_Standard") return PhiHandling::CLAS12_FD_El_Standard;

    return PhiHandling::None;
}

}

MomCorrConfig::MomCorrConfig(const std::string& jsonFile) {
    LoadJSON(jsonFile);
}

void MomCorrConfig::LoadJSON(const std::string& filename) {
    std::ifstream f(filename);

    if (!f.is_open()) {
        throw std::runtime_error("MomCorrConfig: cannot open config file " + filename);
    }

    json j;
    f >> j;

    storedJSON_ = j;

    const auto& a = j.at("analysis");

    experiment_ = a.at("experiment").get<std::string>();
    specifier_ = a.at("specifier").get<std::string>();
    beamEnergy_ = a.at("beam_energy").get<double>();
    bending_ = a.at("bending").get<std::string>();
    channel_ = a.at("channel").get<std::string>();

    mmLow_ = a.at("missing_mass").at("low").get<double>();
    mmHigh_ = a.at("missing_mass").at("high").get<double>();
    mmBin_ = a.at("missing_mass").at("bin_width").get<double>();

    dpLow_ = a.at("dp").at("low").get<double>();
    dpHigh_ = a.at("dp").at("high").get<double>();
    dpBin_ = a.at("dp").at("bin_width").get<double>();

    defaultMomBin_ = a.at("momentum_bin").get<double>();
}

json MomCorrConfig::GetStoredJSON() const { return storedJSON_; }

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
    std::vector<MomCorrParticle> created_particles;

    const auto& particles = storedJSON_.at("particles");

    for (const auto& p : particles) {
        PhiHandling phiMode =
            PhiHandlingFromString(p.at("phi").at("handling").get<std::string>());

        MomCorrParticle new_particle(
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

        created_particles.push_back(new_particle);
    }

    return created_particles;
}
