#include "MomCorrConfig.h"
#include <string>
#include <vector>
#include <fstream>
#include <MomCorrParticle.h>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

//Constructor
MomCorrConfig::MomCorrConfig(const std::string& jsonFile) {
	LoadJSON(jsonFile);
}	

//Setter
void MomCorrConfig::LoadJSON(const std::string& filename)
{
    std::ifstream f(filename);
    
    if (!f.is_open()) {
        throw std::runtime_error("MomCorrConfig: cannot open config file " + filename);
    }
    
    json j;
    f >> j;

    storedJSON_ = j;

    experiment_ = j.at("analysis").at("experiment").get<std::string>();
    specifier_ = j.at("analysis").at("specifier").get<std::string>();
    beamEnergy_ = j.at("analysis").at("beam_energy").get<double>();
    bending_ = j.at("analysis").at("bending").get<std::string>();
    channel_ = j.at("analysis").at("channel").get<std::string>();

    mmLow_ = j.at("analysis").at("missing_mass").at("low").get<double>();
    mmHigh_ = j.at("analysis").at("missing_mass").at("high").get<double>();
    mmBin_ = j.at("analysis").at("missing_mass").at("bin_width").get<double>();

    dpLow_ = j.at("analysis").at("dp").at("dp_low").get<double>();
    dpHigh_ = j.at("analysis").at("dp").at("dp_high").get<double>();
    dpBin_ = j.at("analysis").at("dp").at("dp_bin_width").get<double>();

    defaultMomBin_ = j.at("defaults").at("momentum_bin").get<double>();
}

//Getters

json MomCorrConfig::GetStoredJSON() const { return storedJSON_;}

double MomCorrConfig::GetBeamEnergy() const { return beamEnergy_;}

double MomCorrConfig::GetMissingMassLow() const { return mmLow_;}
double MomCorrConfig::GetMissingMassHigh() const { return mmHigh_;} 
double MomCorrConfig::GetMissingMassBinWidth() const { return mmBin_;}

double MomCorrConfig::GetDpLow() const { return dpLow_;}
double MomCorrConfig::GetDpHigh() const { return dpHigh_;}
double MomCorrConfig::GetDpBinWidth() const { return dpBin_;}

double MomCorrConfig::GetDefaultMomentumBin() const { return defaultMomBin_;}

std::string MomCorrConfig::GetSpecifier() const { return specifier_;}
std::string MomCorrConfig::GetBending() const { return bending_;}
std::string MomCorrConfig::GetChannel() const { return channel_;}
std::string MomCorrConfig::GetExperiment() const { return experiment_;}

std::vector<MomCorrParticle> MomCorrConfig::BuildParticles() const {

	std::vector<MomCorrParticle> created_particles;
	const auto& particles = storedJSON_.at("particles");

	for(const auto& p: particles){
		MomCorrParticle new_particle(p.at("name"), p.at("mass"), p.at("px"), p.at("py"), p.at("pz"), p.at("sector_branch"), p.at("detector"), p.at("sectors"), p.at("momentum").at("min"), p.at("momentum").at("max"), p.at("momentum").at("bin_width"), p.at("phi").at("handling"), p.at("phi").at("enabled"));
		created_particles.push_back(new_particle);
	}
	return created_particles;
}
