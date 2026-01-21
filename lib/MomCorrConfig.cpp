#include "MomCorrConfig.h"
#include <string>
#include <fstream>

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

    specifier_ = j.at("dataset").at("specifier").get<std::string>();
    beamEnergy_ = j.at("dataset").at("beam_energy").get<double>();
    bending_ = j.at("dataset").at("bending").get<std::string>();
    channel_ = j.at("dataset").at("channel").get<std::string>();

    mmLow_ = j.at("missing_mass").at("low").get<double>();
    mmHigh_ = j.at("missing_mass").at("high").get<double>();
    mmBin_ = j.at("missing_mass").at("bin_width").get<double>();

    dpLow_ = j.at("momentum_correction").at("dp_low").get<double>();
    dpHigh_ = j.at("momentum_correction").at("dp_high").get<double>();
    dpBin_ = j.at("momentum_correction").at("dp_bin_width").get<double>();

    defaultMomBin_ = j.at("defaults").at("momentum_bin").get<double>();
}

//Getters
double MomCorrConfig::GetBeamEnergy() const {return beamEnergy_;}

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
