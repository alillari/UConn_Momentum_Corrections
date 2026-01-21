#include "MomCorrConfig.h"
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <unordered_map>

#include <ROOT/RDataFrame.hxx>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

//Constructor
MomCorrConfig::MomCorrConfig(const std::string& jsonFile) :    

//Getters
double MomCorrConfig::getBeamE();

double MomCorrConfig::getmmLow() const { return }
double MomCorrConfig:: 
double MomCorrConfig::

double MomCorrConfig::
double MomCorrConfig::
double MomCorrConfig::

double MomCorrConfig::

std::string MomCorrConfig::
std::string MomCorrConfig::
std::string MomCorrConfig::

//Functions

void MomCorrConfig::LoadJSON(const std::string& filename)
{
    std::ifstream f(filename);
    json j;
    f >> j;

    specifier_ = j["dataset"]["specifier"];
    beamEnergy_ = j["dataset"]["beam_energy"];
    bending_ = j["dataset"]["bending"];
    channel_ = j["dataset"]["channel"];

    mmLow_ = j["missing_mass"]["low"];
    mmHigh_ = j["missing_mass"]["high"];
    mmBin_ = j["missing_mass"]["bin_width"];

    dpLow_ = j["momentum_correction"]["dp_low"];
    dpHigh_ = j["momentum_correction"]["dp_high"];
    dpBin_ = j["momentum_correction"]["dp_bin_width"];

    defaultMomBin_ = j["defaults"]["momentum_bin"];
}
