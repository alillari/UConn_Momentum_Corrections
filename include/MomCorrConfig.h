#ifndef MOMCORRCONFIG_H
#define MOMCORRCONFIG_H

#include "MomCorrParticle.h"

#include <nlohmann/json.hpp>

#include <string>
#include <vector>

using json = nlohmann::json;

class MomCorrConfig {
public:
    explicit MomCorrConfig(const std::string& jsonFile);

    std::string GetSpecifier() const;
    double GetBeamEnergy() const;
    std::string GetBending() const;
    std::string GetChannel() const;
    std::string GetExperiment() const;

    double GetMissingMassLow() const;
    double GetMissingMassHigh() const;
    double GetMissingMassBinWidth() const;

    double GetDpLow() const;
    double GetDpHigh() const;
    double GetDpBinWidth() const;

    double GetDefaultMomentumBin() const;

    json GetStoredJSON() const;

    std::vector<MomCorrParticle> BuildParticles() const;

private:
    void LoadJSON(const std::string& filename);

    json storedJSON_;

    std::string specifier_;
    std::string bending_;
    std::string channel_;
    std::string experiment_;

    double beamEnergy_;

    double mmLow_;
    double mmHigh_;
    double mmBin_;

    double dpLow_;
    double dpHigh_;
    double dpBin_;

    double defaultMomBin_;
};

#endif
