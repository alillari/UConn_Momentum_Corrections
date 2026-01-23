#ifndef MOMCORRCONFIG_H
#define MOMCORRCONFIG_H

#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include <MomCorrParticle.h>

using json = nlohmann::json;

class MomCorrConfig {

public: 
    explicit MomCorrConfig(const std::string& jsonFile);

    // Dataset info
    std::string GetSpecifier() const;
    double GetBeamEnergy() const;
    std::string GetBending() const;
    std::string GetChannel() const;
    std::string GetExperiment() const;

    // Missing mass
    double GetMissingMassLow() const;
    double GetMissingMassHigh() const;
    double GetMissingMassBinWidth() const;

    // Momentum correction
    double GetDpLow() const;
    double GetDpHigh() const;
    double GetDpBinWidth() const;

    double GetDefaultMomentumBin() const;

    json GetStoredJSON() const;

    std::vector<MomCorrParticle> BuildParticles() const;

private:

    void LoadJSON(const std::string& filename);

    // Stored values
    json storedJSON_;
    
    std::string specifier_;
    std::string bending_;
    std::string channel_;
    std::string experiment_;

    double beamEnergy_;
    double mmLow_, mmHigh_, mmBin_;
    double dpLow_, dpHigh_, dpBin_;
    double defaultMomBin_;

};

#endif //MOMCORRCONFIG_H
