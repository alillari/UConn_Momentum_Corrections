#ifndef MOMCORRCONFIG_H
#define MOMCORRCONFIG_H

#include <string>
#include <nlohmann/json.hpp>

using json = nlohmann::json;

class MomCorrConfig {

public: 
	explicit MomCorrConfig(const std::string& jsonFile);

    // Dataset info
    std::string GetSpecifier() const;
    double GetBeamEnergy() const;
    std::string GetBending() const;
    std::string GetChannelName() const;

    // Missing mass
    double GetMissingMassLow() const;
    double GetMissingMassHigh() const;
    double GetMissingMassBinWidth() const;

    // Momentum correction
    double GetDpLow() const;
    double GetDpHigh() const;
    double GetDpBinWidth() const;

    double GetDefaultMomentumBin() const;

private:

    void LoadJSON(const std::string& filename);

    // Stored values
    std::string specifier_;
    std::string bending_;
    std::string channel_;

    double beamEnergy_;
    double mmLow_, mmHigh_, mmBin_;
    double dpLow_, dpHigh_, dpBin_;
    double defaultMomBin_;

};

#endif //MOMCORRCONFIG_H
