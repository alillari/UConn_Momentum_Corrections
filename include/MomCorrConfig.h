#ifndef MOMCORRCONFIG_H
#define MOMCORRCONFIG_H

#include "MomCorrParticle.h"

#include <nlohmann/json.hpp>

#include <string>
#include <unordered_map>
#include <vector>

using json = nlohmann::json;

struct AxisConfig {
    int bins = 0;
    double low = 0.0;
    double high = 0.0;
};

struct ObservableConfig {
    std::string name;
    std::string type;
    std::string output;
    std::vector<std::string> detectedRoles;
    std::unordered_map<std::string, std::string> roles;
    int missingPid = 0;
    double missingMass = 0.0;
};

struct HistogramConfig {
    std::string name;
    std::string particleRole;
    std::string x;
    std::string y;
    AxisConfig xAxis;
    AxisConfig yAxis;
    bool splitSector = true;
    bool splitPhi = true;
};

class MomCorrConfig {
public:
    explicit MomCorrConfig(const std::string& jsonFile);

    int GetSchemaVersion() const;
    std::string GetTreeName() const;
    std::string GetOutputDirectory() const;
    std::string GetDetectorAdapter() const;

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
    const MomCorrParticle& ParticleByRole(const std::vector<MomCorrParticle>& particles,
                                          const std::string& role) const;
    const std::vector<ObservableConfig>& GetObservables() const;
    const std::vector<HistogramConfig>& GetHistograms() const;

private:
    void LoadJSON(const std::string& filename);
    void ValidateKnownKeys(const json& object,
                           const std::vector<std::string>& allowed,
                           const std::string& path) const;
    AxisConfig ParseAxis(const json& j, const std::string& path) const;
    void ValidateAxis(const AxisConfig& axis, const std::string& path) const;

    json storedJSON_;

    int schemaVersion_;
    std::string treeName_;
    std::string outputDirectory_;
    std::string detectorAdapter_;

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

    std::vector<ObservableConfig> observables_;
    std::vector<HistogramConfig> histograms_;
};

#endif
