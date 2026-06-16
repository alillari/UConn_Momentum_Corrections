#ifndef MOMCORRPARTICLE_H
#define MOMCORRPARTICLE_H

#include <ROOT/RDataFrame.hxx>

#include <string>
#include <vector>
#include <unordered_map>

enum class PhiHandling {
    None,
    CLAS12_CD_Standard,
    CLAS12_FD_El_Standard
};

class MomCorrParticle {
public:
    MomCorrParticle(
        const std::string& name,
        int pid,
        double mass,
        const std::string& pxBranch,
        const std::string& pyBranch,
        const std::string& pzBranch,
        const std::string& sectorBranch,
        const std::string& detector,
        const std::vector<int>& sectors,
        double pMin,
        double pMax,
        double binWidth,
        PhiHandling phiHandling,
        bool usePhiBinning
    );

    // Physics identity
    std::string GetName() const;
    int GetPID() const;
    double GetMass() const;

    // Input tree branches from JSON
    std::string GetInputPxBranch() const;
    std::string GetInputPyBranch() const;
    std::string GetInputPzBranch() const;
    std::string GetInputSectorBranch() const;
    std::vector<std::string> GetInputPBranches() const;

    // Backward-compatible aliases
    std::string GetPxBranch() const;
    std::string GetPyBranch() const;
    std::string GetPzBranch() const;
    std::string GetSectorBranch() const;
    std::vector<std::string> GetPBranches() const;

    // Derived RDataFrame branches
    std::string GetDerivedSectorBranch() const;
    std::string GetMagBranch() const;
    std::string GetThetaBranch() const;
    std::string GetPhiBranch() const;
    std::string GetLocalPhiBranch() const;
    std::string GetPhiBinBranch() const;

    // Detector/config
    std::string GetDetector() const;
    int GetDetectorIDNum() const;
    const std::vector<int>& GetSectors() const;

    double GetMomentumMin() const;
    double GetMomentumMax() const;
    double GetMomentumBinWidth() const;
    int GetBins() const;

    PhiHandling GetPhiHandling() const;
    bool IsPhiBinningEnabled() const;

    double ComputeLocalPhi(double phi, double p, int sector) const;
    int PhiBin(double localPhi) const;
    std::unordered_map<int, std::string> GetPhiBinningLabels() const;

    ROOT::RDF::RNode AddBranches(ROOT::RDF::RNode df) const;

private:
    std::string name_;
    int pid_;
    double mass_;

    std::string pxBranch_;
    std::string pyBranch_;
    std::string pzBranch_;
    std::string sectorBranch_;

    std::string detector_;
    std::vector<int> sectors_;

    double pMin_;
    double pMax_;
    double binWidth_;

    PhiHandling phiHandling_;
    bool usePhiBinning_;
};

#endif
