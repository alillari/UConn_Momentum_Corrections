#ifndef MOMCORRPARTICLE_H
#define MOMCORRPARTICLE_H

#include <string>
#include <vector>
#include <unordered_map>

#include <ROOT/RDataFrame.hxx>

//Default PhiHandlings available
enum class PhiHandling {
    None,
    CLAS12_CD_Standard,
    CLAS12_FD_El_Standard
};

class MomCorrParticle {
public:
    // Constructor
    MomCorrParticle(const std::string& name,
                    const double mass,
                    const std::string& pxBranch,
                    const std::string& pyBranch,
                    const std::string& pzBranch,
                    const std::string& sectorBranch,
                    const int detector,
                    const std::vector<int>& sectors,
                    double pMin, double pMax, double binWidth,
                    PhiHandling phiHandling,
                    bool usePhiBinning);

    // Getters
    std::string GetName() const;
    double GetMass() const;
    std::string GetPxBranch() const;
    std::string GetPyBranch() const;
    std::string GetPzBranch() const;
    std::vector<std::string> GetPBranches() const;
    std::string GetSectorBranch() const;
    int GetDetector() const;
    const std::vector<int>& GetSectors() const;

    double GetMomentumMin() const;
    double GetMomentumMax() const;
    double GetMomentumBinWidth() const;
    int GetBins() const;

    PhiHandling GetPhiHandling() const;
    bool IsPhiBinningEnabled() const;

    // Phi-related helpers
    double ComputeLocalPhi(double phi, double p, int sector) const;
    int PhiBin(double localPhi) const;
    std::unordered_map<int, std::string> GetPhiBinningLabels() const;

    // RDataFrame interface
    ROOT::RDF::RNode AddBranches(ROOT::RDF::RNode df) const;

private:
    // Identity
    std::string name_;
    double mass_;

    // Input branches
    std::string pxBranch_;
    std::string pyBranch_;
    std::string pzBranch_;
    std::string sectorBranch_;

    // Configuration
    int detector_;
    std::vector<int> sectors_;
    double pMin_;
    double pMax_;
    double binWidth_;

    // Phi handling
    PhiHandling phiHandling_;
    bool usePhiBinning_;
};

#endif // MOMCORRPARTICLE_H
