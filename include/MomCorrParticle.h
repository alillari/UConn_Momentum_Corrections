#ifndef MOMCORRPARTICLE_H
#define MOMCORRPARTICLE_H

#include <string>
#include <vector>
#include <functional>
#include <unordered_map>

#include <ROOT/RDataFrame.hxx>

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
		    bool usePhiBinning,
		    std::function<double(double, double, int)> phiShiftFunc,
		    std::function<int(double)> phiBinningFunction,
		    std::unordered_map<int, std::string> phiBinningLabels
		    );

    // Getters
    std::string GetName() const;
    double GetMass() const;
    std::string GetPxBranch() const;
    std::string GetPyBranch() const;
    std::string GetPzBranch() const;
    std::string GetSectorBranch() const;
    int GetDetector() const;
    const std::vector<int>& GetSectors() const;
    double GetMomentumMin() const;
    double GetMomentumMax() const;
    double GetBinWidth() const;
    bool IsPhiBinningEnabled() const;
    std::function<double(double, double, int)> GetPhiShiftFunction() const;
    std::function<int(double)> GetPhiBinningFunction() const;
    std::unordered_map<int, std::string> GetPhiBinningLabels() const;

    ROOT::RDF::RNode AddBranches(ROOT::RDF::RNode df) const;

private:
    std::string name_;
    double mass_;
    std::string pxBranch_, pyBranch_, pzBranch_, sectorBranch_;
    int detector_;
    std::vector<int> sectors_;
    double pMin_, pMax_, binWidth_;
    bool usePhiBinning_;
    std::function<double(double, double, int)> phiShiftFunc_;
    std::function<int(double)> phiBinningFunction_;
    std::unordered_map<int, std::string> phiBinningLabels_;
};

#endif // MOMCORRPARTICLE_H

