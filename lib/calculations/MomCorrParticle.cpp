#include "MomCorrParticle.h"
#include <cmath>
#include <functional>
#include <vector>
#include <string>
#include <unordered_map>

#include <ROOT/RDataFrame.hxx>

// Constructor
MomCorrParticle::MomCorrParticle(const std::string& name,
	       			 const double mass,	
                                 const std::string& pxBranch, 
                                 const std::string& pyBranch, 
                                 const std::string& pzBranch, 
                                 const std::string& sectorBranch, 
                                 const std::vector<int>& sectors,
				 const int detector, 
                                 double pMin, double pMax, double binWidth,
				 bool usePhiBinning,
				 std::function<double(double, double, int)> phiShiftFunc,
				 std::function<int(double)> phiBinningFunc
				 std::unordered_map<int, std::string> phiBinningLabels)
    : name_(name), mass_(mass), pxBranch_(pxBranch), pyBranch_(pyBranch), pzBranch_(pzBranch), 
      sectorBranch_(sectorBranch), sectors_(sectors), detector_(detector), pMin_(pMin), pMax_(pMax), binWidth_(binWidth), 
      usePhiBinning_(usePhiBinning), phiShiftFunc_(phiShiftFunc), phiBinningFunc_(phiBinningFunc), phiBinningLabels_(phiBinningLabels) {}

// Getters
std::string MomCorrParticle::GetName() const { return name_; }
double MomCorrParticle::GetMass() const { return mass_; }
std::string MomCorrParticle::GetPxBranch() const { return pxBranch_; }
std::string MomCorrParticle::GetPyBranch() const { return pyBranch_; }
std::string MomCorrParticle::GetPzBranch() const { return pzBranch_; }
std::string MomCorrParticle::GetSectorBranch() const { return sectorBranch_; }
const std::vector<int>& MomCorrParticle::GetSectors() const { return sectors_; }
int MomCorrParticle::GetDetector() const { return detector_; }
double MomCorrParticle::GetMomentumMin() const { return pMin_; }
double MomCorrParticle::GetMomentumMax() const { return pMax_; }
double MomCorrParticle::GetBinWidth() const { return binWidth_; }
bool MomCorrParticle::IsPhiBinningEnabled() const { return usePhiBinning_; }
std::function<double(double, double, int)> MomCorrParticle::GetPhiShiftFunction() const { return phiShiftFunc_; }
std::function<int(double)> MomCorrParticle::GetPhiBinningFunction() const { return phiBinningFunc_; }
std::unordered_map<int, std::string> MomCorrParticle::GetPhiBinningLabels() const { return phiBinningLabels_; }

ROOT::RDF::RNode MomCorrParticle::AddBranches(ROOT::RDF::RNode df) const {
    std::string pBranch = name_ + "_mag";
    std::string thetaBranch = name_ + "_theta";
    std::string phiBranch = name_ + "_phi";
    std::string localPhiSBranch = name_ + "_localphi";
    std::string phiBin = name_ + "_phiBin";

    // Define momentum magnitude
    df = df.Define(pBranch,
                   [this](double px, double py, double pz) {
                       return std::sqrt(px * px + py * py + pz * pz);
                   }, {pxBranch_, pyBranch_, pzBranch_});

    // Define theta (polar angle)
    df = df.Define(thetaBranch,
                   [this](double p, double pz) {
                       return (p > 0) ? std::acos(pz / p) : -999.0;
                   }, {pBranch, pzBranch_});

    // Define phi (azimuthal angle)
    df = df.Define(phiBranch,
                   [this](double px, double py) {
                       return (180.0 / M_PI)*std::atan2(py, px);
                   }, {pxBranch_, pyBranch_});

    df = df.Define(localPhiSBranch,
                   [this](double phi, double p, int sector) { return phiShiftFunc_(phi, p, sector); },
                   {phiBranch, pBranch, sectorBranch_});

    df = df.Define(phiBin,
		    [this](double localPhiS) { return phiBinningFunc_(localPhiS);},
		    {localPhiSBranch});

    return df;
}
