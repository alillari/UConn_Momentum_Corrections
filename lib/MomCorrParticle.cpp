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
                                 const int detector,
				 const std::vector<int>& sectors,
                                 double pMin, double pMax, double binWidth,
				 PhiHandling phiHandling,
				 bool usePhiBinning)
    : name_(name), mass_(mass), pxBranch_(pxBranch), pyBranch_(pyBranch), pzBranch_(pzBranch), 
      sectorBranch_(sectorBranch), detector_(detector), sectors_(sectors), pMin_(pMin), pMax_(pMax), binWidth_(binWidth), 
      phiHandling_(phiHandling), usePhiBinning_(usePhiBinning) {}

// Getters
std::string MomCorrParticle::GetName() const { return name_; }
double MomCorrParticle::GetMass() const { return mass_; }
std::string MomCorrParticle::GetPxBranch() const { return pxBranch_; }
std::string MomCorrParticle::GetPyBranch() const { return pyBranch_; }
std::string MomCorrParticle::GetPzBranch() const { return pzBranch_; }
std::string MomCorrParticle::GetSectorBranch() const { return name_ + "_sec"; }
const std::vector<int>& MomCorrParticle::GetSectors() const { return sectors_; }
int MomCorrParticle::GetDetector() const { return detector_; }
double MomCorrParticle::GetMomentumMin() const { return pMin_; }
double MomCorrParticle::GetMomentumMax() const { return pMax_; }
double MomCorrParticle::GetMomentumBinWidth() const { return binWidth_; }
int MomCorrParticle::GetBins() const { return static_cast<int>(std::round((pMax_ - pMin_) / binWidth_)); }
PhiHandling MomCorrParticle::GetPhiHandling() const { return phiHandling_; }
bool MomCorrParticle::IsPhiBinningEnabled() const { return usePhiBinning_; }

double MomCorrParticle::ComputeLocalPhi(double phi, double p, int sector) const
{
    switch (phiHandling_) {
        case PhiHandling::CLAS12_FD_Standard: {
            double local = phi - (sector - 1) * 60.0;
            return local - (30.0 / p);
        }
        case PhiHandling::None:
        default:
            return phi;
    }
}

std::unordered_map<int, std::string> MomCorrParticle::GetPhiBinningLabels() const
{
    if (!usePhiBinning_) {
        return {};
    }

    switch (phiHandling_) {

	//Inherited from Richard, 1 is less than -5 deg, 2 is between -5 and 5 deg, and 3 is larger than 5 deg
        case PhiHandling::CLAS12_FD_Standard:
            return {
                {1, "negative"},
                {2, "neutral"},
                {3, "positive"}
            };

        case PhiHandling::CLAS12_CD_Standard:
            return {
	        {0, "no phi bin"}
	    };

        case PhiHandling::None:
        default:
            return {
	        {0, "no phi bin"}
	    };
    }
}

int MomCorrParticle::PhiBin(double localPhi) const {
    switch (phiHandling_) {
    	case PhiHandling::CLAS12_FD_Standard:
	    if(localPhi > 5){ return 3;};
	    if(localPhi <= -5){ return 1;}
	    else{return 2;};
	case PhiHandling::CLAS12_CD_Standard:
	    return 0; 

	case PhiHandling::None:
	default:
	    return 0;  
    }
}

ROOT::RDF::RNode MomCorrParticle::AddBranches(ROOT::RDF::RNode df) const {
    std::string secInt = name_ + "_sec"; 
    df = df.Define(secInt, 
		   [this](float esec) {
		   	return static_cast<int>(esec);
		   }, {sectorBranch_});
	
    std::string pBranch = name_ + "_mag";
    std::string thetaBranch = name_ + "_theta";
    std::string phiBranch = name_ + "_phi";
    std::string localPhiSBranch = name_ + "_localphi";
    std::string phiBin = name_ + "_phiBin";

    // Define momentum magnitude
    df = df.Define(pBranch,
                   [this](float px, float py, float pz) {
                       return std::sqrt(px * px + py * py + pz * pz);
                   }, {pxBranch_, pyBranch_, pzBranch_});

    // Define theta (polar angle)
    df = df.Define(thetaBranch,
                   [this](float p, float pz) {
                       return (p > 0) ? std::acos(pz / p) : -999.0;
                   }, {pBranch, pzBranch_});

    // Define phi (azimuthal angle)
    df = df.Define(phiBranch,
                   [this](float px, float py) {
                       return (180.0 / M_PI)*std::atan2(py, px);
                   }, {pxBranch_, pyBranch_});

    df = df.Define(localPhiSBranch,
                   [this](double phi, float p, int sector) { return ComputeLocalPhi(phi, p, sector); },
                   {phiBranch, pBranch, secInt});

    df = df.Define(phiBin,
		    [this](double localPhiS) { return PhiBin(localPhiS);},
		    {localPhiSBranch});

    return df;
}
