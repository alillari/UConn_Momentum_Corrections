#include "MomCorrParticle.h"

#include <cmath>

MomCorrParticle::MomCorrParticle(
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
)
    : name_(name),
      pid_(pid),
      mass_(mass),
      pxBranch_(pxBranch),
      pyBranch_(pyBranch),
      pzBranch_(pzBranch),
      sectorBranch_(sectorBranch),
      detector_(detector),
      sectors_(sectors),
      pMin_(pMin),
      pMax_(pMax),
      binWidth_(binWidth),
      phiHandling_(phiHandling),
      usePhiBinning_(usePhiBinning) {}

std::string MomCorrParticle::GetName() const { return name_; }
int MomCorrParticle::GetPID() const { return pid_; }
double MomCorrParticle::GetMass() const { return mass_; }

std::string MomCorrParticle::GetInputPxBranch() const { return pxBranch_; }
std::string MomCorrParticle::GetInputPyBranch() const { return pyBranch_; }
std::string MomCorrParticle::GetInputPzBranch() const { return pzBranch_; }
std::string MomCorrParticle::GetInputSectorBranch() const { return sectorBranch_; }

std::vector<std::string> MomCorrParticle::GetInputPBranches() const {
    return {pxBranch_, pyBranch_, pzBranch_};
}

std::string MomCorrParticle::GetPxBranch() const { return GetInputPxBranch(); }
std::string MomCorrParticle::GetPyBranch() const { return GetInputPyBranch(); }
std::string MomCorrParticle::GetPzBranch() const { return GetInputPzBranch(); }
std::string MomCorrParticle::GetSectorBranch() const { return GetDerivedSectorBranch(); }
std::vector<std::string> MomCorrParticle::GetPBranches() const { return GetInputPBranches(); }

std::string MomCorrParticle::GetDerivedSectorBranch() const { return name_ + "_sec"; }
std::string MomCorrParticle::GetMagBranch() const { return name_ + "_mag"; }
std::string MomCorrParticle::GetThetaBranch() const { return name_ + "_theta"; }
std::string MomCorrParticle::GetPhiBranch() const { return name_ + "_phi"; }
std::string MomCorrParticle::GetLocalPhiBranch() const { return name_ + "_localphi"; }
std::string MomCorrParticle::GetPhiBinBranch() const { return name_ + "_phiBin"; }

std::string MomCorrParticle::GetDetector() const { return detector_; }

int MomCorrParticle::GetDetectorIDNum() const {
    if (detector_ == "FD") return 2;
    if (detector_ == "CD") return 3;
    return -1;
}

const std::vector<int>& MomCorrParticle::GetSectors() const { return sectors_; }

double MomCorrParticle::GetMomentumMin() const { return pMin_; }
double MomCorrParticle::GetMomentumMax() const { return pMax_; }
double MomCorrParticle::GetMomentumBinWidth() const { return binWidth_; }

int MomCorrParticle::GetBins() const {
    return static_cast<int>(std::round((pMax_ - pMin_) / binWidth_));
}

PhiHandling MomCorrParticle::GetPhiHandling() const { return phiHandling_; }
bool MomCorrParticle::IsPhiBinningEnabled() const { return usePhiBinning_; }

double MomCorrParticle::ComputeLocalPhi(double phi, double p, int sector) const {
    switch (phiHandling_) {
        case PhiHandling::CLAS12_FD_El_Standard: {
            if (((sector == 3 || sector == 4) && phi < 0.0) ||
                (sector > 4 && phi < 90.0)) {
                phi += 360.0;
            }

            double local = phi - (sector - 1) * 60.0;
            local -= 30.0 / p;
            return local;
        }

        case PhiHandling::CLAS12_CD_Standard:
        case PhiHandling::None:
        default:
            return phi;
    }
}

std::unordered_map<int, std::string> MomCorrParticle::GetPhiBinningLabels() const {
    if (!usePhiBinning_) {
        return {{0, "no phi bin"}};
    }

    switch (phiHandling_) {
        case PhiHandling::CLAS12_FD_El_Standard:
            return {
                {1, "negative"},
                {2, "neutral"},
                {3, "positive"}
            };

        case PhiHandling::CLAS12_CD_Standard:
        case PhiHandling::None:
        default:
            return {{0, "no phi bin"}};
    }
}

int MomCorrParticle::PhiBin(double localPhi) const {
    if (!usePhiBinning_) return 0;

    switch (phiHandling_) {
        case PhiHandling::CLAS12_FD_El_Standard:
            if (localPhi > 5.0) return 3;
            if (localPhi <= -5.0) return 1;
            return 2;

        case PhiHandling::CLAS12_CD_Standard:
        case PhiHandling::None:
        default:
            return 0;
    }
}

ROOT::RDF::RNode MomCorrParticle::AddBranches(ROOT::RDF::RNode df) const {
    df = df.Define(
        GetDerivedSectorBranch(),
        [](float sec) {
            return static_cast<int>(sec);
        },
        {GetInputSectorBranch()}
    );

    df = df.Define(
        GetMagBranch(),
        [](float px, float py, float pz) -> double {
            return std::sqrt(px * px + py * py + pz * pz);
        },
        {GetInputPxBranch(), GetInputPyBranch(), GetInputPzBranch()}
    );

    df = df.Define(
        GetThetaBranch(),
        [](double p, float pz) -> double {
            return (p > 0.0) ? std::acos(pz / p) : -999.0;
        },
        {GetMagBranch(), GetInputPzBranch()}
    );

    df = df.Define(
        GetPhiBranch(),
        [](float px, float py) -> double {
            return (180.0 / M_PI) * std::atan2(py, px);
        },
        {GetInputPxBranch(), GetInputPyBranch()}
    );

    df = df.Define(
        GetLocalPhiBranch(),
        [this](double phi, double p, int sector) {
            return ComputeLocalPhi(phi, p, sector);
        },
        {GetPhiBranch(), GetMagBranch(), GetDerivedSectorBranch()}
    );

    df = df.Define(
        GetPhiBinBranch(),
        [this](double localPhi) {
            return PhiBin(localPhi);
        },
        {GetLocalPhiBranch()}
    );

    return df;
}
