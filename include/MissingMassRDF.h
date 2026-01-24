#ifndef MISSING_MASS_RDF_H
#define MISSING_MASS_RDF_H

#include <ROOT/RDataFrame.hxx>
#include "MomCorrParticle.h"
#include "Kinematics.h"

inline std::vector<std::string>
ExtractMomentumBranches(const std::vector<MomCorrParticle>& particles) {
    std::vector<std::string> branches;
    for (const auto& p : particles) {
        branches.push_back(p.GetPxBranch());
        branches.push_back(p.GetPyBranch());
        branches.push_back(p.GetPzBranch());
    }
    return branches;
}

#endif
