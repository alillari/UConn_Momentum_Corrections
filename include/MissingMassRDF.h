#ifndef MISSING_MASS_RDF_H
#define MISSING_MASS_RDF_H

#include <ROOT/RDataFrame.hxx>
#include <Math/Vector4D.h>
#include "MomCorrParticle.h"
#include "Kinematics.h"

std::vector<std::string> ExtractMomentumBranches(const std::vector<MomCorrParticle>& particles) {
    std::vector<std::string> branches;
    for (const auto& particle : particles) {
        branches.push_back(particle.GetPxBranch());
        branches.push_back(particle.GetPyBranch());
        branches.push_back(particle.GetPzBranch());
    }
    return branches;
}

template <typename RDF>
RDF DefineDetectedP4s(RDF& df, const std::vector<MomCorrParticle>& particle_list) {
    std::vector<std::string> branches = ExtractMomentumBranches(particle_list);

    if (branches.size() == 3) {
        return df.Define(
            "detected_p4s",
            [=](float px1, float py1, float pz1) {
                std::vector<ROOT::Math::PxPyPzMVector> p4s;
                p4s.emplace_back(px1, py1, pz1, particle_list[0].GetMass());
                return p4s;
            },
            branches
        );
    }
    else if (branches.size() == 6) {
        return df.Define(
            "detected_p4s",
            [=](float px1, float py1, float pz1,
                float px2, float py2, float pz2) {
                std::vector<ROOT::Math::PxPyPzMVector> p4s;
                p4s.emplace_back(px1, py1, pz1, particle_list[0].GetMass());
                p4s.emplace_back(px2, py2, pz2, particle_list[1].GetMass());
                return p4s;
            },
            branches
        );
    }
    else if (branches.size() == 9) {
        return df.Define(
            "detected_p4s",
            [=](float px1, float py1, float pz1,
                float px2, float py2, float pz2,
                float px3, float py3, float pz3) {
                std::vector<ROOT::Math::PxPyPzMVector> p4s;
                p4s.emplace_back(px1, py1, pz1, particle_list[0].GetMass());
                p4s.emplace_back(px2, py2, pz2, particle_list[1].GetMass());
                p4s.emplace_back(px3, py3, pz3, particle_list[2].GetMass());
                return p4s;
            },
            branches
        );
    }
    else if (branches.size() == 12) {
        return df.Define(
            "detected_p4s",
            [=](float px1, float py1, float pz1,
                float px2, float py2, float pz2,
                float px3, float py3, float pz3,
                float px4, float py4, float pz4) {
                std::vector<ROOT::Math::PxPyPzMVector> p4s;
                p4s.emplace_back(px1, py1, pz1, particle_list[0].GetMass());
                p4s.emplace_back(px2, py2, pz2, particle_list[1].GetMass());
                p4s.emplace_back(px3, py3, pz3, particle_list[2].GetMass());
                p4s.emplace_back(px4, py4, pz4, particle_list[3].GetMass());
                return p4s;
            },
            branches
        );
    }
    else {
        throw std::runtime_error(
            "Unsupported number of particles in DefineDetectedP4s"
        );
    }
}

template <typename RDF>
RDF DefineMissingMass(RDF& df, const std::vector<MomCorrParticle>& particles, const double beam_energy) {
	df = DefineDetectedP4s(df, particles);

	df = df.Define("missing_mass", [beam_energy](const std::vector<ROOT::Math::PxPyPzMVector>& detected_p4s) {
        return ComputeMissingMass(detected_p4s, beam_energy);
	}
	,{"detected_p4s"});

	return df;
}

#endif
