#ifndef COMPUTE_MISSING_MASS_HELPER_H
#define COMPUTE_MISSING_MASS_HELPER_H

#include <vector>
#include <string>
#include <tuple>
#include <ROOT/RVec.hxx>
#include <ROOT/RDataFrame.hxx>
#include <Math/Vector4D.h>
#include "MomCorrParticle.h"

inline std::vector<std::string> ExtractMomentumBranches(const std::vector<MomCorrParticle>& particles) {
    std::vector<std::string> branches;
    for (const auto& particle : particles) {
        branches.push_back(particle.GetPxBranch());
        branches.push_back(particle.GetPyBranch());
        branches.push_back(particle.GetPzBranch());
    }
    return branches;
}

// Function to compute the missing mass from a vector of particle momenta
inline double compute_missing_mass(const double beam_energy, const std::vector<std::tuple<double, double, double, double>>& particles) {
    const double El_mass = 0.000511;        // Electron mass in GeV
    const double Pro_mass = 0.938272;       // Proton mass in GeV

    ROOT::Math::PxPyPzMVector beam_vector(0, 0, std::sqrt(beam_energy * beam_energy - El_mass * El_mass), El_mass);
    ROOT::Math::PxPyPzMVector target_vector(0, 0, 0, Pro_mass);

    ROOT::Math::PxPyPzMVector detected_sum(0, 0, 0, 0);

    for (const auto& [px, py, pz, mass] : particles) {
        detected_sum += ROOT::Math::PxPyPzMVector(px, py, pz, mass);
    }

    auto missing_vector = (beam_vector + target_vector) - detected_sum;
    return missing_vector.M();
}

// Helper function to handle input momenta and call compute_missing_mass
inline double compute_missing_mass_helper(const std::vector<ROOT::RVec<double>>& momenta, 
                                          const std::vector<MomCorrParticle>& particle_list) {
    std::vector<std::tuple<double, double, double, double>> particle_momenta;

    size_t n_particles = particle_list.size();
    for (size_t i = 0; i < n_particles; ++i) {
        double mass = particle_list[i].GetMass();
        particle_momenta.emplace_back(
            momenta[3 * i][0],      // px
            momenta[3 * i + 1][0],  // py
            momenta[3 * i + 2][0],  // pz
            mass
        );
    }

    return compute_missing_mass(particle_momenta);
}

// Main function to define "missing_mass" in the RDF
template <typename T>
inline T DefineMissingMass(T& df, const std::vector<MomCorrParticle>& particle_list) {
    return df.Define("missing_mass",
        [=](const std::vector<ROOT::RVec<double>>& momenta) {
            return compute_missing_mass_helper(momenta, particle_list);
        },
        ExtractMomentumBranches(particle_list)
    );
}

#endif  // COMPUTE_MISSING_MASS_HELPER_H

