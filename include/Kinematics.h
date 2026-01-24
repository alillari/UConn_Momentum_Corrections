#ifndef KINEMATICS_H
#define KINEMATICS_H

#include <Math/Vector4D.h>
#include <vector>

// Units: GeV
constexpr double ELECTRON_MASS = 0.000511;
constexpr double MUON_MASS = .10;
constexpr double PROTON_MASS  = 0.938272;

inline double ComputeMissingMass(const std::vector<ROOT::Math::PxPyPzMVector>& detected, double beam_energy) {
    ROOT::Math::PxPyPzMVector beam(0, 0, std::sqrt(beam_energy*beam_energy - ELECTRON_MASS*ELECTRON_MASS), ELECTRON_MASS);
    ROOT::Math::PxPyPzMVector target(0, 0, 0, PROTON_MASS);

    ROOT::Math::PxPyPzMVector sum;
    for (const auto& p : detected)
        sum += p;

    return (beam + target - sum).M();
}

#endif
