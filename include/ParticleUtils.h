#pragma once

#include "MomCorrParticle.h"

#include <vector>
#include <stdexcept>

inline const MomCorrParticle& FindParticleByPID(
    const std::vector<MomCorrParticle>& particles,
    int pid
) {
    for (const auto& particle : particles) {
        if (particle.GetPID() == pid) {
            return particle;
        }
    }

    throw std::runtime_error(
        "FindParticleByPID: no particle found with PID = " + std::to_string(pid)
    );
}
