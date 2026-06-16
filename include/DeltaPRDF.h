#pragma once

#include "DeltaP.h"
#include "MomCorrParticle.h"
#include "ParticleUtils.h"

#include <ROOT/RDataFrame.hxx>

#include <string>
#include <vector>
#include <stdexcept>

namespace momcorr {

inline ROOT::RDF::RNode DefineDeltaP_epipN(
    ROOT::RDF::RNode df,
    const std::vector<MomCorrParticle>& particles,
    double beam_energy
) {
    const auto& ele = FindParticleByPID(particles, 11);
    const auto& pip = FindParticleByPID(particles, 211);

    const std::string beam_branch = "beam_energy_const";

    const std::string pip_deltaP_fourvec =
        pip.GetName() + "_deltaP_fourvec";

    const std::string pip_deltaP_legacy =
        pip.GetName() + "_deltaP_legacy";

    const std::string pip_deltaP_diff =
        pip.GetName() + "_deltaP_diff";

    df = df.Define(beam_branch, [beam_energy]() {
        return beam_energy;
    });

    df = df.Define(
        pip_deltaP_fourvec,
        momcorr::DeltaPip_epipN_FourVector_RDF,
        {
            beam_branch,
            ele.GetInputPxBranch(), ele.GetInputPyBranch(), ele.GetInputPzBranch(),
            pip.GetInputPxBranch(), pip.GetInputPyBranch(), pip.GetInputPzBranch()
        }
    );

    df = df.Define(
        pip_deltaP_legacy,
        momcorr::DeltaPip_epipN_Legacy_RDF,
        {
            beam_branch,
            ele.GetInputPxBranch(), ele.GetInputPyBranch(), ele.GetInputPzBranch(),
            pip.GetInputPxBranch(), pip.GetInputPyBranch(), pip.GetInputPzBranch()
        }
    );

    df = df.Define(
        pip_deltaP_diff,
        [](double fourvec, double legacy) {
            return fourvec - legacy;
        },
        {pip_deltaP_fourvec, pip_deltaP_legacy}
    );

    return df;
}

inline ROOT::RDF::RNode DefineDeltaP(
    ROOT::RDF::RNode df,
    const std::vector<MomCorrParticle>& particles,
    double beam_energy,
    const std::string& channel
) {
    if (channel == "epipN" || channel == "epipX" || channel == "SP") {
        return DefineDeltaP_epipN(df, particles, beam_energy);
    }

    throw std::runtime_error(
        "DefineDeltaP: unsupported channel '" + channel + "'."
    );
}

}
