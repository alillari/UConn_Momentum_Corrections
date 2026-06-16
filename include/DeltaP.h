#pragma once

#include <Math/Vector4D.h>

namespace momcorr {

using LorentzVector = ROOT::Math::PxPyPzMVector;

double PredictPipMomentum_epipN_FourVector(
    const LorentzVector& beam,
    const LorentzVector& target,
    const LorentzVector& electron,
    const LorentzVector& pion_meas
);

double DeltaPip_epipN_FourVector(
    const LorentzVector& beam,
    const LorentzVector& target,
    const LorentzVector& electron,
    const LorentzVector& pion_meas
);

double PredictPipMomentum_epipN_Legacy(
    double beam_E,
    double ele_p,
    double ele_theta,
    double pip_p,
    double pip_theta,
    double angle_ele_pip
);

double DeltaPip_epipN_Legacy(
    double beam_E,
    double ele_p,
    double ele_theta,
    double pip_p,
    double pip_theta,
    double angle_ele_pip
);

double DeltaPip_epipN_FourVector_RDF(
    double beam_E,
    float ele_px, float ele_py, float ele_pz,
    float pip_px, float pip_py, float pip_pz
);

double DeltaPip_epipN_Legacy_RDF(
    double beam_E,
    float ele_px, float ele_py, float ele_pz,
    float pip_px, float pip_py, float pip_pz
);

}
