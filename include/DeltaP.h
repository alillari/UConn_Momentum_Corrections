#pragma once

namespace momcorr {

double PredictPipMomentum_epipN(
    double beam_E,
    double ele_p,
    double ele_theta,
    double pip_p,
    double pip_theta,
    double angle_ele_pip
);

double DeltaPip_epipN(
    double beam_E,
    double ele_p,
    double ele_theta,
    double pip_p,
    double pip_theta,
    double angle_ele_pip
);

}
