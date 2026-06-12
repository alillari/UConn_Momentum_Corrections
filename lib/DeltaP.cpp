#include "DeltaP.h"

#include <cmath>
#include <limits>

namespace momcorr {

double PredictPipMomentum_epipN(
    double Beam_Energy,
    double ele_p,
    double ele_theta,
    double pip_p,
    double pip_theta,
    double angle_ele_pip
) {
    constexpr double neutronM = 0.93956542052;
    constexpr double Proton_M = 0.93827208816;
    constexpr double Pion_C_M = 0.13957039;

    const double neutronM2 = neutronM * neutronM;

    const double termA =
        (neutronM2 - Proton_M*Proton_M - Pion_C_M*Pion_C_M) / 2.0;

    const double termB =
        Proton_M * (Beam_Energy - ele_p)
        - Beam_Energy * ele_p * (1.0 - std::cos(ele_theta));

    const double termC =
        ele_p * std::cos(angle_ele_pip)
        - Beam_Energy * std::cos(pip_theta);

    const double W =
        Proton_M + Beam_Energy - ele_p;

    const double sqrtTerm =
        (termA - termB) * (termA - termB)
        + Pion_C_M * Pion_C_M * (termC * termC - W * W);

    if (sqrtTerm < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double denominator = (W + termC) * (W - termC);

    if (std::abs(denominator) < 1e-12) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double root = std::sqrt(sqrtTerm);

    const double numeratorP =
        (termA - termB) * termC + W * root;

    const double numeratorM =
        (termA - termB) * termC - W * root;

    const double pip_CalculateP = numeratorP / denominator;
    const double pip_CalculateM = numeratorM / denominator;

    if (std::abs(pip_p - pip_CalculateP) <=
        std::abs(pip_p - pip_CalculateM)) {
        return pip_CalculateP;
    }

    return pip_CalculateM;
}

double DeltaPip_epipN(
    double beam_E,
    double ele_p,
    double ele_theta,
    double pip_p,
    double pip_theta,
    double angle_ele_pip
) {
    const double pred = PredictPipMomentum_epipN(
        beam_E, ele_p, ele_theta, pip_p, pip_theta, angle_ele_pip
    );

    if (!std::isfinite(pred)) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return pred - pip_p;
}

}
