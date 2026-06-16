#include "DeltaP.h"

#include <Math/Vector3D.h>
#include <Math/VectorUtil.h>

#include <cmath>
#include <limits>

namespace momcorr {

namespace {
constexpr double MASS_ELECTRON = 0.00051099895;
constexpr double MASS_PROTON   = 0.93827208816;
constexpr double MASS_NEUTRON  = 0.93956542052;
constexpr double MASS_PIP      = 0.13957039;

double NaN() {
    return std::numeric_limits<double>::quiet_NaN();
}

bool Bad(double x) {
    return !std::isfinite(x);
}
}

double PredictPipMomentum_epipN_FourVector(
    const LorentzVector& beam,
    const LorentzVector& target,
    const LorentzVector& electron,
    const LorentzVector& pion_meas
) {
    const auto K = beam + target - electron;

    const double K_E = K.E();
    const auto K_vec = K.Vect();

    const double p_meas = pion_meas.P();
    if (p_meas <= 0.0) return NaN();

    const auto n_hat = pion_meas.Vect().Unit();

    const double A = K.M2() + MASS_PIP*MASS_PIP - MASS_NEUTRON*MASS_NEUTRON;
    const double B = K_vec.Dot(n_hat);

    const double a = K_E*K_E - B*B;
    const double b = -A * B;
    const double c = K_E*K_E * MASS_PIP*MASS_PIP - 0.25*A*A;

    const double disc = b*b - 4.0*a*c;

    if (disc < 0.0) return NaN();
    if (std::abs(a) < 1e-14) return NaN();

    const double root = std::sqrt(disc);

    const double p1 = (-b + root) / (2.0*a);
    const double p2 = (-b - root) / (2.0*a);

    const bool p1_good = std::isfinite(p1) && p1 > 0.0;
    const bool p2_good = std::isfinite(p2) && p2 > 0.0;

    if (!p1_good && !p2_good) return NaN();
    if (p1_good && !p2_good) return p1;
    if (!p1_good && p2_good) return p2;

    return (std::abs(p_meas - p1) <= std::abs(p_meas - p2)) ? p1 : p2;
}

double DeltaPip_epipN_FourVector(
    const LorentzVector& beam,
    const LorentzVector& target,
    const LorentzVector& electron,
    const LorentzVector& pion_meas
) {
    const double p_pred = PredictPipMomentum_epipN_FourVector(
        beam, target, electron, pion_meas
    );

    if (!std::isfinite(p_pred)) return NaN();

    return p_pred - pion_meas.P();
}

double PredictPipMomentum_epipN_Legacy(
    double Beam_Energy,
    double ele_p,
    double ele_theta,
    double pip_p,
    double pip_theta,
    double angle_ele_pip
) {
    if (Bad(Beam_Energy) || Bad(ele_p) || Bad(ele_theta) ||
        Bad(pip_p) || Bad(pip_theta) || Bad(angle_ele_pip)) {
        return NaN();
    }

    const double neutronM2 = MASS_NEUTRON * MASS_NEUTRON;

    const double termA =
        (neutronM2 - MASS_PROTON*MASS_PROTON - MASS_PIP*MASS_PIP) / 2.0;

    const double termB =
        MASS_PROTON * (Beam_Energy - ele_p)
        - Beam_Energy * ele_p * (1.0 - std::cos(ele_theta));

    const double termC =
        ele_p * std::cos(angle_ele_pip)
        - Beam_Energy * std::cos(pip_theta);

    const double W = MASS_PROTON + Beam_Energy - ele_p;

    const double sqrtTerm =
        (termA - termB) * (termA - termB)
        + MASS_PIP * MASS_PIP * (termC * termC - W * W);

    if (sqrtTerm < 0.0) return NaN();

    const double denominator = (W + termC) * (W - termC);
    if (std::abs(denominator) < 1e-14) return NaN();

    const double root = std::sqrt(sqrtTerm);

    const double numeratorP = (termA - termB) * termC + W * root;
    const double numeratorM = (termA - termB) * termC - W * root;

    const double pP = numeratorP / denominator;
    const double pM = numeratorM / denominator;

    const bool pP_good = std::isfinite(pP) && pP > 0.0;
    const bool pM_good = std::isfinite(pM) && pM > 0.0;

    if (!pP_good && !pM_good) return NaN();
    if (pP_good && !pM_good) return pP;
    if (!pP_good && pM_good) return pM;

    return (std::abs(pip_p - pP) <= std::abs(pip_p - pM)) ? pP : pM;
}

double DeltaPip_epipN_Legacy(
    double beam_E,
    double ele_p,
    double ele_theta,
    double pip_p,
    double pip_theta,
    double angle_ele_pip
) {
    const double p_pred = PredictPipMomentum_epipN_Legacy(
        beam_E, ele_p, ele_theta, pip_p, pip_theta, angle_ele_pip
    );

    if (!std::isfinite(p_pred)) return NaN();

    return p_pred - pip_p;
}

double DeltaPip_epipN_FourVector_RDF(
    double beam_E,
    float ele_px, float ele_py, float ele_pz,
    float pip_px, float pip_py, float pip_pz
) {
    const LorentzVector beam(0.0, 0.0, beam_E, MASS_ELECTRON);
    const LorentzVector target(0.0, 0.0, 0.0, MASS_PROTON);

    const LorentzVector electron(ele_px, ele_py, ele_pz, MASS_ELECTRON);
    const LorentzVector pion(pip_px, pip_py, pip_pz, MASS_PIP);

    return DeltaPip_epipN_FourVector(beam, target, electron, pion);
}

double DeltaPip_epipN_Legacy_RDF(
    double beam_E,
    float ele_px, float ele_py, float ele_pz,
    float pip_px, float pip_py, float pip_pz
) {
    ROOT::Math::XYZVector ele(ele_px, ele_py, ele_pz);
    ROOT::Math::XYZVector pip(pip_px, pip_py, pip_pz);

    const double ele_p = ele.R();
    const double pip_p = pip.R();

    if (ele_p <= 0.0 || pip_p <= 0.0) return NaN();

    const double angle = ROOT::Math::VectorUtil::Angle(ele, pip);

    return DeltaPip_epipN_Legacy(
        beam_E,
        ele_p,
        ele.Theta(),
        pip_p,
        pip.Theta(),
        angle
    );
}

}
