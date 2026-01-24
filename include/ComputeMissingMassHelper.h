#ifndef MISSING_MASS_HELPER_H
#define MISSING_MASS_HELPER_H

#include <ROOT/RDataFrame.hxx>
#include <Math/Vector4D.h>
#include <vector>
#include <string>
#include <tuple>
#include "MomCorrParticle.h"

//Units in GeV
constexpr double ELECTRON_MASS = 0.000511;
constexpr double MUON_MASS = 0.105658;
constexpr double PROTON_MASS = 0.938272;
constexpr double PIP_MASS = 0.13957;
constexpr double NEUTRON_MASS = 0.939565;

std::vector<std::string> ExtractMomentumBranches(const std::vector<MomCorrParticle>& particles) {
    std::vector<std::string> branches;
    for (const auto& particle : particles) {
        branches.push_back(particle.GetPxBranch());
        branches.push_back(particle.GetPyBranch());
        branches.push_back(particle.GetPzBranch());
    }
    return branches;
}

inline double compute_missing_mass(
    const std::vector<std::tuple<float, float, float, double>>& particles,
    double beam_energy
) {
    ROOT::Math::PxPyPzMVector beam_vector(0, 0, sqrt(pow(beam_energy, 2) - pow(ELECTRON_MASS, 2)), ELECTRON_MASS);
    ROOT::Math::PxPyPzMVector target_vector(0, 0, 0, PROTON_MASS);

    ROOT::Math::PxPyPzMVector detected_sum(0, 0, 0, 0);
    for (const auto& [px, py, pz, mass] : particles) {
        detected_sum += ROOT::Math::PxPyPzMVector(px, py, pz, mass);
    }

    auto missing_vector = (beam_vector + target_vector) - detected_sum;
    return missing_vector.M();
}

template <typename T>
inline T DefineMissingMass(T& df, const std::vector<MomCorrParticle>& particle_list, double beam_energy) {
    std::vector<std::string> branches = ExtractMomentumBranches(particle_list);

    if (branches.size() == 3) {
        return df.Define("missing_mass",
            [=](float px1, float py1, float pz1) {
                std::vector<std::tuple<float, float, float, double>> particle_momenta;
                particle_momenta.emplace_back(px1, py1, pz1, particle_list[0].GetMass());
                return compute_missing_mass(particle_momenta, beam_energy);
            },
            branches);
    }
    if (branches.size() == 6) {
        return df.Define("missing_mass",
            [=](float px1, float py1, float pz1, float px2, float py2, float pz2) {
                std::vector<std::tuple<float, float, float, double>> particle_momenta;
                particle_momenta.emplace_back(px1, py1, pz1, particle_list[0].GetMass());
                particle_momenta.emplace_back(px2, py2, pz2, particle_list[1].GetMass());
                return compute_missing_mass(particle_momenta, beam_energy);
            },
            branches);
    }
    else if (branches.size() == 9) {
        return df.Define("missing_mass",
            [=](float px1, float py1, float pz1, float px2, float py2, float pz2, float px3, float py3, float pz3) {
                std::vector<std::tuple<float, float, float, double>> particle_momenta;
                particle_momenta.emplace_back(px1, py1, pz1, particle_list[0].GetMass());
                particle_momenta.emplace_back(px2, py2, pz2, particle_list[1].GetMass());
                particle_momenta.emplace_back(px3, py3, pz3, particle_list[2].GetMass());
                return compute_missing_mass(particle_momenta, beam_energy);
            },
            branches);
    }
    else if (branches.size() == 12) {
        return df.Define("missing_mass",
            [=](float px1, float py1, float pz1, float px2, float py2, float pz2, float px3, float py3, float pz3,
                float px4, float py4, float pz4) {
                std::vector<std::tuple<float, float, float, double>> particle_momenta;
                particle_momenta.emplace_back(px1, py1, pz1, particle_list[0].GetMass());
                particle_momenta.emplace_back(px2, py2, pz2, particle_list[1].GetMass());
                particle_momenta.emplace_back(px3, py3, pz3, particle_list[2].GetMass());
                particle_momenta.emplace_back(px4, py4, pz4, particle_list[3].GetMass());
                return compute_missing_mass(particle_momenta, beam_energy);
            },
            branches);
    }
    else {
        throw std::runtime_error("Unsupported number of particles. Adjust DefineMissingMass function.");
    }
}

//Goal of this is to calculate the true momentum magnitude of the pion
double findTruePionMom(const ROOT::Math::PxPyPzMVector &beam,
		const ROOT::Math::PxPyPzMVector &target,
		const ROOT::Math::PxPyPzMVector &scattered,
		const ROOT::Math::PxPyPzMVector &reconPip) {

	ROOT::Math::PxPyPzMVector known = beam + target - scattered;
	double cosTheta = known.Dot(reconPip) / (known.R() + reconPip.R());

	//We assume the direction is along the measured direction, but magnitude is adjustable
	//Calculated using conservation of 4 momentum
	double alpha = known.P()*known.P()*cosTheta*cosTheta - known.E()*known.E();								//Quadratic term
	double beta = (NEUTRON_MASS*NEUTRON_MASS - PIP_MASS*PIP_MASS - known.M2())*cosTheta*known.P();					        //Linear term
	double gamma = .25*std::pow(NEUTRON_MASS*NEUTRON_MASS - PIP_MASS*PIP_MASS - known.M2(),2) - known.E()*known.E()*PIP_MASS*PIP_MASS;	//Constant term

	if(beta*beta - 4*alpha*gamma < 0){
		std::cout<<"We have an event inconsistent with our assumptions. This is not good."<<std::endl;
		return -1.0;
	}

	double soln1 = (-beta + std::sqrt(beta*beta - 4*alpha*gamma))/(2*alpha);
	double soln2 = (-beta - std::sqrt(beta*beta - 4*alpha*gamma))/(2*alpha);

	if(soln1 - reconPip.P() < soln2 - reconPip.P()){
		return soln1;
	} else {
		return soln2;
	}
}
#endif  // MISSING_MASS_HELPER_H

