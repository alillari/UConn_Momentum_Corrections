#ifndef EXCLUSIVE_KINEMATICS_H
#define EXCLUSIVE_KINEMATICS_H

#include <Math/Vector4D.h>

// ep → eπ⁺n



double SolveTruePionMomentum_epToEpiN(const ROOT::Math::PxPyPzMVector &beam, const ROOT::Math::PxPyPzMVector &target, const ROOT::Math::PxPyPzMVector &scattered, const ROOT::Math::PxPyPzMVector &reconPip) {

        ROOT::Math::PxPyPzMVector known = beam + target - scattered;
	double cosTheta = known.Vect().Dot(reconPip.Vect()) / (known.P() * reconPip.P());

        //We assume the direction is along the measured direction, but magnitude is adjustable
        //Calculated using conservation of 4 momentum
	
        double alpha = known.P()*known.P()*cosTheta*cosTheta - known.E()*known.E();                                                             //Quadratic term
        double beta = (NEUTRON_MASS*NEUTRON_MASS - PIP_MASS*PIP_MASS - known.M2())*cosTheta*known.P();                                          //Linear term
        double gamma = .25*std::pow(NEUTRON_MASS*NEUTRON_MASS - PIP_MASS*PIP_MASS - known.M2(),2) - known.E()*known.E()*PIP_MASS*PIP_MASS;      //Constant term

        if(beta*beta - 4*alpha*gamma < 0){
                std::cout<<"We have an event inconsistent with our assumptions. This is not good."<<std::endl;
                return -1.0;
        }

        double soln1 = (-beta + std::sqrt(beta*beta - 4*alpha*gamma))/(2*alpha);
        double soln2 = (-beta - std::sqrt(beta*beta - 4*alpha*gamma))/(2*alpha);

	if (std::abs(soln1 - reconPip.P()) < std::abs(soln2 - reconPip.P()))
    		return soln1;
	else
    		return soln2;

}

#endif
