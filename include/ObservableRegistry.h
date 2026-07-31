#ifndef OBSERVABLE_REGISTRY_H
#define OBSERVABLE_REGISTRY_H

#include "DeltaP.h"
#include "Kinematics.h"
#include "MomCorrConfig.h"

#include <ROOT/RDataFrame.hxx>

#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace momcorr {

inline std::string DoubleLiteral(double value) {
    std::ostringstream ss;
    ss << std::setprecision(17) << value;
    return ss.str();
}

inline std::string P4Expression(const MomCorrParticle& particle) {
    return "ROOT::Math::PxPyPzMVector(" +
           particle.GetInputPxBranch() + "," +
           particle.GetInputPyBranch() + "," +
           particle.GetInputPzBranch() + "," +
           DoubleLiteral(particle.GetMass()) + ")";
}

inline std::string MissingMassExpression(
    const ObservableConfig& observable,
    const std::vector<MomCorrParticle>& particles,
    const MomCorrConfig& config
) {
    if (observable.detectedRoles.empty()) {
        throw std::runtime_error(
            "missing_mass observable '" + observable.name +
            "' requires detected_roles"
        );
    }

    std::string expression =
        "ComputeMissingMass(std::vector<ROOT::Math::PxPyPzMVector>{";

    for (std::size_t i = 0; i < observable.detectedRoles.size(); ++i) {
        const auto& particle = config.ParticleByRole(
            particles, observable.detectedRoles[i]);
        if (i > 0) expression += ",";
        expression += P4Expression(particle);
    }

    expression += "}," + DoubleLiteral(config.GetBeamEnergy()) + ")";
    return expression;
}

inline const std::string& RequiredRole(const ObservableConfig& observable,
                                       const std::string& key) {
    auto it = observable.roles.find(key);
    if (it == observable.roles.end() || it->second.empty()) {
        throw std::runtime_error(
            "observable '" + observable.name + "' requires role key '" + key + "'"
        );
    }
    return it->second;
}

inline ROOT::RDF::RNode DefineObservable(
    ROOT::RDF::RNode df,
    const ObservableConfig& observable,
    const std::vector<MomCorrParticle>& particles,
    const MomCorrConfig& config
) {
    if (observable.type == "missing_mass") {
        return df.Define(
            observable.output,
            MissingMassExpression(observable, particles, config)
        );
    }

    if (observable.type == "deltaP_epipN_fourvec" ||
        observable.type == "deltaP_epipN_legacy") {
        const auto& electron = config.ParticleByRole(
            particles, RequiredRole(observable, "electron"));
        const auto& pion = config.ParticleByRole(
            particles, RequiredRole(observable, "pion"));

        const std::string functionName =
            (observable.type == "deltaP_epipN_fourvec")
                ? "momcorr::DeltaPip_epipN_FourVector_RDF"
                : "momcorr::DeltaPip_epipN_Legacy_RDF";

        const std::string expression =
            functionName + "(" +
            DoubleLiteral(config.GetBeamEnergy()) + "," +
            electron.GetInputPxBranch() + "," +
            electron.GetInputPyBranch() + "," +
            electron.GetInputPzBranch() + "," +
            pion.GetInputPxBranch() + "," +
            pion.GetInputPyBranch() + "," +
            pion.GetInputPzBranch() + ")";

        return df.Define(observable.output, expression);
    }

    if (observable.type == "difference") {
        const auto& lhs = RequiredRole(observable, "lhs");
        const auto& rhs = RequiredRole(observable, "rhs");
        return df.Define(observable.output, lhs + " - " + rhs);
    }

    throw std::runtime_error("Unknown observable type '" + observable.type + "'");
}

inline ROOT::RDF::RNode DefineObservables(
    ROOT::RDF::RNode df,
    const std::vector<ObservableConfig>& observables,
    const std::vector<MomCorrParticle>& particles,
    const MomCorrConfig& config
) {
    for (const auto& observable : observables) {
        df = DefineObservable(df, observable, particles, config);
    }
    return df;
}

} // namespace momcorr

#endif
