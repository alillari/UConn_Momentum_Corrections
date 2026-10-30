"""Electron/pion skim coordinate and observable adapter, independent of fitting."""

import numpy as np
from .binning import shifted_coordinates
from .kinematics import epip_kinematics


def coordinates(data, role, specification):
    shift = specification["shift"]
    return shifted_coordinates(
        *[data[f"{role}_p{k}"] for k in "xyz"],
        data[f"{role}_sector"],
        shift["coefficient_deg_gev"],
        shift["pole_gev"],
    )


def observables(data, conditions):
    finite = np.logical_and.reduce(
        [
            np.isfinite(data[f"{role}_p{k}"])
            for role in ("electron", "pion")
            for k in "xyz"
        ]
    )
    result = {
        key: np.full(len(finite), np.nan)
        for key in ("missing_mass", "delta_pion_p", "electron_p", "pion_p")
    }
    computed = epip_kinematics(
        {key: value[finite] for key, value in data.items()}, conditions
    )
    for key in result:
        result[key][finite] = computed[key]
    return result
