"""Explicit measured topologies and four-vector observables for charged tracks.

The schema describes a selected event/candidate row; it is not a PID algorithm.
All beams are electrons along +z and the target is a proton at rest. Particle
names encode charge. In particular, a negative pion is never a generic pion.
"""

import json
from pathlib import Path

import numpy as np

from .kinematics import signed_sqrt


# Preserve the constants used by the original RG-K adapter, in GeV.
MASSES = {
    "electron": 0.000511,
    "pi_plus": 0.139570,
    "pi_minus": 0.139570,
    "proton": 0.938272,
    "neutron": 0.9395654,
}
CHARGES = {"electron": -1, "pi_plus": 1, "pi_minus": -1, "proton": 1, "neutron": 0}
TRACKS = frozenset({"electron", "pi_plus", "pi_minus", "proton"})
OBSERVABLE_UNITS = {
    "missing_mass": "GeV/c^2",
    "missing_mass2": "GeV^2/c^4",
    "missing_energy": "GeV",
    "missing_px": "GeV/c",
    "missing_py": "GeV/c",
    "missing_pz": "GeV/c",
    "missing_p": "GeV/c",
    "missing_pt": "GeV/c",
}


def exact_keys(value, keys, label):
    if not isinstance(value, dict) or set(value) != set(keys):
        raise ValueError(f"{label} requires exactly {sorted(keys)}")


def validate_channel(channel):
    exact_keys(channel, {"name", "final_state", "measured"}, "channel")
    if not isinstance(channel["name"], str) or not channel["name"].strip():
        raise ValueError("Channel name must be nonempty")
    for key in ("final_state", "measured"):
        values = channel[key]
        if (
            not isinstance(values, list)
            or not values
            or any(not isinstance(v, str) or v not in MASSES for v in values)
            or len(set(values)) != len(values)
        ):
            raise ValueError(f"Invalid or repeated particle identities in {key}")
    final, measured = set(channel["final_state"]), set(channel["measured"])
    if "electron" not in measured or not measured <= final or not measured <= TRACKS:
        raise ValueError(
            "Measured particles must be charged final-state tracks including electron"
        )
    if sum(CHARGES[p] for p in final) != 0:
        raise ValueError(
            "Final-state charge is inconsistent with electron + proton initial state"
        )
    if len(final & {"proton", "neutron"}) != 1:
        raise ValueError("This adapter requires one final-state baryon")
    four_track = {"electron", "proton", "pi_plus", "pi_minus"}
    if "pi_minus" in final and (final != four_track or measured != four_track):
        raise ValueError(
            "pi_minus requires the fully measured ep pi_plus pi_minus channel"
        )
    return channel


def required_fields(channel):
    return {"run", "event"} | {
        f"{particle}_{suffix}"
        for particle in channel["measured"]
        for suffix in ("px", "py", "pz", "sector")
    }


def validate_input_config(config):
    exact_keys(
        config,
        {"input_schema_version", "dataset", "channel", "branches", "conditions"},
        "input config",
    )
    if config["input_schema_version"] != 2:
        raise ValueError("Expected input schema 2")
    validate_channel(config["channel"])
    if (
        not isinstance(config["dataset"], dict)
        or not isinstance(config["dataset"].get("tree"), str)
        or not config["dataset"]["tree"]
    ):
        raise ValueError("Dataset must declare a nonempty tree name")
    exact_keys(config["branches"], required_fields(config["channel"]), "branch mapping")
    branch_names = list(config["branches"].values())
    if any(not isinstance(v, str) or not v for v in branch_names) or len(
        set(branch_names)
    ) != len(branch_names):
        raise ValueError(
            "Every logical field requires a distinct nonempty source branch"
        )
    c = config["conditions"]
    exact_keys(c, {"beam_energy_gev", "target", "particle_masses_gev"}, "conditions")
    if c["target"] != "proton_at_rest":
        raise ValueError("Only a proton-at-rest target is implemented")
    exact_keys(c["particle_masses_gev"], MASSES, "particle masses")
    masses = list(c["particle_masses_gev"].values())
    if not np.isfinite(masses).all() or any(m <= 0 for m in masses):
        raise ValueError("Particle masses must be finite and positive")
    if (
        not np.isfinite(c["beam_energy_gev"])
        or c["beam_energy_gev"] <= c["particle_masses_gev"]["electron"]
    ):
        raise ValueError("Beam total energy must exceed electron mass")
    return config


def load_input_config(path):
    """Load schema 2 or retain the existing, explicitly pi+ schema-1 adapter."""
    config = json.loads(Path(path).read_text())
    if config.get("input_schema_version") == 2:
        return validate_input_config(config)
    # Keep the old diagnostics interface available to its standalone tools.
    # The cut CLI validates legacy input using that module before reaching here.
    if config.get("diagnostics_schema_version") == 1:
        return config
    raise ValueError("Unknown input configuration schema")


def validate_observable(spec, channel=None):
    exact_keys(spec, {"kind", "omit"}, "observable")
    if spec["kind"] not in OBSERVABLE_UNITS:
        raise ValueError("Unknown four-vector observable")
    omit = spec["omit"]
    if (
        not isinstance(omit, list)
        or any(not isinstance(p, str) or p not in TRACKS for p in omit)
        or len(set(omit)) != len(omit)
    ):
        raise ValueError("Omitted tracks must have unique charged-particle identities")
    if channel is not None and not set(omit) <= set(channel["measured"]):
        raise ValueError(
            "Observable cannot omit a track absent from the measured topology"
        )
    return spec


def residual_four_vector(data, config, omit=()):
    """Return initial minus included tracks as (energy, px, py, pz).

    Omitting a track defines a leave-one-out observable, not a missing-track input
    topology. Every declared measured momentum must remain finite, even if omitted.
    """
    validate_input_config(config)
    validate_observable(
        {"kind": "missing_energy", "omit": list(omit)}, config["channel"]
    )
    n = len(data["run"])
    missing_fields = required_fields(config["channel"]) - set(data)
    if missing_fields:
        raise ValueError(f"Missing measured-track fields: {sorted(missing_fields)}")
    if any(
        np.asarray(data[k]).shape != (n,) for k in required_fields(config["channel"])
    ):
        raise ValueError("Input fields must be aligned one-dimensional arrays")
    conditions = config["conditions"]
    masses = conditions["particle_masses_gev"]
    beam_e = conditions["beam_energy_gev"]
    residual = np.zeros((n, 4), dtype=np.float64)
    residual[:, 0] = beam_e + masses["proton"]
    residual[:, 3] = np.sqrt(beam_e**2 - masses["electron"] ** 2)
    finite = np.ones(n, dtype=bool)
    for particle in config["channel"]["measured"]:
        momentum = np.column_stack([data[f"{particle}_p{k}"] for k in "xyz"]).astype(
            np.float64
        )
        finite &= np.isfinite(momentum).all(axis=1) & (np.sum(momentum**2, axis=1) > 0)
        if particle not in omit:
            residual[:, 0] -= np.sqrt(
                np.sum(momentum**2, axis=1) + masses[particle] ** 2
            )
            residual[:, 1:] -= momentum
    residual[~finite] = np.nan
    return residual


def observable_values(data, config, recipe):
    if config.get("input_schema_version") != 2:
        from .skim import observables

        return observables(data, config["conditions"])[recipe["observable"]]
    spec = validate_observable(recipe["observable"], config["channel"])
    residual = residual_four_vector(data, config, spec["omit"])
    energy, momentum = residual[:, 0], residual[:, 1:]
    kind = spec["kind"]
    if kind in {"missing_mass", "missing_mass2"}:
        mass2 = energy**2 - np.sum(momentum**2, axis=1)
        return mass2 if kind == "missing_mass2" else signed_sqrt(mass2)
    if kind == "missing_energy":
        return energy
    if kind in {"missing_px", "missing_py", "missing_pz"}:
        return momentum[:, {"missing_px": 0, "missing_py": 1, "missing_pz": 2}[kind]]
    return np.linalg.norm(momentum[:, :2] if kind == "missing_pt" else momentum, axis=1)


def validate_workflow(config, binning, recipe):
    """Reject mismatched schemas/topologies before reading or fitting any data."""
    modern = config.get("input_schema_version") == 2
    if modern != (recipe["recipe_schema_version"] == 2) or modern != (
        binning["binning_schema_version"] == 2
    ):
        raise ValueError("Input, binning, and recipe must use the same workflow schema")
    if not modern:
        return
    validate_input_config(config)
    validate_observable(recipe["observable"], config["channel"])
    measured = set(config["channel"]["measured"])
    if not set(binning["particles"]) <= measured or not set(recipe["roles"]) <= set(
        binning["particles"]
    ):
        raise ValueError(
            "Every calibration role requires a measured track and an explicit grid"
        )
