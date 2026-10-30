"""Conservation-law fixtures, not simulated detector data or validated cuts."""

import copy
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]


def configs(name):
    directory = ROOT / "examples/channels" / name
    return tuple(
        json.loads((directory / file).read_text())
        for file in ("input.json", "binning.json", "cuts.json")
    )


def boost(vector, beta):
    energy, momentum = vector[0], vector[1:]
    beta2 = beta @ beta
    if beta2 == 0:
        return vector.copy()
    gamma = 1 / np.sqrt(1 - beta2)
    dot = beta @ momentum
    return np.r_[
        gamma * (energy + dot),
        momentum + ((gamma - 1) * dot / beta2 + gamma * energy) * beta,
    ]


def decay(parent, mass1, mass2, direction):
    mass = np.sqrt(parent[0] ** 2 - parent[1:] @ parent[1:])
    p = np.sqrt(
        (mass**2 - (mass1 + mass2) ** 2) * (mass**2 - (mass1 - mass2) ** 2)
    ) / (2 * mass)
    direction = np.asarray(direction, dtype=float)
    direction /= np.linalg.norm(direction)
    beta = parent[1:] / parent[0]
    return (
        boost(np.r_[np.hypot(p, mass1), p * direction], beta),
        boost(np.r_[np.hypot(p, mass2), -p * direction], beta),
    )


def exact_event(config):
    c = config["conditions"]
    masses = c["particle_masses_gev"]
    initial = np.array(
        [
            c["beam_energy_gev"] + masses["proton"],
            0,
            0,
            np.sqrt(c["beam_energy_gev"] ** 2 - masses["electron"] ** 2),
        ]
    )
    final = set(config["channel"]["final_state"])
    if final == {"electron", "proton"}:
        electron, proton = decay(
            initial, masses["electron"], masses["proton"], [0.5, 0.1, 0.86]
        )
        vectors = {"electron": electron, "proton": proton}
    else:
        electron_p = np.array([2 * np.sin(0.35), 0.0, 2 * np.cos(0.35)])
        electron = np.r_[
            np.sqrt(electron_p @ electron_p + masses["electron"] ** 2), electron_p
        ]
        hadronic = initial - electron
        if "neutron" in final:
            pion, neutron = decay(
                hadronic, masses["pi_plus"], masses["neutron"], [1.0, 0.1, 0.2]
            )
            vectors = {"electron": electron, "pi_plus": pion, "neutron": neutron}
        else:
            proton, pair = decay(hadronic, masses["proton"], 0.65, [0.3, 0.1, 1.0])
            plus, minus = decay(
                pair, masses["pi_plus"], masses["pi_minus"], [0.1, 1.0, 0.4]
            )
            vectors = {
                "electron": electron,
                "proton": proton,
                "pi_plus": plus,
                "pi_minus": minus,
            }
    np.testing.assert_allclose(sum(vectors.values()), initial, atol=1e-13)
    return vectors


def logical_rows(config, vectors, n=1, run=5875, residuals=None):
    data = {
        "run": np.full(n, run, dtype=np.int32),
        "event": np.arange(n, dtype=np.int32),
    }
    for particle in config["channel"]["measured"]:
        for i, axis in enumerate("xyz"):
            data[f"{particle}_p{axis}"] = np.full(n, vectors[particle][i + 1])
        data[f"{particle}_sector"] = np.ones(n, dtype=np.int32)
    if residuals is not None:
        electron = vectors["electron"]
        mass = config["conditions"]["particle_masses_gev"]["electron"]
        scale = np.sqrt((electron[0] - residuals) ** 2 - mass**2) / np.linalg.norm(
            electron[1:]
        )
        for axis in "xyz":
            data[f"electron_p{axis}"] *= scale
    return data


def one_cell(binning):
    binning = copy.deepcopy(binning)
    for specification in binning["particles"].values():
        specification["selected_scheme"] = "fixture"
        specification["schemes"] = {
            "fixture": {
                "sectors": [1],
                "edges": {
                    "p": [0.05, 10.0],
                    "theta": [0.0, 180.0],
                    "shifted_phi": [-1000.0, 1000.0],
                },
            }
        }
    return binning
