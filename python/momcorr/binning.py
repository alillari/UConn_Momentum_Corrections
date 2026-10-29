"""Explicit sector-separated 3D calibration cells in p, theta, shifted phi.

Units are GeV/c and degrees. The shifted azimuth is NOT wrapped after applying
the inverse-momentum shift: doing so would change the requested coordinate.
"""
from dataclasses import dataclass

import numpy as np


AXES = ("p", "theta", "shifted_phi")
REASONS = ("inside", "invalid_momentum", "invalid_sector", "shift_pole",
           "outside_p", "outside_theta", "outside_shifted_phi")


def edges_checked(values, label):
    edges = np.asarray(values, dtype=np.float64)
    if edges.ndim != 1 or len(edges) < 2 or not np.isfinite(edges).all() or not (np.diff(edges) > 0).all():
        raise ValueError(f"{label}: edges must be finite and strictly increasing")
    return edges


def shifted_coordinates(px, py, pz, sector, coefficient, pole):
    """Richard's sector unwrap, then phi_local + coefficient/(p-pole).

    Sector labels must be integral 1..6. Invalid or zero momenta and momenta at
    or below the inverse-shift pole receive an explicit invalid reason.
    """
    px, py, pz, sector = np.broadcast_arrays(*[np.asarray(x, dtype=np.float64) for x in (px, py, pz, sector)])
    p = np.sqrt(px * px + py * py + pz * pz)
    theta = np.degrees(np.arctan2(np.hypot(px, py), pz))
    phi = np.degrees(np.arctan2(py, px))
    unwrap = (((sector == 3) | (sector == 4)) & (phi < 0)) | ((sector > 4) & (phi < 90))
    local = np.where(unwrap, phi + 360, phi) - (sector - 1) * 60
    reasons = np.zeros(p.shape, dtype=np.uint8)
    bad_momentum = ~np.isfinite(p) | (p <= 0)
    reasons[bad_momentum] = 1
    bad_sector = ~np.isfinite(sector) | (sector < 1) | (sector > 6) | (sector != np.floor(sector))
    reasons[(reasons == 0) & bad_sector] = 2
    reasons[(reasons == 0) & (p <= pole)] = 3
    shifted = np.full(p.shape, np.nan)
    good = reasons == 0
    shifted[good] = local[good] + coefficient / (p[good] - pole)
    return {"p": p, "theta": theta, "phi": phi, "local_phi": local,
            "shifted_phi": shifted, "sector": sector, "reason": reasons}


def bin_indices(values, edges):
    """Half-open bins, except the final right edge is included, like NumPy."""
    values = np.asarray(values)
    indices = np.searchsorted(edges, values, side="right") - 1
    indices = np.where(values == edges[-1], len(edges) - 2, indices)
    valid = np.isfinite(values) & (values >= edges[0]) & (values <= edges[-1])
    return indices, valid


@dataclass
class Grid3D:
    sectors: tuple
    edges: dict

    def __post_init__(self):
        if not self.sectors or len(set(self.sectors)) != len(self.sectors) or any(type(s) is not int or s not in range(1, 7) for s in self.sectors):
            raise ValueError("Sectors must be unique integers 1..6")
        if set(self.edges) != set(AXES):
            raise ValueError("Grid must have exactly p, theta, shifted_phi axes")
        self.edges = {key: edges_checked(self.edges[key], key) for key in AXES}

    @property
    def shape(self):
        return (len(self.sectors),) + tuple(len(self.edges[key]) - 1 for key in AXES)

    @property
    def size(self):
        return int(np.prod(self.shape))

    def assign(self, coordinates):
        reasons = coordinates["reason"].copy()
        sector_indices = np.full(reasons.shape, -1, dtype=np.int64)
        for i, sector in enumerate(self.sectors):
            sector_indices[coordinates["sector"] == sector] = i
        reasons[(reasons == 0) & (sector_indices < 0)] = 2
        indices = [sector_indices]
        for code, key in enumerate(AXES, 4):
            ix, valid = bin_indices(coordinates[key], self.edges[key])
            reasons[(reasons == 0) & ~valid] = code
            indices.append(ix)
        inside = reasons == 0
        cells = np.full(reasons.shape, -1, dtype=np.int64)
        if inside.any():
            cells[inside] = np.ravel_multi_index(tuple(ix[inside] for ix in indices), self.shape)
        return cells, reasons

    def bounds(self, cell):
        si, *indices = np.unravel_index(cell, self.shape)
        return {"sector": self.sectors[si], **{key: [float(self.edges[key][i]), float(self.edges[key][i + 1])] for key, i in zip(AXES, indices)}}


def load_binning(config):
    required = {"binning_schema_version", "units", "coordinate_convention", "reference_band_gev", "mass_edges_gev", "occupancy_thresholds", "particles"}
    if set(config) != required or config["binning_schema_version"] != 1:
        raise ValueError("Unexpected 3D-binning configuration schema or keys")
    if config["units"] != {"p": "GeV/c", "theta": "deg", "shifted_phi": "deg"}:
        raise ValueError("Binning requires momentum in GeV/c and angles in degrees")
    if config["coordinate_convention"] != "richard_sector_local_inverse_momentum_shift":
        raise ValueError("Unknown shifted-phi convention")
    reference = edges_checked(config["reference_band_gev"], "reference_band_gev")
    if len(reference) != 2:
        raise ValueError("Reference band requires two edges")
    mass = edges_checked(config["mass_edges_gev"], "mass_edges_gev")
    if reference[0] < mass[0] or reference[-1] > mass[-1]:
        raise ValueError("Reference band must lie within the stored spectrum range")
    if any(type(v) is not int or v < 1 for v in config["occupancy_thresholds"]):
        raise ValueError("Occupancy thresholds must be positive integers")
    if set(config["particles"]) != {"electron", "pion"}:
        raise ValueError("Expected electron and pion binning specifications")
    result = {}
    for role, spec in config["particles"].items():
        if set(spec) != {"shift", "selected_scheme", "schemes"}:
            raise ValueError(f"Unexpected specification keys for {role}")
        if set(spec["shift"]) != {"coefficient_deg_gev", "pole_gev"}:
            raise ValueError("A shift requires its coefficient and pole")
        if not np.isfinite(list(spec["shift"].values())).all():
            raise ValueError("Shift coefficients must be finite")
        if not spec["schemes"] or spec["selected_scheme"] not in spec["schemes"]:
            raise ValueError("Selected scheme must be defined")
        result[role] = {}
        for name, scheme in spec["schemes"].items():
            if not name.replace("_", "").isalnum() or set(scheme) != {"sectors", "edges"}:
                raise ValueError("Invalid scheme name or keys")
            result[role][name] = Grid3D(tuple(scheme["sectors"]), scheme["edges"])
    return result
