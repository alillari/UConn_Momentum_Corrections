"""NumPy physics adapter for an electron/pion skim on a proton target."""

import numpy as np


def signed_sqrt(x):
    return np.copysign(np.sqrt(np.abs(x)), x)


def epip_kinematics(a, conditions):
    """All computations use float64; energy is total energy, masses in GeV."""
    e = np.column_stack([a[f"electron_p{k}"] for k in "xyz"]).astype(np.float64)
    p = np.column_stack([a[f"pion_p{k}"] for k in "xyz"]).astype(np.float64)
    pe, pp = np.linalg.norm(e, axis=1), np.linalg.norm(p, axis=1)
    me, mpi, mp, mn = [
        conditions[k]
        for k in (
            "electron_mass_gev",
            "pion_mass_gev",
            "target_mass_gev",
            "neutron_mass_gev",
        )
    ]
    eb = conditions["beam_energy_gev"]
    pb = np.sqrt(eb * eb - me * me)
    ee, ep = np.hypot(pe, me), np.hypot(pp, mpi)
    qvec = -e.copy()
    qvec[:, 2] += pb
    nu = eb - ee
    xvec = qvec - p
    xe = nu + mp - ep
    mm2 = xe * xe - np.sum(xvec * xvec, axis=1)
    w2 = (nu + mp) ** 2 - np.sum(qvec * qvec, axis=1)
    q2 = np.sum(qvec * qvec, axis=1) - nu * nu
    result = {
        "electron_p": pe,
        "pion_p": pp,
        "missing_mass2": mm2,
        "missing_mass": signed_sqrt(mm2),
        "q2": q2,
        "w": signed_sqrt(w2),
    }
    for name, v in (("electron", e), ("pion", p)):
        result[f"{name}_theta"] = np.degrees(
            np.arctan2(np.hypot(v[:, 0], v[:, 1]), v[:, 2])
        )
        result[f"{name}_phi"] = np.degrees(np.arctan2(v[:, 1], v[:, 0]))

    # Kinematic pion-momentum estimator with an unobserved on-shell neutron.
    # Squaring introduces candidate roots: check the original energy equation.
    # Multiple physical roots remain ambiguous and are explicitly NOT resolved
    # by choosing whichever is closest to the measured pion momentum.
    with np.errstate(invalid="ignore", divide="ignore"):
        nhat = p / pp[:, None]
        ke = nu + mp
        kval2 = np.sum(qvec * qvec, axis=1)
        aa = w2 + mpi * mpi - mn * mn
        bb = np.sum(qvec * nhat, axis=1)
        qa, qb, qc = ke * ke - bb * bb, -aa * bb, ke * ke * mpi * mpi - aa * aa / 4
        disc = qb * qb - 4 * qa * qc
        root = np.sqrt(np.maximum(disc, 0))
        roots = [(-qb + root) / (2 * qa), (-qb - root) / (2 * qa)]
        physical = []
        for r in roots:
            neutron_e = ke - np.hypot(r, mpi)
            residual = ke * np.hypot(r, mpi) - bb * r - aa / 2
            ok = (
                (disc >= 0)
                & (np.abs(qa) > 1e-14)
                & (pp > 0)
                & (r > 0)
                & np.isfinite(r)
                & (neutron_e > 0)
                & (np.abs(residual) < 1e-8 * np.maximum(1, np.abs(aa)))
            )
            physical.append(ok)
        same = np.isclose(roots[0], roots[1], rtol=1e-9, atol=1e-10)
        unique = (physical[0] ^ physical[1]) | (physical[0] & physical[1] & same)
        predicted = np.where(physical[0], roots[0], roots[1])
        result["delta_pion_p"] = np.where(unique, predicted - pp, np.nan)
        result["pion_root_status"] = np.where(
            unique, 1, np.where(physical[0] & physical[1], 2, 0)
        )
    result["missing_energy"] = xe
    return result
