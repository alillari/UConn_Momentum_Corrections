"""Deterministic Poisson fits and model/range stability checks for 1D counts.

No run, particle, sector, histogram name, or cell-specific branch is used here.
Local Fisher errors describe one model only; fit alternatives are correlated
and are never treated as independent measurements to shrink uncertainty.
"""

from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import minimize
from scipy.special import xlogy

from .models import PeakBackground


@dataclass
class Histogram:
    edges: np.ndarray
    counts: np.ndarray

    def __post_init__(self):
        self.edges = np.asarray(self.edges, dtype=float)
        self.counts = np.asarray(self.counts, dtype=float)
        if (
            self.edges.ndim != 1
            or self.counts.ndim != 1
            or len(self.edges) != len(self.counts) + 1
            or not np.isfinite(self.edges).all()
            or not (np.diff(self.edges) > 0).all()
            or not np.isfinite(self.counts).all()
            or (self.counts < 0).any()
            or not np.allclose(self.counts, np.rint(self.counts), rtol=0, atol=1e-8)
        ):
            raise ValueError(
                "Poisson fitting requires ordered edges and finite nonnegative integer counts"
            )

    def window(self, low, high):
        tolerance = (
            16
            * np.finfo(float).eps
            * max(float(np.abs(self.edges).max()), float(np.ptp(self.edges)))
        )
        starts = np.flatnonzero(np.isclose(self.edges, low, rtol=0, atol=tolerance))
        stops = np.flatnonzero(np.isclose(self.edges, high, rtol=0, atol=tolerance))
        if len(starts) != 1 or len(stops) != 1 or starts[0] >= stops[0]:
            raise ValueError("Fit ranges must align with saved histogram edges")
        a, b = starts[0], stops[0]
        return Histogram(self.edges[a : b + 1], self.counts[a:b])


def poisson_deviance(counts, expected):
    expected = np.maximum(expected, 1e-12)
    return float(2 * np.sum(expected - counts + xlogy(counts, counts / expected)))


def fit_one(histogram, model_spec, window, settings):
    h = histogram.window(*window)
    model = PeakBackground(h.edges, **model_spec)
    counts = h.counts
    total = float(counts.sum())
    result = {
        "model": model_spec,
        "window": list(window),
        "entries": int(total),
        "accepted": False,
        "issues": [],
    }
    if total < settings["min_fit_entries"] or len(counts) <= model.parameter_count + 1:
        result["issues"] = ["insufficient_fit_entries"]
        return result
    center_low, center_high = settings["center_range"]
    width_low, width_high = settings["width_range"]
    if not window[0] < center_low < center_high < window[1]:
        raise ValueError("Center search range must be strictly inside every fit range")
    centers = (h.edges[:-1] + h.edges[1:]) / 2
    search = (centers >= center_low) & (centers <= center_high)
    smoothed = gaussian_filter1d(counts, 1)
    peak = float(centers[search][np.argmax(smoothed[search])])
    side = ~search
    baseline = max(
        0.01, float(np.median(counts[side])) if side.any() else float(counts.mean()) / 2
    )
    initial_background = max(1.0, baseline * len(counts))
    initial_signal = max(1.0, total - initial_background)
    upper_log = np.log(max(10.0, 10 * total))
    bounds = [
        (np.log(0.01), upper_log),
        ((center_low - window[0]) / model.span, (center_high - window[0]) / model.span),
        (np.log(width_low / model.span), np.log(width_high / model.span)),
    ] + [(np.log(0.01), upper_log)] * (model.degree + 1)

    def objective(p):
        signal, background, jac = model.components(p)
        expected = np.maximum(signal + background, 1e-12)
        nll = np.sum(expected - xlogy(counts, expected)) / total
        gradient = jac.T @ (1 - counts / expected) / total
        return float(nll), gradient

    best = None
    for seed in settings["width_seeds"]:
        p0 = np.array(
            [
                np.log(initial_signal),
                (peak - window[0]) / model.span,
                np.log(seed / model.span),
            ]
            + [np.log(initial_background / (model.degree + 1))] * (model.degree + 1)
        )
        solution = minimize(
            objective,
            p0,
            jac=True,
            method="L-BFGS-B",
            bounds=bounds,
            options={
                "maxiter": settings["max_iterations"],
                "ftol": 1e-12,
                "gtol": 1e-7,
                "maxls": 30,
            },
        )
        if best is None or solution.fun < best.fun:
            best = solution
    p = best.x
    signal, background, jac = model.components(p)
    expected = np.maximum(signal + background, 1e-12)
    center, width = window[0] + model.span * p[1], model.span * np.exp(p[2])
    deviance = max(0.0, poisson_deviance(counts, expected))
    ndf = len(counts) - model.parameter_count

    # Background-only model is a convex fit in its nonnegative yields. The
    # likelihood gain is a local contrast diagnostic, not a discovery Z-value.
    def null_objective(log_yields):
        components = model.basis * np.exp(log_yields)
        lam = np.maximum(components.sum(axis=1), 1e-12)
        return (
            float(np.sum(lam - xlogy(counts, lam)) / total),
            components.T @ (1 - counts / lam) / total,
        )

    null = minimize(
        null_objective,
        np.full(model.degree + 1, np.log(total / (model.degree + 1))),
        jac=True,
        method="L-BFGS-B",
        bounds=bounds[3:],
        options={"maxiter": settings["max_iterations"], "ftol": 1e-12, "gtol": 1e-7},
    )
    null_expected = model.basis @ np.exp(null.x)
    gain = max(0.0, poisson_deviance(counts, null_expected) - deviance)
    fisher = jac.T @ (jac / expected[:, None])
    center_error = width_error = signal_error = None
    try:
        condition = float(np.linalg.cond(fisher))
        covariance = np.linalg.inv(fisher)
        if (
            condition > settings["max_covariance_condition"]
            or (np.diag(covariance) <= 0).any()
            or not np.isfinite(covariance).all()
        ):
            raise np.linalg.LinAlgError("Ill-conditioned local covariance")
        center_error = float(model.span * np.sqrt(covariance[1, 1]))
        width_error = float(width * np.sqrt(covariance[2, 2]))
        signal_error = float(np.exp(p[0]) * np.sqrt(covariance[0, 0]))
    except np.linalg.LinAlgError:
        result["issues"].append("unreliable_covariance")
    if not best.success:
        result["issues"].append("optimizer_not_converged")
    if not null.success:
        result["issues"].append("background_only_not_converged")
    if deviance / ndf > settings["max_deviance_ndf"]:
        result["issues"].append("structured_residuals")
    if np.sqrt(gain) < settings["min_contrast"]:
        result["issues"].append("weak_peak")
    if (
        min(center - center_low, center_high - center)
        < settings["boundary_center_fraction"] * (center_high - center_low)
        or min(np.log(width / width_low), np.log(width_high / width))
        < settings["boundary_log_width_margin"]
    ):
        result["issues"].append("peak_at_search_boundary")
    if (
        center_error is not None
        and center_error / width > settings["max_center_error_over_width"]
    ):
        result["issues"].append("uncertain_center")
    if (
        width_error is not None
        and width_error / width > settings["max_relative_width_error"]
    ):
        result["issues"].append("uncertain_width")
    k, n = model.parameter_count, len(counts)
    result.update(
        {
            "center": float(center),
            "width": float(width),
            "center_error": center_error,
            "width_error": width_error,
            "signal_yield": float(np.exp(p[0])),
            "signal_error": signal_error,
            "parameters": p.tolist(),
            "deviance": deviance,
            "ndf": ndf,
            "deviance_ndf": deviance / ndf,
            "contrast": float(np.sqrt(gain)),
            "aicc": deviance + 2 * k + 2 * k * (k + 1) / (n - k - 1),
            "optimizer_message": str(best.message),
            "accepted": not result["issues"],
        }
    )
    return result


def fit_ensemble(histogram, recipe):
    """Compare AICc within identical fit domains; range changes are a separate
    sensitivity check. Select one actual fit; do not average correlated errors.
    """
    low, high = recipe["reference_band"]
    band = histogram.window(low, high)
    base = {
        "status": "insufficient",
        "issues": [],
        "reference_entries": int(band.counts.sum()),
        "fits": [],
        "chosen_fit": None,
    }
    if base["reference_entries"] < recipe["min_reference_entries"]:
        base["issues"] = ["insufficient_reference_entries"]
        return base
    winners, plausible = [], []
    for window_id, window in enumerate(recipe["fit_windows"]):
        indices = []
        for model in recipe["models"]:
            fit = fit_one(histogram, model, window, recipe["fit"])
            fit["window_id"] = window_id
            base["fits"].append(fit)
            if fit["accepted"]:
                indices.append(len(base["fits"]) - 1)
        if indices:
            winner = min(indices, key=lambda i: base["fits"][i]["aicc"])
            winners.append(winner)
            plausible.extend(
                i
                for i in indices
                if base["fits"][i]["aicc"]
                <= base["fits"][winner]["aicc"]
                + recipe["stability"]["max_aicc_difference"]
            )
    base["status"] = "review"
    if len(winners) < recipe["stability"]["min_windows"]:
        base["issues"] = ["too_few_acceptable_windows"]
        return base
    fits = [base["fits"][i] for i in plausible]
    centers = np.array([f["center"] for f in fits])
    widths = np.array([f["width"] for f in fits])
    base["center_range"] = [float(centers.min()), float(centers.max())]
    base["width_range"] = [float(widths.min()), float(widths.max())]
    base["winning_fits"] = winners
    base["plausible_fits"] = plausible
    stable = recipe["stability"]
    if np.ptp(centers) > min(
        stable["max_center_span"],
        stable["max_center_span_over_width"] * np.median(widths),
    ):
        base["issues"].append("range_sensitive_center")
    if widths.max() / widths.min() > stable["max_width_ratio"]:
        base["issues"].append("range_sensitive_width")
    preferred = recipe["preferred_window"]
    chosen = next(
        (i for i in winners if base["fits"][i]["window_id"] == preferred),
        min(
            winners,
            key=lambda i: abs(base["fits"][i]["center"] - np.median(centers))
            / np.median(widths)
            + abs(np.log(base["fits"][i]["width"] / np.median(widths))),
        ),
    )
    base["chosen_fit"] = chosen
    if not base["issues"]:
        base["status"] = "ready"
    return base
