"""Serializable cell-local windows and vectorized selection, without file I/O."""

import numpy as np


DECISION_REASONS = (
    "pass",
    "outside_grid",
    "nonfinite_observable",
    "unsupported_cell",
    "outside_window",
)


def window_from_fit(result, lower_widths, upper_widths):
    """Create a candidate window from one actual fit, preserving review status."""
    output = {
        "status": result["status"],
        "low": None,
        "high": None,
        "center": None,
        "width": None,
    }
    if result["chosen_fit"] is None:
        return output
    fit = result["fits"][result["chosen_fit"]]
    center, width = fit["center"], fit["width"]
    low, high = center - lower_widths * width, center + upper_widths * width
    output.update({"low": low, "high": high, "center": center, "width": width})
    if low < fit["window"][0] or high > fit["window"][1]:
        output["status"] = "review"
        result["issues"].append("cut_extends_beyond_fit_domain")
    return output


class WindowTable:
    """Unsupported cells have no applicable cut. No fallback or interpolation.

    A 'ready' status means the configured algorithmic gates passed, not that a
    physics analysis has validated this candidate selection.
    """

    def __init__(self, cells, size):
        self.low = np.full(size, np.nan)
        self.high = np.full(size, np.nan)
        self.center = np.full(size, np.nan)
        self.width = np.full(size, np.nan)
        seen = set()
        for cell in cells:
            i = cell["cell_id"]
            if type(i) is not int or not 0 <= i < size or i in seen:
                raise ValueError("Duplicate or invalid cell ID")
            seen.add(i)
            window = cell["window"]
            if window["status"] not in {"ready", "review", "insufficient"}:
                raise ValueError("Unknown window status")
            if window["status"] == "ready":
                values = [window[k] for k in ("low", "high", "center", "width")]
                if (
                    not np.isfinite(values).all()
                    or not values[0] < values[2] < values[1]
                    or values[3] <= 0
                ):
                    raise ValueError("Invalid ready window")
                self.low[i], self.high[i], self.center[i], self.width[i] = values
        if len(seen) != size:
            raise ValueError("Window table must explicitly account for every grid cell")

    def apply(self, cell_ids, observable, width_scale=1.0):
        cell_ids = np.asarray(cell_ids)
        observable = np.asarray(observable)
        if cell_ids.shape != observable.shape or not np.issubdtype(
            cell_ids.dtype, np.integer
        ):
            raise ValueError(
                "Cell IDs must be an integer array matching the observable shape"
            )
        if not np.isfinite(width_scale) or width_scale <= 0:
            raise ValueError("Width scale must be positive")
        reason = np.full(cell_ids.shape, 1, dtype=np.uint8)
        inside = (cell_ids >= 0) & (cell_ids < len(self.low))
        reason[inside] = 2
        finite = inside & np.isfinite(observable)
        reason[finite] = 3
        positions = np.flatnonzero(finite)
        ids = cell_ids[finite]
        supported = np.isfinite(self.low[ids])
        positions, ids = positions[supported], ids[supported]
        low = self.center[ids] + width_scale * (self.low[ids] - self.center[ids])
        high = self.center[ids] + width_scale * (self.high[ids] - self.center[ids])
        reason.flat[positions] = 4
        passed = (observable.flat[positions] >= low) & (
            observable.flat[positions] <= high
        )
        reason.flat[positions[passed]] = 0
        return reason == 0, reason
