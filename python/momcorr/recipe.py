"""Strict, detector-independent fit/cut recipe validation."""

import json

import numpy as np
from .models import SIGNALS, BACKGROUNDS
from .channels import TRACKS, OBSERVABLE_UNITS, validate_observable


def exact_keys(value, expected, label):
    if not isinstance(value, dict) or set(value) != set(expected):
        raise ValueError(
            f"Unexpected or missing keys in {label}: expected {sorted(expected)}"
        )


def load_recipe(path):
    r = json.loads(path.read_text())
    modern = r.get("recipe_schema_version") == 2
    exact_keys(
        r,
        {
            "recipe_schema_version",
            "name",
            "observable",
            "roles",
            "training_runs",
            "validation_runs",
            "ensemble",
            "cut",
            "annotations",
        }
        | ({"spectrum"} if modern else set()),
        "recipe",
    )
    if not isinstance(r["annotations"], dict) or any(
        not isinstance(k, str) or not isinstance(v, str)
        for k, v in r["annotations"].items()
    ):
        raise ValueError("Recipe annotations must be a string-to-string mapping")
    if modern:
        validate_observable(r["observable"])
        spectrum = r["spectrum"]
        exact_keys(spectrum, {"edges", "label", "unit", "reference_value"}, "spectrum")
        edges = np.asarray(spectrum["edges"], dtype=float)
        if (
            edges.ndim != 1
            or len(edges) < 3
            or not np.isfinite(edges).all()
            or not (np.diff(edges) > 0).all()
        ):
            raise ValueError("Spectrum requires finite increasing bin edges")
        if not isinstance(spectrum["label"], str) or not spectrum["label"].strip():
            raise ValueError("Spectrum label is required")
        if spectrum["unit"] != OBSERVABLE_UNITS[r["observable"]["kind"]]:
            raise ValueError("Spectrum unit does not match observable")
        if spectrum["reference_value"] is not None and not np.isfinite(
            spectrum["reference_value"]
        ):
            raise ValueError("Reference value must be finite or null")
    elif r["recipe_schema_version"] != 1 or r["observable"] != "missing_mass":
        raise ValueError("Unsupported recipe schema/observable adapter")
    if (
        not r["roles"]
        or len(set(r["roles"])) != len(r["roles"])
        or set(r["roles"]) - (TRACKS if modern else {"electron", "pion"})
    ):
        raise ValueError("Invalid or duplicated particle roles")
    for partition in ("training_runs", "validation_runs"):
        if (
            not r[partition]
            or len(set(r[partition])) != len(r[partition])
            or any(type(v) is not int or v < 1 for v in r[partition])
        ):
            raise ValueError("Run partitions require unique positive run IDs")
    if set(r["training_runs"]) & set(r["validation_runs"]):
        raise ValueError("Training and validation runs overlap")
    e = r["ensemble"]
    exact_keys(
        e,
        {
            "reference_band",
            "min_reference_entries",
            "fit_windows",
            "models",
            "fit",
            "stability",
            "preferred_window",
        },
        "ensemble",
    )
    f = e["fit"]
    exact_keys(
        f,
        {
            "center_range",
            "width_range",
            "width_seeds",
            "min_fit_entries",
            "max_iterations",
            "max_deviance_ndf",
            "min_contrast",
            "max_center_error_over_width",
            "max_relative_width_error",
            "max_covariance_condition",
            "boundary_center_fraction",
            "boundary_log_width_margin",
        },
        "fit settings",
    )
    st = e["stability"]
    exact_keys(
        st,
        {
            "min_windows",
            "max_center_span",
            "max_center_span_over_width",
            "max_width_ratio",
            "max_aicc_difference",
        },
        "stability",
    )
    exact_keys(
        r["cut"], {"lower_widths", "upper_widths", "sensitivity_scales"}, "cut policy"
    )
    for pair in [e["reference_band"], f["center_range"], f["width_range"]] + e[
        "fit_windows"
    ]:
        if len(pair) != 2 or not np.isfinite(pair).all() or pair[0] >= pair[1]:
            raise ValueError("Invalid range")
    if modern:
        for low, high in [e["reference_band"]] + e["fit_windows"]:
            if low < edges[0] or high > edges[-1]:
                raise ValueError("Fit/reference ranges must lie in the stored spectrum")
    if (
        f["width_range"][0] <= 0
        or not f["width_seeds"]
        or any(
            not f["width_range"][0] < x < f["width_range"][1] for x in f["width_seeds"]
        )
    ):
        raise ValueError("Width seeds must lie inside positive width bounds")
    for key in ("min_fit_entries", "max_iterations"):
        if type(f[key]) is not int or f[key] < 1:
            raise ValueError("Invalid integer fit setting")
    if type(e["min_reference_entries"]) is not int or e["min_reference_entries"] < 1:
        raise ValueError("Invalid minimum reference entries")
    for v in (
        [
            f[k]
            for k in (
                "max_deviance_ndf",
                "min_contrast",
                "max_center_error_over_width",
                "max_relative_width_error",
                "max_covariance_condition",
                "boundary_center_fraction",
                "boundary_log_width_margin",
            )
        ]
        + [
            st[k]
            for k in (
                "max_center_span",
                "max_center_span_over_width",
                "max_width_ratio",
            )
        ]
        + [r["cut"][k] for k in ("lower_widths", "upper_widths")]
        + r["cut"]["sensitivity_scales"]
    ):
        if not np.isfinite(v) or v <= 0:
            raise ValueError("Quality/cut settings must be finite and positive")
    if not e["models"] or len(e["fit_windows"]) != len(
        {tuple(w) for w in e["fit_windows"]}
    ):
        raise ValueError("Empty model family or duplicate fit ranges")
    if type(e["preferred_window"]) is not int or not 0 <= e["preferred_window"] < len(
        e["fit_windows"]
    ):
        raise ValueError("Invalid preferred fit-window index")
    if type(st["min_windows"]) is not int or not 1 <= st["min_windows"] <= len(
        e["fit_windows"]
    ):
        raise ValueError("Invalid minimum acceptable windows")
    for model in e["models"]:
        exact_keys(model, {"signal", "background", "degree"}, "model")
        if (
            model["signal"] not in SIGNALS
            or model["background"] not in BACKGROUNDS
            or type(model["degree"]) is not int
            or not 0 <= model["degree"] <= 4
        ):
            raise ValueError("Unknown registered model or invalid degree")
    if (
        not np.isfinite(st["max_aicc_difference"])
        or st["max_aicc_difference"] < 0
        or st["max_width_ratio"] < 1
    ):
        raise ValueError("Invalid model-stability setting")
    if not r["cut"]["sensitivity_scales"] or 1.0 not in r["cut"]["sensitivity_scales"]:
        raise ValueError("Sensitivity scales must include the nominal scale 1")
    return r
