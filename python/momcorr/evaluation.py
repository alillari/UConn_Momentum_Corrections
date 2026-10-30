"""Apply cell windows, export ROOT friend masks, and measure held-out behavior."""

from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path

import numpy as np
import uproot

from .binning import load_binning
from .skim import coordinates, observables
from .cuts import WindowTable, DECISION_REASONS
from .provenance import sha256_file, code_state
from .root_skim import root_files, read_file, require_separate_output
from .channels import observable_values, validate_workflow


def apply_selection(input_dir, selection, output, selection_path, write_masks=True):
    if selection.get("selection_schema_version") not in (1, 2):
        raise ValueError("Unknown selection schema")
    if selection != json.loads(Path(selection_path).read_text()):
        raise ValueError("Selection object differs from its provenance file")
    output = require_separate_output(input_dir, output)
    roles = list(selection["cells"])
    grids = load_binning(selection["binning"])
    selected = {
        role: grids[role][selection["binning"]["particles"][role]["selected_scheme"]]
        for role in roles
    }
    tables = {
        role: WindowTable(selection["cells"][role], selected[role].size)
        for role in roles
    }
    diagnostics = selection["diagnostics"]
    recipe = selection["recipe"]
    validate_workflow(diagnostics, selection["binning"], recipe)
    modern = recipe["recipe_schema_version"] == 2
    if selection["selection_schema_version"] != recipe["recipe_schema_version"]:
        raise ValueError("Selection and recipe schemas differ")
    digest = sha256_file(selection_path)
    spectrum_edges = np.array(
        recipe["spectrum"]["edges"]
        if modern
        else selection["binning"]["mass_edges_gev"]
    )
    delta_edges = np.linspace(-0.5, 0.5, 251)
    primary = "observable" if modern else "mass"
    axes = (
        {primary: spectrum_edges}
        if modern
        else {"mass": spectrum_edges, "delta": delta_edges}
    )
    hist = {key + "_edges": edges for key, edges in axes.items()}
    losses, counters, by_run, sensitivity = (
        {},
        {r: {g: Counter() for g in ("training", "validation", "other")} for r in roles},
        {},
        {r: {} for r in roles},
    )
    manifest = []
    friend_dir = output / "masks"
    if write_masks:
        friend_dir.mkdir(exist_ok=True)

    def fill(key, values, edges):
        finite = np.isfinite(values)
        counts = np.histogram(values[finite], edges)[0].astype(np.uint64)
        if key not in hist:
            hist[key] = np.zeros(len(edges) - 1, dtype=np.uint64)
            losses[key] = Counter()
        hist[key] += counts
        losses[key].update(
            {
                "input": len(values),
                "nonfinite": int((~finite).sum()),
                "outside": int(finite.sum() - counts.sum()),
            }
        )

    files = root_files(input_dir)
    for number, path in enumerate(files, 1):
        before = path.stat()
        rows = 0
        target = friend_dir / (path.stem + ".cuts.root")
        temporary = target.with_suffix(".root.part")
        friend = None
        if write_masks:
            friend = uproot.recreate(temporary)
            schema = {"entry": np.int64, "run": np.int32, "evn": np.int32}
            for role in roles:
                schema.update(
                    {
                        role + "_cell": np.int64,
                        role + "_pass": np.bool_,
                        role + "_reason": np.uint8,
                        role + "_grid_reason": np.uint8,
                    }
                )
            friend.mktree("cuts", schema)
        try:
            for start, data in read_file(
                path, diagnostics["dataset"]["tree"], diagnostics["branches"]
            ):
                n = len(data["run"])
                rows += n
                if modern:
                    values = observable_values(data, diagnostics, recipe)
                    variables = {"observable": values}
                else:
                    derived = observables(data, diagnostics["conditions"])
                    values = derived[recipe["observable"]]
                    variables = {"mass": values, "delta": derived["delta_pion_p"]}
                group = {
                    "training": np.isin(data["run"], recipe["training_runs"]),
                    "validation": np.isin(data["run"], recipe["validation_runs"]),
                }
                group["other"] = ~(group["training"] | group["validation"])
                branches = {
                    "entry": np.arange(start, start + n, dtype=np.int64),
                    "run": data["run"].astype(np.int32),
                    "evn": data["event"].astype(np.int32),
                }
                for role, grid in selected.items():
                    ids, grid_reasons = grid.assign(
                        coordinates(data, role, selection["binning"]["particles"][role])
                    )
                    passed, reasons = tables[role].apply(ids, values)
                    branches.update(
                        {
                            role + "_cell": ids,
                            role + "_pass": passed,
                            role + "_reason": reasons,
                            role + "_grid_reason": grid_reasons,
                        }
                    )
                    for name, mask in group.items():
                        c = counters[role][name]
                        c["input"] += int(mask.sum())
                        for code, label in enumerate(DECISION_REASONS):
                            c[label] += int((mask & (reasons == code)).sum())
                        for variable, variable_values in variables.items():
                            for stage, active in (
                                ("before", mask),
                                ("after", mask & passed),
                            ):
                                fill(
                                    f"{role}_{name}_{stage}_{variable}",
                                    variable_values[active],
                                    axes[variable],
                                )
                    for run in np.unique(data["run"]):
                        key = f"{role}_{int(run)}"
                        c = by_run.setdefault(key, Counter())
                        mask = data["run"] == run
                        c["input"] += int(mask.sum())
                        c["supported"] += int((mask & np.isin(reasons, [0, 4])).sum())
                        c["pass"] += int((mask & passed).sum())
                    for scale in recipe["cut"]["sensitivity_scales"]:
                        active, _ = tables[role].apply(ids, values, width_scale=scale)
                        key = f"{scale:.8g}"
                        total = sensitivity[role].setdefault(key, Counter())
                        total["validation_pass"] += int(
                            (active & group["validation"]).sum()
                        )
                        fill(
                            f"{role}_sensitivity_{key}_{primary if modern else 'delta'}",
                            variables[primary if modern else "delta"][
                                active & group["validation"]
                            ],
                            axes[primary if modern else "delta"],
                        )
                if friend is not None:
                    friend["cuts"].extend(branches)
            source_hash = sha256_file(path)
            after = path.stat()
            if (before.st_size, before.st_mtime_ns) != (
                after.st_size,
                after.st_mtime_ns,
            ):
                raise ValueError(f"Input changed during selection: {path}")
            if friend is not None:
                friend["provenance"] = json.dumps(
                    {
                        "source": str(path),
                        "source_sha256": source_hash,
                        "selection_sha256": digest,
                        "entries": rows,
                        "role_combination": "Each role is evaluated independently; there is no implicit AND/OR.",
                        "decision_codes": DECISION_REASONS,
                    }
                )
        finally:
            if friend is not None:
                friend.close()
        if write_masks:
            temporary.replace(target)
        manifest.append(
            {
                "source": str(path),
                "source_sha256": source_hash,
                "entries": rows,
                "mask": str(target) if write_masks else None,
            }
        )
        if number % 50 == 0 or number == len(files):
            print(f"Applied windows to {number}/{len(files)} files", flush=True)
    for groups in counters.values():
        for c in groups.values():
            assert c["input"] == sum(c[k] for k in DECISION_REASONS)
    for key, c in losses.items():
        assert c["input"] == int(hist[key].sum()) + c["nonfinite"] + c["outside"]
    summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "selection_sha256": digest,
        "selection_status": selection["status"],
        "roles_independent": True,
        "decisions": counters,
        "by_run": by_run,
        "width_sensitivity": sensitivity,
        "histogram_accounting": losses,
        "code": code_state(),
        "files": len(files),
        "entries": sum(r["entries"] for r in manifest),
        "interpretation": "Pass fractions are retention, not signal efficiency/purity. Before/after spectra of the selected observable are selection diagnostics, not independent physics validation.",
    }
    (output / "evaluation.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n"
    )
    (output / "mask_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    np.savez_compressed(output / "evaluation_histograms.npz", **hist)
    return summary
