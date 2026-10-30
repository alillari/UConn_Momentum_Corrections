"""Cached run-separated spectra and per-cell fit orchestration.

The cache separates data access from fitting. Recipes can change without
rereading the skim, and the fit engine sees only training spectra.
"""

from concurrent.futures import ProcessPoolExecutor
from collections import Counter
from datetime import datetime, timezone
import gzip
import json
import csv
import multiprocessing
from pathlib import Path

import numpy as np

from .binning import load_binning, bin_indices
from .cuts import window_from_fit
from .skim import coordinates, observables
from .provenance import code_state, sha256_file
from .root_skim import root_files, read_file, require_separate_output
from .channels import observable_values, validate_workflow


def build_cache(input_dir, diagnostics, binning, recipe, output):
    validate_workflow(diagnostics, binning, recipe)
    modern = recipe["recipe_schema_version"] == 2
    grids = load_binning(binning)
    output = require_separate_output(input_dir, output)
    files = root_files(input_dir)
    edges = np.asarray(
        recipe["spectrum"]["edges"] if modern else binning["mass_edges_gev"]
    )
    selected = {
        role: grids[role][binning["particles"][role]["selected_scheme"]]
        for role in recipe["roles"]
    }
    edge_key = "observable_edges" if modern else "mass_edges"
    outside_name = "observable_outside" if modern else "mass_outside"
    outside_accounting = (
        "invalid_or_outside_observable" if modern else "invalid_or_outside_mass"
    )
    arrays = {edge_key: edges}
    accounting = {
        role: {partition: Counter() for partition in ("training", "validation")}
        for role in selected
    }
    for role, grid in selected.items():
        for partition in ("training", "validation"):
            arrays[f"{role}_{partition}_spectra"] = np.zeros(
                (grid.size, len(edges) - 1), dtype=np.uint64
            )
            arrays[f"{role}_{partition}_counts"] = np.zeros(grid.size, dtype=np.uint64)
            arrays[f"{role}_{partition}_{outside_name}"] = np.zeros(
                grid.size, dtype=np.uint64
            )
    run_counts, manifest = Counter(), []
    expected_runs = set(recipe["training_runs"] + recipe["validation_runs"])
    for number, path in enumerate(files, 1):
        before = path.stat()
        file_entries = 0
        for start, data in read_file(
            path, diagnostics["dataset"]["tree"], diagnostics["branches"]
        ):
            file_entries += len(data["run"])
            runs, counts = np.unique(data["run"], return_counts=True)
            if set(runs) - expected_runs:
                raise ValueError(
                    "Input contains runs not assigned to a training/validation partition"
                )
            run_counts.update({str(int(r)): int(n) for r, n in zip(runs, counts)})
            values = observable_values(data, diagnostics, recipe)
            observable_bin, observable_ok = bin_indices(values, edges)
            for role, grid in selected.items():
                ids, reasons = grid.assign(
                    coordinates(data, role, binning["particles"][role])
                )
                for partition in ("training", "validation"):
                    group = np.isin(data["run"], recipe[f"{partition}_runs"])
                    inside = group & (ids >= 0)
                    good = inside & observable_ok
                    np.add.at(arrays[f"{role}_{partition}_counts"], ids[inside], 1)
                    np.add.at(
                        arrays[f"{role}_{partition}_spectra"],
                        (ids[good], observable_bin[good]),
                        1,
                    )
                    np.add.at(
                        arrays[f"{role}_{partition}_{outside_name}"],
                        ids[inside & ~observable_ok],
                        1,
                    )
                    accounting[role][partition].update(
                        {
                            "input": int(group.sum()),
                            "inside_grid": int(inside.sum()),
                            "outside_grid": int((group & ~inside).sum()),
                            "spectrum": int(good.sum()),
                            outside_accounting: int((inside & ~observable_ok).sum()),
                        }
                    )
        digest = sha256_file(path)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"Input changed during caching: {path}")
        manifest.append(
            {
                "path": str(path),
                "bytes": after.st_size,
                "entries": file_entries,
                "sha256": digest,
            }
        )
        if number % 50 == 0 or number == len(files):
            print(
                f"Cached {number}/{len(files)} files; {sum(run_counts.values()):,} entries",
                flush=True,
            )
    for role in selected:
        for partition in ("training", "validation"):
            a = accounting[role][partition]
            assert a["input"] == a["inside_grid"] + a["outside_grid"]
            assert a["inside_grid"] == a["spectrum"] + a[outside_accounting]
            assert np.array_equal(
                arrays[f"{role}_{partition}_spectra"].sum(axis=1)
                + arrays[f"{role}_{partition}_{outside_name}"],
                arrays[f"{role}_{partition}_counts"],
            )
    sources = [
        {"name": p.name, "sha256": sha256_file(p)}
        for p in sorted(Path(input_dir).iterdir())
        if p.is_file() and p.suffix in {".groovy", ".sh", ".dat", ".cshrc"}
    ]
    metadata = {
        "cache_schema_version": 2 if modern else 1,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "files": len(files),
        "entries": sum(run_counts.values()),
        "run_entries": dict(run_counts),
        "training_runs": recipe["training_runs"],
        "validation_runs": recipe["validation_runs"],
        "roles": recipe["roles"],
        "diagnostics": diagnostics,
        "binning": binning,
        "accounting": accounting,
        "code": code_state(),
        "preparation_files": sources,
        "completeness": {
            "snapshot_files": len(files),
            "expected_total": "Not independently established upstream",
        },
    }
    if modern:
        metadata.update(observable=recipe["observable"], spectrum=recipe["spectrum"])
    (output / "cache.json").write_text(
        json.dumps(metadata, indent=2, allow_nan=False) + "\n"
    )
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    np.savez_compressed(output / "histograms.npz", **arrays)
    return metadata


def _fit_task(task):
    from .peakfit import Histogram, fit_ensemble

    role, cell, edges, counts, ensemble, cut = task
    result = fit_ensemble(Histogram(edges, counts), ensemble)
    window = window_from_fit(result, cut["lower_widths"], cut["upper_widths"])
    result["status"] = window["status"]
    return role, cell, result, window


def derive(cache_dir, recipe, output, workers=1):
    if type(workers) is not int or workers < 1:
        raise ValueError("Worker count must be a positive integer")
    cache_dir, output = Path(cache_dir), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    meta = json.loads((cache_dir / "cache.json").read_text())
    validate_workflow(meta["diagnostics"], meta["binning"], recipe)
    modern = recipe["recipe_schema_version"] == 2
    if modern and (
        meta.get("observable") != recipe["observable"]
        or meta.get("spectrum") != recipe["spectrum"]
    ):
        raise ValueError("Observable/spectrum differs from cache; rebuild cache")
    for key in ("training_runs", "validation_runs", "roles"):
        if recipe[key] != meta[key]:
            raise ValueError("Recipe partition/roles differ from cache; rebuild cache")
    grids = load_binning(meta["binning"])
    tasks = []
    with np.load(cache_dir / "histograms.npz") as cache:
        edges = cache["observable_edges" if modern else "mass_edges"]
        for role in recipe["roles"]:
            h = cache[f"{role}_training_spectra"]
            tasks.extend(
                (role, i, edges, row, recipe["ensemble"], recipe["cut"])
                for i, row in enumerate(h)
            )
    records = {role: [] for role in recipe["roles"]}
    cells = {role: [] for role in recipe["roles"]}
    executor = None
    if workers > 1:
        executor = ProcessPoolExecutor(
            max_workers=workers, mp_context=multiprocessing.get_context("spawn")
        )
        iterator = executor.map(_fit_task, tasks, chunksize=16)
    else:
        iterator = map(_fit_task, tasks)
    try:
        for n, (role, cell, result, window) in enumerate(iterator, 1):
            grid = grids[role][meta["binning"]["particles"][role]["selected_scheme"]]
            records[role].append({"cell_id": cell, **result})
            chosen = (
                result["fits"][result["chosen_fit"]]
                if result["chosen_fit"] is not None
                else None
            )
            cells[role].append(
                {
                    "cell_id": cell,
                    "bounds": grid.bounds(cell),
                    "window": window,
                    "issues": result["issues"],
                    "training_reference_entries": result["reference_entries"],
                    "chosen_model": chosen["model"] if chosen else None,
                    "center_error": chosen["center_error"] if chosen else None,
                    "width_error": chosen["width_error"] if chosen else None,
                    "center_range": result.get("center_range"),
                    "width_range": result.get("width_range"),
                }
            )
            if n % 1000 == 0 or n == len(tasks):
                print(f"Processed {n:,}/{len(tasks):,} cell spectra", flush=True)
    finally:
        if executor is not None:
            executor.shutdown()
    selection = {
        "selection_schema_version": 2 if modern else 1,
        "status": "candidate_not_physics_validated",
        "observable": recipe["observable"],
        "recipe": recipe,
        "binning": meta["binning"],
        "diagnostics": meta["diagnostics"],
        "cells": cells,
        "source_cache_sha256": sha256_file(cache_dir / "cache.json"),
        "source_histograms_sha256": sha256_file(cache_dir / "histograms.npz"),
        "code": code_state(),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    (output / "cuts.json").write_text(
        json.dumps(selection, indent=2, allow_nan=False) + "\n"
    )
    for role, rows in cells.items():
        with (output / f"{role}_windows.csv").open("w", newline="") as stream:
            writer = csv.DictWriter(
                stream,
                fieldnames=[
                    "cell_id",
                    "sector",
                    "status",
                    "low",
                    "high",
                    "center",
                    "width",
                    "center_error",
                    "width_error",
                    "training_reference_entries",
                    "background_degree",
                    "issues",
                ],
            )
            writer.writeheader()
            for row in rows:
                writer.writerow(
                    {
                        "cell_id": row["cell_id"],
                        "sector": row["bounds"]["sector"],
                        **row["window"],
                        "center_error": row["center_error"],
                        "width_error": row["width_error"],
                        "training_reference_entries": row["training_reference_entries"],
                        "background_degree": (
                            row["chosen_model"]["degree"]
                            if row["chosen_model"]
                            else None
                        ),
                        "issues": ";".join(row["issues"]),
                    }
                )
    with gzip.open(output / "fits.json.gz", "wt") as stream:
        json.dump(records, stream, allow_nan=False)
    summary = {
        role: dict(Counter(c["window"]["status"] for c in rows))
        for role, rows in cells.items()
    }
    (output / "fit_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print("Candidate windows:", summary, flush=True)
    return selection
