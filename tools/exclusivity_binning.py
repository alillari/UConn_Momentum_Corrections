#!/usr/bin/env python3
"""Explore, count, and export sector-separated 3D exclusivity-study cells.

No fits or exclusivity cuts are applied. The reference mass band is only an
occupancy diagnostic, contains background, and is not a signal/purity estimate.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from momcorr.binning import AXES, REASONS, bin_indices, load_binning, shifted_coordinates
import exclusivity_diagnostics as diag


def coordinate_chunk(a, role, spec):
    return shifted_coordinates(*[a[f"{role}_p{k}"] for k in "xyz"], a[f"{role}_sector"],
        spec["shift"]["coefficient_deg_gev"], spec["shift"]["pole_gev"])


def build(input_dir, diagnostics_config, config, output):
    grids = load_binning(config)
    files = sorted(Path(input_dir).resolve().glob("*.root"))
    if not files:
        raise ValueError("No ROOT files")
    out = Path(output).resolve()
    if out == Path(input_dir).resolve() or Path(input_dir).resolve() in out.parents:
        raise ValueError("Output must be outside the input directory")
    out.mkdir(parents=True, exist_ok=True)
    mass_edges = np.array(config["mass_edges_gev"])
    n_mass = len(mass_edges) - 1
    data, accounting, runs = {}, {}, set()
    for role, schemes in grids.items():
        for name, grid in schemes.items():
            key = f"{role}_{name}"
            data[f"{key}_counts"] = np.zeros(grid.size, dtype=np.uint64)
            data[f"{key}_reference"] = np.zeros(grid.size, dtype=np.uint64)
            accounting[key] = Counter()
        selected = config["particles"][role]["selected_scheme"]
        grid = schemes[selected]
        data[f"{role}_spectra"] = np.zeros((grid.size, n_mass), dtype=np.uint64)
        data[f"{role}_mass_outside"] = np.zeros(grid.size, dtype=np.uint64)
        for coordinate in ("local_phi", "shifted_phi"):
            data[f"{role}_{coordinate}_projection"] = np.zeros((140, 240), dtype=np.uint64)
    data["projection_p_edges"] = np.linspace(0, 7, 141)
    data["projection_phi_edges"] = np.linspace(-120, 120, 241)
    data["mass_edges"] = mass_edges
    manifest, per_run, projection_losses, extrema = [], Counter(), Counter(), {}
    extreme_records = {}
    reference_low, reference_high = config["reference_band_gev"]
    for number, path in enumerate(files, 1):
        before = path.stat()
        digest = diag.sha256_file(path)
        file_entries = 0
        with uproot.open(path, handler=uproot.source.file.MemmapSource,
                         decompression_executor=uproot.source.futures.TrivialExecutor(),
                         interpretation_executor=uproot.source.futures.TrivialExecutor()) as f:
            tree = f[diagnostics_config["dataset"]["tree"]]
            mapping = diagnostics_config["branches"]
            if set(mapping.values()) - set(tree.keys()):
                raise ValueError(f"Missing input branches: {path}")
            for raw in tree.iterate(list(mapping.values()), step_size="40 MB", library="np"):
                a = {k: raw[v] for k, v in mapping.items()}
                file_entries += len(a["run"])
                finite = np.logical_and.reduce([np.isfinite(a[k]) for k in diag.FIELDS[2:8]])
                mm = np.full(len(a["run"]), np.nan)
                mm[finite] = diag.kinematics({k: v[finite] for k, v in a.items()}, diagnostics_config["conditions"])["missing_mass"]
                reference = np.isfinite(mm) & (mm >= reference_low) & (mm < reference_high)
                mass_bin, mass_ok = bin_indices(mm, mass_edges)
                unique_runs, rc = np.unique(a["run"], return_counts=True)
                runs.update(int(x) for x in unique_runs)
                per_run.update({str(int(r)): int(n) for r, n in zip(unique_runs, rc)})
                for role, schemes in grids.items():
                    coords = coordinate_chunk(a, role, config["particles"][role])
                    for axis in AXES:
                        valid = coords["reason"] == 0
                        values = coords[axis][valid]
                        if len(values):
                            key = f"{role}_{axis}"
                            previous = extrema.setdefault(key, [float("inf"), -float("inf")])
                            records = extreme_records.setdefault(key, {})
                            indices = np.flatnonzero(valid)
                            for side, offset, better in (("minimum", 0, np.less), ("maximum", 1, np.greater)):
                                position = int(indices[np.argmin(values) if offset == 0 else np.argmax(values)])
                                value = float(coords[axis][position])
                                if better(value, previous[offset]):
                                    previous[offset] = value
                                    records[side] = {"value": value, "file": str(path), "run": int(a["run"][position]),
                                        "event": int(a["event"][position]), "sector": float(coords["sector"][position]),
                                        "px_py_pz": [float(a[f"{role}_p{k}"][position]) for k in "xyz"]}
                    for coordinate in ("local_phi", "shifted_phi"):
                        key = f"{role}_{coordinate}_projection"
                        ok = (coords["reason"] == 0) & np.isfinite(coords[coordinate])
                        h = np.histogram2d(coords["p"][ok], coords[coordinate][ok],
                            bins=(data["projection_p_edges"], data["projection_phi_edges"]))[0].astype(np.uint64)
                        data[key] += h
                        projection_losses[key] += int(len(mm) - h.sum())
                    for name, grid in schemes.items():
                        key = f"{role}_{name}"
                        cells, reasons = grid.assign(coords)
                        for code, label in enumerate(REASONS):
                            accounting[key][label] += int((reasons == code).sum())
                        inside = cells >= 0
                        np.add.at(data[f"{key}_counts"], cells[inside], 1)
                        np.add.at(data[f"{key}_reference"], cells[inside & reference], 1)
                        if name == config["particles"][role]["selected_scheme"]:
                            good = inside & mass_ok
                            np.add.at(data[f"{role}_spectra"], (cells[good], mass_bin[good]), 1)
                            np.add.at(data[f"{role}_mass_outside"], cells[inside & ~mass_ok], 1)
                            for run in unique_runs:
                                for suffix, mask in (("counts", inside), ("reference", inside & reference)):
                                    rkey = f"{role}_run_{int(run)}_{suffix}"
                                    if rkey not in data:
                                        data[rkey] = np.zeros(grid.size, dtype=np.uint64)
                                    active = mask & (a["run"] == run)
                                    np.add.at(data[rkey], cells[active], 1)
        after = path.stat()
        if (before.st_size, before.st_mtime_ns) != (after.st_size, after.st_mtime_ns):
            raise ValueError(f"Input changed during processing: {path}")
        manifest.append({"path": str(path), "bytes": before.st_size, "sha256": digest, "entries": file_entries})
        if number % 50 == 0 or number == len(files):
            print(f"Binned {number}/{len(files)} files; {sum(per_run.values()):,} entries", flush=True)
    entries = sum(per_run.values())
    summaries = {}
    for role, schemes in grids.items():
        summaries[role] = {}
        for name, grid in schemes.items():
            key = f"{role}_{name}"
            counts, reference = data[f"{key}_counts"], data[f"{key}_reference"]
            occupied = counts > 0
            inside = int(counts.sum())
            reference_total = int(reference.sum())
            assert sum(accounting[key].values()) == entries
            assert inside == accounting[key]["inside"]
            thresholds = {}
            for threshold in config["occupancy_thresholds"]:
                enough = reference >= threshold
                thresholds[str(threshold)] = {
                    "cells": int(enough.sum()),
                    "fraction_of_inside_entries": float(counts[enough].sum() / inside) if inside else 0.,
                    "fraction_of_inside_reference_entries": float(reference[enough].sum() / reference_total) if reference_total else 0.}
            summaries[role][name] = {"shape": grid.shape, "total_cells": grid.size,
                "occupied_cells": int(occupied.sum()), "inside": inside,
                "inside_fraction": inside / entries if entries else 0., "reference_entries": reference_total,
                "reference_count_quantiles_in_occupied_cells": np.percentile(reference[occupied], [0, 25, 50, 75, 100]).tolist() if occupied.any() else [],
                "reason_counts": dict(accounting[key]), "support": thresholds}
        selected = config["particles"][role]["selected_scheme"]
        grid = schemes[selected]
        counts, reference = data[f"{role}_{selected}_counts"], data[f"{role}_{selected}_reference"]
        assert np.array_equal(data[f"{role}_spectra"].sum(axis=1) + data[f"{role}_mass_outside"], counts)
        with (out / f"{role}_cells.csv").open("w", newline="") as stream:
            fields = ["cell_id", "sector"] + [f"{axis}_{side}" for axis in AXES for side in ("low", "high")] + ["entries", "reference_band_entries", "spectrum_entries", "spectrum_outside"] + [f"run_{r}_entries" for r in sorted(runs)]
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            for cell in range(grid.size):
                bounds = grid.bounds(cell)
                row = {"cell_id": cell, "sector": bounds["sector"], "entries": int(counts[cell]),
                       "reference_band_entries": int(reference[cell]), "spectrum_entries": int(data[f"{role}_spectra"][cell].sum()),
                       "spectrum_outside": int(data[f"{role}_mass_outside"][cell])}
                for axis in AXES:
                    row[f"{axis}_low"], row[f"{axis}_high"] = bounds[axis]
                for run in sorted(runs):
                    row[f"run_{run}_entries"] = int(data[f"{role}_run_{run}_counts"][cell])
                writer.writerow(row)
    repo = Path(__file__).resolve().parents[1]
    code = {"commit": subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip(),
        "working_tree_status": subprocess.check_output(["git", "-C", str(repo), "status", "--porcelain"], text=True),
        "frontend_sha256": diag.sha256_file(__file__), "binning_sha256": diag.sha256_file(repo / "python/momcorr/binning.py")}
    summary = {"generated_at_utc": datetime.now(timezone.utc).isoformat(), "entries": entries, "files": len(files),
        "run_entries": dict(per_run), "schemes": summaries, "coordinate_extrema": extrema,
        "coordinate_extreme_examples": extreme_records,
        "projection_excluded_entries": dict(projection_losses), "binning_config": config,
        "diagnostics_config": diagnostics_config, "code": code, "versions": {"python": sys.version, "numpy": np.__version__, "uproot": uproot.__version__},
        "definitions": {"reference_band": "Diagnostic count including background; not an exclusivity cut or a signal-yield estimate.",
            "axis_order": ["sector", "p", "theta", "shifted_phi"], "cell_id": "C-order flat index; use exact saved edges and sector order.",
            "edges": "Half-open intervals; final right edge included. No wrapping after the inverse-momentum phi shift.",
            "selection": "No physics selection, fits, or corrections applied. Rows outside a grid are explicitly counted, not selected away."}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (out / "binning.json").write_text(json.dumps(config, indent=2) + "\n")
    np.savez_compressed(out / "cells.npz", **data)
    return summary, data, grids


def report(summary, data, grids, output):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LogNorm
    import textwrap

    out = Path(output)
    config = summary["binning_config"]
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "pdf.fonttype": 42})
    page = 0
    def figure(title, rows=1, columns=1, note="Sector-separated 3D cells; no fits or exclusivity cuts."):
        fig, ax = plt.subplots(rows, columns, figsize=(11.7, 8.3), squeeze=False)
        fig.suptitle(title, x=.065, y=.97, ha="left", fontsize=17, weight="bold", color="#183d52")
        fig.text(.065, .925, note, va="top", fontsize=9, color="#526571")
        fig.subplots_adjust(left=.08, right=.91, top=.85, bottom=.10, hspace=.44, wspace=.30)
        return fig, ax
    def save(pdf, fig, name):
        nonlocal page
        page += 1
        fig.text(.95, .035, str(page), ha="right", color="#526571")
        pdf.savefig(fig)
        fig.savefig(out / f"{page:02d}_{name}.png", dpi=140)
        plt.close(fig)
    def heat(ax, x, y, counts, vmax=None):
        vmax = max(2, int(counts.max())) if vmax is None else vmax
        return ax.pcolormesh(x, y, np.ma.masked_equal(counts.T, 0), cmap="viridis", norm=LogNorm(1, vmax), rasterized=True)
    def occupancy_page(pdf, role, by_sector=False):
        name = config["particles"][role]["selected_scheme"]
        grid = grids[role][name]
        h = data[f"{role}_{name}_counts"].reshape(grid.shape)
        reference = data[f"{role}_{name}_reference"].reshape(grid.shape)
        if by_sector:
            number = len(grid.sectors)
            title = f"{role.capitalize()}: reference-band population by sector"
            note = "Sum over shifted phi; reference mass band 0.85-1.05 GeV contains background and is not a cut."
            slices = [reference[i].sum(axis=-1) for i in range(number)]
            labels = [f"Sector {s}" for s in grid.sectors]
        else:
            number = grid.shape[-1]
            title = f"{role.capitalize()}: actual 3D cells across shifted-phi slices"
            note = "Each panel is a shifted-phi slice; displayed p/theta cells are summed over sectors only for this overview."
            slices = [h[..., i].sum(axis=0) for i in range(number)]
            labels = [f"{grid.edges['shifted_phi'][i]:g} <= phi shift {'<=' if i == number - 1 else '<'} {grid.edges['shifted_phi'][i+1]:g} deg" for i in range(number)]
        columns = 3 if number <= 6 else 4
        rows = int(np.ceil(number / columns))
        fig, axes = figure(title, rows, columns, note)
        fig.subplots_adjust(left=.07, right=.87, wspace=.32)
        vmax = max(2, max(int(s.max()) for s in slices))
        for i, ax in enumerate(axes.flat):
            if i >= number:
                ax.axis("off")
                continue
            m = heat(ax, grid.edges["p"], grid.edges["theta"], slices[i], vmax)
            ax.set_title(labels[i], fontsize=9)
            ax.set_xlabel("|p| [GeV/c]")
            if i % columns == 0:
                ax.set_ylabel("theta [deg]")
            for edge in grid.edges["p"]:
                ax.axvline(edge, color="white", alpha=.12, lw=.3)
            for edge in grid.edges["theta"]:
                ax.axhline(edge, color="white", alpha=.12, lw=.3)
        fig.colorbar(m, cax=fig.add_axes([.91, .10, .015, .75]), label="Entries / displayed cell (shared scale)")
        save(pdf, fig, role + ("_sectors" if by_sector else "_phi_slices"))

    with PdfPages(out / "binning_report.pdf") as pdf:
        fig, axes = figure("3D binning for the RG-K exclusivity study", note="Coordinates and bin populations established before fitting peaks or defining cuts.")
        ax = axes[0, 0]
        ax.axis("off")
        lines = [f"{summary['files']:,} files; {summary['entries']:,} entries; sectors remain separate.",
                 "Axes: |p| [GeV/c], theta [deg], shifted phi [deg].", "",
                 "phi_local = sector-unwrapped atan2(py, px) - (sector - 1) * 60 deg",
                 "electron: phi_shift = phi_local - 30/p",
                 "pi+:      phi_shift = phi_local + 32/(p - 0.05)",
                 "No wrapping after the shift; momenta are preserved. User-confirmed conventions.", ""]
        for role in ("electron", "pion"):
            name = config["particles"][role]["selected_scheme"]
            s = summary["schemes"][role][name]
            grid = grids[role][name]
            lines += [f"{role.capitalize()} / {name}: {s['occupied_cells']:,} occupied of {s['total_cells']:,} cells.",
                      f"Grid coverage: {100*s['inside_fraction']:.4f}% of input entries."]
            for axis in AXES:
                edges = ", ".join(f"{v:g}" for v in grid.edges[axis])
                lines += textwrap.wrap(f"  {axis}: {edges}", 104, subsequent_indent="    ")
            lines.append("")
        ax.text(0, 1, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=10, family="monospace", linespacing=1.4)
        save(pdf, fig, "scheme")

        fig, axes = figure("Check the shifted-phi coordinates", 2, 2,
            "All input entries with valid coordinates; display range: 0-7 GeV/c and -120 to +120 deg.")
        for row, role in enumerate(("electron", "pion")):
            for column, coord in enumerate(("local_phi", "shifted_phi")):
                ax = axes[row, column]
                key = f"{role}_{coord}_projection"
                m = heat(ax, data["projection_p_edges"], data["projection_phi_edges"], data[key])
                ax.figure.colorbar(m, ax=ax, label="Entries / bin", pad=.02)
                ax.set_xlabel("|p| [GeV/c]")
                ax.set_ylabel(coord.replace("_", " ") + " [deg]")
                ax.set_title(f"{role.capitalize()} | {summary['projection_excluded_entries'][key]:,} outside display or invalid")
                if coord == "shifted_phi":
                    selected = config["particles"][role]["selected_scheme"]
                    edges = grids[role][selected].edges["shifted_phi"]
                    for edge in edges:
                        ax.axhline(edge, color="white", lw=.5, alpha=.5)
        save(pdf, fig, "coordinates")

        fig, axes = figure("Resolution versus per-cell statistics", 2, 2,
            "Statistics refer to the diagnostic mass band, including background; these thresholds are not fit-quality decisions.")
        colors = ["#176b96", "#d28530", "#368264"]
        for row, role in enumerate(("electron", "pion")):
            for color, name in zip(colors, grids[role]):
                reference = data[f"{role}_{name}_reference"]
                counts = data[f"{role}_{name}_counts"]
                occupied = counts > 0
                axes[row, 0].hist(np.log10(np.maximum(reference[occupied], .5)), bins=np.linspace(-.4, 5, 55), histtype="step", label=name, color=color)
                stats = summary["schemes"][role][name]
                thresholds = config["occupancy_thresholds"]
                coverage = [100*stats["support"][str(t)]["fraction_of_inside_reference_entries"] for t in thresholds]
                axes[row, 1].plot(thresholds, coverage, "o-", label=name, color=color)
            axes[row, 0].set_title(role.capitalize())
            axes[row, 0].set_xlabel("log10 reference-band entries / occupied cell")
            axes[row, 0].set_ylabel("Occupied cells")
            axes[row, 0].text(.02, .96, "Zero counts shown at -0.30", transform=axes[row, 0].transAxes,
                                va="top", fontsize=8)
            axes[row, 1].set_xlabel("Minimum reference-band entries per cell")
            axes[row, 1].set_ylabel("Inside-grid reference entries supported [%]")
            axes[row, 1].set_ylim(0, 101)
            for ax in axes[row]:
                ax.legend(fontsize=8)
                ax.grid(alpha=.15)
        save(pdf, fig, "statistics")

        for role in ("electron", "pion"):
            occupancy_page(pdf, role)
            occupancy_page(pdf, role, by_sector=True)

        fig, axes = figure("Missing-mass spectra from individual 3D cells", 2, 3,
            "Examples at 10/50/90th percentiles among cells with at least 100 reference-band entries. No peak fits.")
        centers = (data["mass_edges"][:-1] + data["mass_edges"][1:]) / 2
        for row, role in enumerate(("electron", "pion")):
            name = config["particles"][role]["selected_scheme"]
            ref = data[f"{role}_{name}_reference"]
            populated = np.flatnonzero(ref >= 100)
            ordered = populated[np.argsort(ref[populated])]
            if not len(ordered):
                for ax in axes[row]:
                    ax.text(.5, .5, "No cells with >=100 reference entries", ha="center", transform=ax.transAxes)
                continue
            for column, quantile in enumerate((.1, .5, .9)):
                cell = int(ordered[int(quantile * (len(ordered) - 1))])
                ax = axes[row, column]
                bounds = grids[role][name].bounds(cell)
                ax.stairs(data[f"{role}_spectra"][cell], data["mass_edges"], color="#176b96")
                ax.axvline(summary["diagnostics_config"]["conditions"]["neutron_mass_gev"], color="#d28530", linestyle="--", lw=1)
                label = (f"{role} cell {cell}, sector {bounds['sector']}\n"
                         f"p {bounds['p'][0]:g}-{bounds['p'][1]:g}; theta {bounds['theta'][0]:g}-{bounds['theta'][1]:g}\n"
                         f"phi shift {bounds['shifted_phi'][0]:g}-{bounds['shifted_phi'][1]:g}; band N={ref[cell]:,}")
                ax.set_title(label, fontsize=9)
                ax.set_xlabel("Missing mass [GeV/c^2]")
                ax.set_ylabel("Entries / mass bin")
        save(pdf, fig, "cell_spectra")

        fig, axes = figure("Run coverage inside the selected 3D grids", 1, 2,
            "Fraction of each run falling into each p bin; no luminosity normalization and no exclusivity selection.")
        for ax, role in zip(axes.flat, ("electron", "pion")):
            name = config["particles"][role]["selected_scheme"]
            grid = grids[role][name]
            for run, total in sorted(summary["run_entries"].items(), key=lambda x: int(x[0])):
                h = data[f"{role}_run_{run}_counts"].reshape(grid.shape).sum(axis=(0, 2, 3))
                ax.stairs(h / total, grid.edges["p"], label=run)
            ax.set_title(role.capitalize())
            ax.set_xlabel("|p| [GeV/c]")
            ax.set_ylabel("Fraction of input run entries / p bin")
            ax.legend()
        save(pdf, fig, "run_coverage")

        fig, axes = figure("Accounting and use of the saved cells", note="The grid is a calibration domain, not an event-selection decision.")
        ax = axes[0, 0]
        ax.axis("off")
        lines = []
        for role in ("electron", "pion"):
            name = config["particles"][role]["selected_scheme"]
            stats = summary["schemes"][role][name]
            lines.append(f"{role.capitalize()} / {name}:")
            for reason, count in stats["reason_counts"].items():
                lines.append(f"  {reason:24s} {count:>12,d}")
            lines.append("")
        lines += ["Cell IDs are C-order flat indices with axis order (sector, p, theta, shifted_phi).",
                  "Intervals are half-open, except the final right edge is included.",
                  "CSV files record every cell's bounds, counts, reference counts, and per-run populations.",
                  "cells.npz stores coarse/balanced/fine populations and selected-grid missing-mass spectra.",
                  "Invalid/out-of-range missing masses are stored per cell; all in-grid rows are accounted for.",
                  "Coordinate extrema include file/run/event examples in summary.json for tracing outliers.",
                  "The reference band includes background: its count is neither purity nor signal efficiency.",
                  "No fits, energy-loss corrections, momentum corrections, or exclusivity cuts are applied."]
        ax.text(0, 1, "\n".join(lines), transform=ax.transAxes, va="top", fontsize=10, family="monospace", linespacing=1.4)
        save(pdf, fig, "accounting")
    (out / "report_metadata.json").write_text(json.dumps({"pages": page, "renderer_sha256": diag.sha256_file(__file__),
        "summary_sha256": diag.sha256_file(out / "summary.json"), "matplotlib": matplotlib.__version__}, indent=2) + "\n")
    print("Report:", out / "binning_report.pdf", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--diagnostics-config", required=True, type=Path)
    parser.add_argument("--binning-config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--render-only", action="store_true")
    args = parser.parse_args()
    config = json.loads(args.binning_config.read_text())
    grids = load_binning(config)
    diagnostics_config = diag.load_config(args.diagnostics_config)
    if args.render_only:
        summary = json.loads((args.output / "summary.json").read_text())
        if config != summary["binning_config"] or diagnostics_config != summary["diagnostics_config"]:
            raise ValueError("Configs changed; rerun the analysis")
        with np.load(args.output / "cells.npz") as saved:
            data = {k: saved[k] for k in saved.files}
    else:
        summary, data, grids = build(args.input, diagnostics_config, config, args.output)
    report(summary, data, grids, args.output)


if __name__ == "__main__":
    main()
