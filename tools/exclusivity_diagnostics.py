#!/usr/bin/env python3
"""Stream a two-track ROOT skim into auditable, uncut ep -> e pi+ X diagnostics.

This is an optional Python analysis frontend; it does not change the C++ v2
histogram schema. No correction, exclusivity selection, or event deduplication
is applied. Histogram ranges affect displays only and their losses are counted.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
import sys

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from momcorr.kinematics import signed_sqrt, epip_kinematics as kinematics
from momcorr.provenance import sha256_file


FIELDS = ("run", "event", "electron_px", "electron_py", "electron_pz",
          "pion_px", "pion_py", "pion_pz", "electron_sector", "pion_sector")


def load_config(path):
    config = json.loads(Path(path).read_text())
    required = {"diagnostics_schema_version", "dataset", "branches", "conditions", "axes"}
    if set(config) != required or config["diagnostics_schema_version"] != 1:
        raise ValueError("Expected diagnostics schema 1 and exactly: " + ", ".join(sorted(required)))
    if set(config["branches"]) != set(FIELDS):
        raise ValueError("Branch mappings must contain exactly: " + ", ".join(FIELDS))
    c = config["conditions"]
    for name in ("beam_energy_gev", "electron_mass_gev", "pion_mass_gev", "target_mass_gev", "neutron_mass_gev"):
        if not np.isfinite(c[name]) or c[name] <= 0:
            raise ValueError(f"{name} must be finite and positive")
    if c["beam_energy_gev"] <= c["electron_mass_gev"]:
        raise ValueError("Beam total energy must exceed electron mass")
    for name, axis in config["axes"].items():
        if set(axis) != {"bins", "low", "high", "label"}:
            raise ValueError(f"Invalid axis fields: {name}")
        if type(axis["bins"]) is not int or axis["bins"] <= 0:
            raise ValueError(f"Axis bins must be a positive integer: {name}")
        if not np.isfinite([axis["low"], axis["high"]]).all() or axis["low"] >= axis["high"]:
            raise ValueError(f"Invalid axis range: {name}")
    required_axes = {"electron_p", "pion_p", "electron_theta", "pion_theta", "electron_phi", "pion_phi", "missing_mass", "missing_mass2", "q2", "w", "delta_pion_p"}
    if set(config["axes"]) != required_axes:
        raise ValueError("Unexpected or missing diagnostic axis")
    return config


class Histograms:
    def __init__(self, axes):
        self.axes = axes
        self.edges = {k: np.linspace(v["low"], v["high"], v["bins"] + 1) for k, v in axes.items()}
        self.counts, self.losses, self.spec = {}, {}, {}

    def add(self, key, variables, data):
        variables = tuple(variables)
        self.spec[key] = variables
        finite = np.logical_and.reduce([np.isfinite(x) for x in data])
        inside = finite.copy()
        for name, x in zip(variables, data):
            edges = self.edges[name]
            inside &= (x >= edges[0]) & (x <= edges[-1])
        loss = self.losses.setdefault(key, Counter())
        loss.update({"input": len(finite), "nonfinite": int((~finite).sum()), "outside_display": int((finite & ~inside).sum())})
        if key not in self.counts:
            self.counts[key] = np.zeros(tuple(len(self.edges[n]) - 1 for n in variables), dtype=np.int64)
        if len(variables) == 1:
            h = np.histogram(data[0][finite], self.edges[variables[0]])[0]
        else:
            h = np.histogram2d(data[0][finite], data[1][finite], [self.edges[n] for n in variables])[0]
        self.counts[key] += h.astype(np.int64)


def analyze(input_dir, config, output_dir):
    files = sorted(Path(input_dir).resolve().glob("*.root"))
    if not files:
        raise ValueError("No .root files found")
    out = Path(output_dir).resolve()
    if out == Path(input_dir).resolve() or Path(input_dir).resolve() in out.parents:
        raise ValueError("Output directory must be outside the input data directory")
    out.mkdir(parents=True, exist_ok=True)
    hist = Histograms(config["axes"])
    manifest, schemas, run_counts, sectors = [], {}, Counter(), {"electron": Counter(), "pion": Counter()}
    invalid, numerical, root_counts = Counter(), Counter(), Counter()
    pairs = []
    branchmap = config["branches"]
    tree_name = config["dataset"]["tree"]
    for index, path in enumerate(files):
        before = path.stat()
        digest = sha256_file(path)
        # Explicit local, synchronous I/O also works in restricted runtimes
        # where the default asynchronous source/executor can stall.
        with uproot.open(path, handler=uproot.source.file.MemmapSource,
                         decompression_executor=uproot.source.futures.TrivialExecutor(),
                         interpretation_executor=uproot.source.futures.TrivialExecutor()) as f:
            tree = f[tree_name]
            missing = set(branchmap.values()) - set(tree.keys())
            if missing:
                raise ValueError(f"Missing branches in {path.name}: {sorted(missing)}")
            signature = json.dumps(tree.typenames(), sort_keys=True)
            schemas[signature] = schemas.get(signature, 0) + 1
            file_runs, file_entries = Counter(), 0
            match = re.search(r"rec_clas_(\d+)", path.name)
            filename_run = int(match.group(1)) if match else None
            for raw in tree.iterate(list(branchmap.values()), step_size="40 MB", library="np"):
                a = {k: raw[v] for k, v in branchmap.items()}
                n = len(a["run"])
                file_entries += n
                unique, counts = np.unique(a["run"], return_counts=True)
                rc = {int(k): int(v) for k, v in zip(unique, counts)}
                run_counts.update(rc)
                file_runs.update(rc)
                good_ids = (a["run"] >= 0) & (a["event"] >= 0)
                invalid["negative_run_or_event_id"] += int((~good_ids).sum())
                pairs.append((a["run"][good_ids].astype(np.uint64) << np.uint64(32)) | a["event"][good_ids].astype(np.uint64))
                finite = np.logical_and.reduce([np.isfinite(a[k]) for k in FIELDS[2:8]])
                invalid["nonfinite_momentum_entries"] += int((~finite).sum())
                for name in FIELDS:
                    invalid[f"nonfinite_{name}"] += int((~np.isfinite(a[name])).sum())
                clean = {k: v[finite] for k, v in a.items()}
                d = kinematics(clean, config["conditions"])
                numerical["finite_momenta_entries"] += len(d["electron_p"])
                numerical["zero_electron_momentum"] += int((d["electron_p"] == 0).sum())
                numerical["zero_pion_momentum"] += int((d["pion_p"] == 0).sum())
                numerical["negative_missing_mass2"] += int((d["missing_mass2"] < 0).sum())
                numerical["negative_missing_energy"] += int((d["missing_energy"] < 0).sum())
                numerical["electron_p_above_beam_energy"] += int((d["electron_p"] > config["conditions"]["beam_energy_gev"]).sum())
                root_counts.update({int(k): int(v) for k, v in zip(*np.unique(d["pion_root_status"], return_counts=True))})
                for name in config["axes"]:
                    hist.add(name, (name,), (d[name],))
                for name in ("electron_p", "pion_p", "electron_theta", "pion_theta", "electron_phi", "pion_phi", "q2", "w"):
                    hist.add(f"missing_mass2_vs_{name}", (name, "missing_mass2"), (d[name], d["missing_mass2"]))
                for species in ("electron", "pion"):
                    hist.add(f"{species}_p_theta", (f"{species}_p", f"{species}_theta"), (d[f"{species}_p"], d[f"{species}_theta"]))
                    values, scounts = np.unique(a[f"{species}_sector"], return_counts=True)
                    sectors[species].update({str(float(k)): int(v) for k, v in zip(values, scounts)})
                    for sector in range(1, 7):
                        mask = clean[f"{species}_sector"] == sector
                        hist.add(f"mm2_{species}_sector_{sector}", ("missing_mass2",), (d["missing_mass2"][mask],))
                        hist.add(f"mm2_vs_p_{species}_sector_{sector}", (f"{species}_p", "missing_mass2"), (d[f"{species}_p"][mask], d["missing_mass2"][mask]))
                for run in rc:
                    mask = clean["run"] == run
                    for variable in ("missing_mass2", "electron_p", "pion_p", "delta_pion_p"):
                        hist.add(f"run_{run}_{variable}", (variable,), (d[variable][mask],))
        after = path.stat()
        if before.st_size != after.st_size or before.st_mtime_ns != after.st_mtime_ns:
            raise ValueError(f"Input changed during read: {path}")
        if filename_run is not None and set(file_runs) != {filename_run}:
            raise ValueError(f"Run mismatch between filename and branches: {path.name} vs {dict(file_runs)}")
        manifest.append({"path": str(path), "bytes": before.st_size, "sha256": digest, "entries": file_entries, "runs": dict(file_runs)})
        if (index + 1) % 25 == 0 or index + 1 == len(files):
            print(f"Read {index + 1}/{len(files)} files; {sum(run_counts.values()):,} entries", flush=True)
    identifiers = np.concatenate(pairs)
    _, multiplicities = np.unique(identifiers, return_counts=True)
    identity = {"unique_run_event_pairs": len(multiplicities), "repeated_pairs": int((multiplicities > 1).sum()),
                "entries_beyond_first_per_pair": int(np.sum(multiplicities - 1)),
                "maximum_entries_per_pair": int(multiplicities.max()) if len(multiplicities) else 0,
                "interpretation": "Repeated pairs may be multiple candidates or duplicated input; no entries were removed."}
    del pairs, identifiers
    try:
        repo = Path(__file__).resolve().parents[1]
        commit = subprocess.check_output(["git", "-C", str(repo), "rev-parse", "HEAD"], text=True).strip()
        git_diff = subprocess.check_output(["git", "-C", str(repo), "status", "--porcelain"], text=True)
    except (OSError, subprocess.CalledProcessError):
        commit, git_diff = None, "unavailable"
    summary = {"generated_at_utc": datetime.now(timezone.utc).isoformat(), "files": len(files), "entries": sum(run_counts.values()),
        "bytes": sum(x["bytes"] for x in manifest), "run_entries": dict(run_counts), "schema_variants": len(schemas),
        "schemas": [{"branches": json.loads(k), "files": v} for k, v in schemas.items()],
        "sector_counts": sectors, "invalid": dict(invalid), "numerical_checks": dict(numerical),
        "pion_root_counts": dict(root_counts), "identity": identity, "config": config,
        "selection": "No physics cuts; no corrections; repeated run/event pairs retained. Nonfinite momenta excluded from kinematic plots only.",
        "histogram_accounting": {k: dict(v, displayed=int(hist.counts[k].sum())) for k, v in hist.losses.items()},
        "histogram_variables": hist.spec, "code": {"base_commit": commit, "working_tree_status": git_diff,
            "script_sha256": sha256_file(__file__)},
        "invocation": {"arguments": sys.argv, "working_directory": str(Path.cwd())},
        "versions": {"python": sys.version, "numpy": np.__version__, "uproot": uproot.__version__}}
    (out / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n")
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    np.savez_compressed(out / "histograms.npz", **hist.counts, **{f"edges_{k}": v for k, v in hist.edges.items()})
    return summary, hist


def make_report(summary, hist, output_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.colors import LogNorm
    import textwrap

    out = Path(output_dir)
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
                         "figure.facecolor": "white", "savefig.facecolor": "white", "pdf.fonttype": 42})
    colors = ["#176b96", "#e08b2c", "#33845b", "#a65278", "#6b63a5", "#6f777a"]
    neutron2 = summary["config"]["conditions"]["neutron_mass_gev"]**2
    pages = []

    def figure(title, nrows=1, ncols=1, note="No added physics cuts or corrections. Entries are candidate rows, not necessarily unique events."):
        fig, ax = plt.subplots(nrows, ncols, figsize=(11.7, 8.3), squeeze=False)
        fig.suptitle(title, x=.065, ha="left", y=.965, fontsize=18, weight="bold", color="#183d52")
        fig.text(.065, .925, note, fontsize=9, color="#4a5964", va="top")
        fig.subplots_adjust(left=.08, right=.95, top=.85, bottom=.10, hspace=.38,
                            wspace=.43 if ncols == 3 else .28)
        return fig, ax

    def save(pdf, fig, name):
        number = len(pages) + 1
        fig.text(.95, .035, str(number), ha="right", fontsize=9, color="#627480")
        pdf.savefig(fig)
        fig.savefig(out / f"{name}.png", dpi=140)
        pages.append((name, fig))

    def line(ax, key, label=None, normalize=False, **kwargs):
        variable = hist.spec[key][0]
        y = hist.counts[key].astype(float)
        if normalize and y.sum():
            y /= y.sum()
        ax.stairs(y, hist.edges[variable], label=label, **kwargs)
        ax.set_xlabel(hist.axes[variable]["label"])
        ax.set_ylabel("Fraction of displayed entries" if normalize else "Entries / bin")

    def heat(ax, key, reference=False, colorbar=True, vmax=None):
        x, y = hist.spec[key]
        h = hist.counts[key]
        vmax = max(2, int(h.max())) if vmax is None else vmax
        m = ax.pcolormesh(hist.edges[x], hist.edges[y], np.ma.masked_equal(h.T, 0), norm=LogNorm(1, vmax), cmap="viridis", rasterized=True)
        ax.set_xlabel(hist.axes[x]["label"])
        ax.set_ylabel(hist.axes[y]["label"])
        if colorbar:
            ax.figure.colorbar(m, ax=ax, label="Entries / bin", pad=.02)
        if reference:
            ax.axhline(neutron2, color="#ffb45a", linewidth=1, linestyle="--")
        return m

    with PdfPages(out / "diagnostics.pdf") as pdf:
        fig, axes = figure("RG-K single-pion input audit", note="First diagnostic pass | ep -> e pi+ X | supplied momenta preserved")
        ax = axes[0, 0]
        ax.axis("off")
        c = summary["config"]["conditions"]
        lines = [f"{summary['files']:,} ROOT files   |   {summary['entries']:,} candidate entries   |   {summary['bytes']/1024**2:.1f} MiB",
                 f"{summary['identity']['unique_run_event_pairs']:,} unique (run, event) pairs; {summary['identity']['repeated_pairs']:,} pairs occur more than once.",
                 f"Tree: {summary['config']['dataset']['tree']}   |   {summary['schema_variants']} schema variant(s)",
                 "", f"Nominal beam total energy: {c['beam_energy_gev']:.5f} GeV; target: proton at rest (LH2).",
                 "Dataset identity: RG-K fall 2018, 6535 MeV, Pass 2 (user/Alex identification).",
                 "No energy-loss corrections in input (user confirmation). Other upstream selection details remain unknown.",
                 "No corrections are applied by this diagnostic frontend.", "",
                 "Branches: run, evn, ex, ey, ez, pipx, pipy, pipz, esec, pipsec.",
                 "PID, vertex, timing, detector status, and track-quality fields are absent.",
                 "Sector values are input labels; no FD/CD classification is inferred from them.", "",
                 "Run          Files              Candidate entries"]
        for r, n in sorted(summary["run_entries"].items(), key=lambda v: int(v[0])):
            count = sum(str(r) in [str(k) for k in row["runs"]] for row in json.loads((out / "manifest.json").read_text()))
            lines.append(f"{int(r):5d}       {count:5d}              {n:12,d}")
        ax.text(0, 1, "\n".join(lines), transform=ax.transAxes, va="top", linespacing=1.5, fontsize=11, family="monospace")
        save(pdf, fig, "01_input_audit")

        fig, axes = figure("Supplied track kinematics", 2, 3)
        for ax, key in zip(axes.flat, ("electron_p", "electron_theta", "electron_phi", "pion_p", "pion_theta", "pion_phi")):
            line(ax, key, color=colors[0])
            ax.set_title(key.replace("_", " "))
            ax.grid(alpha=.15)
        save(pdf, fig, "02_track_kinematics")

        fig, axes = figure("Momentum and angle coverage", 1, 2)
        for ax, name in zip(axes.flat, ("electron", "pion")):
            heat(ax, f"{name}_p_theta")
            ax.set_title(name.capitalize())
        save(pdf, fig, "03_coverage")

        fig, axes = figure("Missing-mass spectrum and inclusive kinematics", 2, 2,
            "Nominal dataset beam energy; no per-run RCDB energy retrieval. Dashed lines mark neutron mass or mass squared.")
        for ax, key in zip(axes.flat, ("missing_mass", "missing_mass2", "q2", "w")):
            line(ax, key, color=colors[0])
            ax.set_yscale("log")
            ax.set_ylim(bottom=.8)
            ax.grid(alpha=.15)
            ax.set_title(key.replace("_", " "))
            if key.startswith("missing_mass"):
                ax.axvline(neutron2 if key.endswith("2") else c["neutron_mass_gev"], color=colors[1], linestyle="--")
        save(pdf, fig, "04_missing_mass")

        fig, axes = figure("Missing mass squared versus momentum and angle", 2, 2,
            "Orange line: neutron mass squared. Full input sample; histogram ranges affect displays only.")
        for ax, variable in zip(axes.flat, ("electron_p", "pion_p", "electron_theta", "pion_theta")):
            heat(ax, f"missing_mass2_vs_{variable}", reference=True)
        save(pdf, fig, "05_missing_mass_correlations")

        for species in ("electron", "pion"):
            fig, axes = figure(f"Missing mass squared by {species} sector label", 2, 3,
                "Input sector labels 1-6; orange line: neutron mass squared. No sector-dependent selection.")
            fig.subplots_adjust(right=.88, wspace=.27)
            vmax = max(2, max(int(hist.counts[f"mm2_vs_p_{species}_sector_{s}"].max()) for s in range(1, 7)))
            for sector, ax in enumerate(axes.flat, 1):
                key = f"mm2_vs_p_{species}_sector_{sector}"
                m = heat(ax, key, reference=True, colorbar=False, vmax=vmax)
                if (sector - 1) % 3:
                    ax.set_ylabel("")
                ax.set_title(f"Sector {sector} | {hist.losses[key]['input']:,} entries")
            fig.colorbar(m, cax=fig.add_axes([.91, .10, .018, .75]), label="Entries / bin (shared scale)")
            save(pdf, fig, f"06_{species}_sectors")

        fig, axes = figure("Run comparisons and input sector labels", 2, 2,
            "Run curves normalized within displayed ranges; these are shape comparisons, not luminosity-normalized yields.")
        for color, run in zip(colors, sorted(summary["run_entries"], key=int)):
            for ax, variable in zip(axes.flat, ("missing_mass2", "electron_p", "pion_p")):
                line(ax, f"run_{run}_{variable}", str(run), normalize=True, color=color)
        for ax in list(axes.flat)[:3]:
            ax.legend(fontsize=8)
            ax.grid(alpha=.15)
        ax = axes[1, 1]
        labels = sorted(set(summary["sector_counts"]["electron"]) | set(summary["sector_counts"]["pion"]), key=float)
        positions = np.arange(len(labels))
        for i, species in enumerate(("electron", "pion")):
            values = [summary["sector_counts"][species].get(k, 0) for k in labels]
            ax.bar(positions + (i - .5) * .36, values, .36, label=species, color=colors[i])
        ax.set_xticks(positions, labels)
        ax.set_xlabel("Input sector label")
        ax.set_ylabel("Candidate entries")
        ax.legend()
        save(pdf, fig, "08_run_comparisons")

        fig, axes = figure("Pion momentum inferred from the neutron hypothesis", 1, 2,
            "Diagnostic only: assumes measured electron and pion direction. No inferred momentum is applied to the data.")
        line(axes[0, 0], "delta_pion_p", color=colors[0])
        axes[0, 0].axvline(0, color=colors[1], linestyle="--")
        axes[0, 0].set_yscale("log")
        axes[0, 0].set_ylim(bottom=.8)
        for color, run in zip(colors, sorted(summary["run_entries"], key=int)):
            line(axes[0, 1], f"run_{run}_delta_pion_p", str(run), normalize=True, color=color)
        axes[0, 1].legend()
        save(pdf, fig, "09_pion_constraint")

        fig, axes = figure("Accounting, limitations, and the next selection study", note="No exclusivity cuts have been chosen or applied in this first pass.")
        ax = axes[0, 0]
        ax.axis("off")
        checks = summary["numerical_checks"]
        roots = summary["pion_root_counts"]
        text = [f"Nonfinite momentum entries: {summary['invalid']['nonfinite_momentum_entries']:,}",
                f"Negative missing mass squared: {checks['negative_missing_mass2']:,}; negative missing energy: {checks['negative_missing_energy']:,}.",
                f"Zero electron / pion momentum: {checks['zero_electron_momentum']:,} / {checks['zero_pion_momentum']:,}.",
                f"Electron momentum above nominal beam energy: {checks['electron_p_above_beam_energy']:,}.",
                "", "Pion constraint roots: " + "; ".join(f"{label}: {roots.get(code, roots.get(str(code), 0)):,}" for code, label in [(0, "none"), (1, "unique"), (2, "ambiguous")]),
                "The delta-p plot includes only a unique positive root satisfying the unsquared energy equation.",
                "Ambiguous/no-root candidates remain in all other diagnostics; no nearest-measurement root choice.", "",
                "Displayed / outside-range / nonfinite entries:"]
        for key in ("missing_mass", "missing_mass2", "q2", "w", "delta_pion_p"):
            h = summary["histogram_accounting"][key]
            text.append(f"  {key:16s} {h['displayed']:>12,d} / {h['outside_display']:>10,d} / {h['nonfinite']:>10,d}")
        text += ["", "Next: inspect the neutron band by sector and kinematics, then fit candidate windows.",
                 "Check selection-induced shifts in delta-p, vary windows, and validate on held-out runs.",
                 "Split by run/event when validating future selections; the present input has unique run/event IDs.",
                 "Missing PID/vertex/timing observables cannot be recovered from this skim.",
                 "Nominal run-group conditions are documented; per-run RCDB values were not retrieved."]
        ax.text(0, 1, "\n".join(text), transform=ax.transAxes, va="top", fontsize=10, linespacing=1.5, family="monospace")
        save(pdf, fig, "10_accounting")

        fig, axes = figure("Definitions and provenance", note="Machine-readable configuration, input hashes, histogram counts, and accounting accompany this report.")
        ax = axes[0, 0]
        ax.axis("off")
        statements = [
            "Four-vectors use GeV, with the beam along +z and a proton target at rest. The beam energy is total energy; its momentum includes the electron mass.",
            "Missing four-vector = beam + target - electron - pion. Missing mass squared is retained even when negative. The displayed mass is sign(M2)*sqrt(abs(M2)); negative values are not discarded.",
            "Q2 = -(beam-electron)^2. W is the signed square root of (beam+target-electron)^2. No W, Q2, missing-mass, or kinematic cut is imposed.",
            "All input rows are counted. Only nonfinite momenta are excluded from derived plots. Histogram under/overflow is recorded separately; displayed ranges are not selections.",
            "No energy-loss corrections in the input is user-confirmed. Other prior corrections and upstream selections are not independently reconstructed from the ten branches.",
            "The RG-K analysis note describes data taking and pass-1 processing. Its run-period mapping and target/field conditions are used here, not its selection cuts or pass-1 QA certification for these pass-2 files.",
            "The nominal 6.53536 GeV energy is from the JeffersonLab clas12-config RG-K low-energy configuration. It is a dataset nominal, not a run-by-run beam-energy measurement."
        ]
        y = 1
        for statement in statements:
            wrapped = textwrap.fill(statement, 115)
            ax.text(0, y, wrapped, va="top", transform=ax.transAxes, fontsize=10, linespacing=1.35)
            y -= .055 * (wrapped.count("\n") + 1) + .02
        y -= .01
        for source in c["sources"]:
            wrapped = textwrap.fill(source["url"], 110)
            ax.text(0, y, wrapped, va="top", transform=ax.transAxes, fontsize=8, color="#176b96")
            y -= .04 * (wrapped.count("\n") + 1) + .01
        save(pdf, fig, "11_provenance")
    for _, fig in pages:
        plt.close(fig)
    (out / "report_metadata.json").write_text(json.dumps({
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "renderer_sha256": sha256_file(__file__), "matplotlib": matplotlib.__version__,
        "summary_sha256": sha256_file(out / "summary.json"),
        "histograms_sha256": sha256_file(out / "histograms.npz"),
        "pages": len(pages)
    }, indent=2) + "\n")
    print(f"Report: {out / 'diagnostics.pdf'}", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--render-only", action="store_true", help="Render an existing summary/histogram bundle without rereading ROOT inputs")
    args = parser.parse_args()
    config = load_config(args.config)
    if args.render_only:
        summary = json.loads((args.output / "summary.json").read_text())
        if config != summary["config"]:
            raise ValueError("Configuration differs from saved analysis; rerun without --render-only")
        hist = Histograms(summary["config"]["axes"])
        with np.load(args.output / "histograms.npz") as data:
            hist.counts = {k: data[k] for k in data.files if not k.startswith("edges_")}
        hist.spec = summary["histogram_variables"]
        hist.losses = summary["histogram_accounting"]
    else:
        summary, hist = analyze(args.input, config, args.output)
    make_report(summary, hist, args.output)


if __name__ == "__main__":
    main()
