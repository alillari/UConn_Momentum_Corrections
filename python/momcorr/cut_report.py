"""Render saved fit and held-out selection results; contains no fitting policy."""

from collections import Counter
import gzip
import json
from pathlib import Path
import textwrap

import numpy as np

from .models import PeakBackground
from .peakfit import Histogram
from .provenance import sha256_file


def create_report(cache_dir, model_dir, evaluation_dir, output):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from scipy.special import xlogy

    cache_dir, model_dir, evaluation_dir, output = map(
        Path, (cache_dir, model_dir, evaluation_dir, output)
    )
    output.mkdir(parents=True, exist_ok=True)
    selection = json.loads((model_dir / "cuts.json").read_text())
    evaluation = json.loads((evaluation_dir / "evaluation.json").read_text())
    metadata = json.loads((cache_dir / "cache.json").read_text())
    if selection["source_cache_sha256"] != sha256_file(
        cache_dir / "cache.json"
    ) or selection["source_histograms_sha256"] != sha256_file(
        cache_dir / "histograms.npz"
    ):
        raise ValueError("Fit selection does not match the supplied cache")
    if evaluation["selection_sha256"] != sha256_file(model_dir / "cuts.json"):
        raise ValueError("Evaluation does not match the supplied selection")
    with gzip.open(model_dir / "fits.json.gz", "rt") as stream:
        results = json.load(stream)
    results = {role: {r["cell_id"]: r for r in rows} for role, rows in results.items()}
    with np.load(cache_dir / "histograms.npz") as archive:
        cached = {k: archive[k] for k in archive.files}
    with np.load(evaluation_dir / "evaluation_histograms.npz") as archive:
        hist = {k: archive[k] for k in archive.files}
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 10, "pdf.fonttype": 42})
    roles = list(selection["cells"])
    modern = selection["recipe"]["recipe_schema_version"] == 2
    spectrum_spec = selection["recipe"].get("spectrum", {})
    observable_key = "observable" if modern else "mass"
    observable_label = spectrum_spec.get("label", "Missing mass [GeV/c^2]")
    parameter_unit = spectrum_spec.get("unit", "GeV")
    observable_edges = cached["observable_edges" if modern else "mass_edges"]
    reference_value = (
        spectrum_spec["reference_value"]
        if modern
        else selection["diagnostics"]["conditions"]["neutron_mass_gev"]
    )
    role_groups = [roles[i : i + 2] for i in range(0, len(roles), 2)]
    colors = {"ready": "#277e62", "review": "#d48a2c", "insufficient": "#8697a2"}
    page = 0

    def figure(
        title,
        rows=1,
        cols=1,
        note="Candidate selection; algorithmic fit checks do not establish physics validity.",
    ):
        fig, axes = plt.subplots(rows, cols, figsize=(11.7, 8.3), squeeze=False)
        fig.suptitle(
            title,
            x=0.065,
            y=0.97,
            ha="left",
            fontsize=17,
            weight="bold",
            color="#183d52",
        )
        fig.text(0.065, 0.925, note, va="top", fontsize=9, color="#526571")
        fig.subplots_adjust(
            left=0.08, right=0.94, top=0.84, bottom=0.10, wspace=0.36, hspace=0.45
        )
        return fig, axes

    def save(pdf, fig, name):
        nonlocal page
        page += 1
        fig.text(0.95, 0.035, str(page), ha="right", color="#526571")
        pdf.savefig(fig)
        fig.savefig(output / f"{page:02d}_{name}.png", dpi=140)
        plt.close(fig)

    def text_page(pdf, title, paragraphs, name):
        lines = []
        for p in paragraphs:
            lines += textwrap.wrap(p, 112) + [""]
        for start in range(0, len(lines), 30):
            fig, axes = figure(title + (" (continued)" if start else ""))
            ax = axes[0, 0]
            ax.axis("off")
            ax.text(
                0,
                1,
                "\n".join(lines[start : start + 30]),
                transform=ax.transAxes,
                va="top",
                fontsize=10,
                linespacing=1.4,
            )
            save(pdf, fig, name + f"_{start//30+1}")

    def spectrum(ax, key, edges, label, color, normalize=False):
        counts = hist[key].astype(float)
        if normalize and counts.sum():
            counts /= counts.sum()
        ax.stairs(counts, edges, label=label, color=color)

    with PdfPages(output / "exclusivity_report.pdf") as pdf:
        paragraphs = [
            f"Input: {metadata['files']:,} files and {metadata['entries']:,} entries. Data snapshot only; upstream expected total not independently established.",
            "Training runs: "
            + ", ".join(map(str, selection["recipe"]["training_runs"]))
            + ". Held-out validation runs: "
            + ", ".join(map(str, selection["recipe"]["validation_runs"]))
            + ".",
            "Each particle's sector-separated 3D cell has a Gaussian-peak fit with positive polynomial background alternatives. Fits compare identical domains by AICc; domain changes and plausible-model differences are separate stability checks.",
            "Ready = configured algorithmic gates passed; review = ambiguous or failed checks; insufficient = too few reference-band entries. Unsupported cells receive no inferred or neighboring cut.",
        ]
        if modern:
            channel = metadata["diagnostics"]["channel"]
            paragraphs.insert(
                1,
                f"Channel: {channel['name']}. Measured tracks: {', '.join(channel['measured'])}. Selected observable: {observable_label}; omitted measured tracks: {selection['observable']['omit']}.",
            )
            paragraphs.insert(
                2,
                "Dataset identity: "
                + metadata["diagnostics"]["dataset"].get("identity", "Not supplied"),
            )
        for role, rows in selection["cells"].items():
            statuses = Counter(c["window"]["status"] for c in rows)
            v = evaluation["decisions"][role]["validation"]
            supported = v["pass"] + v["outside_window"]
            paragraphs.append(
                f"{role.capitalize()}: {statuses['ready']:,} ready, {statuses['review']:,} review, {statuses['insufficient']:,} insufficient cells. Validation support covers {100*supported/max(1,v['input']):.2f}% of input rows; {v['pass']:,}/{v['input']:,} rows pass ({100*v['pass']/max(1,v['input']):.2f}%). Conditional retention in supported cells: {100*v['pass']/max(1,supported):.2f}%."
            )
        paragraphs += [
            "Roles are evaluated independently. No implicit combination of particle selections is used. Pass fractions describe retained data, not signal efficiency or purity.",
            "Input momenta are preserved and no energy-loss correction is applied. Nominal beam energy, coordinate conventions, bin edges, and all fit/cut settings are saved with the candidate.",
        ]
        text_page(pdf, "First candidate exclusivity windows", paragraphs, "overview")

        for group in role_groups:
            fig, axes = figure(
                "Fit status and supported validation coverage", 1, len(group)
            )
            for ax, role in zip(axes.flat, group):
                counts = Counter(
                    c["window"]["status"] for c in selection["cells"][role]
                )
                labels = ["ready", "review", "insufficient"]
                ax.bar(
                    labels,
                    [counts[k] for k in labels],
                    color=[colors[k] for k in labels],
                )
                ax.set_title(role.capitalize())
                ax.set_ylabel("3D cells, including empty cells")
            save(pdf, fig, "fit_status")

        for role in roles:
            for status in ("ready", "review"):
                rows = [
                    c
                    for c in selection["cells"][role]
                    if c["window"]["status"] == status
                ]
                rows.sort(key=lambda c: c["training_reference_entries"])
                if status == "review":
                    rows = [
                        c
                        for c in rows
                        if any(
                            "parameters" in f
                            for f in results[role][c["cell_id"]]["fits"]
                        )
                    ]
                fig, axes = figure(
                    f"{role.capitalize()}: {status} cell fit examples",
                    2,
                    3,
                    "Training data only. Upper: counts, peak, background and sum. Lower: signed Poisson-deviance residuals.",
                )
                if not rows:
                    for ax in axes.flat:
                        ax.axis("off")
                    axes[0, 0].text(
                        0,
                        0.8,
                        "No cells in this category",
                        transform=axes[0, 0].transAxes,
                    )
                    save(pdf, fig, f"{role}_{status}_examples")
                    continue
                for col, q in enumerate((0.1, 0.5, 0.9)):
                    row = rows[int(q * (len(rows) - 1))]
                    i = row["cell_id"]
                    record = results[role][i]
                    if record["chosen_fit"] is not None:
                        fit = record["fits"][record["chosen_fit"]]
                    else:
                        available = [f for f in record["fits"] if "parameters" in f]
                        preferred = selection["recipe"]["ensemble"]["preferred_window"]
                        window_id = (
                            preferred
                            if any(f["window_id"] == preferred for f in available)
                            else min(f["window_id"] for f in available)
                        )
                        fit = min(
                            (f for f in available if f["window_id"] == window_id),
                            key=lambda f: f["aicc"],
                        )
                    h = Histogram(
                        observable_edges, cached[f"{role}_training_spectra"][i]
                    ).window(*fit["window"])
                    model = PeakBackground(h.edges, **fit["model"])
                    signal, background, _ = model.components(
                        np.array(fit["parameters"])
                    )
                    expected = np.maximum(signal + background, 1e-12)
                    centers = (h.edges[:-1] + h.edges[1:]) / 2
                    ax = axes[0, col]
                    ax.stairs(h.counts, h.edges, color="#253e50", label="data")
                    ax.plot(centers, expected, color="#176b96", label="sum")
                    ax.plot(centers, signal, color="#277e62", label="peak")
                    ax.plot(centers, background, color="#d48a2c", label="background")
                    window = row["window"]
                    if window["low"] is not None:
                        for edge in (window["low"], window["high"]):
                            ax.axvline(edge, color="#aa5276", linestyle="--", lw=1)
                    issues = ", ".join(row["issues"]) or "stable"
                    ax.set_title(
                        f"Cell {i}, S{row['bounds']['sector']}, N={row['training_reference_entries']}, degree {fit['model']['degree']}\n"
                        + textwrap.fill(issues, 38),
                        fontsize=8,
                    )
                    ax.set_xlabel(observable_label)
                    ax.set_ylabel("Counts / observable bin")
                    ax.legend(fontsize=7)
                    residual = np.sign(h.counts - expected) * np.sqrt(
                        np.maximum(
                            0,
                            2
                            * (
                                expected
                                - h.counts
                                + xlogy(h.counts, h.counts / expected)
                            ),
                        )
                    )
                    axes[1, col].plot(centers, residual, ".", ms=3, color="#176b96")
                    axes[1, col].axhline(0, color="#586976", lw=0.7)
                    axes[1, col].set_xlabel(observable_label)
                    axes[1, col].set_ylabel("Signed deviance residual")
                save(pdf, fig, f"{role}_{status}_examples")

        for group in role_groups:
            fig, axes = figure(
                "Accepted peak parameters across the 3D cells",
                len(group),
                2,
                "One actual fit defines each candidate; alternative models/ranges are not averaged to shrink errors.",
            )
            for row, role in enumerate(group):
                ready = [
                    c
                    for c in selection["cells"][role]
                    if c["window"]["status"] == "ready"
                ]
                p = np.array([np.mean(c["bounds"]["p"]) for c in ready])
                sector = np.array([c["bounds"]["sector"] for c in ready])
                for col, key in enumerate(("center", "width")):
                    values = np.array([c["window"][key] for c in ready])
                    axes[row, col].scatter(
                        p,
                        values,
                        c=sector,
                        cmap="tab10",
                        s=7,
                        alpha=0.6,
                        vmin=1,
                        vmax=6,
                    )
                    axes[row, col].set_xlabel("Cell-center |p| [GeV/c]")
                    axes[row, col].set_ylabel(f"{role} peak {key} [{parameter_unit}]")
                if reference_value is not None:
                    axes[row, 0].axhline(
                        reference_value, color="#d48a2c", linestyle="--", lw=1
                    )
            save(pdf, fig, "parameters")

        variables = [(observable_key, hist[observable_key + "_edges"])]
        if not modern:
            variables.append(("delta", hist["delta_edges"]))
        for variable, edges in variables:
            note = (
                "Held-out runs only. This is the selected observable; this plot does not independently validate the cut."
                if variable == observable_key
                else "Held-out runs only. Delta-p is correlated with missing mass: this is a sculpting diagnostic, not truth closure."
            )
            for group in role_groups:
                fig, axes = figure(
                    "Held-out "
                    + (
                        "selected-observable spectra"
                        if variable == observable_key
                        else "pion momentum-residual spectra"
                    ),
                    len(group),
                    2,
                    note,
                )
                for row, role in enumerate(group):
                    for col, normalize in enumerate((False, True)):
                        for stage, color in (
                            ("before", "#176b96"),
                            ("after", "#d48a2c"),
                        ):
                            spectrum(
                                axes[row, col],
                                f"{role}_validation_{stage}_{variable}",
                                edges,
                                stage,
                                color,
                                normalize,
                            )
                        axes[row, col].set_title(role + "-based selection")
                        axes[row, col].set_xlabel(
                            observable_label
                            if variable == observable_key
                            else "Inferred minus measured pion p [GeV/c]"
                        )
                        axes[row, col].set_ylabel(
                            "Fraction of displayed entries"
                            if normalize
                            else "Counts / bin"
                        )
                        if not normalize:
                            axes[row, col].set_yscale("log")
                            axes[row, col].set_ylim(bottom=0.8)
                        axes[row, col].legend()
                save(pdf, fig, "validation_" + variable)

        for group in role_groups:
            fig, axes = figure(
                "Cut-width sensitivity on held-out data",
                len(group),
                2,
                "Same training-derived cells and fits; vary only the window width. These correlated alternatives are not repeated experiments.",
            )
            for row, role in enumerate(group):
                labels, retained = [], []
                total = evaluation["decisions"][role]["validation"]["input"]
                palette = plt.get_cmap("tab10")
                for index, (scale, counts) in enumerate(
                    evaluation["width_sensitivity"][role].items()
                ):
                    color = palette(index % 10)
                    label = f"{float(scale):.3g}x nominal width"
                    spectrum(
                        axes[row, 0],
                        f"{role}_sensitivity_{scale}_{observable_key if modern else 'delta'}",
                        (
                            hist[observable_key + "_edges"]
                            if modern
                            else hist["delta_edges"]
                        ),
                        label,
                        color,
                        True,
                    )
                    labels.append(f"{float(scale):.3g}x")
                    retained.append(100 * counts["validation_pass"] / max(1, total))
                axes[row, 0].set_xlabel(
                    observable_label
                    if modern
                    else "Inferred minus measured pion p [GeV/c]"
                )
                axes[row, 0].set_ylabel(role + ": fraction displayed")
                axes[row, 0].legend(fontsize=8)
                axes[row, 1].bar(
                    labels,
                    retained,
                    color=[palette(i % 10) for i in range(len(labels))],
                )
                axes[row, 1].set_xlabel("Window width scale")
                axes[row, 1].set_ylabel("All validation rows retained [%]")
            save(pdf, fig, "width_sensitivity")

        for group in role_groups:
            fig, axes = figure(
                "Support and retention by run",
                1,
                len(group),
                "Support loss and rejection by the observable window are separate quantities. T/V denotes training/validation.",
            )
            for ax, role in zip(axes.flat, group):
                runs = sorted(
                    int(key.split("_")[-1])
                    for key in evaluation["by_run"]
                    if key.startswith(role + "_")
                )
                support = []
                passed = []
                for run in runs:
                    c = evaluation["by_run"][f"{role}_{run}"]
                    support.append(100 * c["supported"] / c["input"])
                    passed.append(100 * c["pass"] / c["input"])
                x = np.arange(len(runs))
                ax.bar(x - 0.18, support, 0.36, label="supported", color="#176b96")
                ax.bar(x + 0.18, passed, 0.36, label="pass", color="#d48a2c")
                ax.set_xticks(
                    x,
                    [
                        f"{r}\n{'T' if r in selection['recipe']['training_runs'] else 'V'}"
                        for r in runs
                    ],
                )
                ax.set_title(role.capitalize())
                ax.set_ylabel("Input rows [%]")
                ax.legend()
            save(pdf, fig, "runs")

        paragraphs = [
            "Likelihood: independent unweighted Poisson bin counts. Gaussian peak and positive Bernstein-polynomial backgrounds are integrated over each bin. Empty bins participate. Analytic derivatives are used for optimization and local Fisher covariance.",
            "Model comparison: AICc is compared only on identical fit ranges. Plausible alternatives within the configured AICc difference, and winners from different ranges, must satisfy center/width stability gates. One actual fit supplies the cut; statistical errors are never divided by the number of alternatives.",
            "The likelihood-gain contrast is an algorithmic local diagnostic, not a discovery significance. Covariance errors describe the selected model locally and do not establish coverage under an unknown true background.",
            "Cut policy: center minus/plus the configured multiples of fitted width. Boundaries extending beyond the fitted domain require review. Sparse, failed, boundary-limited, poor-covariance, weak, or model/range-sensitive cases have no applied window.",
            "Validation: run groups are disjoint. No fit is made to validation spectra. Before/after spectra directly reflect selection on that same variable. Independent signal-efficiency, purity, or unbiased correction claims require further physics validation.",
            "Artifacts: cuts.json contains every cell and its status; fits.json.gz preserves all alternatives and failures. ROOT friend masks retain the original row order and run/event IDs, with separate pass, reason, cell, and grid-reason fields per particle role. Input file hashes and selection hashes are recorded.",
            "Input provenance note: "
            + metadata["diagnostics"]["dataset"].get(
                "upstream_selection", "Not documented in the input metadata."
            )
            + " Preparation-file hashes, when present, are saved in cache.json. Snapshot size: "
            + str(metadata["files"])
            + " ROOT files; upstream expected total not independently established.",
        ]
        paragraphs += list(selection["recipe"].get("annotations", {}).values())
        text_page(
            pdf, "Method, provenance, and remaining validation", paragraphs, "method"
        )
    (output / "report_metadata.json").write_text(
        json.dumps(
            {
                "pages": page,
                "renderer_sha256": sha256_file(__file__),
                "selection_sha256": sha256_file(model_dir / "cuts.json"),
                "evaluation_sha256": sha256_file(evaluation_dir / "evaluation.json"),
                "matplotlib": matplotlib.__version__,
            },
            indent=2,
        )
        + "\n"
    )
    print("Report:", output / "exclusivity_report.pdf", flush=True)
