#!/usr/bin/env python3
"""CLI for reusable cached-spectrum fitting and candidate cut application."""
import argparse
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "python"))
from momcorr.binning import load_binning
from momcorr.calibration import build_cache, derive
from momcorr.evaluation import apply_selection
from momcorr.recipe import load_recipe
from momcorr.channels import load_input_config
import exclusivity_diagnostics as diagnostics


def cache_args(parser):
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument(
        "--dataset-config",
        "--diagnostics-config",
        dest="diagnostics_config",
        required=True,
        type=Path,
        help="Explicit channel/input config (schema 2) or legacy RG-K diagnostics config",
    )
    parser.add_argument("--binning-config", required=True, type=Path)
    parser.add_argument("--recipe", required=True, type=Path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)
    cache = sub.add_parser(
        "cache", help="Read ROOT once into training/validation cell spectra"
    )
    cache_args(cache)
    cache.add_argument("--output", required=True, type=Path)
    fit = sub.add_parser(
        "derive",
        help="Fit only cached training spectra; save candidate windows and review results",
    )
    fit.add_argument("--cache", required=True, type=Path)
    fit.add_argument("--recipe", required=True, type=Path)
    fit.add_argument("--output", required=True, type=Path)
    fit.add_argument("--workers", type=int, default=4)
    apply = sub.add_parser(
        "apply", help="Apply saved windows independently by particle role"
    )
    apply.add_argument("--input", required=True, type=Path)
    apply.add_argument("--selection", required=True, type=Path)
    apply.add_argument("--output", required=True, type=Path)
    apply.add_argument("--no-masks", action="store_true")
    report = sub.add_parser(
        "report", help="Render fit/review and held-out diagnostics from saved artifacts"
    )
    report.add_argument("--cache", required=True, type=Path)
    report.add_argument("--model", required=True, type=Path)
    report.add_argument("--evaluation", required=True, type=Path)
    report.add_argument("--output", required=True, type=Path)
    run = sub.add_parser("run", help="Execute cache, derive, apply, and report stages")
    cache_args(run)
    run.add_argument("--output", required=True, type=Path)
    run.add_argument("--workers", type=int, default=4)
    run.add_argument("--no-masks", action="store_true")
    args = parser.parse_args()
    if args.command in {"run", "cache"}:
        config = load_input_config(args.diagnostics_config)
        if config.get("diagnostics_schema_version") == 1:
            config = diagnostics.load_config(args.diagnostics_config)
        binning = json.loads(args.binning_config.read_text())
        load_binning(binning)
        recipe = load_recipe(args.recipe)
        directory = args.output / "cache" if args.command == "run" else args.output
        build_cache(args.input, config, binning, recipe, directory)
    if args.command in {"run", "derive"}:
        recipe = load_recipe(args.recipe)
        cache = args.output / "cache" if args.command == "run" else args.cache
        model = args.output / "model" if args.command == "run" else args.output
        derive(cache, recipe, model, args.workers)
    if args.command in {"run", "apply"}:
        path = (
            args.output / "model/cuts.json" if args.command == "run" else args.selection
        )
        selection = json.loads(path.read_text())
        directory = args.output / "evaluation" if args.command == "run" else args.output
        apply_selection(
            args.input, selection, directory, path, write_masks=not args.no_masks
        )
    if args.command in {"run", "report"}:
        from momcorr.cut_report import create_report

        cache = args.output / "cache" if args.command == "run" else args.cache
        model = args.output / "model" if args.command == "run" else args.model
        evaluation = (
            args.output / "evaluation" if args.command == "run" else args.evaluation
        )
        create_report(cache, model, evaluation, args.output)


if __name__ == "__main__":
    main()
