"""Run-split cache isolation and one-to-one ROOT friend-mask integration."""

import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
import uproot

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
sys.path.insert(0, str(ROOT / "tools"))
import exclusivity_diagnostics as diag
from momcorr.binning import load_binning
from momcorr.calibration import build_cache, observables, derive
from momcorr.evaluation import apply_selection
from momcorr.provenance import sha256_file
from momcorr.recipe import load_recipe
from momcorr.models import PeakBackground


class WorkflowChecks(unittest.TestCase):
    def test_changing_validation_spectra_does_not_change_fitted_windows(self):
        diagnostics = diag.load_config(
            ROOT / "CLAS12/JSON_configs/rgk_epip_diagnostics.json"
        )
        binning = json.loads(
            (ROOT / "CLAS12/JSON_configs/rgk_epip_binning.json").read_text()
        )
        for role in binning["particles"]:
            binning["particles"][role]["schemes"] = {
                "balanced": {
                    "sectors": [1],
                    "edges": {
                        "p": [0, 7],
                        "theta": [0, 90],
                        "shifted_phi": [-150, 150],
                    },
                }
            }
        recipe = load_recipe(ROOT / "CLAS12/JSON_configs/rgk_epip_cuts.json")
        recipe["roles"] = ["electron"]
        edges = np.array(binning["mass_edges_gev"])
        model = PeakBackground(edges, 1)
        parameters = np.array(
            [
                np.log(3000),
                (0.955 - 0.65) / 0.7,
                np.log(0.018 / 0.7),
                np.log(1600),
                np.log(3500),
            ]
        )
        training = np.random.default_rng(12).poisson(model.expected(parameters))[
            None, :
        ]
        meta = {
            "roles": recipe["roles"],
            "training_runs": recipe["training_runs"],
            "validation_runs": recipe["validation_runs"],
            "binning": binning,
            "diagnostics": diagnostics,
        }
        with tempfile.TemporaryDirectory() as temporary:
            selections = []
            for i in range(2):
                p = Path(temporary) / str(i)
                p.mkdir()
                (p / "cache.json").write_text(json.dumps(meta))
                np.savez_compressed(
                    p / "histograms.npz",
                    mass_edges=edges,
                    electron_training_spectra=training,
                    electron_validation_spectra=np.full_like(training, i * 100000),
                )
                selections.append(derive(p, recipe, p / "model", workers=1))
            self.assertEqual(selections[0]["cells"], selections[1]["cells"])
            self.assertEqual(
                selections[0]["cells"]["electron"][0]["window"]["status"], "ready"
            )

    def test_run_isolation_mask_alignment_and_input_preservation(self):
        diagnostics = diag.load_config(
            ROOT / "CLAS12/JSON_configs/rgk_epip_diagnostics.json"
        )
        binning = json.loads(
            (ROOT / "CLAS12/JSON_configs/rgk_epip_binning.json").read_text()
        )
        recipe = load_recipe(ROOT / "CLAS12/JSON_configs/rgk_epip_cuts.json")
        recipe["training_runs"], recipe["validation_runs"] = [1], [2]
        for role in binning["particles"]:
            binning["particles"][role]["schemes"] = {
                "balanced": {
                    "sectors": [1],
                    "edges": {
                        "p": [0, 7],
                        "theta": [0, 90],
                        "shifted_phi": [-150, 150],
                    },
                }
            }
        with tempfile.TemporaryDirectory() as temporary:
            p = Path(temporary)
            source = p / "input"
            source.mkdir()
            for run in (1, 2):
                arrays = {
                    "run": np.full(2, run, dtype=np.int32),
                    "evn": np.array([10, 11], dtype=np.int32),
                    "ex": np.array([0.6, 0.6], dtype=np.float32),
                    "ey": np.zeros(2, dtype=np.float32),
                    "ez": np.array([2.8, 2.8], dtype=np.float32),
                    "pipx": np.array([0.3, 0.3], dtype=np.float32),
                    "pipy": np.zeros(2, dtype=np.float32),
                    "pipz": np.array([1.0, 1.4], dtype=np.float32),
                    "esec": np.ones(2, dtype=np.float32),
                    "pipsec": np.ones(2, dtype=np.float32),
                }
                with uproot.recreate(source / f"run_{run}.root") as f:
                    f.mktree("h22", {k: v.dtype for k, v in arrays.items()})
                    f["h22"].extend(arrays)
            hashes = {f.name: sha256_file(f) for f in source.glob("*.root")}
            meta = build_cache(source, diagnostics, binning, recipe, p / "cache")
            self.assertEqual(meta["entries"], 4)
            for role in recipe["roles"]:
                self.assertEqual(meta["accounting"][role]["training"]["input"], 2)
                self.assertEqual(meta["accounting"][role]["validation"]["input"], 2)
            logical = {k: arrays[v] for k, v in diagnostics["branches"].items()}
            mass = observables(logical, diagnostics["conditions"])["missing_mass"]
            self.assertGreater(abs(mass[0] - mass[1]), 0.01)
            cells = {
                role: [
                    {
                        "cell_id": 0,
                        "window": {
                            "status": "ready",
                            "low": float(mass[0] - 0.001),
                            "high": float(mass[0] + 0.001),
                            "center": float(mass[0]),
                            "width": 0.0003,
                        },
                    }
                ]
                for role in recipe["roles"]
            }
            selection = {
                "selection_schema_version": 1,
                "status": "test_candidate",
                "observable": "missing_mass",
                "recipe": recipe,
                "diagnostics": diagnostics,
                "binning": binning,
                "cells": cells,
            }
            path = p / "selection.json"
            path.write_text(json.dumps(selection))
            result = apply_selection(source, selection, p / "evaluation", path)
            for role in recipe["roles"]:
                self.assertEqual(result["decisions"][role]["validation"]["pass"], 1)
                self.assertEqual(
                    result["decisions"][role]["validation"]["outside_window"], 1
                )
            masks = list((p / "evaluation/masks").glob("*.root"))
            self.assertEqual(len(masks), 2)
            for mask in masks:
                with uproot.open(
                    mask,
                    handler=uproot.source.file.MemmapSource,
                    decompression_executor=uproot.source.futures.TrivialExecutor(),
                    interpretation_executor=uproot.source.futures.TrivialExecutor(),
                ) as f:
                    a = f["cuts"].arrays(library="np")
                    np.testing.assert_array_equal(a["entry"], [0, 1])
                    np.testing.assert_array_equal(a["evn"], [10, 11])
                    np.testing.assert_array_equal(a["electron_pass"], [True, False])
            self.assertEqual(
                hashes, {f.name: sha256_file(f) for f in source.glob("*.root")}
            )


if __name__ == "__main__":
    unittest.main()
