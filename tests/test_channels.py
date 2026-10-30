"""Topology, charged identity, four-vector closure, and multi-track workflow."""

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
from channel_fixtures import configs, exact_event, logical_rows, one_cell
from momcorr.channels import (
    validate_input_config,
    validate_channel,
    observable_values,
    residual_four_vector,
    validate_workflow,
    load_input_config,
)
from momcorr.binning import load_binning
from momcorr.recipe import load_recipe
from momcorr.calibration import build_cache, derive
from momcorr.evaluation import apply_selection
from momcorr.kinematics import epip_kinematics
from momcorr.skim import coordinates


class ChannelChecks(unittest.TestCase):
    def test_example_schemas_and_mismatched_topology(self):
        for name in ("epipn", "ep_elastic", "eppippim"):
            config, binning, recipe = configs(name)
            self.assertEqual(
                load_input_config(ROOT / "examples/channels" / name / "input.json"),
                config,
            )
            self.assertEqual(
                load_recipe(ROOT / "examples/channels" / name / "cuts.json"), recipe
            )
            load_binning(binning)
            validate_workflow(config, binning, recipe)
        config, binning, recipe = configs("ep_elastic")
        recipe["roles"].append("pi_minus")
        with self.assertRaisesRegex(ValueError, "measured track"):
            validate_workflow(config, binning, recipe)

    def test_negative_pion_requires_all_four_measured_tracks(self):
        channel = configs("eppippim")[0]["channel"]
        for particle in channel["measured"]:
            bad = copy.deepcopy(channel)
            bad["measured"].remove(particle)
            with self.assertRaises(ValueError):
                validate_channel(bad)
        bad = {
            "name": "wrong_charge",
            "final_state": ["electron", "pi_minus", "neutron"],
            "measured": ["electron", "pi_minus"],
        }
        with self.assertRaises(ValueError):
            validate_channel(bad)

    def test_charge_names_and_branch_mapping_cannot_alias_tracks(self):
        config, _, recipe = configs("eppippim")
        config["branches"]["pi_minus_px"] = config["branches"]["pi_plus_px"]
        with self.assertRaisesRegex(ValueError, "distinct"):
            validate_input_config(config)
        config, binning, recipe = configs("eppippim")
        binning["particles"]["pion"] = binning["particles"].pop("pi_minus")
        with self.assertRaises(ValueError):
            load_binning(binning)
        config, binning, recipe = configs("ep_elastic")
        binning["particles"]["proton"]["detector"] = "CD"
        with self.assertRaisesRegex(ValueError, "FD"):
            load_binning(binning)

    def test_elastic_and_double_pion_exact_closure_and_leave_one_out_masses(self):
        for name in ("ep_elastic", "eppippim"):
            config, _, recipe = configs(name)
            vectors = exact_event(config)
            data = logical_rows(config, vectors)
            np.testing.assert_allclose(
                residual_four_vector(data, config), 0, atol=2e-14
            )
            np.testing.assert_allclose(
                observable_values(data, config, recipe), 0, atol=2e-14
            )
            for particle in config["channel"]["measured"]:
                recipe["observable"] = {"kind": "missing_mass2", "omit": [particle]}
                expected = config["conditions"]["particle_masses_gev"][particle] ** 2
                np.testing.assert_allclose(
                    observable_values(data, config, recipe), expected, atol=2e-14
                )
                np.testing.assert_allclose(
                    residual_four_vector(data, config, [particle])[0],
                    vectors[particle],
                    atol=2e-14,
                )

    def test_new_pi_plus_adapter_matches_legacy_and_neutron_closure(self):
        config, _, recipe = configs("epipn")
        data = logical_rows(config, exact_event(config), n=3)
        legacy = {
            key.replace("pi_plus_", "pion_"): value for key, value in data.items()
        }
        old_conditions = json.loads(
            (ROOT / "CLAS12/JSON_configs/rgk_epip_diagnostics.json").read_text()
        )["conditions"]
        values = observable_values(data, config, recipe)
        np.testing.assert_allclose(
            values, config["conditions"]["particle_masses_gev"]["neutron"], atol=1e-14
        )
        np.testing.assert_allclose(
            values, epip_kinematics(legacy, old_conditions)["missing_mass"], atol=1e-14
        )

    def test_each_closure_observable_has_expected_sign_and_nonfinite_policy(self):
        config, _, recipe = configs("eppippim")
        data = logical_rows(config, exact_event(config), n=2)
        data["pi_minus_px"][0] += 0.1
        data["pi_minus_px"][1] = np.nan
        residual = residual_four_vector(data, config)
        self.assertAlmostEqual(residual[0, 1], -0.1, places=13)
        expected = {
            "missing_px": residual[:, 1],
            "missing_py": residual[:, 2],
            "missing_pz": residual[:, 3],
            "missing_energy": residual[:, 0],
            "missing_p": np.linalg.norm(residual[:, 1:], axis=1),
            "missing_pt": np.linalg.norm(residual[:, 1:3], axis=1),
        }
        for kind, values in expected.items():
            recipe["observable"] = {"kind": kind, "omit": []}
            np.testing.assert_allclose(
                observable_values(data, config, recipe), values, equal_nan=True
            )
        recipe["observable"] = {"kind": "missing_mass", "omit": ["pi_minus"]}
        self.assertTrue(np.isnan(observable_values(data, config, recipe)[1]))
        for axis in "xyz":
            data[f"pi_minus_p{axis}"][0] = 0.0
        self.assertTrue(np.isnan(observable_values(data, config, recipe)[0]))
        del data["pi_minus_pz"]
        with self.assertRaisesRegex(ValueError, "Missing measured-track"):
            residual_four_vector(data, config)

    def test_distinct_particle_shifts_are_used_without_charge_substitution(self):
        config, binning, _ = configs("eppippim")
        data = logical_rows(config, exact_event(config))
        for particle, coefficient in (
            ("pi_plus", 17.0),
            ("pi_minus", -11.0),
            ("proton", 3.0),
        ):
            spec = copy.deepcopy(binning["particles"][particle])
            spec["shift"] = {"coefficient_deg_gev": coefficient, "pole_gev": 0.0}
            coords = coordinates(data, particle, spec)
            np.testing.assert_allclose(
                coords["shifted_phi"] - coords["local_phi"], coefficient / coords["p"]
            )

    def test_spectrum_units_and_omission_must_match(self):
        config, binning, recipe = configs("ep_elastic")
        recipe["observable"]["omit"] = ["pi_minus"]
        with self.assertRaisesRegex(ValueError, "absent"):
            validate_workflow(config, binning, recipe)
        recipe["observable"]["omit"] = []
        recipe["spectrum"]["unit"] = "GeV/c"
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "recipe.json"
            path.write_text(json.dumps(recipe))
            with self.assertRaisesRegex(ValueError, "unit"):
                load_recipe(path)


class ChannelWorkflowChecks(unittest.TestCase):
    def test_elastic_and_four_track_pipeline_keeps_all_role_masks_aligned(self):
        for name in ("ep_elastic", "eppippim"):
            with self.subTest(channel=name), tempfile.TemporaryDirectory() as temporary:
                directory = Path(temporary)
                config, binning, recipe = configs(name)
                binning = one_cell(binning)
                vectors = exact_event(config)
                source = directory / "input"
                source.mkdir()
                expected = {}
                rng = np.random.default_rng(243)
                for run, n in ((5875, 6000), (5880, 3000)):
                    residuals = np.r_[
                        rng.normal(0.005, 0.018, int(n * 0.7)),
                        rng.uniform(-0.4, 0.4, n - int(n * 0.7)),
                    ]
                    rng.shuffle(residuals)
                    logical = logical_rows(config, vectors, n, run, residuals)
                    expected[run] = logical
                    raw = {
                        config["branches"][key]: value for key, value in logical.items()
                    }
                    with uproot.recreate(source / f"run_{run}.root") as file:
                        file.mktree(
                            "h22", {key: value.dtype for key, value in raw.items()}
                        )
                        file["h22"].extend(raw)
                metadata = build_cache(
                    source, config, binning, recipe, directory / "cache"
                )
                self.assertEqual(metadata["entries"], 9000)
                for particle in config["channel"]["measured"]:
                    self.assertEqual(
                        metadata["accounting"][particle]["training"]["spectrum"], 6000
                    )
                selection = derive(
                    directory / "cache", recipe, directory / "model", workers=1
                )
                for particle in config["channel"]["measured"]:
                    window = selection["cells"][particle][0]["window"]
                    self.assertEqual(window["status"], "ready")
                    self.assertAlmostEqual(window["center"], 0.005, delta=0.002)
                    self.assertAlmostEqual(window["width"], 0.018, delta=0.002)
                summary = apply_selection(
                    source,
                    selection,
                    directory / "evaluation",
                    directory / "model/cuts.json",
                )
                self.assertEqual(summary["entries"], 9000)
                manifest = json.loads(
                    (directory / "evaluation/mask_manifest.json").read_text()
                )
                options = dict(
                    handler=uproot.source.file.MemmapSource,
                    decompression_executor=uproot.TrivialExecutor(),
                    interpretation_executor=uproot.TrivialExecutor(),
                )
                for record in manifest:
                    with uproot.open(record["mask"], **options) as file:
                        rows = file["cuts"].arrays(library="np")
                    run = int(rows["run"][0])
                    logical = expected[run]
                    np.testing.assert_array_equal(
                        rows["entry"], np.arange(len(logical["run"]))
                    )
                    np.testing.assert_array_equal(rows["evn"], logical["event"])
                    values = observable_values(logical, config, recipe)
                    for particle in config["channel"]["measured"]:
                        window = selection["cells"][particle][0]["window"]
                        np.testing.assert_array_equal(
                            rows[particle + "_pass"],
                            (values >= window["low"]) & (values <= window["high"]),
                        )
                        self.assertEqual(
                            summary["decisions"][particle]["validation"]["input"], 3000
                        )
                changed = copy.deepcopy(recipe)
                changed["observable"] = {"kind": "missing_px", "omit": []}
                with self.assertRaisesRegex(ValueError, "rebuild cache"):
                    derive(
                        directory / "cache", changed, directory / "changed", workers=1
                    )
                changed = copy.deepcopy(recipe)
                changed["spectrum"]["edges"][1] += 0.0001
                with self.assertRaisesRegex(ValueError, "rebuild cache"):
                    derive(
                        directory / "cache", changed, directory / "changed", workers=1
                    )

    def test_missing_negative_pion_branch_fails_before_creating_a_cache(self):
        config, binning, recipe = configs("eppippim")
        data = logical_rows(config, exact_event(config))
        raw = {
            config["branches"][key]: value
            for key, value in data.items()
            if key != "pi_minus_px"
        }
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "input"
            source.mkdir()
            with uproot.recreate(source / "missing.root") as file:
                file.mktree("h22", {key: value.dtype for key, value in raw.items()})
                file["h22"].extend(raw)
            with self.assertRaisesRegex(ValueError, "Missing branches.*pimx"):
                build_cache(
                    source, config, one_cell(binning), recipe, directory / "cache"
                )
            self.assertFalse((directory / "cache/cache.json").exists())


if __name__ == "__main__":
    unittest.main()
