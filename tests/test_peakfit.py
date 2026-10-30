"""Injected-peak recovery, count likelihood, model derivatives, and cut semantics."""

import copy
from pathlib import Path
import sys
import unittest

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from momcorr.models import PeakBackground
from momcorr.peakfit import Histogram, fit_one, fit_ensemble, poisson_deviance
from momcorr.recipe import load_recipe
from momcorr.cuts import WindowTable, window_from_fit


class PeakFitChecks(unittest.TestCase):
    def setUp(self):
        self.recipe = load_recipe(ROOT / "CLAS12/JSON_configs/rgk_epip_cuts.json")[
            "ensemble"
        ]
        self.edges = np.linspace(0.65, 1.35, 351)

    def sample(self, seed=12):
        model = PeakBackground(self.edges, 1)
        center, width = 0.955, 0.018
        p = np.array(
            [
                np.log(3000),
                (center - 0.65) / 0.7,
                np.log(width / 0.7),
                np.log(1600),
                np.log(3500),
            ]
        )
        return Histogram(
            self.edges, np.random.default_rng(seed).poisson(model.expected(p))
        )

    def test_model_integrals_and_analytic_derivatives(self):
        m = PeakBackground(self.edges, 2)
        p = np.array(
            [
                np.log(1000),
                0.43,
                np.log(0.025 / 0.7),
                np.log(100),
                np.log(200),
                np.log(300),
            ]
        )
        signal, background, jac = m.components(p)
        self.assertAlmostEqual(signal.sum(), 1000, places=6)
        self.assertAlmostEqual(background.sum(), 600, places=10)
        for i in range(len(p)):
            step = np.zeros(len(p))
            step[i] = 1e-6
            numeric = (m.expected(p + step) - m.expected(p - step)) / (2e-6)
            np.testing.assert_allclose(jac[:, i], numeric, rtol=2e-5, atol=1e-6)

    def test_recovers_injected_center_and_width_with_sloped_background(self):
        result = fit_one(
            self.sample(),
            {"signal": "gaussian", "background": "bernstein", "degree": 1},
            [0.78, 1.16],
            self.recipe["fit"],
        )
        self.assertTrue(result["accepted"], result)
        self.assertLess(abs(result["center"] - 0.955), 0.002)
        self.assertLess(abs(result["width"] - 0.018), 0.002)
        self.assertGreater(result["center_error"], 0)

    def test_range_and_background_ensemble_produces_candidate(self):
        result = fit_ensemble(self.sample(), self.recipe)
        self.assertEqual(result["status"], "ready", result)
        self.assertGreaterEqual(len(result["winning_fits"]), 2)
        window = window_from_fit(result, 3, 3)
        self.assertEqual(window["status"], "ready")
        self.assertLess(window["low"], 0.955)
        self.assertGreater(window["high"], 0.955)

    def test_fit_is_invariant_to_observable_unit_scale(self):
        h = self.sample()
        model = {"signal": "gaussian", "background": "bernstein", "degree": 1}
        first = fit_one(h, model, [0.78, 1.16], self.recipe["fit"])
        settings = copy.deepcopy(self.recipe["fit"])
        for key in ("center_range", "width_range", "width_seeds"):
            settings[key] = [v * 1000 for v in settings[key]]
        scaled = fit_one(
            Histogram(h.edges * 1000, h.counts), model, [780, 1160], settings
        )
        self.assertTrue(scaled["accepted"])
        self.assertLess(abs(scaled["center"] / 1000 - first["center"]), 1e-7)
        self.assertLess(abs(scaled["width"] / 1000 - first["width"]), 1e-7)

    def test_empty_histogram_does_not_get_a_window(self):
        result = fit_ensemble(Histogram(self.edges, np.zeros(350)), self.recipe)
        self.assertEqual(result["status"], "insufficient")
        self.assertIsNone(window_from_fit(result, 3, 3)["low"])

    def test_background_only_is_not_accepted_as_peak(self):
        h = Histogram(
            self.edges, np.random.default_rng(21).poisson(np.linspace(10, 30, 350))
        )
        result = fit_ensemble(h, self.recipe)
        self.assertNotEqual(result["status"], "ready")

    def test_poisson_empty_bins_and_invalid_histograms(self):
        self.assertAlmostEqual(poisson_deviance(np.array([0.0]), np.array([2.0])), 4)
        with self.assertRaises(ValueError):
            Histogram([0, 1, 2], [1, -1])
        with self.assertRaises(ValueError):
            Histogram([0, 1, 2], [0.5, 1])
        with self.assertRaises(ValueError):
            self.sample().window(0.781, 1.16)


class WindowChecks(unittest.TestCase):
    def table(self):
        return WindowTable(
            [
                {
                    "cell_id": 0,
                    "window": {
                        "status": "ready",
                        "low": 0.9,
                        "high": 1.0,
                        "center": 0.95,
                        "width": 0.02,
                    },
                },
                {
                    "cell_id": 1,
                    "window": {
                        "status": "review",
                        "low": 0.9,
                        "high": 1.0,
                        "center": 0.95,
                        "width": 0.02,
                    },
                },
            ],
            2,
        )

    def test_explicit_reasons_no_negative_id_fallback(self):
        passed, reason = self.table().apply(
            np.array([0, 0, 1, -1, 0]), np.array([0.95, 1.1, 0.95, 0.95, np.nan])
        )
        np.testing.assert_array_equal(passed, [True, False, False, False, False])
        np.testing.assert_array_equal(reason, [0, 4, 3, 1, 2])

    def test_endpoints_included_and_width_sensitivity(self):
        passed, _ = self.table().apply(np.array([0, 0]), np.array([0.9, 1.0]))
        self.assertTrue(passed.all())
        passed, _ = self.table().apply(np.array([0]), np.array([1.02]), width_scale=2)
        self.assertTrue(passed[0])

    def test_missing_or_duplicate_table_entries_rejected(self):
        with self.assertRaises(ValueError):
            WindowTable([], 1)
        with self.assertRaises(ValueError):
            WindowTable(
                [
                    {
                        "cell_id": 0,
                        "window": {
                            "status": "ready",
                            "low": 1.0,
                            "high": 0.9,
                            "center": 0.95,
                            "width": 0.02,
                        },
                    }
                ],
                1,
            )

    def test_overlapping_training_and_validation_runs_rejected(self):
        import json, tempfile

        r = load_recipe(ROOT / "CLAS12/JSON_configs/rgk_epip_cuts.json")
        r["validation_runs"] = r["training_runs"]
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "bad.json"
            path.write_text(json.dumps(r))
            with self.assertRaises(ValueError):
                load_recipe(path)


if __name__ == "__main__":
    unittest.main()
