"""Coordinate, edge, sector, and ROOT accounting checks for 3D cells."""
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
from momcorr.binning import Grid3D, load_binning, shifted_coordinates
import exclusivity_binning as frontend
import exclusivity_diagnostics as diag


def vector(p, theta, phi):
    theta, phi = np.radians(theta), np.radians(phi)
    return p*np.sin(theta)*np.cos(phi), p*np.sin(theta)*np.sin(phi), p*np.cos(theta)


class BinningChecks(unittest.TestCase):
    def test_electron_sector_unwrap_and_inverse_shift(self):
        p = np.array([3., 3., 3.])
        px, py, pz = vector(p, 15., np.array([10., -170., -50.]))
        c = shifted_coordinates(px, py, pz, np.array([1, 4, 6]), -30., 0.)
        np.testing.assert_allclose(c["local_phi"], [10, 10, 10], atol=1e-12)
        np.testing.assert_allclose(c["shifted_phi"], 0, atol=1e-12)
        np.testing.assert_allclose(c["theta"], 15, atol=1e-12)
        self.assertTrue((c["reason"] == 0).all())

    def test_pion_shift_is_not_wrapped_after_transform(self):
        px, py, pz = vector(np.array([1.05, .1]), 30., np.array([-32., 0.]))
        c = shifted_coordinates(px, py, pz, np.array([1, 1]), 32., .05)
        self.assertAlmostEqual(c["shifted_phi"][0], 0, places=12)
        self.assertAlmostEqual(c["shifted_phi"][1], 640, places=10)

    def test_bad_sector_and_shift_pole_are_explicit(self):
        c = shifted_coordinates([1, 1, 0, .04], [0, 0, 0, 0], [0, 0, 0, 0], [1.5, 0, 1, 1], 32., .05)
        np.testing.assert_array_equal(c["reason"], [2, 2, 1, 3])
        self.assertTrue(np.isnan(c["shifted_phi"]).all())

    def test_exact_edges_all_three_axes_and_sectors(self):
        grid = Grid3D((1, 2), {"p": [0, 1, 2], "theta": [0, 10, 20], "shifted_phi": [-10, 0, 10]})
        coordinates = {"p": np.array([0, 1, 2, 1]), "theta": np.array([0, 10, 20, 10]),
            "shifted_phi": np.array([-10, 0, 10, 0]), "sector": np.array([1, 1, 1, 2]), "reason": np.zeros(4, dtype=np.uint8)}
        ids, reasons = grid.assign(coordinates)
        np.testing.assert_array_equal(ids, [0, 7, 7, 15])
        np.testing.assert_array_equal(reasons, 0)
        self.assertEqual(grid.bounds(15), {"sector": 2, "p": [1., 2.], "theta": [10., 20.], "shifted_phi": [0., 10.]})

    def test_outside_dimensions_have_disjoint_accounting(self):
        grid = Grid3D((1,), {"p": [0, 2], "theta": [0, 20], "shifted_phi": [-10, 10]})
        coordinates = {"p": np.array([3, 1, 1, 1]), "theta": np.array([0, 30, 10, 10]),
            "shifted_phi": np.array([0, 0, 20, 0]), "sector": np.array([1, 1, 1, 2]), "reason": np.zeros(4, dtype=np.uint8)}
        ids, reasons = grid.assign(coordinates)
        np.testing.assert_array_equal(ids, -1)
        np.testing.assert_array_equal(reasons, [4, 5, 6, 2])

    def test_schema_rejects_units_and_duplicate_edges(self):
        c = json.loads((ROOT / "CLAS12/JSON_configs/rgk_epip_binning.json").read_text())
        load_binning(c)
        c["units"]["theta"] = "rad"
        with self.assertRaises(ValueError):
            load_binning(c)
        with self.assertRaises(ValueError):
            Grid3D((1,), {"p": [0, 1, 1], "theta": [0, 20], "shifted_phi": [-10, 10]})

    def test_root_cells_and_spectrum_accounting(self):
        c = json.loads((ROOT / "CLAS12/JSON_configs/rgk_epip_binning.json").read_text())
        for role in c["particles"]:
            c["particles"][role]["schemes"] = {"balanced": {"sectors": [1], "edges": {"p": [0, 7], "theta": [0, 90], "shifted_phi": [-100, 100]}}}
        diagnostic = diag.load_config(ROOT / "CLAS12/JSON_configs/rgk_epip_diagnostics.json")
        arrays = {"run": np.array([5875, 5875], dtype=np.int32), "evn": np.array([10, 11], dtype=np.int32)}
        for role, p, theta, phi in (("electron", 3, 15, 10), ("pion", 1.05, 30, -32)):
            v = vector(p, theta, phi)
            for k, value in zip("xyz", v):
                arrays[diagnostic["branches"][f"{role}_p{k}"]] = np.array([value, value], dtype=np.float32)
            arrays[diagnostic["branches"][f"{role}_sector"]] = np.ones(2, dtype=np.float32)
        with tempfile.TemporaryDirectory() as temporary:
            directory = Path(temporary)
            source = directory / "input"
            source.mkdir()
            with uproot.recreate(source / "fixture.root") as f:
                f.mktree("h22", {k: v.dtype for k, v in arrays.items()})
                f["h22"].extend(arrays)
            s, data, _ = frontend.build(source, diagnostic, c, directory / "out")
            self.assertEqual(s["entries"], 2)
            for role in c["particles"]:
                self.assertEqual(s["schemes"][role]["balanced"]["inside"], 2)
                self.assertEqual(int(data[f"{role}_spectra"].sum() + data[f"{role}_mass_outside"].sum()), 2)
                self.assertEqual(int(data[f"{role}_run_5875_counts"].sum()), 2)


if __name__ == "__main__":
    unittest.main()
