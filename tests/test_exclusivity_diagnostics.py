"""Physics/accounting checks for the optional diagnostic frontend."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
import uproot

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("diagnostics", ROOT / "tools/exclusivity_diagnostics.py")
d = importlib.util.module_from_spec(spec)
spec.loader.exec_module(d)


class DiagnosticChecks(unittest.TestCase):
    def setUp(self):
        self.config = d.load_config(ROOT / "CLAS12/JSON_configs/rgk_epip_diagnostics.json")

    def test_exact_two_body_final_state_closes_on_neutron(self):
        """Generate pi+n in their CM, boost to lab; recover an on-shell neutron."""
        c = self.config["conditions"]
        me, mpi, mn, mp = [c[k] for k in ("electron_mass_gev", "pion_mass_gev", "neutron_mass_gev", "target_mass_gev")]
        eb = c["beam_energy_gev"]
        electron = np.array([2 * np.sin(.35), 0, 2 * np.cos(.35)])
        ke = eb + mp - np.sqrt(4 + me**2)
        kv = np.array([0., 0., np.sqrt(eb**2 - me**2)]) - electron
        w = np.sqrt(ke**2 - kv @ kv)
        pcm = np.sqrt((w**2 - (mpi + mn)**2) * (w**2 - (mpi - mn)**2)) / (2 * w)
        pstar = np.array([pcm, 0., 0.])
        estar = np.hypot(pcm, mpi)
        beta, gamma = kv / ke, ke / w
        pion = pstar + (((gamma - 1) * (beta @ pstar) / (beta @ beta)) + gamma * estar) * beta
        a = {f"{name}_p{k}": np.array([vec[i]]) for name, vec in (("electron", electron), ("pion", pion)) for i, k in enumerate("xyz")}
        result = d.kinematics(a, c)
        self.assertAlmostEqual(result["missing_mass2"][0], mn**2, places=12)
        self.assertAlmostEqual(result["w"][0], w, places=12)
        self.assertEqual(result["pion_root_status"][0], 1)
        self.assertAlmostEqual(result["delta_pion_p"][0], 0, places=10)

    def test_negative_mass_squared_is_retained(self):
        self.assertTrue(np.array_equal(d.signed_sqrt(np.array([-4., 0., 9.])), [-2, 0, 3]))
        a = {f"{name}_p{k}": np.array([v]) for name, vec in (("electron", (0, 0, 5.)), ("pion", (4., 0, 0))) for k, v in zip("xyz", vec)}
        r = d.kinematics(a, self.config["conditions"])
        self.assertLess(r["missing_mass2"][0], 0)
        self.assertLess(r["missing_mass"][0], 0)

    def test_histograms_account_for_all_rows_and_right_endpoint(self):
        h = d.Histograms({"x": {"bins": 2, "low": 0, "high": 2, "label": "x"}})
        h.add("x", ("x",), (np.array([-1., 0., .5, 2., 3., np.nan]),))
        self.assertEqual(h.counts["x"].sum(), 3)
        self.assertEqual(h.losses["x"]["outside_display"], 2)
        self.assertEqual(h.losses["x"]["nonfinite"], 1)

    def test_root_integration_retains_multiple_candidates(self):
        with tempfile.TemporaryDirectory() as temporary:
            p = Path(temporary)
            source = p / "input"
            source.mkdir()
            arrays = {"run": np.array([5875, 5875], dtype=np.int32), "evn": np.array([10, 10], dtype=np.int32)}
            for branch in list(self.config["branches"].values())[2:]:
                arrays[branch] = np.array([1., 1.], dtype=np.float32)
            with uproot.recreate(source / "rec_clas_005875.fixture.root") as f:
                f.mktree("h22", {k: v.dtype for k, v in arrays.items()})
                f["h22"].extend(arrays)
            summary, _ = d.analyze(source, self.config, p / "output")
            self.assertEqual(summary["entries"], 2)
            self.assertEqual(summary["identity"]["unique_run_event_pairs"], 1)
            self.assertEqual(summary["identity"]["entries_beyond_first_per_pair"], 1)
            for h in summary["histogram_accounting"].values():
                self.assertEqual(h["input"], h["displayed"] + h["outside_display"] + h["nonfinite"])


if __name__ == "__main__":
    unittest.main()
