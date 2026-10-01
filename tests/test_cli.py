import csv
import json
import os
import tempfile
import unittest

from molecular_geometry.cli import run

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXAMPLES = os.path.join(ROOT, "examples")


class CliTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.out = os.path.join(self.tmp.name, "out")

    def args(self, **over):
        base = {
            "topology": os.path.join(EXAMPLES, "topology.json"),
            "coordinates": os.path.join(EXAMPLES, "coordinates.csv"),
            "unit": "angstrom",
            "config": os.path.join(EXAMPLES, "config.json"),
            "output_dir": self.out,
        }
        base.update(over)
        return [f"--{k.replace('_', '-')}={v}" for k, v in base.items()]

    def test_examples_end_to_end(self):
        self.assertEqual(run(self.args()), 0)
        with open(os.path.join(self.out, "analysis.json"),
                  encoding="utf-8") as fh:
            report = json.load(fh)
        confs = report["conformers"]
        self.assertEqual(list(confs), ["linear", "overlap", "rotated"])
        self.assertAlmostEqual(
            confs["rotated"]["reference_deviation"]["rmsd"], 0.0
        )
        self.assertEqual(
            confs["linear"]["summary"]["counts"]["torsion"]["undefined"], 1
        )
        overlap = confs["overlap"]
        self.assertEqual(overlap["summary"]["counts"]["angle"]["valid"], 0)
        self.assertIsNone(overlap["summary"]["angle"])
        self.assertEqual(overlap["summary"]["bond_length"]["min"], 0.0)
        with open(os.path.join(self.out, "geometry.csv"),
                  encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        undef = [r for r in rows if r["value"] == ""]
        self.assertTrue(undef)
        self.assertTrue(all(r["reason"] for r in undef))
        for row in rows:
            if row["value"]:
                self.assertNotIn("nan", row["value"].lower())

    def test_input_error_no_output(self):
        bad = os.path.join(self.tmp.name, "bad.json")
        with open(bad, "w", encoding="utf-8") as fh:
            json.dump({"atoms": [], "bonds": []}, fh)
        self.assertEqual(run(self.args(topology=bad)), 1)
        self.assertFalse(os.path.exists(self.out))

    def test_unknown_reference_no_output(self):
        cfg = os.path.join(self.tmp.name, "cfg.json")
        with open(cfg, "w", encoding="utf-8") as fh:
            json.dump({"reference_conformer": "nope"}, fh)
        self.assertEqual(run(self.args(config=cfg)), 1)
        self.assertFalse(os.path.exists(self.out))


if __name__ == "__main__":
    unittest.main()
