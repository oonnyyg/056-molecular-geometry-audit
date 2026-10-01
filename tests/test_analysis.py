import math
import unittest

from molecular_geometry import analysis
from molecular_geometry.io import Topology


def make_topology():
    return Topology(
        atom_ids=("A1", "A2", "A3", "A4"),
        masses={"A1": 1.0, "A2": 2.0, "A3": 3.0, "A4": 4.0},
        bonds=(("A1", "A2"), ("A2", "A3"), ("A3", "A4")),
    )


class EnumerateTest(unittest.TestCase):
    def test_chain_measurements(self):
        bonds, angles, torsions = analysis.enumerate_measurements(
            make_topology()
        )
        self.assertEqual(len(bonds), 3)
        self.assertEqual(angles, [("A1", "A2", "A3"), ("A2", "A3", "A4")])
        self.assertEqual(torsions, [("A1", "A2", "A3", "A4")])

    def test_star_has_no_torsion(self):
        topo = Topology(
            atom_ids=("C", "X", "Y", "Z"),
            masses={"C": 1.0, "X": 1.0, "Y": 1.0, "Z": 1.0},
            bonds=(("C", "X"), ("C", "Y"), ("C", "Z")),
        )
        bonds, angles, torsions = analysis.enumerate_measurements(topo)
        self.assertEqual(len(angles), 3)
        self.assertEqual(torsions, [])


class SizeTest(unittest.TestCase):
    def test_com_and_gyration_mass_weighted(self):
        coords = {"A": (0.0, 0.0, 0.0), "B": (2.0, 0.0, 0.0)}
        masses = {"A": 1.0, "B": 3.0}
        self.assertEqual(
            analysis.center_of_mass(coords, masses), (1.5, 0.0, 0.0)
        )
        self.assertAlmostEqual(
            analysis.radius_of_gyration(coords, masses), math.sqrt(0.75)
        )


class DeviationTest(unittest.TestCase):
    def test_rigid_transform_invariant(self):
        ref = {"A": (0, 0, 0), "B": (1, 0, 0), "C": (0, 1, 0)}
        moved = {"A": (5, 5, 5), "B": (5, 6, 5), "C": (5, 5, 4)}
        got = analysis.pairwise_distance_deviation(moved, ref, ("A", "B", "C"))
        self.assertAlmostEqual(got["rmsd"], 0.0)
        self.assertAlmostEqual(got["max_abs_diff"], 0.0)

    def test_max_pair(self):
        ref = {"A": (0, 0, 0), "B": (1, 0, 0), "C": (0, 2, 0)}
        other = {"A": (0, 0, 0), "B": (1, 0, 0), "C": (0, 4, 0)}
        got = analysis.pairwise_distance_deviation(other, ref, ("A", "B", "C"))
        self.assertEqual(got["max_pair"], ["A", "C"])
        self.assertAlmostEqual(got["max_abs_diff"], 2.0)

    def test_single_atom(self):
        got = analysis.pairwise_distance_deviation(
            {"A": (0, 0, 0)}, {"A": (1, 1, 1)}, ("A",)
        )
        self.assertIsNone(got["rmsd"])
        self.assertIsNotNone(got["reason"])


class SummarizeTest(unittest.TestCase):
    def test_counts_and_null_stats(self):
        results = [
            {"kind": "bond", "atoms": ("A", "B"), "value": 1.0,
             "reason": None},
            {"kind": "angle", "atoms": ("A", "B", "C"), "value": None,
             "reason": "zero-length bond"},
        ]
        summary = analysis.summarize(results)
        self.assertEqual(summary["counts"]["angle"]["undefined"], 1)
        self.assertIsNone(summary["angle"])
        self.assertEqual(summary["bond_length"]["mean"], 1.0)


if __name__ == "__main__":
    unittest.main()
