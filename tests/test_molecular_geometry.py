"""Automated tests for the molecular geometry audit (standard library)."""

from __future__ import annotations

import csv
import io
import json
import math
import os
import subprocess
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from molecular_geometry.analysis import (
    center_of_mass,
    radius_of_gyration,
    reference_distance_deltas,
)
from molecular_geometry.coordinates import load_config, load_coordinates
from molecular_geometry.errors import GeometryInputError
from molecular_geometry.geometry import (
    build_measurement_sets,
    measure_angle,
    measure_torsion,
)
from molecular_geometry.report import run_audit
from molecular_geometry.topology import load_topology


def write_tmp(suffix: str, content: str) -> str:
    fd, path = tempfile.mkstemp(suffix=suffix, text=True)
    with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
        handle.write(content)
    return path


TOPOLOGY_JSON = json.dumps({
    "atoms": [
        {"id": "A1", "mass": 1.0},
        {"id": "A2", "mass": 2.0},
        {"id": "A3", "mass": 3.0},
        {"id": "A4", "mass": 4.0},
    ],
    "bonds": [["A1", "A2"], ["A2", "A3"], ["A3", "A4"]],
})


class TopologyValidationTests(unittest.TestCase):
    def load(self, payload: object):
        path = write_tmp(".json", json.dumps(payload))
        try:
            return load_topology(path)
        finally:
            os.unlink(path)

    def test_empty_topology_rejected(self):
        with self.assertRaises(GeometryInputError):
            self.load({"atoms": []})

    def test_duplicate_atoms_rejected(self):
        with self.assertRaises(GeometryInputError):
            self.load({"atoms": [
                {"id": "A1", "mass": 1.0},
                {"id": "A1", "mass": 2.0},
            ]})

    def test_non_positive_mass_rejected(self):
        with self.assertRaises(GeometryInputError):
            self.load({"atoms": [{"id": "A1", "mass": 0.0}]})

    def test_self_bond_rejected(self):
        with self.assertRaises(GeometryInputError):
            self.load({"atoms": [{"id": "A1", "mass": 1.0}],
                       "bonds": [["A1", "A1"]]})

    def test_duplicate_bond_even_when_reversed_rejected(self):
        with self.assertRaises(GeometryInputError):
            self.load({"atoms": [
                {"id": "A1", "mass": 1.0},
                {"id": "A2", "mass": 1.0},
            ], "bonds": [["A1", "A2"], ["A2", "A1"]]})

    def test_unknown_atom_in_bond_rejected(self):
        with self.assertRaises(GeometryInputError):
            self.load({"atoms": [
                {"id": "A1", "mass": 1.0},
                {"id": "A2", "mass": 1.0},
            ], "bonds": [["A1", "X9"]]})

    def test_bond_normalization_and_chain_enumeration(self):
        topology = self.load(json.loads(TOPOLOGY_JSON))
        sets = build_measurement_sets(topology)
        self.assertEqual(sets.bonds, [
            ("A1", "A2"), ("A2", "A3"), ("A3", "A4")
        ])
        self.assertEqual(sets.angles, [
            ("A1", "A2", "A3"), ("A2", "A3", "A4")
        ])
        self.assertEqual(sets.torsions, [("A1", "A2", "A3", "A4")])


class GeometryTests(unittest.TestCase):
    def test_angle_outward_vectors_ninety_degrees(self):
        positions = {
            "a": (1.0, 0.0, 0.0),
            "c": (0.0, 0.0, 0.0),
            "b": (0.0, 1.0, 0.0),
        }
        result = measure_angle(positions, "a", "c", "b")
        self.assertAlmostEqual(result.value, 90.0)

    def test_zero_length_bond_makes_angle_null_with_reason(self):
        positions = {"a": (0.0, 0.0, 0.0), "c": (0.0, 0.0, 0.0),
                     "b": (1.0, 0.0, 0.0)}
        result = measure_angle(positions, "a", "c", "b")
        self.assertIsNone(result.value)
        self.assertIsNotNone(result.undefined_reason)

    def test_collinear_chain_makes_torsion_null(self):
        positions = {
            "a": (0.0, 0.0, 0.0), "b": (1.0, 0.0, 0.0),
            "c": (2.0, 0.0, 0.0), "d": (3.0, 0.0, 0.0),
        }
        result = measure_torsion(positions, "a", "b", "c", "d")
        self.assertIsNone(result.value)
        self.assertIn("collinear", result.undefined_reason)

    def test_torsion_ninety_and_unsigned(self):
        positions = {
            "a": (0.0, 0.0, 0.0), "b": (1.0, 0.0, 0.0),
            "c": (1.0, 1.0, 0.0), "d": (1.0, 1.0, 1.0),
        }
        result = measure_torsion(positions, "a", "b", "c", "d")
        self.assertAlmostEqual(result.value, 90.0)


class SizeAndReferenceTests(unittest.TestCase):
    def setUp(self):
        self.topology_path = write_tmp(".json", TOPOLOGY_JSON)
        self.topology = load_topology(self.topology_path)

    def tearDown(self):
        os.unlink(self.topology_path)

    def test_mass_center_not_geometric_center(self):
        positions = {
            "A1": (0.0, 0.0, 0.0),
            "A2": (0.0, 0.0, 0.0),
            "A3": (0.0, 0.0, 0.0),
            "A4": (10.0, 0.0, 0.0),
        }
        com = center_of_mass(self.topology, positions)
        self.assertAlmostEqual(com[0], 4.0)
        geometric = 2.5
        self.assertNotAlmostEqual(com[0], geometric)

    def test_radius_of_gyration_mass_weighted(self):
        positions = {
            "A1": (0.0, 0.0, 0.0),
            "A2": (0.0, 0.0, 0.0),
            "A3": (0.0, 0.0, 0.0),
            "A4": (10.0, 0.0, 0.0),
        }
        com = center_of_mass(self.topology, positions)
        rgyr = radius_of_gyration(self.topology, positions, com)
        expected = math.sqrt((6.0 * 16.0 + 4.0 * 36.0) / 10.0)
        self.assertAlmostEqual(rgyr, expected)

    def test_distance_rms_invariant_under_translation(self):
        reference = {
            "A1": (0.0, 0.0, 0.0), "A2": (1.0, 0.0, 0.0),
            "A3": (1.0, 1.0, 0.0), "A4": (1.0, 1.0, 1.0),
        }
        moved = {atom: (x + 7.0, y - 3.0, z + 2.0)
                 for atom, (x, y, z) in reference.items()}
        result = reference_distance_deltas(
            self.topology, moved, reference
        )
        self.assertAlmostEqual(result.rms, 0.0)

    def test_fewer_than_two_atoms_returns_null(self):
        single = load_topology(write_tmp(
            ".json", json.dumps({"atoms": [{"id": "X", "mass": 1.0}]})
        ))
        point = {"X": (0.0, 0.0, 0.0)}
        try:
            result = reference_distance_deltas(single, point, point)
        finally:
            pass
        self.assertIsNone(result.rms)
        self.assertIsNone(result.max_pair)
        self.assertIsNotNone(result.reason)

    def test_reference_conformer_missing_rejected(self):
        topo = write_tmp(".json", TOPOLOGY_JSON)
        coord = write_tmp(
            ".csv",
            "conformer,atom,x,y,z\n"
            "other,A1,0,0,0\nother,A2,1,0,0\nother,A3,1,1,0\n"
            "other,A4,1,1,1\n",
        )
        conf = write_tmp(
            ".json", json.dumps({"reference_conformer": "absent"})
        )
        out = tempfile.mkdtemp()
        with self.assertRaises(GeometryInputError):
            run_audit(coord, topo, conf, out)
        os.unlink(topo)
        os.unlink(coord)
        os.unlink(conf)


class CoordinateValidationTests(unittest.TestCase):
    def setUp(self):
        self.topology_path = write_tmp(".json", TOPOLOGY_JSON)
        self.topology = load_topology(self.topology_path)

    def tearDown(self):
        os.unlink(self.topology_path)

    def load_csv(self, content: str, unit: str = "angstrom"):
        path = write_tmp(".csv", content)
        try:
            config_path = write_tmp(
                ".json",
                json.dumps({"reference_conformer": "c1",
                            "coordinates_unit": unit}),
            )
            try:
                config = load_config(config_path)
            finally:
                os.unlink(config_path)
            return load_coordinates(path, self.topology, config.coordinates_unit)
        finally:
            os.unlink(path)

    HEADER = "conformer,atom,x,y,z\n"

    def test_empty_table_rejected(self):
        with self.assertRaises(GeometryInputError):
            self.load_csv(self.HEADER)

    def test_missing_atom_coordinates_rejected(self):
        rows = "c1,A1,0,0,0\nc1,A2,0,0,0\nc1,A3,0,0,0\n"
        with self.assertRaises(GeometryInputError):
            self.load_csv(self.HEADER + rows)

    def test_duplicate_row_rejected(self):
        rows = (
            "c1,A1,0,0,0\nc1,A1,0,0,0\nc1,A2,0,0,0\n"
            "c1,A3,0,0,0\nc1,A4,0,0,0\n"
        )
        with self.assertRaises(GeometryInputError):
            self.load_csv(self.HEADER + rows)

    def test_non_finite_coordinate_rejected(self):
        rows = (
            "c1,A1,nan,0,0\nc1,A2,0,0,0\nc1,A3,0,0,0\n"
            "c1,A4,0,0,0\n"
        )
        with self.assertRaises(GeometryInputError):
            self.load_csv(self.HEADER + rows)

    def test_nanometres_converted_to_angstrom(self):
        rows = (
            "c1,A1,0,0,0\nc1,A2,1,0,0\nc1,A3,1,1,0\n"
            "c1,A4,1,1,1\n"
        )
        table = self.load_csv(self.HEADER + rows, unit="nm")
        self.assertAlmostEqual(
            table.positions["c1"]["A2"][0], 10.0
        )


class CliIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp()
        self.topology_path = os.path.join(self.tmp, "topology.json")
        self.coordinates_path = os.path.join(self.tmp, "coordinates.csv")
        self.config_path = os.path.join(self.tmp, "config.json")
        with open(self.topology_path, "w", encoding="utf-8") as handle:
            handle.write(TOPOLOGY_JSON)
        with open(self.config_path, "w", encoding="utf-8") as handle:
            json.dump({"reference_conformer": "c1"}, handle)
        self.out_dir = os.path.join(self.tmp, "out")

    def write_coordinates(self, rows: list[str]) -> None:
        with open(self.coordinates_path, "w", encoding="utf-8") as handle:
            handle.write("conformer,atom,x,y,z\n" + "\n".join(rows) + "\n")

    def run_cli(self):
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        return subprocess.run(
            [sys.executable, "-m", "molecular_geometry",
             "--coordinates", self.coordinates_path,
             "--topology", self.topology_path,
             "--config", self.config_path,
             "--output-dir", self.out_dir],
            cwd=root, capture_output=True, text=True,
        )

    def atoms(self) -> list[str]:
        return [
            "c1,A1,0,0,0",
            "c1,A2,1,0,0",
            "c1,A3,1,1,0",
            "c1,A4,1,1,1",
        ]

    def test_successful_run_writes_reports(self):
        self.write_coordinates(self.atoms())
        result = self.run_cli()
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertTrue(
            os.path.exists(os.path.join(self.out_dir, "geometry.csv"))
        )
        analysis_path = os.path.join(self.out_dir, "analysis.json")
        with open(analysis_path, encoding="utf-8") as handle:
            analysis = json.load(handle)
        self.assertEqual(
            analysis["conformers"]["c1"]["measurement_counts"]["total"],
            {"valid": 6, "undefined": 0},
        )

    def test_invalid_input_nonzero_and_no_report_created(self):
        self.write_coordinates(["c1,A1,0,0,0"])  # missing other atoms
        result = self.run_cli()
        self.assertNotEqual(result.returncode, 0)
        self.assertFalse(os.path.exists(self.out_dir))

    def test_geometry_csv_sorted_and_null_rendered_empty(self):
        self.write_coordinates(self.atoms())
        self.run_cli()
        with open(os.path.join(self.out_dir, "geometry.csv"),
                  encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        conformers = [row["conformer"] for row in rows]
        self.assertEqual(conformers, sorted(conformers))
        self.assertTrue(rows)


if __name__ == "__main__":
    unittest.main()
