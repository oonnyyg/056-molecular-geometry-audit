import math
import unittest

from molecular_geometry.geometry import (
    angle_between,
    bond_angle,
    distance,
    torsion_angle,
)


class GeometryTest(unittest.TestCase):
    def test_distance(self):
        self.assertAlmostEqual(distance((0, 0, 0), (3, 4, 0)), 5.0)
        self.assertEqual(distance((1, 1, 1), (1, 1, 1)), 0.0)

    def test_angle(self):
        self.assertAlmostEqual(
            bond_angle((1, 0, 0), (0, 0, 0), (0, 1, 0)), 90.0
        )
        self.assertAlmostEqual(
            bond_angle((1, 0, 0), (0, 0, 0), (-2, 0, 0)), 180.0
        )
        self.assertIsNone(bond_angle((0, 0, 0), (0, 0, 0), (1, 0, 0)))
        self.assertIsNone(angle_between((0, 0, 0), (1, 0, 0)))

    def test_torsion_planar(self):
        # cis planar chain: normals coincide -> 0 degrees
        value, reason = torsion_angle(
            (0, 1, 0), (0, 0, 0), (1, 0, 0), (1, 1, 0)
        )
        self.assertIsNone(reason)
        self.assertAlmostEqual(value, 0.0)
        # trans planar chain: normals opposite -> 180 degrees
        value, _ = torsion_angle((0, -1, 0), (0, 0, 0), (1, 0, 0), (1, 1, 0))
        self.assertAlmostEqual(value, 180.0)

    def test_torsion_right_angle(self):
        value, reason = torsion_angle(
            (0, 0, 1), (0, 0, 0), (1, 0, 0), (1, 1, 0)
        )
        self.assertIsNone(reason)
        self.assertAlmostEqual(value, 90.0)

    def test_torsion_undefined(self):
        # collinear middle triple
        value, reason = torsion_angle(
            (0, 0, 0), (1, 0, 0), (2, 0, 0), (3, 0, 0)
        )
        self.assertIsNone(value)
        self.assertIn("collinear", reason)
        # zero-length bond
        value, reason = torsion_angle(
            (0, 0, 0), (1, 0, 0), (1, 0, 0), (2, 1, 0)
        )
        self.assertIsNone(value)
        self.assertIn("zero-length", reason)


if __name__ == "__main__":
    unittest.main()
