import json
import os
import tempfile
import unittest

from molecular_geometry import io


def write(tmp, name, text):
    path = os.path.join(tmp, name)
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return path


TOPO = {
    "atoms": [
        {"id": "A", "mass": 1.0},
        {"id": "B", "mass": 2.0},
    ],
    "bonds": [["A", "B"]],
}


class TopologyTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)

    def topo(self, data):
        return io.load_topology(
            write(self.tmp.name, "t.json", json.dumps(data))
        )

    def test_ok_and_bond_normalized(self):
        t = self.topo({**TOPO, "bonds": [["B", "A"]]})
        self.assertEqual(t.bonds, (("A", "B"),))
        self.assertEqual(t.atom_ids, ("A", "B"))

    def test_rejections(self):
        bad = [
            {"atoms": [], "bonds": []},
            {"atoms": [{"id": "A", "mass": 1}, {"id": "A", "mass": 1}],
             "bonds": []},
            {"atoms": [{"id": "A", "mass": 0}], "bonds": []},
            {"atoms": [{"id": "A", "mass": float("nan")}], "bonds": []},
            {"atoms": [{"id": "", "mass": 1}], "bonds": []},
            {**TOPO, "bonds": [["A", "A"]]},
            {**TOPO, "bonds": [["A", "B"], ["B", "A"]]},
            {**TOPO, "bonds": [["A", "C"]]},
        ]
        for data in bad:
            with self.assertRaises(io.InputError, msg=str(data)):
                self.topo(data)


class CoordinateTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.topo = io.load_topology(
            write(self.tmp.name, "t.json", json.dumps(TOPO))
        )

    def coords(self, text, unit="angstrom"):
        path = write(self.tmp.name, "c.csv", text)
        return io.load_coordinates(path, self.topo, unit)

    def test_ok_and_nm_scaling(self):
        text = "conformer,atom,x,y,z\nc1,A,0,0,0\nc1,B,1,0,0\n"
        self.assertEqual(self.coords(text)["c1"]["B"], (1.0, 0.0, 0.0))
        self.assertEqual(
            self.coords(text, unit="nm")["c1"]["B"], (10.0, 0.0, 0.0)
        )

    def test_rejections(self):
        bad = [
            "conformer,atom,x,y,z\n",
            "conformer,atom,x,y,z\nc1,A,0,0,0\n",
            "conformer,atom,x,y,z\nc1,A,0,0,0\nc1,A,1,0,0\nc1,B,0,0,0\n",
            "conformer,atom,x,y,z\nc1,A,0,0,0\nc1,B,nan,0,0\n",
            "conformer,atom,x,y,z\nc1,A,0,0,0\nc1,B,inf,0,0\n",
            "conformer,atom,x,y,z\nc1,A,0,0,0\nc1,C,0,0,0\nc1,B,0,0,0\n",
        ]
        for text in bad:
            with self.assertRaises(io.InputError, msg=text):
                self.coords(text)

    def test_bad_unit(self):
        with self.assertRaises(io.InputError):
            self.coords(
                "conformer,atom,x,y,z\nc1,A,0,0,0\nc1,B,1,0,0\n", unit="pm"
            )


if __name__ == "__main__":
    unittest.main()
