"""Core geometric measurements: vectors, bonds, angles and torsions."""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .topology import Topology

Vec3 = tuple[float, float, float]


def _sub(a: Vec3, b: Vec3) -> Vec3:
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def _dot(a: Vec3, b: Vec3) -> float:
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _cross(a: Vec3, b: Vec3) -> Vec3:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(a: Vec3) -> float:
    return math.sqrt(_dot(a, a))


def _angle_between(a: Vec3, b: Vec3) -> float:
    """Angle in degrees between two vectors (0..180)."""
    length_a = _norm(a)
    length_b = _norm(b)
    cosine = _dot(a, b) / (length_a * length_b)
    cosine = min(1.0, max(-1.0, cosine))
    return math.degrees(math.acos(cosine))


@dataclass(frozen=True)
class Measurement:
    """One local geometric measurement for one conformer."""

    kind: str  # "bond_length" | "bond_angle" | "torsion"
    atoms: tuple[str, ...]
    value: float | None
    undefined_reason: str | None = None


@dataclass
class MeasurementSets:
    """Measurement definitions shared by every conformer."""

    bonds: list[tuple[str, str]] = field(default_factory=list)
    angles: list[tuple[str, str, str]] = field(default_factory=list)
    torsions: list[tuple[str, str, str, str]] = field(default_factory=list)


def build_measurement_sets(topology: Topology) -> MeasurementSets:
    """Enumerate bonds, adjacent bond pairs (angles) and three-bond chains.

    * Each undirected bond is stored with its lexicographically smaller atom
      first.
    * An angle is formed by every pair of bonds sharing a centre atom.
    * A torsion chain a-b-c-d requires four distinct atoms and the bonds
      a-b, b-c, c-d; the reversed chain d-c-b-a is the same torsion, kept
      once using the lexicographically smaller atom-id sequence.
    """
    bonds = sorted(topology.bonds)

    angles: set[tuple[str, str, str]] = set()
    for center, neighbors in topology.adjacency.items():
        n = len(neighbors)
        for i in range(n):
            for j in range(i + 1, n):
                first, second = sorted((neighbors[i], neighbors[j]))
                angles.add((first, center, second))

    torsions: set[tuple[str, str, str, str]] = set()
    for middle_a, middle_b in bonds:
        # Bond middle_a--middle_b is the central bond of the chain.
        for outer_a in topology.adjacency[middle_a]:
            if outer_a == middle_b:
                continue
            for outer_b in topology.adjacency[middle_b]:
                if outer_b == middle_a or outer_b == outer_a:
                    continue
                chain = (outer_a, middle_a, middle_b, outer_b)
                reversed_chain = (outer_b, middle_b, middle_a, outer_a)
                torsions.add(min(chain, reversed_chain))

    return MeasurementSets(
        bonds=bonds,
        angles=sorted(angles),
        torsions=sorted(torsions),
    )


def measure_bond(positions: dict[str, Vec3], a: str, b: str) -> Measurement:
    """Distance between two bonded atoms; zero length is a valid value."""
    distance = _norm(_sub(positions[a], positions[b]))
    return Measurement("bond_length", (a, b), distance)


def measure_angle(
    positions: dict[str, Vec3], a: str, center: str, b: str
) -> Measurement:
    """Angle between outward vectors center->a and center->b, in degrees."""
    vec_a = _sub(positions[a], positions[center])
    vec_b = _sub(positions[b], positions[center])
    if _norm(vec_a) == 0.0:
        return Measurement(
            "bond_angle", (a, center, b), None,
            f"zero-length bond {center}-{a}",
        )
    if _norm(vec_b) == 0.0:
        return Measurement(
            "bond_angle", (a, center, b), None,
            f"zero-length bond {center}-{b}",
        )
    return Measurement("bond_angle", (a, center, b), _angle_between(vec_a, vec_b))


def measure_torsion(
    positions: dict[str, Vec3], a: str, b: str, c: str, d: str
) -> Measurement:
    """Unsigned torsion for chain a-b-c-d, in the range 0..180 degrees.

    The angle is measured between the normals of the two planes a-b-c and
    b-c-d (cross products of consecutive bond vectors along the chain).
    Degenerate cases return null with an explanatory reason:

    * any bond of the chain has zero length;
    * the three consecutive atoms a-b-c or b-c-d are collinear, which makes
      the corresponding plane normal undefined.
    """
    atoms = (a, b, c, d)
    bond_vectors = [
        _sub(positions[b], positions[a]),
        _sub(positions[c], positions[b]),
        _sub(positions[d], positions[c]),
    ]
    for (p, q), vector in zip(((a, b), (b, c), (c, d)), bond_vectors):
        if _norm(vector) == 0.0:
            return Measurement(
                "torsion", atoms, None, f"zero-length bond {p}-{q}"
            )

    normal_a = _cross(bond_vectors[0], bond_vectors[1])
    normal_b = _cross(bond_vectors[1], bond_vectors[2])
    if _norm(normal_a) == 0.0:
        return Measurement(
            "torsion", atoms, None, f"atoms {a}-{b}-{c} are collinear"
        )
    if _norm(normal_b) == 0.0:
        return Measurement(
            "torsion", atoms, None, f"atoms {b}-{c}-{d} are collinear"
        )
    return Measurement("torsion", atoms, _angle_between(normal_a, normal_b))


def measure_conformer(
    topology: Topology,
    sets: MeasurementSets,
    positions: dict[str, Vec3],
) -> list[Measurement]:
    """Compute every local measurement for one conformer."""
    results: list[Measurement] = []
    results.extend(measure_bond(positions, a, b) for a, b in sets.bonds)
    results.extend(
        measure_angle(positions, a, center, b)
        for a, center, b in sets.angles
    )
    results.extend(
        measure_torsion(positions, a, b, c, d)
        for a, b, c, d in sets.torsions
    )
    return results
