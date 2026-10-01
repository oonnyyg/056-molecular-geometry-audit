"""Per-conformer analysis: global size, summaries and reference distance RMS."""

from __future__ import annotations

import math
from dataclasses import dataclass

from .coordinates import Vec3
from .geometry import Measurement
from .topology import Topology


@dataclass(frozen=True)
class PairDistanceDelta:
    """Result of the reference distance-difference analysis."""

    rms: float | None
    max_pair: tuple[str, str] | None
    max_abs_delta: float | None
    reason: str | None = None


def center_of_mass(topology: Topology, positions: dict[str, Vec3]) -> Vec3:
    """Mass-weighted centre (not the geometric centroid)."""
    total_mass = 0.0
    weighted: list[float] = [0.0, 0.0, 0.0]
    for atom_id in topology.atom_ids:
        mass = topology.masses[atom_id]
        total_mass += mass
        point = positions[atom_id]
        for axis in range(3):
            weighted[axis] += mass * point[axis]
    return (
        weighted[0] / total_mass,
        weighted[1] / total_mass,
        weighted[2] / total_mass,
    )


def radius_of_gyration(
    topology: Topology, positions: dict[str, Vec3], com: Vec3
) -> float:
    """Mass-weighted radius of gyration in angstrom."""
    total_mass = sum(topology.masses.values())
    accumulated = 0.0
    for atom_id in topology.atom_ids:
        point = positions[atom_id]
        dx = point[0] - com[0]
        dy = point[1] - com[1]
        dz = point[2] - com[2]
        accumulated += topology.masses[atom_id] * (
            dx * dx + dy * dy + dz * dz
        )
    return math.sqrt(accumulated / total_mass)


def summarize(values: list[float]) -> dict[str, float | None]:
    """Minimum, maximum and mean, or all null when there is no valid value."""
    if not values:
        return {"min": None, "max": None, "mean": None}
    return {
        "min": min(values),
        "max": max(values),
        "mean": sum(values) / len(values),
    }


def measurement_counts(
    measurements: list[Measurement],
) -> dict[str, dict[str, int]]:
    """Count valid vs. undefined measurements per measurement kind."""
    counts: dict[str, dict[str, int]] = {
        kind: {"valid": 0, "undefined": 0}
        for kind in ("bond_length", "bond_angle", "torsion")
    }
    for measurement in measurements:
        if measurement.value is None:
            counts[measurement.kind]["undefined"] += 1
        else:
            counts[measurement.kind]["valid"] += 1
    counts["total"] = {
        "valid": sum(bucket["valid"] for bucket in counts.values()),
        "undefined": sum(bucket["undefined"] for bucket in counts.values()),
    }
    return counts


def distance(point_a: Vec3, point_b: Vec3) -> float:
    dx = point_a[0] - point_b[0]
    dy = point_a[1] - point_b[1]
    dz = point_a[2] - point_b[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)


def reference_distance_deltas(
    topology: Topology,
    positions: dict[str, Vec3],
    reference: dict[str, Vec3],
) -> PairDistanceDelta:
    """RMS of all unordered inter-atom distance differences vs. reference.

    The metric compares internal pairwise distances, never raw coordinates,
    so it is invariant under any translation or rotation of the conformer.
    With fewer than two atoms no atom pair exists and null is returned.
    The pair with the largest absolute difference is reported; ties are
    resolved lexicographically by atom id (smaller pair sequence wins).
    """
    atom_ids = sorted(topology.atom_ids)
    if len(atom_ids) < 2:
        return PairDistanceDelta(
            None, None, None,
            "fewer than two atoms: no atom pair exists",
        )

    squared_sum = 0.0
    pair_count = 0
    max_pair: tuple[str, str] | None = None
    max_abs_delta = -1.0
    for i, atom_a in enumerate(atom_ids):
        for atom_b in atom_ids[i + 1 :]:
            observed = distance(positions[atom_a], positions[atom_b])
            reference_distance_value = distance(
                reference[atom_a], reference[atom_b]
            )
            delta = observed - reference_distance_value
            abs_delta = abs(delta)
            squared_sum += delta * delta
            pair_count += 1
            if abs_delta > max_abs_delta:
                max_abs_delta = abs_delta
                max_pair = (atom_a, atom_b)
    return PairDistanceDelta(
        math.sqrt(squared_sum / pair_count),
        max_pair,
        max_abs_delta,
    )
