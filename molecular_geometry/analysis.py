"""Per-conformer measurements, summaries and reference comparison."""
import math
from itertools import combinations

from .geometry import bond_angle, distance, torsion_angle


def enumerate_measurements(topology):
    """Build the static measurement list from the topology."""
    neighbors = {a: set() for a in topology.atom_ids}
    for a, b in topology.bonds:
        neighbors[a].add(b)
        neighbors[b].add(a)
    bonds = list(topology.bonds)
    angles = set()
    for center, nbrs in neighbors.items():
        for a, b in combinations(sorted(nbrs), 2):
            angles.add((a, center, b))
    torsions = set()
    for b, c in topology.bonds:
        for a in neighbors[b] - {c}:
            for d in neighbors[c] - {b}:
                if a == d:
                    continue
                torsions.add(min((a, b, c, d), (d, c, b, a)))
    return bonds, sorted(angles), sorted(torsions)


def measure_conformer(coords, measurements):
    """Compute all measurements for one conformer (list of dicts)."""
    bonds, angles, torsions = measurements
    results = []
    for a, b in bonds:
        results.append({
            "kind": "bond", "atoms": (a, b),
            "value": distance(coords[a], coords[b]), "reason": None,
        })
    for a, c, b in angles:
        value = bond_angle(coords[a], coords[c], coords[b])
        results.append({
            "kind": "angle", "atoms": (a, c, b), "value": value,
            "reason": None if value is not None else "zero-length bond",
        })
    for chain in torsions:
        value, reason = torsion_angle(*[coords[a] for a in chain])
        results.append({
            "kind": "torsion", "atoms": chain,
            "value": value, "reason": reason,
        })
    return results


def center_of_mass(coords, masses):
    total = sum(masses.values())
    return tuple(
        sum(masses[a] * coords[a][i] for a in masses) / total for i in range(3)
    )


def radius_of_gyration(coords, masses):
    com = center_of_mass(coords, masses)
    total = sum(masses.values())
    acc = sum(
        masses[a] * sum((coords[a][i] - com[i]) ** 2 for i in range(3))
        for a in masses
    )
    return math.sqrt(acc / total)


def summarize(results):
    """Counts plus min/max/mean stats for bonds and valid angles."""
    summary = {"counts": {}, "bond_length": None, "angle": None}
    for kind in ("bond", "angle", "torsion"):
        items = [r for r in results if r["kind"] == kind]
        valid = [r["value"] for r in items if r["value"] is not None]
        summary["counts"][kind] = {
            "total": len(items),
            "valid": len(valid),
            "undefined": len(items) - len(valid),
        }
    for kind, key in (("bond", "bond_length"), ("angle", "angle")):
        valid = [
            r["value"]
            for r in results
            if r["kind"] == kind and r["value"] is not None
        ]
        if valid:
            summary[key] = {
                "min": min(valid),
                "max": max(valid),
                "mean": sum(valid) / len(valid),
            }
    return summary


def pairwise_distance_deviation(coords, ref_coords, atom_ids):
    """RMS of unordered atom-pair distance differences vs reference.

    Invariant to translation/rotation. Returns dict with rmsd and the
    pair (lexicographically first on ties) of maximal absolute
    difference, or a dict with reason when fewer than two atoms.
    """
    if len(atom_ids) < 2:
        return {
            "rmsd": None,
            "max_pair": None,
            "max_abs_diff": None,
            "reason": "fewer than two atoms; no pairs to compare",
        }
    sq_sum = 0.0
    count = 0
    best_pair = None
    best_diff = -1.0
    for a, b in combinations(sorted(atom_ids), 2):
        diff = abs(
            distance(coords[a], coords[b])
            - distance(ref_coords[a], ref_coords[b])
        )
        sq_sum += diff * diff
        count += 1
        if diff > best_diff:  # pairs iterated in lexicographic order
            best_diff = diff
            best_pair = (a, b)
    return {
        "rmsd": math.sqrt(sq_sum / count),
        "max_pair": list(best_pair),
        "max_abs_diff": best_diff,
        "reason": None,
    }
