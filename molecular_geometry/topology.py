"""Topology loading and validation.

The topology is a JSON object describing atom identifiers, positive finite
masses, and undirected chemical bonds.  Bond orders are deliberately not
inferred or required: each bond is an unordered pair of distinct atom ids.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass

from .errors import GeometryInputError


@dataclass(frozen=True)
class Topology:
    """Validated molecular topology.

    Attributes:
        atom_ids: Atom identifiers in the order declared by the input.
        masses: Mapping of atom id to mass.
        bonds: Frozenset of undirected bonds, each stored as a sorted pair.
        adjacency: Atom id to lexicographically sorted neighbor ids.
    """

    atom_ids: tuple[str, ...]
    masses: dict[str, float]
    bonds: frozenset[tuple[str, str]]
    adjacency: dict[str, tuple[str, ...]]


def _load_json(path: str) -> object:
    try:
        with open(path, "r", encoding="utf-8") as handle:
            return json.load(handle)
    except FileNotFoundError as exc:
        raise GeometryInputError(f"topology file not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GeometryInputError(
            f"topology file is not valid JSON ({path}): {exc.msg}"
        ) from exc
    except OSError as exc:
        raise GeometryInputError(f"cannot read topology file {path}: {exc}") from exc


def _coerce_atom_entry(entry: object, index: int) -> tuple[str, float]:
    if isinstance(entry, dict):
        atom_id = entry.get("id", entry.get("name"))
        mass = entry.get("mass")
    elif isinstance(entry, (list, tuple)) and len(entry) == 2:
        atom_id, mass = entry
    else:
        raise GeometryInputError(
            f"atom entry #{index} must be an object with 'id' and 'mass'"
        )
    if not isinstance(atom_id, str) or not atom_id:
        raise GeometryInputError(f"atom entry #{index} has a missing or empty id")
    return atom_id, _coerce_mass(mass, atom_id)


def _coerce_mass(mass: object, atom_id: str) -> float:
    if isinstance(mass, bool) or not isinstance(mass, (int, float)):
        raise GeometryInputError(f"atom {atom_id!r} mass must be a number")
    value = float(mass)
    if not math.isfinite(value) or value <= 0.0:
        raise GeometryInputError(
            f"atom {atom_id!r} mass must be a positive finite number, got {mass!r}"
        )
    return value


def _coerce_bond(entry: object, index: int) -> tuple[str, str]:
    if isinstance(entry, dict):
        endpoints = entry.get("atoms", entry.get("ends"))
    else:
        endpoints = entry
    if not isinstance(endpoints, (list, tuple)) or len(endpoints) != 2:
        raise GeometryInputError(
            f"bond #{index} must list exactly two atom endpoints"
        )
    a, b = endpoints
    if not isinstance(a, str) or not isinstance(b, str):
        raise GeometryInputError(f"bond #{index} endpoints must be atom id strings")
    return tuple(sorted((a, b)))  # type: ignore[return-value]


def load_topology(path: str) -> Topology:
    """Load and validate a topology JSON file."""
    data = _load_json(path)
    if not isinstance(data, dict):
        raise GeometryInputError("topology JSON root must be an object")

    raw_atoms = data.get("atoms")
    if not isinstance(raw_atoms, list) or not raw_atoms:
        raise GeometryInputError(
            "topology must contain a non-empty 'atoms' list"
        )

    atom_ids: list[str] = []
    masses: dict[str, float] = {}
    seen: set[str] = set()
    for index, entry in enumerate(raw_atoms):
        atom_id, mass = _coerce_atom_entry(entry, index)
        if atom_id in seen:
            raise GeometryInputError(f"duplicate atom id {atom_id!r}")
        seen.add(atom_id)
        atom_ids.append(atom_id)
        masses[atom_id] = mass

    raw_bonds = data.get("bonds", [])
    if not isinstance(raw_bonds, list):
        raise GeometryInputError("topology 'bonds' must be a list")

    bonds: set[tuple[str, str]] = set()
    adjacency: dict[str, list[str]] = {atom_id: [] for atom_id in atom_ids}
    for index, entry in enumerate(raw_bonds):
        a, b = _coerce_bond(entry, index)
        if a == b:
            raise GeometryInputError(f"bond #{index} connects {a!r} to itself")
        if a not in seen or b not in seen:
            unknown = a if a not in seen else b
            raise GeometryInputError(
                f"bond #{index} references unknown atom {unknown!r}"
            )
        if (a, b) in bonds:
            raise GeometryInputError(
                f"duplicate bond between {a!r} and {b!r}"
                " (reversed endpoints count as the same bond)"
            )
        bonds.add((a, b))
        adjacency[a].append(b)
        adjacency[b].append(a)

    sorted_adjacency = {
        atom_id: tuple(sorted(neighbors)) for atom_id, neighbors in adjacency.items()
    }
    return Topology(
        atom_ids=tuple(atom_ids),
        masses=masses,
        bonds=frozenset(bonds),
        adjacency=sorted_adjacency,
    )
