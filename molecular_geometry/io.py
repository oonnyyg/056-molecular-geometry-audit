"""Input loading and validation for topology JSON and coordinate CSV."""
from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass


class InputError(Exception):
    """Raised for any invalid user input; CLI exits non-zero."""


def _require_id(value, where):
    if not isinstance(value, str) or not value.strip():
        raise InputError(f"{where}: id must be a non-empty string")
    return value


@dataclass(frozen=True)
class Topology:
    atom_ids: tuple
    masses: dict
    bonds: tuple  # endpoint pairs sorted lexicographically


def load_topology(path):
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise InputError("topology JSON must be an object")
    atoms = data.get("atoms")
    bonds = data.get("bonds")
    if not isinstance(atoms, list) or not atoms:
        raise InputError("topology must contain a non-empty 'atoms' list")
    if not isinstance(bonds, list):
        raise InputError("topology must contain a 'bonds' list")

    masses = {}
    for i, entry in enumerate(atoms):
        if not isinstance(entry, dict):
            raise InputError(f"atoms[{i}] must be an object with id and mass")
        atom_id = _require_id(entry.get("id"), f"atoms[{i}].id")
        mass = entry.get("mass")
        if not isinstance(mass, (int, float)) or isinstance(mass, bool):
            raise InputError(f"atom {atom_id!r}: mass must be a number")
        mass = float(mass)
        if not math.isfinite(mass) or mass <= 0.0:
            raise InputError(f"atom {atom_id!r}: mass must be positive and finite")
        if atom_id in masses:
            raise InputError(f"duplicate atom id {atom_id!r}")
        masses[atom_id] = mass

    seen = set()
    norm_bonds = []
    for i, bond in enumerate(bonds):
        if (
            not isinstance(bond, (list, tuple))
            or len(bond) != 2
            or not all(isinstance(x, str) for x in bond)
        ):
            raise InputError(f"bonds[{i}] must be a pair of atom id strings")
        a, b = bond
        if a == b:
            raise InputError(f"bonds[{i}]: self bond {a!r}-{b!r} is not allowed")
        for atom in (a, b):
            if atom not in masses:
                raise InputError(f"bonds[{i}]: unknown atom {atom!r}")
        key = (a, b) if a < b else (b, a)
        if key in seen:
            raise InputError(f"duplicate bond {a!r}-{b!r} (reversed counts as same)")
        seen.add(key)
        norm_bonds.append(key)

    return Topology(
        atom_ids=tuple(sorted(masses)),
        masses=masses,
        bonds=tuple(sorted(norm_bonds)),
    )


UNIT_TO_ANGSTROM = {"angstrom": 1.0, "nm": 10.0}


def load_coordinates(path, topology, unit):
    """Return {conformer_id: {atom_id: (x, y, z)}} in angstrom."""
    if unit not in UNIT_TO_ANGSTROM:
        raise InputError(
            f"unit must be one of {sorted(UNIT_TO_ANGSTROM)}, got {unit!r}"
        )
    scale = UNIT_TO_ANGSTROM[unit]
    conformers = {}
    with open(path, "r", encoding="utf-8", newline="") as fh:
        reader = csv.DictReader(fh)
        required = {"conformer", "atom", "x", "y", "z"}
        if reader.fieldnames is None or not required.issubset(
            set(reader.fieldnames)
        ):
            raise InputError(
                "coordinate CSV must have columns: conformer, atom, x, y, z"
            )
        for lineno, row in enumerate(reader, start=2):
            conf = _require_id(row.get("conformer"), f"line {lineno} conformer")
            atom = _require_id(row.get("atom"), f"line {lineno} atom")
            if atom not in topology.masses:
                raise InputError(f"line {lineno}: unknown atom {atom!r}")
            coords = []
            for axis in ("x", "y", "z"):
                raw = row.get(axis)
                try:
                    value = float(raw)
                except (TypeError, ValueError):
                    raise InputError(
                        f"line {lineno}: coordinate {axis}={raw!r} is not a number"
                    ) from None
                if not math.isfinite(value):
                    raise InputError(
                        f"line {lineno}: coordinate {axis} is not finite"
                    )
                coords.append(value * scale)
            slot = conformers.setdefault(conf, {})
            if atom in slot:
                raise InputError(
                    f"line {lineno}: duplicate atom {atom!r} in conformer {conf!r}"
                )
            slot[atom] = tuple(coords)
    if not conformers:
        raise InputError("coordinate CSV contains no conformers")
    for conf in sorted(conformers):
        missing = set(topology.masses) - set(conformers[conf])
        if missing:
            raise InputError(
                f"conformer {conf!r} misses atoms: {sorted(missing)}"
            )
    return conformers
