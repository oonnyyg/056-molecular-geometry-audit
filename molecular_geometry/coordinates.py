"""Coordinate table and configuration loading.

Coordinate CSV files list one atom position per row::

    conformer,atom,x,y,z
    c1,C,0.0,0.0,0.0

Header names are matched case-insensitively and may use ``conf``/``id`` for
the conformer column.  Coordinates are given globally in angstrom (default)
or nanometres, selected via the configuration; nanometres are converted to
angstrom so every reported number uses angstrom.
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import dataclass

from .errors import GeometryInputError
from .topology import Topology

Vec3 = tuple[float, float, float]

_CONFORMER_HEADERS = ("conformer", "conformer_id", "conf", "conf_id", "id")
_ATOM_HEADERS = ("atom", "atom_id", "name")
_UNIT_ALIASES = {"a": "angstrom", "angstrom": "angstrom", "angstroms": "angstrom",
                 "nm": "nanometer", "nanometer": "nanometer",
                 "nanometers": "nanometer", "nanometre": "nanometer",
                 "nanometres": "nanometer"}


@dataclass(frozen=True)
class Config:
    """Run configuration."""

    reference_conformer: str
    coordinates_unit: str


@dataclass(frozen=True)
class CoordinateTable:
    """Coordinates grouped by conformer, ordered by lexicographic conformer id."""

    conformers: tuple[str, ...]
    positions: dict[str, dict[str, Vec3]]
    input_unit: str


def load_config(path: str) -> Config:
    """Load the run configuration JSON.

    Required key: ``reference_conformer`` (or ``reference``).
    Optional key: ``coordinates_unit`` (``"angstrom"``/``"nm"``).
    """
    try:
        with open(path, "r", encoding="utf-8") as handle:
            data = json.load(handle)
    except FileNotFoundError as exc:
        raise GeometryInputError(f"config file not found: {path}") from exc
    except json.JSONDecodeError as exc:
        raise GeometryInputError(
            f"config file is not valid JSON ({path}): {exc.msg}"
        ) from exc
    except OSError as exc:
        raise GeometryInputError(f"cannot read config file {path}: {exc}") from exc

    if not isinstance(data, dict):
        raise GeometryInputError("config JSON root must be an object")
    reference = data.get("reference_conformer", data.get("reference"))
    if not isinstance(reference, str) or not reference:
        raise GeometryInputError(
            "config must set non-empty 'reference_conformer'"
        )
    raw_unit = data.get("coordinates_unit", data.get("unit", "angstrom"))
    if not isinstance(raw_unit, str):
        raise GeometryInputError("'coordinates_unit' must be a string")
    unit = _UNIT_ALIASES.get(raw_unit.strip().lower())
    if unit is None:
        raise GeometryInputError(
            f"unsupported coordinates_unit {raw_unit!r}; use 'angstrom' or 'nm'"
        )
    return Config(reference_conformer=reference, coordinates_unit=unit)


def _resolve_columns(fieldnames: object) -> tuple[int, int, int, int, int]:
    if not fieldnames:
        raise GeometryInputError("coordinate CSV is empty (no header row)")
    lowered = [name.strip().lower() if name is not None else ""
               for name in fieldnames]

    def find(candidates: tuple[str, ...], label: str) -> int:
        for candidate in candidates:
            if candidate in lowered:
                return lowered.index(candidate)
        raise GeometryInputError(f"coordinate CSV header missing {label} column")

    return (
        find(_CONFORMER_HEADERS, "conformer id"),
        find(_ATOM_HEADERS, "atom id"),
        find(("x",), "x"),
        find(("y",), "y"),
        find(("z",), "z"),
    )


def _parse_float(token: str, conformer: str, atom: str, axis: str) -> float:
    try:
        value = float(token.strip())
    except (ValueError, AttributeError) as exc:
        raise GeometryInputError(
            f"non-numeric {axis} coordinate for conformer {conformer!r}, "
            f"atom {atom!r}: {token!r}"
        ) from exc
    if not math.isfinite(value):
        raise GeometryInputError(
            f"non-finite {axis} coordinate for conformer {conformer!r}, "
            f"atom {atom!r}"
        )
    return value


def load_coordinates(path: str, topology: Topology, unit: str) -> CoordinateTable:
    """Load, validate and unit-convert a coordinate CSV.

    Every conformer must cover exactly the atoms declared in the topology:
    an empty table, duplicate rows, missing atoms or non-finite coordinates
    are rejected.
    """
    try:
        handle = open(path, "r", encoding="utf-8", newline="")
    except FileNotFoundError as exc:
        raise GeometryInputError(f"coordinate file not found: {path}") from exc
    except OSError as exc:
        raise GeometryInputError(
            f"cannot read coordinate file {path}: {exc}"
        ) from exc

    with handle:
        reader = csv.reader(handle)
        try:
            header = next(reader)
        except StopIteration as exc:
            raise GeometryInputError("coordinate CSV is empty") from exc
        i_conf, i_atom, i_x, i_y, i_z = _resolve_columns(header)

        rows: dict[str, dict[str, Vec3]] = {}
        seen_pairs: set[tuple[str, str]] = set()
        row_no = 1
        for row in reader:
            row_no += 1
            if not row or all(not cell.strip() for cell in row):
                continue  # tolerate blank lines
            if len(row) < len(header):
                raise GeometryInputError(
                    f"coordinate CSV row {row_no} has too few columns"
                )
            conformer = row[i_conf].strip()
            atom = row[i_atom].strip()
            if not conformer:
                raise GeometryInputError(
                    f"coordinate CSV row {row_no} has an empty conformer id"
                )
            if not atom:
                raise GeometryInputError(
                    f"coordinate CSV row {row_no} has an empty atom id"
                )
            if atom not in topology.masses:
                raise GeometryInputError(
                    f"coordinate CSV row {row_no} references unknown atom {atom!r}"
                )
            pair = (conformer, atom)
            if pair in seen_pairs:
                raise GeometryInputError(
                    f"duplicate coordinates for conformer {conformer!r}, "
                    f"atom {atom!r} (row {row_no})"
                )
            seen_pairs.add(pair)
            point = (
                _parse_float(row[i_x], conformer, atom, "x"),
                _parse_float(row[i_y], conformer, atom, "y"),
                _parse_float(row[i_z], conformer, atom, "z"),
            )
            rows.setdefault(conformer, {})[atom] = point

    if not rows:
        raise GeometryInputError("coordinate CSV contains no data rows")

    factor = 10.0 if unit == "nanometer" else 1.0
    expected = set(topology.atom_ids)
    for conformer, positions in rows.items():
        present = set(positions)
        if present != expected:
            missing = sorted(expected - present)
            if missing:
                raise GeometryInputError(
                    f"conformer {conformer!r} is missing coordinates for "
                    f"atom(s): {', '.join(missing)}"
                )
    converted = {
        conformer: {
            atom: (x * factor, y * factor, z * factor)
            for atom, (x, y, z) in positions.items()
        }
        for conformer, positions in rows.items()
    }
    return CoordinateTable(
        conformers=tuple(sorted(rows)),
        positions=converted,
        input_unit=unit,
    )
