"""Orchestrate the audit and write the geometry CSV and analysis JSON."""

from __future__ import annotations

import csv
import json
import os

from .analysis import (
    center_of_mass,
    measurement_counts,
    radius_of_gyration,
    reference_distance_deltas,
    summarize,
)
from .coordinates import (
    CoordinateTable,
    load_config,
    load_coordinates,
)
from .errors import GeometryInputError
from .geometry import (
    build_measurement_sets,
    measure_conformer,
)
from .topology import load_topology

GEOMETRY_CSV_NAME = "geometry.csv"
ANALYSIS_JSON_NAME = "analysis.json"


def _measurement_sort_key(measurement) -> tuple:
    return (measurement.atoms, measurement.kind)


def _prepare_output_dir(output_dir: str) -> None:
    if os.path.exists(output_dir) and not os.path.isdir(output_dir):
        raise GeometryInputError(
            f"output path is not a directory: {output_dir}"
        )
    os.makedirs(output_dir, exist_ok=True)


def run_audit(
    coordinates_path: str,
    topology_path: str,
    config_path: str,
    output_dir: str,
) -> None:
    """Validate all inputs, compute the audit, and atomically write reports.

    All validation and computation happens before any file is created, so an
    invalid run never creates or partially overwrites report files.
    """
    config = load_config(config_path)
    topology = load_topology(topology_path)
    table = load_coordinates(
        coordinates_path, topology, config.coordinates_unit
    )

    if config.reference_conformer not in table.positions:
        raise GeometryInputError(
            f"reference conformer {config.reference_conformer!r} is not "
            "present in the coordinate table"
        )

    sets = build_measurement_sets(topology)
    reference_positions = table.positions[config.reference_conformer]

    geometry_rows: list[dict] = []
    conformer_reports: dict[str, dict] = {}

    for conformer in table.conformers:
        positions = table.positions[conformer]
        measurements = measure_conformer(topology, sets, positions)
        measurements.sort(key=_measurement_sort_key)

        com = center_of_mass(topology, positions)
        rgyr = radius_of_gyration(topology, positions, com)

        bond_values = [
            m.value for m in measurements
            if m.kind == "bond_length" and m.value is not None
        ]
        angle_values = [
            m.value for m in measurements
            if m.kind == "bond_angle" and m.value is not None
        ]
        counts = measurement_counts(measurements)
        deltas = reference_distance_deltas(
            topology, positions, reference_positions
        )

        conformer_reports[conformer] = {
            "center_of_mass_angstrom": list(com),
            "radius_of_gyration_angstrom": rgyr,
            "measurement_counts": counts,
            "bond_length_angstrom": summarize(bond_values),
            "bond_angle_degrees": summarize(angle_values),
            "distance_difference_rms_angstrom": deltas.rms,
            "max_distance_difference": {
                "pair": list(deltas.max_pair) if deltas.max_pair else None,
                "absolute_difference_angstrom": deltas.max_abs_delta,
                "reason": deltas.reason,
            },
        }

        for measurement in measurements:
            geometry_rows.append({
                "conformer": conformer,
                "measurement": measurement.kind,
                "atoms": json.dumps(list(measurement.atoms)),
                "value": "" if measurement.value is None
                else repr(measurement.value),
                "unit": _unit_for(measurement.kind),
                "undefined_reason": measurement.undefined_reason or "",
            })

    analysis = {
        "schema_version": 1,
        "units": {
            "length": "angstrom",
            "angle": "degree",
            "input_coordinates": table.input_unit,
        },
        "geometry_conventions": {
            "bond_length": "Euclidean distance between bonded atoms",
            "bond_angle": (
                "angle between outward vectors from the centre atom, "
                "0..180 degrees"
            ),
            "torsion": (
                "unsigned angle between the a-b-c and b-c-d plane normals "
                "of the three-bond chain a-b-c-d, 0..180 degrees; the "
                "reversed chain is reported once"
            ),
            "radius_of_gyration": (
                "mass-weighted RMS distance to the mass centre"
            ),
            "distance_difference_rms": (
                "RMS over all unordered atom pairs of the difference "
                "between within-conformer and reference pairwise "
                "distances; invariant to translation and rotation; not a "
                "chirality test or an energy score"
            ),
        },
        "reference_conformer": config.reference_conformer,
        "conformer_order": list(table.conformers),
        "conformers": conformer_reports,
    }

    _prepare_output_dir(output_dir)
    _atomic_write_geometry_csv(
        os.path.join(output_dir, GEOMETRY_CSV_NAME), geometry_rows
    )
    _atomic_write_analysis_json(
        os.path.join(output_dir, ANALYSIS_JSON_NAME), analysis
    )


def _unit_for(kind: str) -> str:
    if kind == "bond_length":
        return "angstrom"
    return "degree"


def _atomic_write_bytes(path: str, payload: bytes) -> None:
    tmp_path = path + ".tmp"
    with open(tmp_path, "wb") as handle:
        handle.write(payload)
    os.replace(tmp_path, path)


def _atomic_write_geometry_csv(path: str, rows: list[dict]) -> None:
    from io import StringIO

    buffer = StringIO(newline="")
    fieldnames = [
        "conformer", "measurement", "atoms", "value",
        "unit", "undefined_reason",
    ]
    writer = csv.DictWriter(buffer, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)
    _atomic_write_bytes(path, buffer.getvalue().encode("utf-8"))


def _atomic_write_analysis_json(path: str, analysis: dict) -> None:
    payload = json.dumps(analysis, indent=2, ensure_ascii=False, sort_keys=False)
    _atomic_write_bytes(path, (payload + "\n").encode("utf-8"))
