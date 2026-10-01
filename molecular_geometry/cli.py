"""Command-line interface: python -m molecular_geometry."""
import argparse
import csv
import json
import os
import sys

from . import analysis, io


def build_parser():
    parser = argparse.ArgumentParser(
        prog="python -m molecular_geometry",
        description=(
            "Audit molecular conformation geometry: bond lengths, angles, "
            "torsions, size metrics and reference-conformer deviation."
        ),
    )
    parser.add_argument("--topology", required=True, help="topology JSON path")
    parser.add_argument("--coordinates", required=True, help="coordinate CSV path")
    parser.add_argument(
        "--unit",
        required=True,
        choices=sorted(io.UNIT_TO_ANGSTROM),
        help="unit of the input coordinates (converted to angstrom)",
    )
    parser.add_argument("--config", required=True, help="config JSON path")
    parser.add_argument("--output-dir", required=True, help="output directory")
    return parser


def load_config(path, conformers):
    with open(path, "r", encoding="utf-8") as fh:
        data = json.load(fh)
    if not isinstance(data, dict):
        raise io.InputError("config JSON must be an object")
    ref = data.get("reference_conformer")
    if not isinstance(ref, str) or not ref.strip():
        raise io.InputError("config must set 'reference_conformer' (string)")
    if ref not in conformers:
        raise io.InputError(
            f"reference conformer {ref!r} not present in coordinates"
        )
    return ref


def run(argv=None):
    args = build_parser().parse_args(argv)
    try:
        topology = io.load_topology(args.topology)
        conformers = io.load_coordinates(args.coordinates, topology, args.unit)
        reference = load_config(args.config, conformers)
    except (io.InputError, OSError, json.JSONDecodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1

    measurements = analysis.enumerate_measurements(topology)
    geometry_rows = []
    report = {
        "unit": "angstrom",
        "angle_unit": "degree",
        "reference_conformer": reference,
        "conformers": {},
    }
    for conf in sorted(conformers):
        coords = conformers[conf]
        results = analysis.measure_conformer(coords, measurements)
        for r in results:
            geometry_rows.append([
                conf,
                r["kind"],
                "-".join(r["atoms"]),
                "" if r["value"] is None else repr(r["value"]),
                "" if r["reason"] is None else r["reason"],
            ])
        deviation = analysis.pairwise_distance_deviation(
            coords, conformers[reference], topology.atom_ids
        )
        report["conformers"][conf] = {
            "center_of_mass": analysis.center_of_mass(coords, topology.masses),
            "radius_of_gyration": analysis.radius_of_gyration(
                coords, topology.masses
            ),
            "summary": analysis.summarize(results),
            "reference_deviation": deviation,
        }

    os.makedirs(args.output_dir, exist_ok=True)
    geometry_path = os.path.join(args.output_dir, "geometry.csv")
    with open(geometry_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["conformer", "kind", "atoms", "value", "reason"])
        writer.writerows(geometry_rows)
    analysis_path = os.path.join(args.output_dir, "analysis.json")
    with open(analysis_path, "w", encoding="utf-8") as fh:
        json.dump(report, fh, indent=2, ensure_ascii=False)
        fh.write("\n")
    return 0
