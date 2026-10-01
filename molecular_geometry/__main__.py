"""Command line entry point: ``python -m molecular_geometry ...``"""

from __future__ import annotations

import argparse
import sys

from .errors import GeometryInputError
from .report import ANALYSIS_JSON_NAME, GEOMETRY_CSV_NAME, run_audit


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="python -m molecular_geometry",
        description=(
            "Audit molecular conformer geometry: compute bond lengths, bond "
            "angles, torsions, mass-weighted size and reference distance "
            "differences, then write geometry.csv and analysis.json."
        ),
    )
    parser.add_argument(
        "--coordinates", required=True,
        help="conformer coordinate CSV (columns: conformer, atom, x, y, z)",
    )
    parser.add_argument(
        "--topology", required=True,
        help="topology JSON with non-empty 'atoms' and optional 'bonds'",
    )
    parser.add_argument(
        "--config", required=True,
        help=(
            "config JSON with 'reference_conformer' and optional "
            "'coordinates_unit' ('angstrom' or 'nm')"
        ),
    )
    parser.add_argument(
        "--output-dir", required=True,
        help="directory receiving geometry.csv and analysis.json",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        run_audit(
            coordinates_path=args.coordinates,
            topology_path=args.topology,
            config_path=args.config,
            output_dir=args.output_dir,
        )
    except GeometryInputError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(
        f"wrote {GEOMETRY_CSV_NAME} and {ANALYSIS_JSON_NAME} "
        f"to {args.output_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
