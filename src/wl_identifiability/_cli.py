"""Command-line entry point for quick single-molecule identifiability checks.

Usage
-----
    pixi run wl-check "CCCN"
    pixi run wl-check "CCCN" "c1ccccc1" "CC(=O)O"
"""

from __future__ import annotations

import argparse
import sys

from .experiments import is_smi_identifiable


def main() -> None:
    parser = argparse.ArgumentParser(
        prog="wl-check",
        description="Check whether one or more molecules are identifiable by 1-WL.",
        epilog='Example: wl-check "CCCN" "c1ccccc1"',
    )
    parser.add_argument(
        "smiles",
        nargs="+",
        metavar="SMILES",
        help="One or more SMILES strings to check.",
    )
    args = parser.parse_args()

    exit_code = 0
    for smi in args.smiles:
        try:
            result = is_smi_identifiable(smi)
            label = "identifiable" if result else "NOT identifiable"
            print(f"{smi}: {label}")
        except ValueError as exc:
            print(f"{smi}: ERROR — {exc}", file=sys.stderr)
            exit_code = 1

    sys.exit(exit_code)
