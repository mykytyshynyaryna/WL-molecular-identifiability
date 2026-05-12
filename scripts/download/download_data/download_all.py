"""
Run all dataset download scripts in sequence.

Each individual script is called as a subprocess so that failures in one
dataset do not abort the remaining downloads.  Exit codes from each script
are collected and a summary is printed at the end.

Note: download_zinc_from_uri.py is not included here because it requires a
local URI list file. Run it separately when you have that file available.

Usage:
    python scripts/download/download_all.py
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

SCRIPTS_DIR = Path(__file__).resolve().parent

DOWNLOAD_SCRIPTS = [
    SCRIPTS_DIR / "download_mutag.py",
    SCRIPTS_DIR / "download_nci1.py",
    SCRIPTS_DIR / "download_nci109.py",
    SCRIPTS_DIR / "download_dd.py",
    SCRIPTS_DIR / "download_enzymes.py",
    SCRIPTS_DIR / "download_zinc.py",
]


def run_script(script: Path) -> int:
    """
    Execute script with the current Python interpreter.

    Returns the script's exit code (0 = success).
    """
    print(f"\n{'=' * 60}")
    print(f"Running: {script.name}")
    print("=" * 60)

    result = subprocess.run(
        [sys.executable, str(script)],
        check=False,
    )
    return result.returncode


def main() -> int:
    """Run all download scripts and print a summary."""
    print("=== download_all.py: starting all dataset downloads ===\n")

    results: dict[str, int] = {}

    for script in DOWNLOAD_SCRIPTS:
        if not script.exists():
            print(f"[WARN] Script not found, skipping: {script}")
            results[script.name] = -1
            continue

        code = run_script(script)
        results[script.name] = code

    print(f"\n{'=' * 60}")
    print("SUMMARY")
    print("=" * 60)

    n_ok = sum(1 for c in results.values() if c == 0)
    n_fail = len(results) - n_ok

    for name, code in results.items():
        status = "OK" if code == 0 else f"FAILED (exit {code})"
        print(f"  {name:<35} {status}")

    print(f"\nTotal: {n_ok} succeeded, {n_fail} failed.")

    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
