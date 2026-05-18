from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PIPELINE = Path(__file__).resolve().parent / "run_pipeline.py"

FIXED_DATASETS = [
    "data/processed/MUTAG/mutag_smiles.smi",
    "data/processed/NCI1/nci1_smiles.smi",
    "data/processed/NCI109/nci109_smiles.smi",
    "data/raw/ZINC/zinc250k.smi",
]


def collect_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for rel in FIXED_DATASETS:
        p = root / rel
        if p.exists():
            files.append(p)
        else:
            print(f"  [SKIP] {rel} not found — run download/parse steps first")
    zinc20_dir = root / "data/raw/ZINC20"
    if zinc20_dir.is_dir():
        zinc20_files = sorted(zinc20_dir.glob("*.smi"))
        if zinc20_files:
            files.extend(zinc20_files)
        else:
            print(f"  [SKIP] No .smi files found in {zinc20_dir}")
    return files


def main() -> int:
    p = argparse.ArgumentParser(description="Run the reproducibility pipeline over all datasets.")
    p.add_argument("--workers", type=int, default=8, help="Worker processes per dataset (default: 8)")
    args = p.parse_args()

    files = collect_files(ROOT)
    if not files:
        print("No dataset files found. Run the download and parse steps first (README sections 2-3).")
        return 1

    print(f"Processing {len(files)} dataset file(s) with {args.workers} workers each.")
    failed: list[Path] = []
    for i, f in enumerate(files, 1):
        print(f"\n[{i}/{len(files)}] {f.relative_to(ROOT)}")
        result = subprocess.run(
            [sys.executable, str(PIPELINE), "--data", str(f), "--workers", str(args.workers)],
            cwd=str(ROOT),
        )
        if result.returncode != 0:
            print(f"  [FAIL] exited with code {result.returncode}")
            failed.append(f)

    if failed:
        print(f"\n{len(failed)} file(s) failed:")
        for f in failed:
            print(f"  {f.relative_to(ROOT)}")
        return 1

    print(f"\nAll {len(files)} file(s) processed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
