"""
Download the NCI1 chemical graph dataset.

Primary source  : https://www.chrsmrrs.com/graphkerneldatasets/NCI1.zip  (TU Dortmund)
Target          : data/raw/NCI1/

The archive extracts into NCI1/ with plain-text files:
    NCI1_A.txt              - edge list
    NCI1_graph_indicator.txt
    NCI1_graph_labels.txt
    NCI1_node_labels.txt
    README.txt

NCI1 contains ~4110 chemical compounds screened for activity against
non-small-cell lung cancer.  Node labels encode atom type.

Note on SDF format
------------------
The original NCI compound-set SDF files were distributed by the NCI DTP program.
# TODO: Verify and add a direct, stable NCI SDF download URL:
#   https://dtp.cancer.gov/databases_tools/bulk_data.htm
#   (access may require navigating the NCI DTP download portal)
#
# Alternative mirror (if available):
#   https://raw.githubusercontent.com/CHEMPHY/datasets/main/NCI1/NCI1.sdf
# TODO: Replace with a confirmed stable URL before using in production.

Usage:
    python scripts/download/download_nci1.py
"""
from __future__ import annotations

import sys
import urllib.request
import zipfile
from pathlib import Path

DATASET_NAME = "NCI1"
SOURCE_URL = "https://www.chrsmrrs.com/graphkerneldatasets/NCI1.zip"

ROOT_DIR = Path(__file__).resolve().parents[2]
TARGET_DIR = ROOT_DIR / "data" / "raw" / DATASET_NAME



def download_file(url: str, dest: Path) -> None:
    """Download url to dest; skip silently if dest already exists."""
    if dest.exists():
        print(f"[SKIP] Already present: {dest.name}")
        return

    print(f"[DOWNLOAD] {url}")
    print(f"        -> {dest}")
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
        with urllib.request.urlopen(req, timeout=120) as response:
            dest.write_bytes(response.read())
        size_kb = dest.stat().st_size / 1024
        print(f"[OK] Saved {size_kb:.1f} KB")
    except Exception as exc:
        print(f"[ERROR] Download failed: {exc}")
        raise


def extract_zip(archive: Path, target_dir: Path) -> None:
    """Extract all contents of archive into target_dir."""
    print(f"[EXTRACT] {archive.name} -> {target_dir}/")
    with zipfile.ZipFile(archive, "r") as zf:
        zf.extractall(target_dir)
    n_files = sum(1 for _ in target_dir.rglob("*") if _.is_file())
    print(f"[OK] Extracted {n_files} file(s)")



def main() -> int:
    """Download and extract the NCI1 dataset."""
    print(f"=== Downloading {DATASET_NAME} ===")
    TARGET_DIR.mkdir(parents=True, exist_ok=True)

    archive = TARGET_DIR / f"{DATASET_NAME}.zip"

    try:
        download_file(SOURCE_URL, archive)
    except Exception:
        print("[ABORT] Could not download archive.")
        return 1

    try:
        extract_zip(archive, TARGET_DIR)
    except Exception as exc:
        print(f"[ERROR] Extraction failed: {exc}")
        return 1

    print(f"[DONE] {DATASET_NAME} is ready in {TARGET_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
