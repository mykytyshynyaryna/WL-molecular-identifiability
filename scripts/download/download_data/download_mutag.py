"""
Download the MUTAG graph dataset from the TU Dortmund benchmark collection.

Source : https://www.chrsmrrs.com/graphkerneldatasets/MUTAG.zip
Target : data/raw/MUTAG/

The archive extracts into a folder named MUTAG/ with plain-text files:
    MUTAG_A.txt             - edge list
    MUTAG_graph_indicator.txt
    MUTAG_graph_labels.txt
    MUTAG_node_labels.txt
    README.txt

MUTAG contains 188 chemical compounds labelled by mutagenicity.
Node labels encode atom type; edge labels encode bond type.

Note on SDF format
------------------
If you need MUTAG in .sdf format, a curated SDF version is available through
the CDK (Chemistry Development Kit) project and various cheminformatics mirrors.
# TODO: Add a direct, stable SDF source URL here if needed.

Usage:
    python scripts/download/download_mutag.py
"""

from __future__ import annotations

import sys
import urllib.request
import zipfile
from pathlib import Path

DATASET_NAME = "MUTAG"
SOURCE_URL = "https://www.chrsmrrs.com/graphkerneldatasets/MUTAG.zip"

ROOT_DIR = Path(__file__).resolve().parents[3]
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
    """Download and extract the MUTAG dataset."""
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
