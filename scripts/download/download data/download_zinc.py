"""
Download a ZINC molecular dataset subset in .smi format.

This script downloads the ZINC250k subset — 250 000 drug-like molecules used
as a standard benchmark in molecular graph learning (e.g., Benchmarking GNNs,
junction tree VAE, and chemical VAE papers).

Primary source  : ZINC15 via the aspuru-guzik-group chemical_vae repository
                  (provides the 250k SMILES strings as a CSV with properties)
Target          : data/raw/ZINC/zinc250k.smi

Output format
-------------
A plain .smi file: one SMILES string per line (no header, no ID column).
The CSV download is automatically converted to .smi by stripping the header
and keeping only the first (SMILES) column.

# TODO: ZINC15 direct .smi endpoints (e.g. tranches) require selecting
#   specific subsets from https://zinc15.docking.org/tranches/home/ and
#   are session-based.  The CSV mirror below is more stable:
#     SOURCE_URL = (
#         "https://raw.githubusercontent.com/"
#         "aspuru-guzik-group/chemical_vae/master/"
#         "models/zinc_properties/250k_rndm_zinc_drugs_clean_3.csv"
#     )
#
# Alternative: download directly from ZINC15 tranches as .smi files by
#   selecting "2D" -> "SMILES" format in the download interface and
#   saving the resulting URI list, then using scripts/download/download_zinc_from_uri.py.

Usage:
    python scripts/download/download_zinc.py
"""
from __future__ import annotations

import csv
import io
import sys
import urllib.request
from pathlib import Path

DATASET_NAME = "ZINC"

SOURCE_URL = (
    "https://raw.githubusercontent.com/"
    "aspuru-guzik-group/chemical_vae/master/"
    "models/zinc_properties/250k_rndm_zinc_drugs_clean_3.csv"
)

ROOT_DIR = Path(__file__).resolve().parents[2]
TARGET_DIR = ROOT_DIR / "data" / "raw" / DATASET_NAME
SMI_FILE = TARGET_DIR / "zinc250k.smi"



def download_csv_as_bytes(url: str) -> bytes:
    """Fetch url and return raw bytes; raises on HTTP error."""
    print(f"[DOWNLOAD] {url}")
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req, timeout=180) as response:
        data = response.read()
    print(f"[OK] Received {len(data) / 1024:.1f} KB")
    return data


def csv_bytes_to_smi(data: bytes) -> str:
    """
    Parse a CSV byte-string and extract the SMILES column (column 0).

    Returns a newline-joined string of SMILES, one per line.
    """
    text = data.decode("utf-8", errors="replace")
    reader = csv.reader(io.StringIO(text))

    smiles_lines = []
    for i, row in enumerate(reader):
        if i == 0:
            if row and not _looks_like_smiles(row[0]):
                continue
        if row:
            smiles_lines.append(row[0].strip())

    return "\n".join(smiles_lines)


def _looks_like_smiles(token: str) -> bool:
    """Heuristic: SMILES strings contain element symbols or bracket characters."""
    return any(c in token for c in ("C", "N", "O", "S", "c", "n", "[", "("))


def save_smi(content: str, dest: Path) -> None:
    """Write the SMILES content to dest."""
    dest.write_text(content, encoding="utf-8")
    n_lines = content.count("\n") + 1
    print(f"[OK] Wrote {n_lines} SMILES to {dest}")



def main() -> int:
    """Download the ZINC250k CSV and convert it to a .smi file."""
    print(f"=== Downloading {DATASET_NAME} (250k subset) ===")

    if SMI_FILE.exists():
        print(f"[SKIP] Already present: {SMI_FILE}")
        print(f"[DONE] {DATASET_NAME} is ready in {TARGET_DIR}")
        return 0

    TARGET_DIR.mkdir(parents=True, exist_ok=True)

    try:
        raw_bytes = download_csv_as_bytes(SOURCE_URL)
    except Exception as exc:
        print(f"[ERROR] Download failed: {exc}")
        print(
            "[HINT] Check your internet connection or update SOURCE_URL in this script."
        )
        return 1

    print("[CONVERT] CSV -> .smi")
    try:
        smi_content = csv_bytes_to_smi(raw_bytes)
    except Exception as exc:
        print(f"[ERROR] CSV parsing failed: {exc}")
        return 1

    try:
        save_smi(smi_content, SMI_FILE)
    except Exception as exc:
        print(f"[ERROR] Could not write file: {exc}")
        return 1

    print(f"[DONE] {DATASET_NAME} is ready in {TARGET_DIR}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
