"""
Download ZINC .smi tranche files from a local URI list file.

This script reads a ZINC URI list file (downloaded from the ZINC20 download
interface) and fetches each individual .smi tranche file from the URLs it
contains.  Downloads run in parallel (--workers, default 8) so 1000+ files
complete in minutes rather than hours.

It is distinct from download_zinc.py, which downloads the ZINC250k benchmark
subset directly from GitHub.

Usage (from project root):
    python scripts/download/download_zinc_from_uri.py
    python scripts/download/download_zinc_from_uri.py --uri ZINC-downloader-2D-smi.uri
    python scripts/download/download_zinc_from_uri.py --uri path/to/file.uri --n 100
    python scripts/download/download_zinc_from_uri.py --workers 16
"""

from __future__ import annotations

import argparse
import sys
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock

ROOT_DIR = Path(__file__).resolve().parents[3]
DEFAULT_URI = ROOT_DIR / "ZINC-downloader-2D-smi.uri"
DEFAULT_OUT = ROOT_DIR / "data" / "raw" / "ZINC20"
DEFAULT_WORKERS = 3
DEFAULT_RETRIES = 5
DEFAULT_BACKOFF = 2.0


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Download all ZINC .smi tranche files listed in a .uri file.")
    p.add_argument(
        "--uri",
        default=str(DEFAULT_URI),
        help=f"Path to the ZINC URI list file (default: {DEFAULT_URI.name} in project root).",
    )
    p.add_argument(
        "--n",
        type=int,
        default=None,
        help="Limit to the first N files (default: all).",
    )
    p.add_argument(
        "--out",
        default=str(DEFAULT_OUT),
        help=f"Output directory (default: {DEFAULT_OUT.relative_to(ROOT_DIR)}).",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=f"Number of parallel download threads (default: {DEFAULT_WORKERS}).",
    )
    p.add_argument(
        "--retries",
        type=int,
        default=DEFAULT_RETRIES,
        help=f"Max attempts per file before giving up (default: {DEFAULT_RETRIES}).",
    )
    p.add_argument(
        "--backoff",
        type=float,
        default=DEFAULT_BACKOFF,
        help=f"Initial retry wait in seconds, doubled each attempt (default: {DEFAULT_BACKOFF}).",
    )
    return p.parse_args()


_print_lock = Lock()


def _download_one(
    url: str,
    out_dir: Path,
    idx: int,
    total: int,
    retries: int = DEFAULT_RETRIES,
    backoff: float = DEFAULT_BACKOFF,
) -> tuple[str, bool, str]:
    """Download a single URL to out_dir with exponential-backoff retries.

    Returns (filename, ok, message).
    """
    filename = url.split("/")[-1]
    dest = out_dir / filename

    if dest.exists():
        with _print_lock:
            print(f"  [{idx:>4}/{total}] SKIP (exists)  {filename}")
        return filename, True, "skip"

    last_exc: str = "unknown"
    wait = backoff

    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(
                url,
                headers={
                    "User-Agent": (
                        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 Chrome/124.0 Safari/537.36"
                    )
                },
            )
            with urllib.request.urlopen(req, timeout=60) as response:
                dest.write_bytes(response.read())
            with _print_lock:
                suffix = f" (attempt {attempt})" if attempt > 1 else ""
                print(f"  [{idx:>4}/{total}] OK            {filename}{suffix}")
            return filename, True, "ok"

        except Exception as exc:
            last_exc = str(exc)
            if attempt < retries:
                with _print_lock:
                    print(f"  [{idx:>4}/{total}] RETRY {attempt}/{retries}  {filename}  ({exc})  — waiting {wait:.0f}s")
                time.sleep(wait)
                wait *= 2
            else:
                with _print_lock:
                    print(f"  [{idx:>4}/{total}] FAIL          {filename}  ({exc})")

    return filename, False, last_exc


def main() -> int:
    args = parse_args()

    uri_path = Path(args.uri)
    if not uri_path.exists():
        print(f"ERROR: URI file not found: {uri_path}")
        print("  Make sure the file exists or pass --uri <path>")
        return 1

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    urls = [
        line.strip().replace("http://", "https://", 1)
        for line in uri_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    if args.n is not None:
        urls = urls[: args.n]

    total = len(urls)
    print(f"URI file : {uri_path}")
    print(f"Output   : {out_dir}")
    print(f"Files    : {total}")
    print(f"Workers  : {args.workers}")
    print(f"Retries  : {args.retries}  (backoff starts at {args.backoff}s, doubles each attempt)")
    print()

    ok_count = 0
    failed: list[str] = []

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {
            pool.submit(_download_one, url, out_dir, i + 1, total, args.retries, args.backoff): url
            for i, url in enumerate(urls)
        }
        for future in as_completed(futures):
            _filename, success, _msg = future.result()
            if success:
                ok_count += 1
            else:
                failed.append(futures[future])

    print(f"\nDone: {ok_count}/{total} succeeded, {len(failed)} failed.")
    if failed:
        print("Failed URLs:")
        for u in failed:
            print(f"  {u}")
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
