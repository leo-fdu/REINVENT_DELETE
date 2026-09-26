#!/usr/bin/env python3
"""Summarise per-target Tanimoto similarity distributions.

Reads every ``results/<target>/<target>_tanimoto_similarity.csv`` produced by
``crystal_actives_similarity.py`` and writes one CSV row per target with the
mean, median, (sample) variance, maximum, minimum, 90th percentile and 75th
percentile of the similarities to the co-crystal ligand.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

try:
    import numpy as np
except ImportError as exc:
    raise SystemExit("This script requires numpy. Install it before running.") from exc

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RESULT_ROOT = SCRIPT_DIR.parent / "results"
DEFAULT_SUMMARY_NAME = "similarity_summary_stats.csv"

SIMILARITY_COLUMN = "tanimoto_similarity"
SUMMARY_COLUMNS = (
    "target",
    "n_actives",
    "mean",
    "median",
    "variance",
    "max",
    "min",
    "p90",
    "p75",
)
DECIMALS = 6


def find_similarity_tables(result_root: Path) -> Dict[str, Path]:
    """Return {target: csv path} for every per-target similarity table."""
    tables = {}
    for child in sorted(result_root.iterdir()):
        if not child.is_dir():
            continue
        csv_path = child / f"{child.name}_tanimoto_similarity.csv"
        if csv_path.is_file():
            tables[child.name] = csv_path
    return tables


def load_similarities(csv_path: Path) -> np.ndarray:
    """Load the similarity column of one target's CSV table."""
    values: List[float] = []
    with csv_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            values.append(float(row[SIMILARITY_COLUMN]))
    if not values:
        raise ValueError(f"no similarities found in {csv_path}")
    return np.asarray(values)


def summarise_target(target: str, values: np.ndarray) -> Tuple:
    """Compute the summary statistics of one target's similarities."""
    variance = float(np.var(values, ddof=1)) if values.size > 1 else 0.0
    return (
        target,
        values.size,
        round(float(np.mean(values)), DECIMALS),
        round(float(np.median(values)), DECIMALS),
        round(variance, DECIMALS),
        round(float(np.max(values)), DECIMALS),
        round(float(np.min(values)), DECIMALS),
        round(float(np.percentile(values, 90)), DECIMALS),
        round(float(np.percentile(values, 75)), DECIMALS),
    )


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Write a per-target summary (mean, median, variance, max, min, "
            "p90, p75) of the Tanimoto similarity between ChEMBL actives and "
            "the co-crystal ligand."
        )
    )
    parser.add_argument(
        "result_root",
        nargs="?",
        type=Path,
        default=DEFAULT_RESULT_ROOT,
        help="directory with per-target result folders (default: similarity_test/results)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="summary CSV path (default: <result_root>/similarity_summary_stats.csv)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the summary command."""
    parser = build_argument_parser()
    arguments = parser.parse_args(argv)
    result_root = arguments.result_root.expanduser().resolve()

    if not result_root.is_dir():
        parser.error(f"result directory does not exist: {result_root}")

    tables = find_similarity_tables(result_root)
    if not tables:
        print(f"No similarity CSV tables found below {result_root}", file=sys.stderr)
        return 1

    rows = []
    failures = []
    for target, csv_path in tables.items():
        try:
            rows.append(summarise_target(target, load_similarities(csv_path)))
        except Exception as exc:
            failures.append((target, str(exc)))
            print(f"[FAILED] {target}: {exc}", file=sys.stderr)

    output_path = arguments.output or (result_root / DEFAULT_SUMMARY_NAME)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(SUMMARY_COLUMNS)
        writer.writerows(rows)
    print(f"Summarised {len(rows)}/{len(tables)} targets -> {output_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
