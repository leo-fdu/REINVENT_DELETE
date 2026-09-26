#!/usr/bin/env python3
"""Plot per-target density of Tanimoto similarity to the co-crystal ligand.

Reads every ``result/<target>/<target>_tanimoto_similarity.csv`` produced by
``crystal_actives_similarity.py`` and draws one Gaussian KDE curve per target.
Each density is normalised to integrate to one, so the y value at a given
similarity is the fraction of that target's ChEMBL actives per unit
similarity.  All curves are overlaid in a single figure saved next to the
CSV tables.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    from scipy import stats
except ImportError as exc:
    raise SystemExit(
        "This script requires matplotlib, numpy and scipy. Install them before running."
    ) from exc

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_RESULT_ROOT = SCRIPT_DIR.parent / "result"
DEFAULT_FIGURE_NAME = "tanimoto_similarity_density.png"

SIMILARITY_COLUMN = "tanimoto_similarity"
GRID_POINTS = 500


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


def kde_density(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Evaluate the Gaussian KDE of ``values`` on ``grid`` (integrates to one)."""
    if values.size == 1:
        density = np.zeros_like(grid)
        density[np.argmin(np.abs(grid - values[0]))] = 1.0
        return density
    return stats.gaussian_kde(values)(grid)


def plot_densities(
    similarities: Dict[str, np.ndarray],
    figure_path: Path,
    dpi: int,
) -> None:
    """Draw all per-target density curves in one figure."""
    grid = np.linspace(0.0, 1.0, GRID_POINTS)
    colors = plt.get_cmap("tab20").colors

    fig, ax = plt.subplots(figsize=(10, 6))
    for index, (target, values) in enumerate(sorted(similarities.items())):
        ax.plot(
            grid,
            kde_density(values, grid),
            color=colors[index % len(colors)],
            linewidth=1.5,
            label=f"{target} (n={values.size})",
        )

    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(bottom=0.0)
    ax.set_xlabel("Tanimoto similarity to co-crystal ligand (Morgan FP, r=2, 2048 bits)")
    ax.set_ylabel("Density (fraction of actives per unit similarity)")
    ax.set_title("ChEMBL actives similarity to the co-crystal ligand, 16 targets")
    ax.legend(fontsize=8, ncol=2, frameon=False)
    ax.grid(alpha=0.3, linewidth=0.5)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=dpi)
    plt.close(fig)


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Overlay per-target KDE density curves of the Tanimoto similarity "
            "between ChEMBL actives and the co-crystal ligand."
        )
    )
    parser.add_argument(
        "result_root",
        nargs="?",
        type=Path,
        default=DEFAULT_RESULT_ROOT,
        help="directory with per-target result folders (default: similarity_test/result)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="figure path (default: <result_root>/tanimoto_similarity_density.png)",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="figure resolution in dpi (default: 300)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the density plotting command."""
    parser = build_argument_parser()
    arguments = parser.parse_args(argv)
    result_root = arguments.result_root.expanduser().resolve()

    if not result_root.is_dir():
        parser.error(f"result directory does not exist: {result_root}")

    tables = find_similarity_tables(result_root)
    if not tables:
        print(f"No similarity CSV tables found below {result_root}", file=sys.stderr)
        return 1

    similarities = {}
    failures = []
    for target, csv_path in tables.items():
        try:
            similarities[target] = load_similarities(csv_path)
        except Exception as exc:
            failures.append((target, str(exc)))
            print(f"[FAILED] {target}: {exc}", file=sys.stderr)

    figure_path = arguments.output or (result_root / DEFAULT_FIGURE_NAME)
    plot_densities(similarities, figure_path, dpi=arguments.dpi)
    print(f"Plotted {len(similarities)}/{len(tables)} targets -> {figure_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
