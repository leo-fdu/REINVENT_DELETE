#!/usr/bin/env python3
"""Draw 2D molecular depictions for MOL2 files in a dataset."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
except ImportError as exc:
    raise SystemExit("This script requires RDKit. Install it before running.") from exc


DEFAULT_DATASET_ROOT = Path(__file__).resolve().parents[1] / "real-world_dataset"
DEFAULT_IMAGE_SIZE = (800, 600)


def positive_integer(value: str) -> int:
    """Parse a strictly positive integer for an argparse option."""
    parsed_value = int(value)
    if parsed_value <= 0:
        raise argparse.ArgumentTypeError("value must be greater than zero")
    return parsed_value


def find_mol2_files(dataset_root: Path) -> List[Path]:
    """Return all MOL2 files below the dataset root in stable order."""
    return sorted(
        path
        for path in dataset_root.rglob("*")
        if path.is_file() and path.suffix.lower() == ".mol2"
    )


def draw_mol2_file(
    mol2_path: Path,
    image_path: Path,
    image_size: Tuple[int, int],
) -> None:
    """Read one MOL2 file and write its 2D depiction as a PNG file."""
    molecule = Chem.MolFromMol2File(
        str(mol2_path),
        sanitize=True,
        removeHs=True,
    )
    if molecule is None:
        molecule = Chem.MolFromMol2File(
            str(mol2_path),
            sanitize=False,
            removeHs=True,
        )
    if molecule is None:
        raise ValueError("RDKit could not parse the MOL2 file")
    if molecule.GetNumAtoms() == 0:
        raise ValueError("the MOL2 file contains no atoms")

    AllChem.Compute2DCoords(molecule)
    Draw.MolToFile(
        molecule,
        str(image_path),
        size=image_size,
        imageType="png",
    )


def render_dataset(
    dataset_root: Path,
    image_size: Tuple[int, int],
) -> Tuple[int, List[Tuple[Path, str]]]:
    """Render every MOL2 file and return successes plus failures."""
    mol2_files = find_mol2_files(dataset_root)
    failures: List[Tuple[Path, str]] = []
    rendered_count = 0

    for mol2_path in mol2_files:
        image_path = mol2_path.with_name(f"{mol2_path.stem}_2d.png")
        try:
            draw_mol2_file(mol2_path, image_path, image_size)
        except Exception as exc:
            failures.append((mol2_path, str(exc)))
            print(f"[FAILED] {mol2_path}: {exc}", file=sys.stderr)
        else:
            rendered_count += 1
            print(f"[OK] {image_path}")

    return rendered_count, failures


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Read all MOL2 files below a dataset directory and save 2D PNG "
            "depictions next to the source files."
        )
    )
    parser.add_argument(
        "dataset_root",
        nargs="?",
        type=Path,
        default=DEFAULT_DATASET_ROOT,
        help="dataset root containing target folders (default: real-world_dataset)",
    )
    parser.add_argument(
        "--width",
        type=positive_integer,
        default=DEFAULT_IMAGE_SIZE[0],
        help=f"PNG width in pixels (default: {DEFAULT_IMAGE_SIZE[0]})",
    )
    parser.add_argument(
        "--height",
        type=positive_integer,
        default=DEFAULT_IMAGE_SIZE[1],
        help=f"PNG height in pixels (default: {DEFAULT_IMAGE_SIZE[1]})",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the MOL2 rendering command."""
    parser = build_argument_parser()
    arguments = parser.parse_args(argv)
    dataset_root = arguments.dataset_root.expanduser().resolve()

    if not dataset_root.is_dir():
        parser.error(f"dataset directory does not exist: {dataset_root}")

    mol2_files = find_mol2_files(dataset_root)
    if not mol2_files:
        print(f"No MOL2 files found below {dataset_root}", file=sys.stderr)
        return 1

    rendered_count, failures = render_dataset(
        dataset_root,
        image_size=(arguments.width, arguments.height),
    )
    print(
        f"Rendered {rendered_count}/{len(mol2_files)} MOL2 files "
        f"under {dataset_root}."
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
