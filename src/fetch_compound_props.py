#!/usr/bin/env python3
"""Download ChEMBL molecule properties for every compound in the raw activity tables."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set

import pandas as pd
import requests

from chembl_api import fetch_molecules

DEFAULT_DATASET_ROOT = Path(__file__).resolve().parents[1] / "real-world_dataset"
DEFAULT_OUTPUT_NAME = "chembl_compound_properties.csv"

PROPERTY_FIELDS = [
    "molecule_chembl_id",
    "pref_name",
    "max_phase",
    "molecule_type",
    "first_approval",
    "usan_stem",
    "standard_inchi_key",
    "canonical_smiles",
    "mw_freebase",
    "full_mwt",
    "alogp",
    "cx_logp",
    "psa",
    "hba",
    "hbd",
    "num_ro5_violations",
    "ro3_pass",
    "acd_most_apka",
    "acd_most_bpka",
    "acd_logp",
    "acd_logd",
]


def collect_molecule_ids(dataset_root: Path) -> List[str]:
    """Return sorted unique parent molecule IDs across all raw activity tables."""
    ids: Set[str] = set()
    for raw_path in sorted(dataset_root.glob("*/chembl_*/activities_raw.csv")):
        frame = pd.read_csv(raw_path, usecols=["parent_molecule_chembl_id"], dtype=str)
        ids.update(frame["parent_molecule_chembl_id"].dropna().unique().tolist())
        print(f"  {raw_path.relative_to(dataset_root)}: cumulative {len(ids)} ids")
    return sorted(ids)


def load_done_ids(partial_path: Path) -> Set[str]:
    """Return molecule IDs already present in a partial output file."""
    if not partial_path.is_file():
        return set()
    frame = pd.read_csv(partial_path, usecols=["molecule_chembl_id"], dtype=str)
    return set(frame["molecule_chembl_id"].dropna())


def flatten_molecule(record: Dict[str, object]) -> Dict[str, object]:
    """Flatten one molecule record, expanding nested property/structure dicts."""
    flat: Dict[str, object] = {
        "molecule_chembl_id": record.get("molecule_chembl_id"),
        "pref_name": record.get("pref_name"),
        "max_phase": record.get("max_phase"),
        "molecule_type": (record.get("molecule_type") or {}).get("molecule_type")
        if isinstance(record.get("molecule_type"), dict)
        else record.get("molecule_type"),
        "first_approval": record.get("first_approval"),
        "usan_stem": record.get("usan_stem"),
    }
    structures = record.get("molecule_structures") or {}
    flat["standard_inchi_key"] = structures.get("standard_inchi_key") if isinstance(structures, dict) else None
    flat["canonical_smiles"] = structures.get("canonical_smiles") if isinstance(structures, dict) else None
    properties = record.get("molecule_properties") or {}
    if isinstance(properties, dict):
        for key in PROPERTY_FIELDS:
            if key in properties:
                flat[key] = properties.get(key)
    for key in PROPERTY_FIELDS:
        flat.setdefault(key, None)
    return flat


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Collect unique parent_molecule_chembl_id values from all "
            "activities_raw.csv tables and download their molecule properties."
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
        "--output-name",
        default=DEFAULT_OUTPUT_NAME,
        help=f"output CSV inside dataset_root (default: {DEFAULT_OUTPUT_NAME})",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="re-download even if the output already exists",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the molecule property download."""
    parser = build_argument_parser()
    arguments = parser.parse_args(argv)
    dataset_root = arguments.dataset_root.expanduser().resolve()
    output_path = dataset_root / arguments.output_name

    if output_path.exists() and not arguments.force:
        print(f"[skip] {output_path} already exists (use --force to re-download)")
        return 0

    if not dataset_root.is_dir():
        parser.error(f"dataset directory does not exist: {dataset_root}")

    print("Collecting molecule IDs from raw activity tables...")
    molecule_ids = collect_molecule_ids(dataset_root)
    if not molecule_ids:
        print("No molecule IDs found; run fetch_chembl_activities.py first.", file=sys.stderr)
        return 1

    partial_path = output_path.with_suffix(".partial.csv")
    done_ids = set() if arguments.force else load_done_ids(partial_path)
    if done_ids:
        print(f"Resuming: {len(done_ids)} molecules already in {partial_path.name}")
    todo_ids = [mid for mid in molecule_ids if mid not in done_ids]
    print(f"Found {len(molecule_ids)} unique parent molecule IDs; {len(todo_ids)} to fetch.")

    if todo_ids:
        write_header = not partial_path.is_file() or arguments.force
        if arguments.force and partial_path.is_file():
            partial_path.unlink()
            write_header = True
        with requests.Session() as session:
            session.headers.update({"User-Agent": "chembl-reference-fetch/1.0"})
            try:
                for record in fetch_molecules(todo_ids, session=session):
                    row = flatten_molecule(record)
                    with partial_path.open("a", newline="", encoding="utf-8") as handle:
                        writer = csv.DictWriter(handle, fieldnames=PROPERTY_FIELDS)
                        if write_header:
                            writer.writeheader()
                            write_header = False
                        writer.writerow(row)
            except Exception as exc:
                print(f"[FAILED] molecule download: {exc}", file=sys.stderr)
                print(f"Progress kept in {partial_path}; re-run to resume.", file=sys.stderr)
                return 1

    if not partial_path.is_file():
        print("Nothing to merge into the final properties table.", file=sys.stderr)
        return 1
    frame = pd.read_csv(partial_path, dtype=str, low_memory=False)
    frame = frame.drop_duplicates(subset=["molecule_chembl_id"]).sort_values("molecule_chembl_id")
    frame.to_csv(output_path, index=False)
    partial_path.unlink(missing_ok=True)
    print(f"[ok] {len(frame)} molecules -> {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
