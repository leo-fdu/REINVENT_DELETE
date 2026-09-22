#!/usr/bin/env python3
"""Download all ChEMBL activity records for the 16 benchmark targets."""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import requests

from chembl_api import get_json, iter_activity_pages

DEFAULT_DATASET_ROOT = Path(__file__).resolve().parents[1] / "real-world_dataset"
DEFAULT_TARGETS_FILE = "chembl_targets.csv"

ACTIVITY_FIELDS = [
    "activity_id",
    "molecule_chembl_id",
    "parent_molecule_chembl_id",
    "canonical_smiles",
    "molecule_pref_name",
    "standard_type",
    "standard_relation",
    "standard_value",
    "standard_units",
    "pchembl_value",
    "assay_chembl_id",
    "assay_type",
    "assay_description",
    "bao_format",
    "bao_label",
    "document_chembl_id",
    "document_journal",
    "document_year",
    "data_validity_comment",
    "potential_duplicate",
    "standard_flag",
    "target_chembl_id",
    "target_pref_name",
    "target_organism",
    "ligand_efficiency",
]


def load_targets(targets_path: Path) -> List[Dict[str, str]]:
    """Read the target metadata table."""
    with targets_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def output_dir(dataset_root: Path, target: Dict[str, str]) -> Path:
    """Return the chembl_<species> output directory for one target row."""
    return dataset_root / target["target_dir"] / f"chembl_{target['species']}"


def flatten_activity(record: Dict[str, object]) -> Dict[str, object]:
    """Flatten one activity record, expanding ligand efficiency into columns."""
    flat = {field: record.get(field) for field in ACTIVITY_FIELDS if field != "ligand_efficiency"}
    efficiency = record.get("ligand_efficiency") or {}
    if isinstance(efficiency, dict):
        for key in ("bei", "le", "lle", "sei"):
            flat[f"le_{key}"] = efficiency.get(key)
    else:
        for key in ("bei", "le", "lle", "sei"):
            flat[f"le_{key}"] = None
    return flat


def count_existing_rows(partial_path: Path) -> int:
    """Return the number of data rows already written to a partial CSV."""
    if not partial_path.is_file():
        return 0
    with partial_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        try:
            next(reader)
        except StopIteration:
            return 0
        return sum(1 for _ in reader)


def fetch_target(
    target: Dict[str, str],
    dataset_root: Path,
    force: bool,
    session: requests.Session,
    chembl_db_version: str,
) -> Optional[Path]:
    """Download one target's activities to CSV. Return the path, or None if skipped."""
    out_dir = output_dir(dataset_root, target)
    out_dir.mkdir(parents=True, exist_ok=True)
    raw_path = out_dir / "activities_raw.csv"
    partial_path = out_dir / "activities_raw.partial.csv"
    meta_path = out_dir / "activities_raw.meta.csv"

    target_id = target["target_chembl_id"]
    label = f"{target['target_dir']}/{target['species']} ({target_id})"

    if raw_path.exists() and not force:
        print(f"[skip] {label}: {raw_path} already exists")
        return raw_path

    fieldnames = [f for f in ACTIVITY_FIELDS if f != "ligand_efficiency"] + [
        "le_bei",
        "le_le",
        "le_lle",
        "le_sei",
    ]
    if force and partial_path.exists():
        partial_path.unlink()
    rows_written = count_existing_rows(partial_path)
    mode = "a" if rows_written else "w"
    if rows_written:
        print(f"[resume] {label}: {rows_written} rows already in {partial_path.name}")

    print(f"[fetch] {label}")
    with partial_path.open(mode, newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        if mode == "w":
            writer.writeheader()
        for record in iter_activity_pages(
            target_id,
            start_offset=rows_written,
            session=session,
        ):
            writer.writerow(flatten_activity(record))
            rows_written += 1

    partial_path.replace(raw_path)
    with meta_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["target_chembl_id", "n_rows", "chembl_db_version", "fetched_at_utc"],
        )
        writer.writeheader()
        writer.writerow(
            {
                "target_chembl_id": target_id,
                "n_rows": rows_written,
                "chembl_db_version": chembl_db_version,
                "fetched_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
        )

    print(f"[ok] {label}: {rows_written} rows -> {raw_path}")
    return raw_path


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Download all ChEMBL activity records for every target listed in "
            "chembl_targets.csv and store them under each target folder."
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
        "--targets-file",
        default=DEFAULT_TARGETS_FILE,
        help=f"target metadata CSV inside dataset_root (default: {DEFAULT_TARGETS_FILE})",
    )
    parser.add_argument(
        "--only",
        default=None,
        help="comma-separated target_dir values to fetch (default: all)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="re-download even if activities_raw.csv already exists",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the ChEMBL activity download."""
    parser = build_argument_parser()
    arguments = parser.parse_args(argv)
    dataset_root = arguments.dataset_root.expanduser().resolve()
    targets_path = dataset_root / arguments.targets_file

    if not targets_path.is_file():
        parser.error(f"target metadata file does not exist: {targets_path}")

    targets = load_targets(targets_path)
    if arguments.only:
        wanted = {name.strip() for name in arguments.only.split(",") if name.strip()}
        targets = [t for t in targets if t["target_dir"] in wanted]
        if not targets:
            parser.error(f"no targets matched --only {arguments.only}")

    failures: List[str] = []
    with requests.Session() as session:
        session.headers.update({"User-Agent": "chembl-reference-fetch/1.0"})
        status = get_json("https://www.ebi.ac.uk/chembl/api/data/status.json", session=session)
        chembl_db_version = status.get("chembl_db_version", "unknown")
        print(f"ChEMBL database version: {chembl_db_version}")
        for target in targets:
            try:
                fetch_target(
                    target,
                    dataset_root,
                    arguments.force,
                    session,
                    chembl_db_version,
                )
            except Exception as exc:
                label = f"{target['target_dir']}/{target['species']}"
                failures.append(label)
                print(f"[FAILED] {label}: {exc}", file=sys.stderr)

    print(f"Fetched {len(targets) - len(failures)}/{len(targets)} target sets.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
