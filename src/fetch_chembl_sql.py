#!/usr/bin/env python3
"""Extract per-target ChEMBL activity tables from a local ChEMBL SQLite database."""

from __future__ import annotations

import argparse
import csv
import sqlite3
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Set, TextIO

DEFAULT_DATASET_ROOT = Path(__file__).resolve().parents[1] / "real-world_dataset"
DEFAULT_TARGETS_FILE = "chembl_targets.csv"
DEFAULT_PROPS_NAME = "chembl_compound_properties.csv"
DEFAULT_DB_PATH = Path.home() / "Desktop" / "data" / "chembl36" / "chembl_36.db"

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
    "assay_confidence_score",
    "assay_relationship_type",
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
    "le_bei",
    "le_le",
    "le_lle",
    "le_sei",
]

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
    "psa",
    "hba",
    "hbd",
    "rtb",
    "aromatic_rings",
    "heavy_atoms",
    "qed_weighted",
    "np_likeness_score",
    "num_ro5_violations",
    "ro3_pass",
    "full_molformula",
]

ACTIVITY_QUERY = """
SELECT
    act.activity_id                          AS activity_id,
    md.chembl_id                             AS molecule_chembl_id,
    pmd.chembl_id                            AS parent_molecule_chembl_id,
    cs.canonical_smiles                      AS canonical_smiles,
    md.pref_name                             AS molecule_pref_name,
    act.standard_type                        AS standard_type,
    act.standard_relation                    AS standard_relation,
    act.standard_value                       AS standard_value,
    act.standard_units                       AS standard_units,
    act.pchembl_value                        AS pchembl_value,
    a.chembl_id                              AS assay_chembl_id,
    a.assay_type                             AS assay_type,
    a.description                            AS assay_description,
    a.confidence_score                       AS assay_confidence_score,
    a.relationship_type                      AS assay_relationship_type,
    a.bao_format                             AS bao_format,
    bao.label                                AS bao_label,
    d.chembl_id                              AS document_chembl_id,
    d.journal                                AS document_journal,
    d.year                                   AS document_year,
    act.data_validity_comment                AS data_validity_comment,
    act.potential_duplicate                  AS potential_duplicate,
    act.standard_flag                        AS standard_flag,
    td.chembl_id                             AS target_chembl_id,
    td.pref_name                             AS target_pref_name,
    td.organism                              AS target_organism,
    le.bei                                   AS le_bei,
    le.le                                    AS le_le,
    le.lle                                   AS le_lle,
    le.sei                                   AS le_sei
FROM activities act
JOIN assays a ON act.assay_id = a.assay_id
JOIN target_dictionary td ON a.tid = td.tid
LEFT JOIN molecule_dictionary md ON act.molregno = md.molregno
LEFT JOIN molecule_hierarchy mh ON act.molregno = mh.molregno
LEFT JOIN molecule_dictionary pmd
       ON pmd.molregno = COALESCE(mh.parent_molregno, act.molregno)
LEFT JOIN compound_structures cs ON act.molregno = cs.molregno
LEFT JOIN bioassay_ontology bao ON a.bao_format = bao.bao_id
LEFT JOIN docs d ON act.doc_id = d.doc_id
LEFT JOIN ligand_eff le ON act.activity_id = le.activity_id
WHERE a.tid IN ({tid_list})
ORDER BY td.chembl_id, act.activity_id
"""

PROPERTY_QUERY = """
SELECT DISTINCT
    pmd.chembl_id                AS molecule_chembl_id,
    pmd.pref_name                AS pref_name,
    pmd.max_phase                AS max_phase,
    pmd.molecule_type            AS molecule_type,
    pmd.first_approval           AS first_approval,
    pmd.usan_stem                AS usan_stem,
    pcs.standard_inchi_key       AS standard_inchi_key,
    pcs.canonical_smiles         AS canonical_smiles,
    pcp.mw_freebase              AS mw_freebase,
    pcp.full_mwt                 AS full_mwt,
    pcp.alogp                    AS alogp,
    pcp.psa                      AS psa,
    pcp.hba                      AS hba,
    pcp.hbd                      AS hbd,
    pcp.rtb                      AS rtb,
    pcp.aromatic_rings           AS aromatic_rings,
    pcp.heavy_atoms              AS heavy_atoms,
    pcp.qed_weighted             AS qed_weighted,
    pcp.np_likeness_score        AS np_likeness_score,
    pcp.num_ro5_violations       AS num_ro5_violations,
    pcp.ro3_pass                 AS ro3_pass,
    pcp.full_molformula          AS full_molformula
FROM activities act
JOIN assays a ON act.assay_id = a.assay_id
LEFT JOIN molecule_hierarchy mh ON act.molregno = mh.molregno
LEFT JOIN molecule_dictionary pmd
       ON pmd.molregno = COALESCE(mh.parent_molregno, act.molregno)
LEFT JOIN compound_structures pcs ON pmd.molregno = pcs.molregno
LEFT JOIN compound_properties pcp ON pmd.molregno = pcp.molregno
WHERE a.tid IN ({tid_list})
"""


def load_targets(targets_path: Path) -> List[Dict[str, str]]:
    """Read the target metadata table."""
    with targets_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def output_dir(dataset_root: Path, target: Dict[str, str]) -> Path:
    """Return the chembl_<species> output directory for one target row."""
    return dataset_root / target["target_dir"] / f"chembl_{target['species']}"


def read_chembl_version(connection: sqlite3.Connection) -> str:
    """Return the ChEMBL release name recorded in the database."""
    row = connection.execute(
        "SELECT name FROM version WHERE name GLOB 'ChEMBL_[0-9]*' ORDER BY name DESC LIMIT 1"
    ).fetchone()
    return row[0] if row else "unknown"


def write_activity_tables(
    connection: sqlite3.Connection,
    targets: Sequence[Dict[str, str]],
    dataset_root: Path,
    chembl_db_version: str,
) -> int:
    """Stream all activities for the requested targets into per-target CSVs."""
    tid_list = ",".join(str(int(target["tid"])) for target in targets)
    query = ACTIVITY_QUERY.format(tid_list=tid_list)

    writers: Dict[str, csv.DictWriter] = {}
    handles: Dict[str, TextIO] = {}
    counts: Dict[str, int] = {}
    for target in targets:
        key = f"{target['target_dir']}/{target['species']}"
        out_dir = output_dir(dataset_root, target)
        out_dir.mkdir(parents=True, exist_ok=True)
        handle = (out_dir / "activities_raw.csv").open("w", newline="", encoding="utf-8")
        writer = csv.DictWriter(handle, fieldnames=ACTIVITY_FIELDS)
        writer.writeheader()
        handles[key] = handle
        writers[key] = writer
        counts[key] = 0

    total = 0
    try:
        cursor = connection.execute(query)
        for row in cursor:
            record = dict(row)
            key = f"{target_dir_from_chembl(targets, record['target_chembl_id'])}"
            writers[key].writerow(record)
            counts[key] += 1
            total += 1
            if total % 50000 == 0:
                print(f"  {total} activity rows written...")
    finally:
        for handle in handles.values():
            handle.close()

    for target in targets:
        key = f"{target['target_dir']}/{target['species']}"
        out_dir = output_dir(dataset_root, target)
        with (out_dir / "activities_raw.meta.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["target_chembl_id", "n_rows", "chembl_db_version", "fetched_at_utc"],
            )
            writer.writeheader()
            writer.writerow(
                {
                    "target_chembl_id": target["target_chembl_id"],
                    "n_rows": counts[key],
                    "chembl_db_version": chembl_db_version,
                    "fetched_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                }
            )
        print(f"[ok] {key}: {counts[key]} rows")
    return total


def target_dir_from_chembl(targets: Sequence[Dict[str, str]], chembl_id: str) -> str:
    """Return '<target_dir>/<species>' for a ChEMBL target ID."""
    for target in targets:
        if target["target_chembl_id"] == chembl_id:
            return f"{target['target_dir']}/{target['species']}"
    raise KeyError(f"unexpected target chembl id: {chembl_id}")


def write_property_table(
    connection: sqlite3.Connection,
    targets: Sequence[Dict[str, str]],
    output_path: Path,
) -> int:
    """Write molecule properties for every parent compound seen in the activities."""
    tid_list = ",".join(str(int(target["tid"])) for target in targets)
    query = PROPERTY_QUERY.format(tid_list=tid_list)
    seen: Set[str] = set()
    written = 0
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PROPERTY_FIELDS)
        writer.writeheader()
        for row in connection.execute(query):
            record = dict(row)
            molecule_id = record.get("molecule_chembl_id")
            if molecule_id is None or molecule_id in seen:
                continue
            seen.add(molecule_id)
            writer.writerow(record)
            written += 1
    return written


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Extract per-target activity tables and molecule properties from a "
            "local ChEMBL SQLite database into real-world_dataset/<target>/chembl_<species>/."
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
        "--db",
        type=Path,
        default=DEFAULT_DB_PATH,
        help=f"path to the ChEMBL SQLite database (default: {DEFAULT_DB_PATH})",
    )
    parser.add_argument(
        "--targets-file",
        default=DEFAULT_TARGETS_FILE,
        help=f"target metadata CSV inside dataset_root (default: {DEFAULT_TARGETS_FILE})",
    )
    parser.add_argument(
        "--props-name",
        default=DEFAULT_PROPS_NAME,
        help=f"molecule properties CSV inside dataset_root (default: {DEFAULT_PROPS_NAME})",
    )
    parser.add_argument(
        "--only",
        default=None,
        help="comma-separated target_dir values to fetch (default: all)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the local ChEMBL SQLite extraction."""
    parser = build_argument_parser()
    arguments = parser.parse_args(argv)
    dataset_root = arguments.dataset_root.expanduser().resolve()
    targets_path = dataset_root / arguments.targets_file
    db_path = arguments.db.expanduser().resolve()

    if not targets_path.is_file():
        parser.error(f"target metadata file does not exist: {targets_path}")
    if not db_path.is_file():
        parser.error(f"ChEMBL SQLite database does not exist: {db_path}")

    targets = load_targets(targets_path)
    if arguments.only:
        wanted = {name.strip() for name in arguments.only.split(",") if name.strip()}
        targets = [t for t in targets if t["target_dir"] in wanted]
        if not targets:
            parser.error(f"no targets matched --only {arguments.only}")

    connection = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    connection.row_factory = sqlite3.Row
    try:
        chembl_db_version = read_chembl_version(connection)
        print(f"ChEMBL database version: {chembl_db_version}")
        print(f"Extracting activities for {len(targets)} target sets...")
        total = write_activity_tables(connection, targets, dataset_root, chembl_db_version)
        print(f"Total activity rows: {total}")

        props_path = dataset_root / arguments.props_name
        print("Extracting molecule properties...")
        n_props = write_property_table(connection, targets, props_path)
        print(f"[ok] {n_props} parent molecules -> {props_path}")
    finally:
        connection.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
