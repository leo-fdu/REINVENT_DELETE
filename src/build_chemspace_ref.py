#!/usr/bin/env python3
"""Build the ChEMBL active-molecule chemical-space reference from raw activity tables."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Descriptors, QED
except ImportError as exc:
    raise SystemExit("This script requires RDKit. Install it before running.") from exc

DEFAULT_DATASET_ROOT = Path(__file__).resolve().parents[1] / "real-world_dataset"
DEFAULT_TARGETS_FILE = "chembl_targets.csv"
DEFAULT_PROPS_FILE = "chembl_compound_properties.csv"
DEFAULT_SUMMARY_NAME = "chembl_chemspace_summary.csv"
DEFAULT_PCHEMBL_THRESHOLD = 6.0

ACTIVE_STANDARD_TYPES = ("IC50", "Ki", "Kd", "EC50", "AC50", "XC50")
ACTIVE_ASSAY_TYPES = ("B", "F")

# Green lane: CYP3A4 is a drug-metabolizing enzyme, so most of its ChEMBL records
# are tagged assay_type=A (ADME / DDI liability). Those measurements report the same
# pocket interaction as B/F records (99.6% of A rows and 99.9% of B/F rows carry
# "CYP3A4/P450/microsome" in assay_description), so A is only a curation label here.
# For the other 15 targets, A-type rows are misclassified pharmacology and stay excluded.
ADME_ASSAY_TARGETS = {"cp3a4"}
PATHOGEN_PATTERN = re.compile(
    r"neisseria|mycobacterium|pneumocystis|plasmodium|bacterial|fungal|candida|cryptosporidium",
    re.IGNORECASE,
)

DESCRIPTOR_COLUMNS = [
    "mw",
    "clogp",
    "tpsa",
    "hbd",
    "hba",
    "rotb",
    "n_rings",
    "n_arom_rings",
    "fsp3",
    "qed",
]


def load_targets(targets_path: Path) -> List[Dict[str, str]]:
    """Read the target metadata table."""
    with targets_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def filter_active_rows(frame: pd.DataFrame, target_dir: str, threshold: float) -> pd.DataFrame:
    """Apply the active-record filter to one raw activity table."""
    filtered = frame.copy()
    filtered["pchembl_value"] = pd.to_numeric(filtered["pchembl_value"], errors="coerce")
    filtered["standard_flag"] = pd.to_numeric(filtered["standard_flag"], errors="coerce")
    filtered["potential_duplicate"] = pd.to_numeric(filtered["potential_duplicate"], errors="coerce")

    assay_types = ACTIVE_ASSAY_TYPES
    if target_dir in ADME_ASSAY_TARGETS:
        assay_types = ACTIVE_ASSAY_TYPES + ("A",)

    mask = (
        filtered["standard_type"].isin(ACTIVE_STANDARD_TYPES)
        & (filtered["standard_relation"] == "=")
        & (filtered["pchembl_value"] >= threshold)
        & (filtered["standard_flag"] == 1)
        & (filtered["potential_duplicate"] == 0)
        & (filtered["assay_type"].isin(assay_types))
        & filtered["data_validity_comment"].isna()
    )
    if target_dir == "dyr":
        descriptions = filtered["assay_description"].fillna("")
        mask &= ~descriptions.str.contains(PATHOGEN_PATTERN)

    return filtered.loc[mask].copy()


def classify_evidence(bao_labels: Sequence[str]) -> str:
    """Summarise assay-format evidence for one compound."""
    unique = {label for label in bao_labels if isinstance(label, str) and label}
    if unique == {"single protein format"}:
        return "single_protein"
    if "single protein format" in unique:
        return "mixed"
    if unique & {"cell-based format", "assay format"}:
        return "cell_based" if "cell-based format" in unique else "other"
    return "other"


def aggregate_actives(filtered: pd.DataFrame, smiles_lookup: Dict[str, str]) -> pd.DataFrame:
    """Collapse activity rows to one row per parent molecule."""
    if filtered.empty:
        return pd.DataFrame()

    rows = []
    for parent_id, group in filtered.groupby("parent_molecule_chembl_id", sort=True):
        pchembl = group["pchembl_value"].astype(float)
        rows.append(
            {
                "molecule_chembl_id": group["molecule_chembl_id"].iloc[0],
                "parent_molecule_chembl_id": parent_id,
                "canonical_smiles": smiles_lookup.get(parent_id)
                or group["canonical_smiles"].dropna().iloc[0]
                if group["canonical_smiles"].notna().any()
                else smiles_lookup.get(parent_id),
                "pref_name": group["molecule_pref_name"].dropna().iloc[0]
                if group["molecule_pref_name"].notna().any()
                else None,
                "pchembl_median": float(pchembl.median()),
                "pchembl_n": int(len(pchembl)),
                "pchembl_min": float(pchembl.min()),
                "pchembl_max": float(pchembl.max()),
                "standard_types": ";".join(sorted(group["standard_type"].dropna().unique())),
                "assay_types": ";".join(sorted(group["assay_type"].dropna().unique())),
                "bao_labels": ";".join(sorted(set(group["bao_label"].dropna()))),
                "evidence_level": classify_evidence(group["bao_label"].dropna().tolist()),
                "assay_confidence_scores": ";".join(
                    sorted({str(int(v)) for v in group["assay_confidence_score"].dropna().unique()})
                )
                if "assay_confidence_score" in group
                else None,
                "assay_relationship_types": ";".join(
                    sorted(set(group["assay_relationship_type"].dropna()))
                )
                if "assay_relationship_type" in group
                else None,
                "document_years": ";".join(
                    sorted({str(int(year)) for year in group["document_year"].dropna().unique()})
                ),
            }
        )
    return pd.DataFrame(rows)


def compute_descriptors(smiles: str) -> Optional[Dict[str, float]]:
    """Compute RDKit physicochemical descriptors for one SMILES string."""
    molecule = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if molecule is None:
        return None
    return {
        "mw": Descriptors.MolWt(molecule),
        "clogp": Descriptors.MolLogP(molecule),
        "tpsa": Descriptors.TPSA(molecule),
        "hbd": Descriptors.NumHDonors(molecule),
        "hba": Descriptors.NumHAcceptors(molecule),
        "rotb": Descriptors.NumRotatableBonds(molecule),
        "n_rings": Descriptors.RingCount(molecule),
        "n_arom_rings": Descriptors.NumAromaticRings(molecule),
        "fsp3": Descriptors.FractionCSP3(molecule),
        "qed": QED.qed(molecule),
    }


def write_sdf(chemspace: pd.DataFrame, sdf_path: Path) -> int:
    """Write 2D-embedded SDF records with pchembl annotations."""
    writer = Chem.SDWriter(str(sdf_path))
    written = 0
    for row in chemspace.itertuples(index=False):
        molecule = Chem.MolFromSmiles(row.canonical_smiles)
        if molecule is None:
            continue
        AllChem.Compute2DCoords(molecule)
        molecule.SetProp("_Name", str(row.parent_molecule_chembl_id))
        molecule.SetProp("molecule_chembl_id", str(row.molecule_chembl_id))
        molecule.SetProp("pchembl_median", f"{row.pchembl_median:.3f}")
        molecule.SetProp("pchembl_n", str(row.pchembl_n))
        molecule.SetProp("evidence_level", str(row.evidence_level))
        molecule.SetProp("qed", f"{row.qed:.3f}")
        molecule.SetProp("mw", f"{row.mw:.2f}")
        molecule.SetProp("clogp", f"{row.clogp:.2f}")
        writer.write(molecule)
        written += 1
    writer.close()
    return written


def summarise(chemspace: pd.DataFrame, target: Dict[str, str]) -> Dict[str, object]:
    """Compute quantiles for one target's chemical space."""
    summary: Dict[str, object] = {
        "target_dir": target["target_dir"],
        "species": target["species"],
        "target_chembl_id": target["target_chembl_id"],
        "n_actives": len(chemspace),
        "n_single_protein": int((chemspace["evidence_level"] == "single_protein").sum()),
    }
    if chemspace.empty:
        for column in DESCRIPTOR_COLUMNS + ["pchembl_median"]:
            for stat in ("min", "p5", "p25", "median", "p75", "p95", "max"):
                summary[f"{column}_{stat}"] = None
        return summary

    for column in DESCRIPTOR_COLUMNS + ["pchembl_median"]:
        quantiles = chemspace[column].quantile([0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0])
        summary[f"{column}_min"] = float(quantiles.loc[0.0])
        summary[f"{column}_p5"] = float(quantiles.loc[0.05])
        summary[f"{column}_p25"] = float(quantiles.loc[0.25])
        summary[f"{column}_median"] = float(quantiles.loc[0.5])
        summary[f"{column}_p75"] = float(quantiles.loc[0.75])
        summary[f"{column}_p95"] = float(quantiles.loc[0.95])
        summary[f"{column}_max"] = float(quantiles.loc[1.0])
    return summary


def process_target(
    target: Dict[str, str],
    dataset_root: Path,
    threshold: float,
    props: Optional[pd.DataFrame],
) -> Tuple[Optional[Path], Dict[str, object]]:
    """Filter, aggregate and annotate one target's actives. Return chemspace path + summary."""
    out_dir = dataset_root / target["target_dir"] / f"chembl_{target['species']}"
    raw_path = out_dir / "activities_raw.csv"
    label = f"{target['target_dir']}/{target['species']}"

    if not raw_path.is_file():
        print(f"[skip] {label}: missing {raw_path}", file=sys.stderr)
        return None, summarise(pd.DataFrame(), target)

    raw = pd.read_csv(raw_path, dtype=str, low_memory=False)
    filtered = filter_active_rows(raw, target["target_dir"], threshold)
    print(f"[filter] {label}: {len(raw)} raw -> {len(filtered)} active rows")

    smiles_lookup: Dict[str, str] = {}
    max_phase_lookup: Dict[str, object] = {}
    if props is not None and not props.empty:
        for row in props.itertuples(index=False):
            if isinstance(row.canonical_smiles, str) and row.canonical_smiles:
                smiles_lookup[row.molecule_chembl_id] = row.canonical_smiles
            max_phase_lookup[row.molecule_chembl_id] = row.max_phase

    actives = aggregate_actives(filtered, smiles_lookup)
    if actives.empty:
        print(f"[warn] {label}: no actives passed the filter", file=sys.stderr)
        return None, summarise(pd.DataFrame(), target)

    actives["max_phase"] = actives["parent_molecule_chembl_id"].map(max_phase_lookup)

    descriptors = actives["canonical_smiles"].map(compute_descriptors)
    invalid = descriptors.isna().sum()
    if invalid:
        print(f"[warn] {label}: {invalid} SMILES failed RDKit parsing and were dropped")
        keep = descriptors.notna()
        actives = actives.loc[keep].copy()
        descriptors = descriptors.loc[keep]
    descriptor_frame = pd.DataFrame(list(descriptors), index=actives.index)
    chemspace = pd.concat([actives.reset_index(drop=True), descriptor_frame.reset_index(drop=True)], axis=1)
    chemspace = chemspace.sort_values("pchembl_median", ascending=False).reset_index(drop=True)

    actives_csv = out_dir / "actives.csv"
    chemspace_csv = out_dir / "chemspace.csv"
    sdf_path = out_dir / "actives.sdf"

    actives.to_csv(actives_csv, index=False)
    chemspace.to_csv(chemspace_csv, index=False)
    n_sdf = write_sdf(chemspace, sdf_path)
    print(f"[ok] {label}: {len(chemspace)} actives -> {chemspace_csv.name}, {sdf_path.name} ({n_sdf} records)")

    return chemspace_csv, summarise(chemspace, target)


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Filter ChEMBL activities into an active-molecule set (pchembl >= threshold), "
            "compute physicochemical descriptors, and write per-target chemical-space "
            "references plus a cross-target summary."
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
        "--props-file",
        default=DEFAULT_PROPS_FILE,
        help=f"molecule properties CSV inside dataset_root (default: {DEFAULT_PROPS_FILE})",
    )
    parser.add_argument(
        "--pchembl-threshold",
        type=float,
        default=DEFAULT_PCHEMBL_THRESHOLD,
        help=f"minimum pchembl_value for actives (default: {DEFAULT_PCHEMBL_THRESHOLD})",
    )
    parser.add_argument(
        "--summary-name",
        default=DEFAULT_SUMMARY_NAME,
        help=f"summary CSV inside dataset_root (default: {DEFAULT_SUMMARY_NAME})",
    )
    parser.add_argument(
        "--only",
        default=None,
        help="comma-separated target_dir values to process (default: all)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the chemical-space reference build."""
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

    props_path = dataset_root / arguments.props_file
    props = pd.read_csv(props_path, dtype=str, low_memory=False) if props_path.is_file() else None
    if props is None:
        print(f"[warn] {props_path} not found; using activity-table SMILES only", file=sys.stderr)

    summaries = []
    failures = 0
    for target in targets:
        try:
            _, summary = process_target(
                target,
                dataset_root,
                arguments.pchembl_threshold,
                props,
            )
            summaries.append(summary)
        except Exception as exc:
            failures += 1
            label = f"{target['target_dir']}/{target['species']}"
            print(f"[FAILED] {label}: {exc}", file=sys.stderr)

    summary_path = dataset_root / arguments.summary_name
    pd.DataFrame(summaries).to_csv(summary_path, index=False)
    print(f"Wrote {len(summaries)} summary rows -> {summary_path}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
