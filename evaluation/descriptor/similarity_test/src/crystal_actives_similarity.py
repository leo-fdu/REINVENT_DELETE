#!/usr/bin/env python3
"""Tanimoto similarity between co-crystal ligands and ChEMBL actives.

For every target in the real-world dataset this script

1. converts the co-crystal ligand (``crystal.mol2``) to a sanitized RDKit
   molecule,
2. computes its Morgan fingerprint (ECFP4-like: radius 2, 2048 bits),
3. computes the Morgan fingerprints of all ChEMBL active molecules listed in
   ``chembl_human/actives.csv`` for the same target, and
4. writes a CSV table with the columns ``chembl_id``, ``smiles`` and
   ``tanimoto_similarity`` to ``results/<target>/``.

MOL2 cleanup
------------
The ``crystal.mol2`` files carry Gasteiger partial charges, no explicit
hydrogens and a few inconsistent aromatic annotations.  RDKit's MOL2 parser
therefore estimates bogus formal charges (e.g. ``[N-2]``) and marks every
atom ``NoImplicit``.  The loader below

* zeroes all formal charges (the heuristics used for H-less MOL2 files are
  unreliable) and re-enables implicit hydrogens,
* repairs non-ring "aromatic" annotations (a carboxyl group in the ppara
  MOL2 is marked aromatic) by converting them to single/double bonds,
* kekulizes the ring systems; if that fails because a pyrrole-type aromatic
  nitrogen lacks its hydrogen (five-membered aromatic rings in adrb1, akt1,
  gria2 and parp1), explicit hydrogens are added to candidate nitrogens until
  kekulization succeeds.  Every attempt runs on a fresh copy of the molecule
  because a failed ``Kekulize`` leaves the molecule partially modified.

The conversion is validated by a SMILES round-trip: the canonical SMILES of
the cleaned molecule is parsed again and must yield an identical Morgan
fingerprint (Tanimoto 1.0).
"""

from __future__ import annotations

import argparse
import csv
import sys
from itertools import combinations
from pathlib import Path
from typing import List, Optional, Sequence, Tuple

try:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import AllChem
except ImportError as exc:
    raise SystemExit("This script requires RDKit. Install it before running.") from exc

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[3]
DEFAULT_DATASET_ROOT = REPO_ROOT / "real-world_dataset"
DEFAULT_OUTPUT_ROOT = SCRIPT_DIR.parent / "results"

CRYSTAL_FILENAME = "crystal.mol2"
ACTIVES_FILENAME = Path("chembl_human") / "actives.csv"
ACTIVES_ID_COLUMN = "molecule_chembl_id"
ACTIVES_SMILES_COLUMN = "canonical_smiles"
RESULT_COLUMNS = ("chembl_id", "smiles", "tanimoto_similarity")


def find_targets(dataset_root: Path) -> List[str]:
    """Return target names that provide a crystal MOL2 and a ChEMBL actives table."""
    targets = []
    for child in sorted(dataset_root.iterdir()):
        if not child.is_dir():
            continue
        if (child / CRYSTAL_FILENAME).is_file() and (child / ACTIVES_FILENAME).is_file():
            targets.append(child.name)
    return targets


def fix_nonring_aromaticity(mol: Chem.Mol) -> None:
    """Repair aromatic annotations on atoms/bonds that are not part of a ring.

    Some MOL2 files mark exocyclic groups (e.g. a carboxyl carbon/oxygens) as
    aromatic.  Such bonds are converted to double bonds for degree-1 O/S and
    to single bonds otherwise; a carbon carrying two exocyclic C=O bonds is
    turned into a carboxylic acid (one bond reduced to single, chosen
    deterministically by atom index).
    """
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.AROMATIC or bond.IsInRing():
            continue
        begin, end = bond.GetBeginAtom(), bond.GetEndAtom()
        exocyclic = None
        if begin.IsInRing() and not end.IsInRing():
            exocyclic = end
        elif end.IsInRing() and not begin.IsInRing():
            exocyclic = begin
        reference = exocyclic if exocyclic is not None else end
        if reference.GetDegree() == 1 and reference.GetSymbol() in ("O", "S"):
            bond.SetBondType(Chem.BondType.DOUBLE)
        else:
            bond.SetBondType(Chem.BondType.SINGLE)
        bond.SetIsAromatic(False)
    for atom in mol.GetAtoms():
        if atom.GetIsAromatic() and not atom.IsInRing():
            atom.SetIsAromatic(False)
    mol.UpdatePropertyCache(strict=False)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() != "C":
            continue
        carbonyl_bonds = sorted(
            (
                bond
                for bond in atom.GetBonds()
                if bond.GetBondType() == Chem.BondType.DOUBLE
                and bond.GetOtherAtom(atom).GetSymbol() == "O"
                and bond.GetOtherAtom(atom).GetDegree() == 1
            ),
            key=lambda bond: bond.GetOtherAtom(atom).GetIdx(),
        )
        for bond in carbonyl_bonds[1:]:
            bond.SetBondType(Chem.BondType.SINGLE)
    mol.UpdatePropertyCache(strict=False)


def try_kekulize(mol: Chem.Mol) -> Optional[Chem.Mol]:
    """Return a kekulized copy of ``mol`` or ``None`` if kekulization fails.

    A fresh copy is used for every attempt because a failed ``Kekulize`` call
    leaves the molecule partially modified.
    """
    candidate = Chem.Mol(mol)
    try:
        Chem.Kekulize(candidate, clearAromaticFlags=False)
    except Exception:
        return None
    return candidate


def kekulize_with_pyrrole_fix(mol: Chem.Mol) -> Chem.Mol:
    """Kekulize ``mol``, adding pyrrole-type hydrogens to aromatic N if needed.

    Aromatic five-membered-ring nitrogens are sometimes annotated ``N.ar``
    without their hydrogen (pyrrole-type N), which makes kekulization
    impossible.  Candidate nitrogens (aromatic, degree 2, neutral, no
    hydrogens) are given one explicit hydrogen - alone or in pairs - until
    kekulization succeeds.
    """
    kekulized = try_kekulize(mol)
    if kekulized is not None:
        return kekulized
    candidates = sorted(
        (
            atom
            for atom in mol.GetAtoms()
            if atom.GetSymbol() == "N"
            and atom.GetIsAromatic()
            and atom.GetDegree() == 2
            and atom.GetFormalCharge() == 0
            and atom.GetTotalNumHs() == 0
        ),
        key=lambda atom: atom.GetIdx(),
    )
    for size in range(1, min(2, len(candidates)) + 1):
        for combo in combinations(candidates, size):
            for atom in combo:
                atom.SetNumExplicitHs(1)
            mol.UpdatePropertyCache(strict=False)
            kekulized = try_kekulize(mol)
            if kekulized is not None:
                return kekulized
            for atom in combo:
                atom.SetNumExplicitHs(0)
            mol.UpdatePropertyCache(strict=False)
    raise ValueError("could not kekulize the molecule")


def load_crystal_ligand(mol2_path: Path) -> Chem.Mol:
    """Read a co-crystal ligand MOL2 file and return a sanitized molecule."""
    mol = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=True)
    if mol is None:
        raise ValueError("RDKit could not parse the MOL2 file")
    for atom in mol.GetAtoms():
        atom.SetFormalCharge(0)
        atom.SetNoImplicit(False)
    mol.UpdatePropertyCache(strict=False)
    Chem.FastFindRings(mol)
    fix_nonring_aromaticity(mol)
    mol = kekulize_with_pyrrole_fix(mol)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistryFrom3D(mol, replaceExistingTags=True)
    return mol


def morgan_fingerprint(mol: Chem.Mol, radius: int, n_bits: int):
    """Compute the Morgan fingerprint (bit vector) of a molecule."""
    return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=n_bits)


def load_actives(actives_path: Path) -> Tuple[List[Tuple[str, str]], int]:
    """Return (chembl_id, canonical SMILES) pairs and the number of bad rows."""
    actives: List[Tuple[str, str]] = []
    skipped = 0
    with actives_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            chembl_id = (row.get(ACTIVES_ID_COLUMN) or "").strip()
            smiles = (row.get(ACTIVES_SMILES_COLUMN) or "").strip()
            mol = Chem.MolFromSmiles(smiles) if smiles else None
            if not chembl_id or mol is None:
                skipped += 1
                continue
            actives.append((chembl_id, Chem.MolToSmiles(mol)))
    return actives, skipped


def compute_similarity_table(
    crystal_mol: Chem.Mol,
    actives: List[Tuple[str, str]],
    radius: int,
    n_bits: int,
) -> List[Tuple[str, str, float]]:
    """Compute the Tanimoto similarity of every active against the crystal ligand."""
    reference_fp = morgan_fingerprint(crystal_mol, radius, n_bits)
    active_mols = [Chem.MolFromSmiles(smiles) for _, smiles in actives]
    active_fps = [morgan_fingerprint(mol, radius, n_bits) for mol in active_mols]
    scores = DataStructs.BulkTanimotoSimilarity(reference_fp, active_fps)
    return [
        (chembl_id, smiles, round(score, 6))
        for (chembl_id, smiles), score in zip(actives, scores)
    ]


def write_similarity_table(
    rows: List[Tuple[str, str, float]],
    output_dir: Path,
    target: str,
) -> Path:
    """Write the similarity CSV table for one target."""
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"{target}_tanimoto_similarity.csv"
    with csv_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(RESULT_COLUMNS)
        writer.writerows(rows)
    return csv_path


def process_target(
    target: str,
    dataset_root: Path,
    output_root: Path,
    radius: int,
    n_bits: int,
) -> Tuple[int, int]:
    """Run the similarity calculation for one target and return row counts."""
    target_dir = dataset_root / target
    crystal_mol = load_crystal_ligand(target_dir / CRYSTAL_FILENAME)
    crystal_smiles = Chem.MolToSmiles(crystal_mol)

    roundtrip = Chem.MolFromSmiles(crystal_smiles)
    if roundtrip is None:
        raise ValueError("canonical SMILES of the crystal ligand does not re-parse")
    same = DataStructs.TanimotoSimilarity(
        morgan_fingerprint(crystal_mol, radius, n_bits),
        morgan_fingerprint(roundtrip, radius, n_bits),
    )
    if same != 1.0:
        raise ValueError(f"crystal ligand SMILES round-trip changed the fingerprint ({same:.4f})")

    actives, skipped = load_actives(target_dir / ACTIVES_FILENAME)
    rows = compute_similarity_table(crystal_mol, actives, radius, n_bits)
    csv_path = write_similarity_table(rows, output_root / target, target)
    print(f"[OK] {target}: crystal SMILES {crystal_smiles}")
    print(f"     {len(rows)} actives -> {csv_path}")
    return len(rows), skipped


def build_argument_parser() -> argparse.ArgumentParser:
    """Build the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Convert each target's co-crystal ligand (crystal.mol2) to SMILES "
            "and compute the Tanimoto similarity of its Morgan fingerprint "
            "against all ChEMBL actives of the same target."
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
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="directory for per-target result folders (default: similarity_test/results)",
    )
    parser.add_argument(
        "--radius",
        type=int,
        default=2,
        help="Morgan fingerprint radius (default: 2)",
    )
    parser.add_argument(
        "--n-bits",
        type=int,
        default=2048,
        help="Morgan fingerprint size in bits (default: 2048)",
    )
    parser.add_argument(
        "--targets",
        nargs="+",
        default=None,
        help="restrict the run to these target names (default: all targets)",
    )
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Run the similarity command."""
    parser = build_argument_parser()
    arguments = parser.parse_args(argv)
    dataset_root = arguments.dataset_root.expanduser().resolve()
    output_root = arguments.output_root.expanduser().resolve()

    if not dataset_root.is_dir():
        parser.error(f"dataset directory does not exist: {dataset_root}")

    RDLogger.DisableLog("rdApp.*")

    targets = arguments.targets or find_targets(dataset_root)
    if not targets:
        print(f"No targets with {CRYSTAL_FILENAME} and {ACTIVES_FILENAME} below {dataset_root}", file=sys.stderr)
        return 1

    failures = []
    total_rows = 0
    for target in targets:
        try:
            rows, skipped = process_target(
                target,
                dataset_root,
                output_root,
                radius=arguments.radius,
                n_bits=arguments.n_bits,
            )
            total_rows += rows
            if skipped:
                print(f"     [WARN] skipped {skipped} actives with invalid SMILES", file=sys.stderr)
        except Exception as exc:
            failures.append((target, str(exc)))
            print(f"[FAILED] {target}: {exc}", file=sys.stderr)

    print(f"Processed {len(targets) - len(failures)}/{len(targets)} targets, {total_rows} actives in total.")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
