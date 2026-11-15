"""Export the three supplied, successfully prepared complexes as explicit-H PDBs.

Requires RDKit. Original PDBs and Maestro outputs are never overwritten.
The source structures already contain the chemical preparation; this exporter
does not infer protonation, build loops, or run a new minimization.
"""

import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import io
import json
from pathlib import Path
import re

from rdkit import Chem, rdBase


TARGETS = {"aa2ar": (5, "ZMA", 401), "akt1": (6, "0XZ", 501),
           "cdk2": (7, "FAP", 1300)}
STANDARD = set("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split())
POLYMER = STANDARD | {"TPO", "ACE"}


def residue(atom):
    info = atom.GetPDBResidueInfo()
    if info is None:
        raise ValueError("Missing PDB residue metadata")
    return (info.GetChainId().strip(), info.GetResidueNumber(),
            info.GetInsertionCode().strip(), info.GetResidueName().strip())


def read_mae(path):
    with gzip.open(path, "rb") as handle:
        data = handle.read()
    if len(re.findall(rb"(?m)^f_m_ct\s*\{", data)) != 1:
        raise ValueError(f"Expected exactly one readable complex: {path}")
    mol = next(Chem.MaeMolSupplier(io.BytesIO(data), sanitize=False, removeHs=False))
    if mol is None:
        raise ValueError(f"Unreadable complex: {path}")
    if not mol.HasProp("b_ppw_prepared") or mol.GetProp("b_ppw_prepared") != "1":
        raise ValueError(f"Source is not marked prepared: {path}")
    Chem.SanitizeMol(Chem.Mol(mol))
    return mol


def ordered_segments(mol):
    """Group hydrogens with residues; order polymers by actual peptide bonds.

    Number sorting would break AA2AR's fusion numbering. Peptide connectivity
    also puts AKT1 TPO308 between residues 307 and 309, despite source ordering.
    """
    groups = defaultdict(list)
    for atom in mol.GetAtoms():
        groups[residue(atom)].append(atom.GetIdx())
    successors, predecessors = {}, {}
    for bond in mol.GetBonds():
        a, b = bond.GetBeginAtom(), bond.GetEndAtom()
        if a.GetPDBResidueInfo().GetName().strip() == "N":
            a, b = b, a
        if (a.GetPDBResidueInfo().GetName().strip(),
                b.GetPDBResidueInfo().GetName().strip()) != ("C", "N"):
            continue
        ra, rb = residue(a), residue(b)
        if ra == rb or ra[3] not in POLYMER or rb[3] not in POLYMER:
            continue
        if ra in successors or rb in predecessors or ra[0] != rb[0]:
            raise ValueError("Unexpected branched or cross-chain peptide")
        successors[ra], predecessors[rb] = rb, ra
    segments, seen = [], set()
    for root in groups:
        if root[3] not in POLYMER or root in predecessors:
            continue
        segment, current = [], root
        while current is not None:
            if current in seen:
                raise ValueError("Cyclic peptide ordering")
            seen.add(current)
            segment.extend(groups[current])
            current = successors.get(current)
        segments.append((segment, True))
    if seen != {r for r in groups if r[3] in POLYMER}:
        raise ValueError("Unordered polymer residues")
    segments.extend((indices, False) for r, indices in groups.items()
                    if r[3] not in POLYMER)
    return segments


def write_pdb(mol, path, target, original):
    _, ligand, number = TARGETS[target]
    lines = [f"TITLE     {target.upper()} PREPARED COMPLEX WITH EXPLICIT HYDROGENS",
             "REMARK 900 DERIVED FROM SUCCESSFUL SCHRODINGER PROTEIN PREPARATION.",
             "REMARK 900 RETAINS PREPARED COORDINATES, FORMAL CHARGES AND BONDS.",
             "REMARK 900 COORDINATES ROUNDED TO PDB PRECISION (0.001 ANGSTROM).",
             f"REMARK 900 EXCLUDE POCKET LIGAND {ligand} A {number} WHEN MAKING GRID."]
    # Residue identity annotations remain valid; old missing-atom annotations
    # and old atom serial-dependent records must not be copied to repaired PDBs.
    lines.extend(l for l in original.read_text().splitlines() if l.startswith("MODRES"))
    serials, order, serial = {}, [], 1
    conformer = mol.GetConformer()
    for segment, polymer in ordered_segments(mol):
        for index in segment:
            atom = mol.GetAtomWithIdx(index)
            info = atom.GetPDBResidueInfo()
            chain, resnum, insertion, resname = residue(atom)
            name, element = info.GetName(), atom.GetSymbol()
            charge = atom.GetFormalCharge()
            if (len(name) != 4 or len(chain) > 1 or len(insertion) > 1
                    or len(resname) > 3 or abs(charge) > 9 or not -999 <= resnum <= 9999):
                raise ValueError("Atom metadata exceeds PDB field widths")
            charge_text = f"{abs(charge)}{'+' if charge > 0 else '-'}" if charge else "  "
            x, y, z = conformer.GetAtomPosition(index)
            record = "ATOM  " if resname in STANDARD else "HETATM"
            line = (f"{record}{serial:5d} {name} {resname:>3} {chain:1}{resnum:4d}{insertion:1}   "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{info.GetOccupancy():6.2f}"
                    f"{info.GetTempFactor():6.2f}          {element:>2}{charge_text}")
            if len(line) != 80:
                raise ValueError("PDB atom record overflow")
            lines.append(line)
            serials[index] = serial
            order.append(index)
            serial += 1
        if polymer:
            chain, resnum, insertion, resname = residue(mol.GetAtomWithIdx(segment[-1]))
            lines.append(f"TER   {serial:5d}      {resname:>3} {chain:1}{resnum:4d}{insertion:1}")
            serial += 1
    # Write every bond in both directions; repeated neighbors encode bond order.
    # This includes S-S, peptide links to TPO/ACE, and all explicit hydrogen bonds.
    neighbors = defaultdict(list)
    for bond in mol.GetBonds():
        bond_order = bond.GetBondTypeAsDouble()
        if bond_order not in (1, 2, 3):
            raise ValueError(f"Unsupported bond order for PDB: {bond_order}")
        a, b = serials[bond.GetBeginAtomIdx()], serials[bond.GetEndAtomIdx()]
        neighbors[a].extend([b] * int(bond_order))
        neighbors[b].extend([a] * int(bond_order))
    for atom_serial, bonded in sorted(neighbors.items()):
        bonded.sort()
        for offset in range(0, len(bonded), 4):
            lines.append("CONECT" + f"{atom_serial:5d}" +
                         "".join(f"{i:5d}" for i in bonded[offset:offset + 4]))
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")
    return Chem.RenumberAtoms(mol, order)


def validate(path, expected):
    # Disable proximity bonding: connectivity must survive through CONECT alone.
    actual = Chem.MolFromPDBFile(str(path), sanitize=False, removeHs=False,
                               proximityBonding=False)
    if actual is None or actual.GetNumAtoms() != expected.GetNumAtoms():
        raise ValueError("PDB readback lost atoms")
    Chem.SanitizeMol(Chem.Mol(actual))
    for a, b in zip(expected.GetAtoms(), actual.GetAtoms()):
        if (a.GetAtomicNum(), a.GetFormalCharge(), residue(a), a.GetPDBResidueInfo().GetName()) != (
                b.GetAtomicNum(), b.GetFormalCharge(), residue(b), b.GetPDBResidueInfo().GetName()):
            raise ValueError("PDB readback changed atom identity or charge")
    def bonds(m):
        return {(min(b.GetBeginAtomIdx(), b.GetEndAtomIdx()),
                 max(b.GetBeginAtomIdx(), b.GetEndAtomIdx()), b.GetBondTypeAsDouble())
                for b in m.GetBonds()}
    if bonds(expected) != bonds(actual):
        raise ValueError("PDB readback changed bond connectivity or order")
    displacement = max((expected.GetConformer().GetAtomPosition(i) -
                        actual.GetConformer().GetAtomPosition(i)).Length()
                       for i in range(actual.GetNumAtoms()))
    if displacement > 0.000867:
        raise ValueError("Coordinate differences exceed PDB rounding")
    keys = [(residue(a), a.GetPDBResidueInfo().GetName()) for a in actual.GetAtoms()]
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate residue atom names")
    return {"atoms": actual.GetNumAtoms(), "bonds": actual.GetNumBonds(),
            "hydrogens": sum(a.GetAtomicNum() == 1 for a in actual.GetAtoms()),
            "formal_charge": sum(a.GetFormalCharge() for a in actual.GetAtoms()),
            "max_coordinate_rounding_angstrom": displacement,
            "rdkit_valence_and_exact_bond_readback": "passed"}


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preparation-root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    args.output_dir.mkdir(parents=True, exist_ok=False)
    report = {"rdkit_version": rdBase.rdkitVersion, "glide_rerun_verified": False,
              "coordinate_policy": "prepared coordinates, rounded to 0.001 A; no further minimization",
              "structures": []}
    for target, (job, ligand, number) in TARGETS.items():
        source = args.preparation_root / f"proteinprep_{job}" / f"proteinprep_{job}-out.maegz"
        original = root / "glide_preparation/real-world_dataset/receptors_for_upload" / f"{target}.pdb"
        dest = args.output_dir / f"{target}_fixed.pdb"
        mol = read_mae(source)
        expected = write_pdb(mol, dest, target, original)
        stats = validate(dest, expected)
        retained = Counter(residue(a) for a in mol.GetAtoms())
        if not retained[("A", number, "", ligand)]:
            raise ValueError("Pocket reference ligand missing")
        report["structures"].append({"target": target, "source": str(source),
            "source_sha256": sha256(source), "original_pdb_sha256": sha256(original),
            "output": dest.name, "output_sha256": sha256(dest),
            "pocket_reference": [ligand, "A", number], **stats})
        print(target, stats)
    (args.output_dir / "validation.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
