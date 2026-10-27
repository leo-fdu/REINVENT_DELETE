"""Assemble traceable Maestro preparation inputs; no chemical preparation.

Run from any directory. Existing output directories are never overwritten.
Only Python's standard library is required. Raw input files are never changed.
"""

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REFERENCES = Path("glide_preparation/real-world_dataset/pdb_complexes")
DEFAULT_OUTPUT = ROOT / "glide_preparation/real-world_dataset/receptors_for_upload"
# target, PDB, source, retained protein chains, ligand chain/name/number,
# additional retained HET residue names, reason
TARGETS = [
    ("aa2ar", "3EML", "original", "A", "A", "ZMA", 401, (), "Retain original receptor and its side-chain adjustments; omit STE, SO4 and waters."),
    ("abl1", "2HZI", "original", "AB", "A", "JIN", 600, (), "Preserve original A/B pocket fragments; use the A-chain ligand site."),
    ("aces", "1E66", "original", "A", "A", "HUX", 803, (), "Retain original receptor."),
    ("adrb1", "2VT4", "original", "B", "B", "P32", 400, (), "Select the B-chain monomer matching the original ligand coordinates; retain engineered construct."),
    ("akt1", "4GV1", "reference", "A", "A", "0XZ", 501, ("TPO",), "Use deposited A chain to restore TPO308; reference alternate conformations are retained for PrepWizard review."),
    ("cdk2", "1H00", "reference", "A", "A", "FAP", 1300, ("ACE",), "Restore deposited protein context and existing ACE cap; no missing loops or side chains are modeled."),
    ("cp3a4", "3NXU", "reference", "A", "A", "RIT", 600, ("HEM",), "Keep A-chain HEM/Fe; omit duplicate B monomer and DMS."),
    ("dpp4", "2I78", "original", "ABCD", "B", "KIQ", 901, (), "Preserve original multichain context; the selected pocket is on B, not A."),
    ("drd3", "3PBL", "original", "A", "A", "ETQ", 1200, (), "Correct pocket reference only: ETQ A/1200. Original design ligand and generated molecules are not replaced."),
    ("dyr", "3NXO", "reference", "A", "A", "D2B", 187, ("NDP",), "Retain NADPH NDP A/193; omit sulfates and waters."),
    ("egfr", "2RGP", "original", "A", "A", "HYZ", 1, (), "Retain original receptor; unresolved loops remain unresolved."),
    ("gria2", "3KGC", "original", "B", "B", "ZK1", 263, (), "Retain original B-chain numbering and adjusted coordinates; use ZK1 pocket, not GLU/LY7."),
    ("jak2", "3LPB", "reference", "A", "A", "NVB", 1133, ("PTR",), "Replace incorrect original JNK3 receptor with JAK2 A; retain PTR1007/PTR1008."),
    ("kif11", "3CJO", "reference", "A", "A", "K30", 1, ("ADP", "MG"), "Retain deposited A-chain ADP/Mg as receptor context; K30 is the pocket reference to exclude from the grid."),
    ("parp1", "3L3M", "original", "A", "A", "A92", 351, (), "Retain original receptor and its engineered construct."),
    ("ppara", "2P54", "original", "AB", "A", "735", 469, (), "Retain original receptor A and coactivator peptide B."),
]


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def residue(line):
    return line[17:20].strip(), line[21], int(line[22:26]), line[26]


def xyz(line):
    return tuple(float(line[start:start + 8]) for start in (30, 38, 46))


def is_polymer(line):
    # Modified amino acids/caps are HETATM records but must not terminate a chain.
    return line.startswith("ATOM  ") or line[17:20].strip() in {"TPO", "PTR", "ACE"}


def heavy_signature(lines):
    return Counter((line[76:78].strip().upper(), *xyz(line)) for line in lines
                   if line[76:78].strip().upper() not in ("H", "D"))


def mol2_signature(path):
    active = False
    values = []
    for line in path.read_text().splitlines():
        if line.startswith("@<TRIPOS>"):
            active = line == "@<TRIPOS>ATOM"
        elif active and line.strip():
            fields = line.split()
            element = fields[5].split(".")[0].upper()
            if element not in ("H", "D"):
                values.append((element, *map(float, fields[2:5])))
    return Counter(values)


def assemble(spec, output):
    target, pdb, kind, chains, ligchain, ligname, ligseq, keephets, reason = spec
    refpath = REFERENCES / f"{target}_{pdb}.pdb"
    srcpath = refpath if kind == "reference" else Path(f"real-world_dataset/{target}/receptor_out.pdb")
    sources = {p: (ROOT / p).read_text().splitlines() for p in (srcpath, refpath)}
    ligand_key = (ligname, ligchain, ligseq, " ")
    ligand = [l for l in sources[refpath] if l.startswith("HETATM") and residue(l) == ligand_key]
    if not ligand:
        raise ValueError(f"Missing pocket reference for {target}")
    design_match = heavy_signature(ligand) == mol2_signature(ROOT / f"real-world_dataset/{target}/crystal.mol2")
    if design_match != (target != "drd3"):
        raise ValueError(f"Unexpected original design/reference coordinate correspondence: {target}")

    records = [(srcpath, l) for l in sources[srcpath]
               if (l.startswith("ATOM  ") and l[21] in chains)
               or (l.startswith("HETATM") and l[21] in chains and l[17:20].strip() in keephets)]
    records += [(refpath, l) for l in ligand]
    ids = [(p, int(l[6:11])) for p, l in records]
    if len(set(ids)) != len(ids):
        raise ValueError(f"Duplicate source atom IDs: {target}")
    present = {residue(l) for _, l in records}
    lines = [f"TITLE     {target.upper()} CURATED COMPLEX FOR MAESTRO PREPARATION",
             "REMARK 900 NOT PREPARED: RUN PROTEIN PREPARATION BEFORE GRID GENERATION.",
             f"REMARK 900 EXCLUDE POCKET LIGAND {ligname} {ligchain} {ligseq} WHEN MAKING GRID.",
             "REMARK 900 COORDINATES UNCHANGED; ATOM SERIAL NUMBERS REASSIGNED."]
    # Retain trustworthy source sequence/missing-residue annotations only;
    # never transplant reference SEQRES onto renumbered/cropped original inputs.
    for line in sources[srcpath]:
        if line.startswith("SEQRES") and line[11] in chains:
            lines.append(line)
        elif line.startswith(("REMARK 465", "REMARK 470")):
            is_data = len(line) > 26 and line[15:18].strip().isalpha() and line[19].strip()
            if not is_data or line[19] in chains:
                lines.append(line)
    for p, source in sources.items():
        for line in source:
            if line.startswith("MODRES"):
                key = (line[12:15].strip(), line[16], int(line[18:22]), line[22])
                if key in present and line not in lines:
                    lines.append(line)
            elif line.startswith("SSBOND"):
                keys = [("CYS", line[15], int(line[17:21]), line[21]),
                        ("CYS", line[29], int(line[31:35]), line[35])]
                if all(k in present for k in keys) and line not in lines:
                    lines.append(line)
            elif line.startswith("LINK  "):
                keys = [(line[17:20].strip(), line[21], int(line[22:26]), line[26]),
                        (line[47:50].strip(), line[51], int(line[52:56]), line[56])]
                if all(k in present for k in keys) and line not in lines:
                    lines.append(line)
    serial = 1
    remap = {}
    previous = None
    for p, line in records:
        if previous and is_polymer(previous) and (not is_polymer(line) or previous[21] != line[21]):
            lines.append(f"TER   {serial:5d}      {previous[17:27]}")
            serial += 1
        remap[p, int(line[6:11])] = serial
        lines.append(line[:6] + f"{serial:5d}" + line[11:])
        serial += 1
        previous = line
    for p, source in sources.items():
        for line in source:
            if not line.startswith("CONECT"):
                continue
            old = [int(line[i:i + 5]) for i in range(6, len(line), 5) if line[i:i + 5].strip()]
            if old and (p, old[0]) in remap:
                neighbors = [remap[p, n] for n in old[1:] if (p, n) in remap]
                if neighbors:
                    lines.append("CONECT" + "".join(f"{n:5d}" for n in [remap[p, old[0]], *neighbors]))
    lines.append("END")
    dest = output / f"{target}.pdb"
    dest.write_text("\n".join(lines) + "\n")
    # Re-read to detect serial, coordinate, atom or component corruption.
    reread = dest.read_text().splitlines()
    atoms = [l for l in reread if l.startswith(("ATOM  ", "HETATM"))]
    if len(atoms) != len(records) or any(a[11:] != b[11:] for a, (_, b) in zip(atoms, records)):
        raise ValueError(f"Atom preservation failure: {target}")
    serials = {int(l[6:11]) for l in atoms}
    if len(serials) != len(atoms):
        raise ValueError(f"Nonunique output atom serials: {target}")
    for line in reread:
        if line.startswith("CONECT"):
            if not all(int(line[i:i + 5]) in serials for i in range(6, len(line), 5) if line[i:i + 5].strip()):
                raise ValueError(f"Dangling connectivity: {target}")
    ligcoords = [xyz(l) for l in ligand if l[76:78].strip() not in ("H", "D")]
    hets = sorted({residue(l)[:3] for l in atoms if l.startswith("HETATM")})
    return dict(target=target, pdb_id=pdb, output_file=dest.name, output_sha256=sha256(dest),
                protein_source=str(srcpath), protein_source_kind=kind,
                source_sha256={str(p): sha256(ROOT / p) for p in sources},
                retained_protein_chains=list(chains), pocket_reference=dict(
                    source=str(refpath), resname=ligname, author_chain=ligchain, residue_number=ligseq,
                    heavy_atom_centroid=[sum(v[i] for v in ligcoords) / len(ligcoords) for i in range(3)],
                    matches_original_design_heavy_atoms_and_coordinates=design_match),
                retained_het_residues=hets, atom_records=len(atoms),
                hydrogen_records=sum(l[76:78].strip() in ("H", "D") for l in atoms),
                altloc_atom_records=sum(bool(l[16].strip()) for l in atoms),
                reason=reason, status="curated_input_not_chemically_prepared_or_docking_validated")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=False)
    entries = [assemble(spec, args.output_dir) for spec in TARGETS]
    manifest = dict(schema_version=1, structures=entries,
                    policy="No coordinate optimization, hydrogen changes, loop modeling, altloc choice or grid creation. No waters retained. Each PDB includes one pocket-reference ligand; exclude it from the grid.")
    (args.output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    sums = [f"{e['output_sha256']}  {e['output_file']}" for e in entries]
    sums.append(f"{sha256(args.output_dir / 'manifest.json')}  manifest.json")
    (args.output_dir / "SHA256SUMS").write_text("\n".join(sums) + "\n")
    print(f"Wrote and validated {len(entries)} curated complexes in {args.output_dir}")
    for e in entries:
        print(e["target"], e["atom_records"], e["protein_source_kind"], e["retained_het_residues"])


if __name__ == "__main__":
    main()
