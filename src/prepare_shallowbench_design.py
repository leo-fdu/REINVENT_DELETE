"""Extract coordinate-preserving design candidates using archived RCSB CCD bonds.

No protein preparation, pocket matching, sampling or ligand optimization is done.
The incomplete 4N4F peptide is a display-only graph with explicit open boundaries.
"""
from pathlib import Path
from argparse import ArgumentParser
from collections import defaultdict
import hashlib
import json
import math
import re
import shlex

from rdkit import Chem, rdBase
from rdkit.Chem import rdCIPLabeler, rdDepictor

ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "shallowbench-7targets"
OUT = BASE / "task_design"
IDS = {"4qyg": 135853, "4u1o": 81720, "2roy": 68990, "4n4f": 9015,
       "3vhw": 115690, "4o5k": 147001, "4pwj": 69315}
NOTES = {
    "4qyg": "3DW A/301：仅选择一个晶体副本供二维设计；ShallowBench 口袋对应关系尚未确认。保留稠合核心，优先考虑侧链或环系之间的非环单键。",
    "4u1ofwf": "FWF A/302：与 KAI 是不同配体/位点的候选，尚未确认条目 81720 对应哪一个。先设计不代表选定正式实验配体。",
    "4u1okai": "KAI A/301：与 FWF 分开记录；环内键不能切，两个端片段可能较小，且切键可能影响手性恢复。",
    "2roy": "P28 A/128：保留原始硝基、羧酸和立体化学；晶体副本与多链口袋仍待确认。",
    "4n4f": "仅展示已解析的 C 链 7–12；两端 * 表示未解析肽链的延续，不是已设计连接点。未补建原子或封端；此项暂不可切割，可记录备注/跳过。",
    "3vhw": "VHW A/1：多手性中心和共轭双键；避开会破坏保留片段立体化学的切键，记录时核对 stereo 标记。",
    "4o5k": "0KX A/405：含磷酸基和核糖结构；保留 CCD 化学状态，尚未选定溶液质子化状态。PLANET 不表示 DNA/金属环境。",
    "4pwj": "30Z A/201：仅选择一个晶体副本供二维设计；两个芳环之间有非环连接链，需检查手性恢复及多链口袋。",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path, prefix):
    """Read the simple atom/bond loops of the archived CCD files, fail closed."""
    lines = path.read_text().splitlines()
    start = next(i for i, line in enumerate(lines) if line.startswith(prefix + "."))
    columns = []
    while start < len(lines) and lines[start].startswith(prefix + "."):
        columns.append(lines[start].strip().split(".", 1)[1])
        start += 1
    rows = []
    while start < len(lines) and not lines[start].startswith(("#", "loop_", "_")):
        fields = shlex.split(lines[start])
        if fields:
            if len(fields) != len(columns):
                raise ValueError(f"Unexpected CCD loop row in {path}: {lines[start]}")
            rows.append(dict(zip(columns, fields)))
        start += 1
    return rows


def template(code):
    path = OUT / "ccd" / f"{code}.cif"
    atoms = table(path, "_chem_comp_atom")
    bonds = table(path, "_chem_comp_bond")
    mol = Chem.RWMol()
    index = {}
    coords = Chem.Conformer(len(atoms))
    for row in atoms:
        atom = Chem.Atom(row["type_symbol"])
        atom.SetFormalCharge(int(row["charge"]))
        atom.SetProp("ccd_name", row["atom_id"])
        i = mol.AddAtom(atom)
        index[row["atom_id"]] = i
        coords.SetAtomPosition(i, [float(row[f"pdbx_model_Cartn_{axis}_ideal"]) for axis in "xyz"])
    orders = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE,
              "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
    for row in bonds:
        mol.AddBond(index[row["atom_id_1"]], index[row["atom_id_2"]], orders[row["value_order"]])
    result = mol.GetMol()
    result.AddConformer(coords)
    Chem.SanitizeMol(result)
    Chem.AssignStereochemistryFrom3D(result)
    result = Chem.RemoveHs(result)
    declared = {row["atom_id"]: row["pdbx_stereo_config"] for row in atoms}
    for atom in result.GetAtoms():
        if declared[atom.GetProp("ccd_name")] not in ("R", "S"):
            atom.SetChiralTag(Chem.ChiralType.CHI_UNSPECIFIED)
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    rdCIPLabeler.AssignCIPLabels(result)
    for atom in result.GetAtoms():
        expected = declared[atom.GetProp("ccd_name")]
        if expected in ("R", "S"):
            if not atom.HasProp("_CIPCode"):
                raise ValueError(f"Cannot assign declared stereo: {code}/{atom.GetProp('ccd_name')}")
            if atom.GetProp("_CIPCode") != expected:
                atom.InvertChirality()
    Chem.AssignStereochemistry(result, cleanIt=True, force=True)
    rdCIPLabeler.AssignCIPLabels(result)
    for atom in result.GetAtoms():
        expected = declared[atom.GetProp("ccd_name")]
        if expected in ("R", "S") and atom.GetProp("_CIPCode") != expected:
            raise ValueError(f"Cannot reproduce CCD stereo: {code}/{atom.GetProp('ccd_name')}")
    return result, {a.GetProp("ccd_name"): a.GetIdx() for a in result.GetAtoms()}


def selected_lines(path, residues):
    groups = defaultdict(list)
    for line in path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and (line[21], int(line[22:26]), line[17:20].strip()) in residues:
            if line[76:78].strip() not in ("H", "D"):
                groups[(line[21], line[22:27], line[17:20])].append(line)
    selected = []
    for lines in groups.values():
        labels = {line[16] for line in lines} - {" "}
        label = min(labels, key=lambda alt: (-sum(float(l[54:60]) for l in lines if l[16] == alt) /
                                           sum(l[16] == alt for l in lines), alt)) if labels else " "
        selected += [line for line in lines if line[16] in (" ", label)]
    return sorted(selected, key=lambda line: int(line[6:11]))


def extract(path, residues, peptide=False):
    lines = selected_lines(path, residues)
    groups = defaultdict(list)
    for line in lines:
        groups[(line[21], int(line[22:26]), line[17:20].strip())].append(line)
    if set(groups) != set(residues):
        raise ValueError("Missing residue")
    result = Chem.RWMol()
    atom_rows = []
    residue_indices = {}
    for residue, records in sorted(groups.items()):
        ref, names = template(residue[2])
        observed = [l[12:16].strip() for l in records]
        if len(set(observed)) != len(observed) or set(observed) - set(names):
            raise ValueError(f"Duplicate or unmatched atom names: {residue}")
        missing = set(names) - set(observed)
        if missing and not (peptide and missing == {"OXT"}):
            raise ValueError(f"Missing heavy atoms: {residue}: {missing}")
        # Renumbering also updates chiral neighbor order before copying atoms.
        order = [names[name] for name in observed] + [i for i in range(ref.GetNumAtoms()) if i not in {names[n] for n in observed}]
        ref = Chem.RenumberAtoms(ref, order)
        start = result.GetNumAtoms()
        residue_indices[residue[1]] = {}
        for i, line in enumerate(records):
            atom = Chem.Atom(ref.GetAtomWithIdx(i))
            if atom.GetSymbol().upper() != line[76:78].strip().upper():
                raise ValueError("Element mismatch")
            atom.SetIntProp("source_id", start + i + 1)
            result.AddAtom(atom)
            name = line[12:16].strip()
            residue_indices[residue[1]][name] = start + i
            atom_rows.append({"source_id": start + i + 1, "pdb_serial": int(line[6:11]),
                              "chain": line[21], "residue": residue[1], "ccd": residue[2],
                              "atom_name": name, "altloc": line[16].strip(),
                              "xyz": [float(line[j:j+8]) for j in (30, 38, 46)]})
        stereo_bonds = []
        for bond in ref.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            if max(a, b) >= len(records):
                continue
            result.AddBond(start+a, start+b, bond.GetBondType())
            new = result.GetBondBetweenAtoms(start+a, start+b)
            new.SetIsAromatic(bond.GetIsAromatic())
            stereo_atoms = list(bond.GetStereoAtoms())
            if stereo_atoms and max(stereo_atoms) < len(records):
                stereo_bonds.append((start+a, start+b, stereo_atoms, bond.GetStereo()))
        for a, b, stereo_atoms, stereo in stereo_bonds:
            new = result.GetBondBetweenAtoms(a, b)
            new.SetStereoAtoms(*(start+i for i in stereo_atoms))
            new.SetStereo(stereo)
    if peptide:
        numbers = sorted(residue_indices)
        for previous, current in zip(numbers, numbers[1:]):
            if current != previous + 1:
                raise ValueError("Non-contiguous observed peptide")
            result.AddBond(residue_indices[previous]["C"], residue_indices[current]["N"], Chem.BondType.SINGLE)
        for number in numbers:
            atom = result.GetAtomWithIdx(residue_indices[number]["N"])
            atom.SetNumExplicitHs(max(0, atom.GetNumExplicitHs()-1))
        for endpoint in (residue_indices[numbers[0]]["N"], residue_indices[numbers[-1]]["C"]):
            dummy = Chem.Atom(0)
            dummy.SetIntProp("source_id", result.GetNumAtoms()+1)
            i = result.AddAtom(dummy)
            result.AddBond(endpoint, i, Chem.BondType.SINGLE)
            # Display boundary only, not a modeled atom or crystallographic coordinate.
            atom_rows.append({"source_id": i+1, "virtual_boundary": True, "attached_to": endpoint+1})
    mol = result.GetMol()
    Chem.SanitizeMol(mol)
    Chem.SetDoubleBondNeighborDirections(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    if len(Chem.GetMolFrags(mol)) != 1 or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("Disconnected or radical extraction")
    conf = Chem.Conformer(mol.GetNumAtoms())
    for row in atom_rows:
        xyz = row.get("xyz")
        if xyz is None:
            p = atom_rows[row["attached_to"]-1]["xyz"]
            xyz = [p[0]+1.3, p[1], p[2]]
        conf.SetAtomPosition(row["source_id"]-1, xyz)
    mol.AddConformer(conf)
    return mol, atom_rows, lines


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true", help="Explicitly regenerate derived inputs; may invalidate existing designs")
    args = parser.parse_args()
    if (OUT/"dataset.json").exists() and not args.rebuild:
        parser.error("Dataset already exists. Review saved designs before using --rebuild.")
    OUT.mkdir(exist_ok=True)
    entries = {}
    specs = [("4qyg", "4QYG", [("A",301,"3DW")]),
             ("4u1ofwf", "4U1O", [("A",302,"FWF")]),
             ("4u1okai", "4U1O", [("A",301,"KAI")]),
             ("2roy", "2ROY", [("A",128,"P28")]),
             ("4n4f", "4n4f", [("C",n,c) for n,c in zip(range(7,13),["GLY","LYS","GLY","LEU","GLY","ALY"])]),
             ("3vhw", "3VHW", [("A",1,"VHW")]),
             ("4o5k", "4O5K", [("A",405,"0KX")]),
             ("4pwj", "4PWJ", [("A",201,"30Z")])]
    listing = (BASE / "785_target").read_text()
    for key,pdb,residues in specs:
        path = BASE / f"{pdb}.pdb"
        mol, rows, records = extract(path, residues, peptide=key=="4n4f")
        folder = OUT / "ligands" / key
        folder.mkdir(parents=True, exist_ok=True)
        mol.SetProp("_Name", key)
        mol.SetProp("provenance", "CCD bonds/stereo; PDB heavy atom coordinates; see atoms.json")
        graph = Chem.Mol(mol)
        for atom in graph.GetAtoms():
            atom.SetAtomMapNum(atom.GetIdx()+1)
        mapped_smiles = Chem.MolToSmiles(graph)
        mol.SetProp("authoritative_mapped_smiles", mapped_smiles)
        with Chem.SDWriter(str(folder / "ligand.sdf")) as writer:
            writer.write(mol)
        depiction = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(depiction)
        with Chem.SDWriter(str(folder / "design_graph.sdf")) as writer:
            writer.write(depiction)
        (folder / "observed.pdb").write_text("\n".join(records)+"\nEND\n")
        (folder / "atoms.json").write_text(json.dumps(rows,indent=2)+"\n")
        identifier = IDS[pdb.lower()]
        block = next(b for b in listing.split("--- Target") if f"  id: {identifier}\n" in b)
        center = list(map(float,re.search(r"values=\[([^]]+)\]",block).group(1).split()))
        coords = [r["xyz"] for r in rows if "xyz" in r]
        centroid = [sum(p[i] for p in coords)/len(coords) for i in range(3)]
        source = folder / "ligand.sdf"
        entries[key] = {"pdb": pdb.upper(), "ccd": "+".join(r[2] for r in residues),
                        "label": f"{pdb.upper()} · {residues[0][2]}" if key!="4n4f" else "4N4F · 肽片段",
                        "source_ligand": str(source.relative_to(ROOT)), "source_sha256": sha(source),
                        "design_graph_sdf": str((folder/"design_graph.sdf").relative_to(ROOT)),
                        "source_pdb": str(path.relative_to(ROOT)), "source_pdb_sha256": sha(path),
                        "residues": residues, "atom_mapping": str((folder/"atoms.json").relative_to(ROOT)),
                        "smiles": Chem.MolToSmiles(mol), "mapped_smiles": mapped_smiles,
                        "heavy_atoms": mol.GetNumHeavyAtoms(),
                        "shallowbench_id": identifier, "shallowbench_center": center,
                        "observed_centroid": centroid, "center_distance_angstrom": math.dist(center,centroid),
                        "editable": key!="4n4f", "warning": NOTES[key], "pocket_mapping_confirmed": False}
        if key == "3vhw":
            # Initial live preview used legacy CIP labels; only stereo changed.
            # Atom IDs and connectivity are identical, so user cuts can migrate.
            entries[key]["draft_compatible_hashes"] = ["3e77ffea8a2dafa6ee94db0aff307a6fe34d9cfb247bfd61ed7506bef1bbb84c"]
        print(key, mol.GetNumHeavyAtoms(), Chem.MolToSmiles(mol))
    data = {"schema_version": 1, "title": "ShallowBench · 7 个结构 / 8 个配体候选",
            "rdkit_version": rdBase.rdkitVersion, "targets": entries,
            "ccd_sources": {p.stem: {"url": f"https://files.rcsb.org/ligands/download/{p.name}",
                                     "sha256": sha(p)} for p in sorted((OUT/"ccd").glob("*.cif"))}}
    (OUT/"dataset.json").write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n")


if __name__ == "__main__":
    main()
