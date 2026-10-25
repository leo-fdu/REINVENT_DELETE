"""Optional, explicitly configured ligand collection for the existing editor."""
from pathlib import Path
import hashlib
import json
import re

from rdkit import Chem
from rdkit.Chem import rdDepictor, rdCIPLabeler


class DatasetAtlas:
    def __init__(self, path, helpers):
        self.helpers = helpers
        self.ROOT = helpers.ROOT
        self.data = json.loads(Path(path).read_text())
        if self.data.get("schema_version") != 1:
            raise ValueError("Unknown ligand dataset schema")
        self.entries = self.data["targets"]
        if not self.entries or any(not re.fullmatch(r"[a-z0-9]+", key) for key in self.entries):
            raise ValueError("Invalid dataset target identifiers")
        self.CODES = {key: info["ccd"] for key, info in self.entries.items()}
        self.PDB = {key: info["pdb"] for key, info in self.entries.items()}
        self.RAW_SHA256 = {key: info["source_sha256"] for key, info in self.entries.items()}
        self.title = self.data.get("title", "手动任务设计")

    def __getattr__(self, name):
        return getattr(self.helpers, name)

    def source_path(self, target):
        path = (self.ROOT / self.entries[target]["source_ligand"]).resolve()
        if not path.is_relative_to(self.ROOT):
            raise ValueError("Ligand must belong to the project")
        return path

    def load_ligand(self, target):
        info = self.entries[target]
        path = self.source_path(target)
        if hashlib.sha256(path.read_bytes()).hexdigest() != self.RAW_SHA256[target]:
            raise ValueError(f"Ligand changed: {target}")
        coordinates = Chem.SDMolSupplier(str(path), sanitize=False, removeHs=False)[0]
        if coordinates is None or not coordinates.HasProp("authoritative_mapped_smiles"):
            raise ValueError(f"Missing authoritative atom-mapped graph: {target}")
        if coordinates.GetProp("authoritative_mapped_smiles") != info["mapped_smiles"]:
            raise ValueError("Dataset graph differs from SDF provenance")
        mol = Chem.MolFromSmiles(info["mapped_smiles"])
        if mol is None:
            raise ValueError("Invalid reference graph")
        mol = Chem.RenumberAtoms(mol, sorted(range(mol.GetNumAtoms()),
                                           key=lambda i: mol.GetAtomWithIdx(i).GetAtomMapNum()))
        if mol.GetNumAtoms() != coordinates.GetNumAtoms():
            raise ValueError("Coordinate/graph atom count mismatch")
        for atom, raw in zip(mol.GetAtoms(), coordinates.GetAtoms()):
            if atom.GetAtomMapNum() != atom.GetIdx()+1 or atom.GetAtomicNum() != raw.GetAtomicNum():
                raise ValueError("Coordinate/graph atom mapping mismatch")
            atom.SetIntProp("source_id", atom.GetIdx()+1)
            atom.SetAtomMapNum(0)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        rdCIPLabeler.AssignCIPLabels(mol)
        # Drawing never changes the archived crystallographic coordinates.
        rdDepictor.Compute2DCoords(mol)
        return mol
