"""Chemical provenance and isolation checks for the additional design dataset."""
from pathlib import Path
import hashlib
import json
import sys

import pytest
from rdkit import Chem

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "real-world_dataset/task_design_manual"))
sys.path.insert(0, str(ROOT / "src"))
import design
from dataset import DatasetAtlas
from prepare_shallowbench_design import table
from export_manual_design import export


@pytest.fixture
def collection(monkeypatch):
    original = design.atlas
    atlas = DatasetAtlas(ROOT / "shallowbench-7targets/task_design/dataset.json", original)
    monkeypatch.setattr(design, "atlas", atlas)
    design.ligand.cache_clear()
    yield atlas
    design.ligand.cache_clear()


def test_source_atoms_coordinates_bonds_and_declared_stereo(collection):
    assert len(collection.entries) == 8
    assert len(set(collection.PDB.values())) == 7
    for code, provenance in collection.data["ccd_sources"].items():
        ccd = ROOT / "shallowbench-7targets/task_design/ccd" / f"{code}.cif"
        assert hashlib.sha256(ccd.read_bytes()).hexdigest() == provenance["sha256"]
    for target, info in collection.entries.items():
        mol = design.ligand(target)
        assert Chem.MolToSmiles(mol) == info["smiles"]
        graph_file = Chem.SDMolSupplier(str(ROOT / info["design_graph_sdf"]))[0]
        assert Chem.MolToSmiles(graph_file) == info["smiles"]
        source = ROOT / info["source_pdb"]
        assert hashlib.sha256(source.read_bytes()).hexdigest() == info["source_pdb_sha256"]
        records = {int(line[6:11]): line for line in source.read_text().splitlines()
                   if line.startswith(("ATOM  ", "HETATM"))}
        raw = Chem.SDMolSupplier(str(collection.source_path(target)), sanitize=False)[0]
        rows = json.loads((ROOT / info["atom_mapping"]).read_text())
        assert len(rows) == mol.GetNumAtoms()
        for row in rows:
            if row.get("virtual_boundary"):
                assert target == "4n4f"
                continue
            line = records[row["pdb_serial"]]
            expected = [float(line[i:i+8]) for i in (30,38,46)]
            assert row["xyz"] == expected
            assert list(raw.GetConformer().GetAtomPosition(row["source_id"]-1)) == expected
            ccd = ROOT / "shallowbench-7targets/task_design/ccd" / (row["ccd"]+".cif")
            atom = next(a for a in table(ccd, "_chem_comp_atom") if a["atom_id"] == row["atom_name"])
            observed = mol.GetAtomWithIdx(row["source_id"]-1)
            if atom["pdbx_stereo_config"] in ("R", "S"):
                assert observed.GetProp("_CIPCode") == atom["pdbx_stereo_config"]
        if target == "3vhw":
            assert sum(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds()) == 3
            names = {row["atom_name"]: row["source_id"]-1 for row in rows}
            for row in table(ROOT / "shallowbench-7targets/task_design/ccd/VHW.cif", "_chem_comp_bond"):
                if row["pdbx_stereo_config"] in ("E", "Z"):
                    bond = mol.GetBondBetweenAtoms(names[row["atom_id_1"]], names[row["atom_id_2"]])
                    assert bond.GetProp("_CIPCode") == row["pdbx_stereo_config"]
        if target == "4o5k":
            assert all(a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED
                       for a in mol.GetAtoms() if a.GetSymbol() == "P")


def test_incomplete_peptide_cannot_be_submitted_as_designed(collection, tmp_path):
    with pytest.raises(ValueError, match="完整设计输入"):
        design.cuts_for("4n4f", "libinvent", [[1,2]], complete=True)
    result = design.target_info("4n4f")
    assert not result["editable"]
    assert 'class="bond-hit"' not in result["svg"]
    folder = design.save_snapshot({"4n4f": {"libinvent": {"skip": True, "note": "暂缓"}}}, tmp_path)
    data = json.loads((folder / "designs.json").read_text())
    assert data["status_counts"] == {"designed": 0, "skipped": 1, "incomplete": 0, "not_started": 15}
    assert data["targets"]["4n4f"]["provenance"]["editable"] is False
    assert "source_ligand" in data["targets"]["4qyg"]
    assert "source_mol2" not in data["targets"]["4qyg"]


def test_actual_prior_validation_and_nonoverwriting_records(collection, tmp_path):
    cuts = [(14,17)]
    parts = design.fragments_for("4qyg", cuts)
    retained = max(parts, key=lambda m: m.GetNumHeavyAtoms())
    spec = {"cuts": [list(c) for c in cuts], "retained": [sorted(design.atlas.ids(retained))]}
    task = design.build_task("4qyg", "libinvent", spec, design.models()["libinvent"])
    assert task["validation"] == {"connectivity_match": True, "stereochemistry_match": True}
    first = design.save_snapshot({"4qyg": {"libinvent": spec}}, tmp_path)
    second = design.save_snapshot({"4qyg": {"libinvent": spec}}, tmp_path)
    assert first != second
    assert len(list((first / "figures").glob("*.svg"))) == 16


def test_tampered_source_rejected(collection, tmp_path):
    # Hash validation remains enforced for new datasets as for the original atlas.
    collection.RAW_SHA256["4qyg"] = "0" * 64
    with pytest.raises(ValueError, match="Ligand changed"):
        collection.load_ligand("4qyg")


def test_generic_ligand_export_preserves_provenance(collection, tmp_path):
    selections = {target: {mode: {"skip": True, "note": "test"} for mode in design.MODES}
                  for target in collection.CODES}
    parts = design.fragments_for("4qyg", [(14,17)])
    kept = max(parts, key=lambda p: p.GetNumHeavyAtoms())
    selections["4qyg"]["libinvent"] = {"cuts": [[14,17]], "retained": [sorted(design.atlas.ids(kept))]}
    record = design.save_snapshot(selections, tmp_path)
    output = tmp_path / "export"
    manifest = export(record/"designs.json", output)
    task = next(t for t in manifest["tasks"] if t["status"] == "designed")
    assert task["source_ligand"].endswith("/4qyg/ligand.sdf")
    assert task["provenance"]["pocket_mapping_confirmed"] is False
    assert (output/task["smiles_file"]).read_text().strip() == task["input_smiles"]
