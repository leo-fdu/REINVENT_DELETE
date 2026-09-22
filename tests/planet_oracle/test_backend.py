import importlib
from pathlib import Path
import pytest
from rdkit import Chem

from planet_oracle.chemistry import checked_pocket, prepare_molecule

PROJECT = Path(__file__).resolve().parents[2]


@pytest.fixture
def upstream(monkeypatch):
    root = PROJECT / "PLANET"
    if not (root / "chemutils.py").exists():
        pytest.skip("Local PLANET source required for upstream graph/pocket checks")
    monkeypatch.syspath_prepend(str(root))
    return importlib.import_module("chemutils")


def pdb_atom(serial, atom="CA", residue="ALA", number=1, x=0.0):
    element = "N" if atom == "N" else "C"
    return (f"ATOM  {serial:5d} {atom:>4s} {residue:3s} A{number:4d}    "
            f"{x:8.3f}{0.:8.3f}{0.:8.3f}{1.:6.2f}{20.:6.2f}          {element:>2s}  \n")


def pocket_cfg(tmp_path, contents):
    path = tmp_path / "protein.pdb"
    path.write_text(contents)
    return {"protein_pdb": str(path), "center": [0., 0., 0.]}


@pytest.mark.parametrize("smiles,error", [("", "empty_molecule"), (" ", "empty_molecule"),
    ("invalid", "invalid_smiles"), ("C*", "dummy_atom"), ("C.C", "disconnected_structure"),
    ("C(C)(C)(C)(C)C", "invalid_smiles"),
    ("[Fe](C)(C)(C)(C)(C)(C)(C)(C)(C)(C)C", "graph_capacity_exceeded")])
def test_invalid_molecules(smiles, error):
    mol, observed = prepare_molecule(smiles, 10)
    assert mol is None and observed == error


def test_add_hydrogens_charge_stereo_and_graph(upstream):
    smi = "[NH3+][C@@H](C)C(=O)[O-]"
    original = Chem.MolFromSmiles(smi)
    mol, error = prepare_molecule(smi, upstream.MAX_NB)
    assert error is None and mol.GetNumAtoms() > original.GetNumAtoms()
    assert [a.GetFormalCharge() for a in mol.GetAtoms()][:original.GetNumAtoms()] == [a.GetFormalCharge() for a in original.GetAtoms()]
    assert Chem.MolToSmiles(Chem.RemoveHs(mol), isomericSmiles=True) == Chem.MolToSmiles(original, isomericSmiles=True)
    graph = upstream.mol_batch_to_graph([mol], auto_detect=False)
    official = upstream.mol_batch_to_graph([Chem.AddHs(Chem.MolFromSmiles(smi))], auto_detect=False)
    for observed, expected in zip(graph[:4], official[:4]):
        assert observed.equal(expected)
    assert graph[4] == [(0, mol.GetNumAtoms())]
    other, _ = prepare_molecule(smi.replace("@@", "@"), upstream.MAX_NB)
    assert not upstream.mol_batch_to_graph([other], auto_detect=False)[0].equal(graph[0])


def test_pocket_selection_and_center(upstream, tmp_path):
    cfg = pocket_cfg(tmp_path, pdb_atom(1) + pdb_atom(2, number=2, x=12.) + pdb_atom(3, number=3, x=12.1))
    pocket, info = checked_pocket(upstream, cfg)
    assert pocket.res_count == 2 and info["residues"] == [["1", "A"], ["2", "A"]]
    assert info["radius_angstrom"] == 12
    mol = Chem.MolFromSmiles("C")
    conformer = Chem.Conformer(mol.GetNumAtoms())
    conformer.SetAtomPosition(0, (0, 0, 0))
    mol.AddConformer(conformer)
    sdf = tmp_path / "reference.sdf"
    with Chem.SDWriter(str(sdf)) as writer:
        writer.write(mol)
    del cfg["center"]
    cfg["ligand_sdf"] = str(sdf)
    assert checked_pocket(upstream, cfg)[1]["center"] == [0., 0., 0.]


@pytest.mark.parametrize("contents,match", [
    (pdb_atom(1, x=30.), "Empty pocket"),
    (pdb_atom(1) + pdb_atom(2, residue="UNK", number=2), "would skip.*UNK"),
    (pdb_atom(1, atom="N"), "C-alpha"),
    (pdb_atom(1) + pdb_atom(2), "C-alpha"),
    (pdb_atom(1, x=float("nan")), "Non-finite"),
])
def test_bad_pockets_fail(upstream, tmp_path, contents, match):
    with pytest.raises(ValueError, match=match):
        checked_pocket(upstream, pocket_cfg(tmp_path, contents))


def test_backend_load_once_device_gradients_and_provenance(upstream, monkeypatch, tmp_path):
    """Real Torch/graph/pocket code, tiny injected network; never loads PLANET weights."""
    import torch
    import PLANET_model
    from planet_oracle.backend import PlanetPredictor

    events = []

    class TinyPlanet(torch.nn.Module):
        def __init__(self, *args, device):
            super().__init__()
            events.append(("construct", args, str(device)))
            self.weight = torch.nn.Parameter(torch.zeros(1))

        def load_state_dict(self, state):
            events.append(("load", state))

        def cal_res_features_helper(self, features, coordinates):
            events.append(("precompute", torch.is_grad_enabled()))
            return features + self.weight

        def cal_res_features(self, features, size):
            return features, size

        def screening(self, features, size, graph):
            assert not torch.is_grad_enabled()
            return torch.full((size,), 7.)

    monkeypatch.setattr(PLANET_model, "PLANET", TinyPlanet)
    cfg = pocket_cfg(tmp_path, pdb_atom(1))
    checkpoint = tmp_path / "weights"
    torch.save({}, checkpoint)
    cfg.update(target_id="t", planet_root=str(PROJECT / "PLANET"), checkpoint=str(checkpoint),
               device="cpu", seed=42, batch_size=2, port=8765)
    predictor = PlanetPredictor(cfg)
    mol, _ = predictor.prepare("CCO")
    assert predictor.predict([mol, mol]) == [7., 7.]
    assert predictor.predict([mol]) == [7.]
    assert [event[0] for event in events] == ["construct", "load", "precompute"]
    assert events[0][1] == (300, 8, 300, 300, 3, 10, 1)
    assert events[2][1] is False
    assert not predictor.model.training and not predictor.model.weight.requires_grad
    assert predictor.device == "cpu" and len(predictor.oracle_id) == 64
    assert predictor.manifest["source_sha256"]["chemutils.py"]
    assert predictor.manifest["input_sha256"]["checkpoint"]
    again = PlanetPredictor(cfg)
    assert again.oracle_id == predictor.oracle_id
    cfg = {**cfg, "target_id": "other"}
    assert PlanetPredictor(cfg).oracle_id != predictor.oracle_id


@pytest.mark.parametrize("device", ["cuda", "cuda:0"])
def test_unavailable_cuda_has_no_fallback(upstream, monkeypatch, device):
    import torch
    from planet_oracle.backend import PlanetPredictor

    monkeypatch.setattr(torch.cuda, "is_available", lambda: False)
    with pytest.raises(ValueError, match="unavailable"):
        PlanetPredictor({"device": device})
