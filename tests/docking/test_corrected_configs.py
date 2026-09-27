"""Tests for corrected JAK2/DRD3 receptor + box preparation (traceable inputs).

These tests use the real versioned inputs (official source PDBs and
configs/manual_planet_rl/target_inputs.json) and need no RDKit/Meeko/Vina.
"""

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "evaluation" / "docking"))

from prepare_receptors import (  # noqa: E402
    discover_targets,
    load_corrections,
    prepare_corrected_receptor_pdb,
    sanitize_receptor_pdb,
    write_config,
)

TARGET_INPUTS = ROOT / "configs" / "manual_planet_rl" / "target_inputs.json"

# Values frozen in configs/manual_planet_rl/README.md (traceable correction).
EXPECTED = {
    "jak2": {
        "center": [114.671, 66.114, 10.555],
        "size": [25.033, 23.929, 26.283],
        "resname": "NVB",
        "residue_id": "1133",
    },
    "drd3": {
        "center": [0.085, -14.828, 10.432],
        "size": [20.796, 25.596, 22.782],
        "resname": "ETQ",
        "residue_id": "1200",
    },
}


@pytest.fixture(scope="module")
def corrections():
    return load_corrections(TARGET_INPUTS)


def test_corrections_cover_exactly_jak2_drd3(corrections):
    assert set(corrections) == {"jak2", "drd3"}
    assert set(corrections) <= set(discover_targets())


@pytest.mark.parametrize("target", sorted(EXPECTED))
def test_corrected_box_matches_frozen_values(target, corrections, tmp_path):
    cfg_path = write_config(target, corrections[target], config_dir=tmp_path)
    cfg = json.loads(cfg_path.read_text())
    expected = EXPECTED[target]
    assert cfg["box"]["center"] == expected["center"]
    assert cfg["box"]["size"] == expected["size"]
    assert cfg["box"]["padding_angstrom"] == 8.0
    assert expected["resname"] in cfg["box"]["definition"]
    assert cfg["receptor_pdb"] == f"evaluation/docking/prepared/{target}/receptor.pdb"
    assert cfg["vina"] == {"scoring": "vina", "exhaustiveness": 16,
                           "num_modes": 1, "seed": 42}


@pytest.mark.parametrize("target", sorted(EXPECTED))
def test_corrected_config_provenance(target, corrections, tmp_path):
    cfg_path = write_config(target, corrections[target], config_dir=tmp_path)
    cfg = json.loads(cfg_path.read_text())
    correction = cfg["input_correction"]
    reference = correction["pocket_reference"]
    frozen = json.loads(TARGET_INPUTS.read_text())["targets"][target]
    assert correction["reason"] == frozen["reason"]
    assert reference["resname"] == EXPECTED[target]["resname"]
    assert reference["residue_id"] == EXPECTED[target]["residue_id"]
    # Paths in the versioned config stay project-relative.
    assert not Path(reference["source_pdb"]).is_absolute()
    assert Path(reference["source_pdb"]).is_relative_to("real-world_dataset")
    assert len(reference["source_sha256"]) == 64
    assert correction["receptor"]["source_sha256"] == frozen["receptor"]["source_sha256"]


@pytest.mark.parametrize("target", sorted(EXPECTED))
def test_corrected_receptor_pdb_single_chain(target, corrections, tmp_path):
    pdb_path = prepare_corrected_receptor_pdb(target, corrections[target], prepared_dir=tmp_path)
    lines = pdb_path.read_text().splitlines()
    atom_lines = [l for l in lines if l.startswith(("ATOM  ", "HETATM"))]
    assert atom_lines, "derived receptor must contain atoms"
    assert {l[21] for l in atom_lines} == {"A"}
    assert any(l.startswith("ATOM  ") for l in atom_lines)


def test_uncorrected_targets_have_no_correction(tmp_path):
    cfg_path = write_config("cdk2", None, config_dir=tmp_path)
    cfg = json.loads(cfg_path.read_text())
    assert "input_correction" not in cfg
    assert cfg["receptor_pdb"] == "real-world_dataset/cdk2/receptor_out.pdb"
    assert cfg["box"]["center"] == [1.97, 27.559, 8.825]  # design crystal.mol2


def _atom(serial, name, resname, resseq):
    return (f"ATOM  {serial:>5} {name:^4} {resname:>3} A{resseq:>5}   "
            "  0.000   0.000   0.000  1.00  0.00           "
            + name.strip()[0].rjust(2) + "  ")


def test_sanitize_relabels_only_backbone_cb_residues(tmp_path):
    src = tmp_path / "in.pdb"
    atoms = [
        _atom(1, "N", "LYS", 9), _atom(2, "CA", "LYS", 9), _atom(3, "C", "LYS", 9),
        _atom(4, "O", "LYS", 9), _atom(5, "CB", "LYS", 9),          # truncated LYS
        _atom(6, "N", "LEU", 10), _atom(7, "CA", "LEU", 10), _atom(8, "C", "LEU", 10),
        _atom(9, "O", "LEU", 10), _atom(10, "CB", "LEU", 10), _atom(11, "CG", "LEU", 10),
        _atom(12, "CD1", "LEU", 10), _atom(13, "CD2", "LEU", 10),   # complete LEU
        _atom(14, "N", "GLY", 11), _atom(15, "CA", "GLY", 11), _atom(16, "C", "GLY", 11),
        _atom(17, "O", "GLY", 11),                                   # backbone-only GLY
        _atom(18, "N", "ALA", 12), _atom(19, "CA", "ALA", 12), _atom(20, "C", "ALA", 12),
        _atom(21, "O", "ALA", 12), _atom(22, "CB", "ALA", 12),       # genuine ALA
    ]
    src.write_text("\n".join(atoms) + "\nTER\nEND\n")
    dst = tmp_path / "out.pdb"
    relabeled = sanitize_receptor_pdb(src, dst, tmp_path / "out.log")
    assert relabeled == ["A:9 LYS->ALA"]
    lines = dst.read_text().splitlines()
    assert [l[17:20] for l in lines if l.startswith("ATOM")][:5] == ["ALA"] * 5
    assert {l[17:20] for l in lines if l.startswith("ATOM")} == {"ALA", "LEU", "GLY"}
    # Coordinates/atom names untouched; only the resname field changed.
    for before, after in zip(src.read_text().splitlines(), lines):
        if before.startswith("ATOM") and before[17:20] == "LYS":
            assert before[:17] == after[:17] and before[20:] == after[20:]
        else:
            assert before == after


def test_sanitize_noop_on_complete_receptor(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text(_atom(1, "N", "LEU", 1) + "\n" + _atom(2, "CA", "LEU", 1) + "\n"
                   + _atom(3, "C", "LEU", 1) + "\n" + _atom(4, "O", "LEU", 1) + "\n"
                   + _atom(5, "CB", "LEU", 1) + "\n" + _atom(6, "CG", "LEU", 1) + "\n"
                   + _atom(7, "CD1", "LEU", 1) + "\n" + _atom(8, "CD2", "LEU", 1) + "\n")
    dst = tmp_path / "out.pdb"
    relabeled = sanitize_receptor_pdb(src, dst, tmp_path / "out.log")
    assert relabeled == []
    assert dst.read_text() == src.read_text()
