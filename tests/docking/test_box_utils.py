"""Unit tests for evaluation/docking/box_utils.py (pure Python, no RDKit)."""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "evaluation" / "docking"))

from box_utils import box_from_coords, box_from_mol2, read_mol2_heavy_atom_coords

DATASET = ROOT / "real-world_dataset"
TARGETS = sorted(p.name for p in DATASET.iterdir() if (p / "crystal.mol2").is_file())


def test_all_sixteen_targets_present():
    assert len(TARGETS) == 16


@pytest.mark.parametrize("target", TARGETS)
def test_mol2_parses_heavy_atoms(target):
    coords = read_mol2_heavy_atom_coords(DATASET / target / "crystal.mol2")
    assert len(coords) >= 5  # real ligands are all larger than this
    for coord in coords:
        assert all(isinstance(v, float) for v in coord)


@pytest.mark.parametrize("target", TARGETS)
def test_box_geometry_matches_definition(target):
    path = DATASET / target / "crystal.mol2"
    coords = read_mol2_heavy_atom_coords(path)
    box = box_from_mol2(path, padding=8.0)
    for axis, center, size in zip(zip(*coords), box["center"], box["size"]):
        assert min(axis) <= center <= max(axis)
        expected = max(axis) - min(axis) + 16.0  # 8 A padding on EACH side
        assert abs(size - expected) <= 0.001 + 1e-9  # box is rounded to 3 decimals
        assert 10.0 <= size <= 80.0  # sane Vina box dimension


def test_hydrogens_are_excluded(tmp_path):
    mol2 = tmp_path / "toy.mol2"
    mol2.write_text(
        "@<TRIPOS>MOLECULE\ntoy\n 3 2 0 0 0\nSMALL\nGASTEIGER\n\n"
        "@<TRIPOS>ATOM\n"
        " 1 C 0.0 0.0 0.0 C.3 1 UNL1 0.0\n"
        " 2 H 1.0 0.0 0.0 H 1 UNL1 0.0\n"
        " 3 O 0.0 2.0 0.0 O.2 1 UNL1 -0.4\n"
        "@<TRIPOS>BOND\n 1 1 2 1\n 2 1 3 2\n"
    )
    coords = read_mol2_heavy_atom_coords(mol2)
    assert coords == [(0.0, 0.0, 0.0), (0.0, 2.0, 0.0)]


def test_box_from_coords_padding_and_centroid():
    coords = [(0.0, 0.0, 0.0), (2.0, 4.0, 6.0)]
    box = box_from_coords(coords, padding=8.0)
    assert box["center"] == [1.0, 2.0, 3.0]
    assert box["size"] == [18.0, 20.0, 22.0]


def test_empty_atom_section_raises(tmp_path):
    bad = tmp_path / "bad.mol2"
    bad.write_text("@<TRIPOS>MOLECULE\nempty\n 0 0 0 0 0\nSMALL\n\n@<TRIPOS>ATOM\n@<TRIPOS>BOND\n")
    with pytest.raises(ValueError, match="no heavy atoms"):
        read_mol2_heavy_atom_coords(bad)
