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


def _atom(serial, name, resname, resseq, xyz=(0.0, 0.0, 0.0), altloc=" ", occ=1.0):
    return (f"ATOM  {serial:5d} {name:^4}{altloc}{resname:>3} A{resseq:4d}    "
            f"{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{occ:6.2f}{0.0:6.2f}          "
            f"{name.strip()[0]:>2}  ")


# Sane toy geometries: every side-chain bond is 1.45 A (within [1.0, 2.0]).
LEU_COORDS = {"N": (0, 0, 0), "CA": (1.45, 0, 0), "C": (2.9, 0, 0), "O": (4.35, 0, 0),
              "CB": (1.45, 1.45, 0), "CG": (2.9, 1.45, 0),
              "CD1": (4.35, 1.45, 0), "CD2": (2.9, 2.9, 0)}
ARG_COORDS = {"N": (0, 0, 0), "CA": (1.45, 0, 0), "C": (2.9, 0, 0), "O": (4.35, 0, 0),
              "CB": (1.45, 1.45, 0), "CG": (2.9, 1.45, 0), "CD": (4.35, 1.45, 0),
              "NE": (5.8, 1.45, 0), "CZ": (7.25, 1.45, 0),
              "NH1": (8.7, 1.45, 0), "NH2": (7.25, 2.9, 0)}
BACKBONE_CB_COORDS = {"N": (0, 0, 0), "CA": (1.45, 0, 0), "C": (2.9, 0, 0),
                      "O": (4.35, 0, 0), "CB": (1.45, 1.45, 0)}


def _residue(resname, resseq, coords, start_serial=1, altloc=" ", occ=1.0):
    return [_atom(start_serial + i, name, resname, resseq, xyz=xyz,
                  altloc=altloc, occ=occ)
            for i, (name, xyz) in enumerate(coords.items())]


def test_sanitize_relabels_only_backbone_cb_residues(tmp_path):
    src = tmp_path / "in.pdb"
    atoms = (_residue("LYS", 9, BACKBONE_CB_COORDS)          # truncated LYS
             + _residue("LEU", 10, LEU_COORDS, 10)          # complete LEU
             + _residue("GLY", 11, {k: v for k, v in BACKBONE_CB_COORDS.items()
                                    if k != "CB"}, 20)      # backbone-only GLY
             + _residue("ALA", 12, BACKBONE_CB_COORDS, 30))  # genuine ALA
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


def test_sanitize_truncates_partial_side_chains(tmp_path):
    src = tmp_path / "in.pdb"
    atoms = _residue("LYS", 857, {k: v for k, v in ARG_COORDS.items()
                                  if k in ("N", "CA", "C", "O", "CB", "CG", "CD")})
    src.write_text("\n".join(atoms) + "\nEND\n")
    dst = tmp_path / "out.pdb"
    relabeled = sanitize_receptor_pdb(src, dst, tmp_path / "out.log")
    assert relabeled == ["A:857 LYS->ALA"]
    lines = [l for l in dst.read_text().splitlines() if l.startswith("ATOM")]
    assert [l[12:16].strip() for l in lines] == ["N", "CA", "C", "O", "CB"]
    assert {l[17:20] for l in lines} == {"ALA"}


def test_sanitize_flags_broken_side_chain_geometry(tmp_path):
    src = tmp_path / "in.pdb"
    broken = dict(ARG_COORDS)
    broken["NH1"] = (7.25, 5.0, 0)  # CZ-NH1 = 3.55 A (scrambled altloc coords)
    atoms = _residue("ARG", 88, broken)
    src.write_text("\n".join(atoms) + "\nEND\n")
    dst = tmp_path / "out.pdb"
    relabeled = sanitize_receptor_pdb(src, dst, tmp_path / "out.log")
    assert relabeled == ["A:88 ARG->ALA"]
    log = (tmp_path / "out.log").read_text()
    assert "broken side-chain bond CZ-NH1" in log
    lines = [l for l in dst.read_text().splitlines() if l.startswith("ATOM")]
    assert [l[12:16].strip() for l in lines] == ["N", "CA", "C", "O", "CB"]


def test_sanitize_broken_ca_cb_becomes_gly(tmp_path):
    src = tmp_path / "in.pdb"
    coords = dict(BACKBONE_CB_COORDS)
    coords["CB"] = (1.45, 4.5, 0)  # CA-CB = 4.5 A
    atoms = _residue("ALA", 536, coords)
    src.write_text("\n".join(atoms) + "\nEND\n")
    dst = tmp_path / "out.pdb"
    relabeled = sanitize_receptor_pdb(src, dst, tmp_path / "out.log")
    assert relabeled == ["A:536 ALA->GLY"]
    lines = [l for l in dst.read_text().splitlines() if l.startswith("ATOM")]
    assert [l[12:16].strip() for l in lines] == ["N", "CA", "C", "O"]
    assert {l[17:20] for l in lines} == {"GLY"}


def test_sanitize_resolves_alternate_locations(tmp_path):
    src = tmp_path / "in.pdb"
    atoms = []
    for label, occ in (("A", 0.5), ("B", 0.5)):
        atoms += _residue("ARG", 88, ARG_COORDS, start_serial=len(atoms) + 1,
                          altloc=label, occ=occ)
    src.write_text("\n".join(atoms) + "\nEND\n")
    dst = tmp_path / "out.pdb"
    relabeled = sanitize_receptor_pdb(src, dst, tmp_path / "out.log")
    assert relabeled == []  # complete after picking one conformation
    lines = [l for l in dst.read_text().splitlines() if l.startswith("ATOM")]
    assert len(lines) == 11
    assert {l[16] for l in lines} == {"A"}  # tie broken lexicographically
    assert {l[17:20] for l in lines} == {"ARG"}
    assert "kept 11 of 22" in (tmp_path / "out.log").read_text()


def test_sanitize_noop_on_complete_receptor(tmp_path):
    src = tmp_path / "in.pdb"
    src.write_text("\n".join(_residue("LEU", 1, LEU_COORDS)) + "\n")
    dst = tmp_path / "out.pdb"
    relabeled = sanitize_receptor_pdb(src, dst, tmp_path / "out.log")
    assert relabeled == []
    assert dst.read_text() == src.read_text()
