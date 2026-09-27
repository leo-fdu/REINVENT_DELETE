#!/usr/bin/env python3
"""Prepare receptors and Vina grid boxes for all real-world targets (one-time).

For every target directory in ``real-world_dataset/`` that contains both
``receptor_out.pdb`` and ``crystal.mol2`` this script:

1. converts ``receptor_out.pdb`` to
   ``evaluation/docking/prepared/<target>/receptor.pdbqt`` using Meeko's
   ``mk_prepare_receptor.py`` (must be on PATH, see README.md), and
2. derives the docking grid box from ``crystal.mol2`` and freezes it, together
   with the shared Vina parameters, into
   ``evaluation/docking/configs/<target>.json``.

The JSON configs are versioned (they are the fairness contract shared by both
models); the prepared PDBQTs are derived files that can be regenerated at any
time. The script is deterministic and safe to re-run.

For JAK2 and DRD3 the original receptors/pockets were found to be wrong (see
``configs/manual_planet_rl/README.md``). When
``configs/manual_planet_rl/target_inputs.json`` is present, those targets are
prepared from the same traceable corrections as the PLANET RL workflow: the
receptor is the extracted author chain written to
``evaluation/docking/prepared/<target>/receptor.pdb``, and the grid box is
derived from the pocket-reference crystallographic ligand instead of the
design ``crystal.mol2``. All other targets are unchanged.

Receptor sanitization (uniform for all targets): several source receptors
contain residues with truncated side chains (e.g. CDK2's G-loop LYS A:9 keeps
only backbone+CB) and alternate locations. Meeko cannot template-match
truncated residues, and dropping them would punch holes in the binding site.
``sanitize_receptor_pdb`` therefore (1) resolves alternate locations with the
same policy as the PLANET RL receptor extraction (highest mean occupancy,
ties -> lexicographic label), and (2) truncates incomplete standard side
chains to backbone+CB and relabels them ALA. Coordinates of kept atoms never
change. The sanitized copy ``prepared/<target>/receptor_sanitized.pdb`` is
what Meeko actually reads; all decisions are logged to
``prepared/<target>/receptor_sanitized.log``.

Usage:
    python evaluation/docking/prepare_receptors.py                 # all targets
    python evaluation/docking/prepare_receptors.py --targets adrb1 cdk2
    python evaluation/docking/prepare_receptors.py --box-only      # skip PDBQT conversion
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

from box_utils import box_from_coords, box_from_mol2

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from planet_target_inputs import (  # noqa: E402
    corrected_inputs,
    reference_coords,
    select_alternates,
)

DATASET_DIR = ROOT / "real-world_dataset"
DOCKING_DIR = ROOT / "evaluation" / "docking"
CONFIG_DIR = DOCKING_DIR / "configs"
PREPARED_DIR = DOCKING_DIR / "prepared"
TARGET_INPUTS = ROOT / "configs" / "manual_planet_rl" / "target_inputs.json"

# Frozen experiment settings (shared by all targets and both models).
BOX_PADDING_ANGSTROM = 8.0  # padding added on EACH side of the ligand bounding box
VINA_PARAMS = {"scoring": "vina", "exhaustiveness": 16, "num_modes": 1, "seed": 42}


def discover_targets(dataset_dir: Path = DATASET_DIR) -> list[str]:
    """Return sorted names of target dirs that have both required raw files."""
    return [
        p.name
        for p in sorted(dataset_dir.iterdir())
        if p.is_dir() and (p / "receptor_out.pdb").is_file() and (p / "crystal.mol2").is_file()
    ]


def _project_relative(path: str) -> str:
    return str(Path(path).resolve().relative_to(ROOT))


def load_corrections(target_inputs: Path = TARGET_INPUTS) -> dict:
    """Load traceable receptor/pocket corrections, validating source hashes."""
    if not target_inputs.is_file():
        return {}
    corrections, _ = corrected_inputs(ROOT, target_inputs)
    return corrections


def write_config(target: str, correction: dict | None = None,
                 config_dir: Path = CONFIG_DIR) -> Path:
    """Compute the grid box and freeze the config.

    Uncorrected targets use the design ``crystal.mol2``; corrected targets use
    the pocket-reference crystallographic ligand and the extracted author-chain
    receptor, identical to the PLANET RL workflow.
    """
    definition = ("center = crystal ligand heavy-atom centroid; "
                  "size = ligand bounding box + padding on each side")
    receptor_pdb = f"real-world_dataset/{target}/receptor_out.pdb"
    extra: dict = {}
    if correction is None:
        box = box_from_mol2(DATASET_DIR / target / "crystal.mol2", padding=BOX_PADDING_ANGSTROM)
    else:
        reference = correction["pocket_reference"]
        coords = reference_coords(reference["source_pdb"], reference["chain"],
                                  reference["resname"], reference["residue_id"])
        box = box_from_coords(coords, padding=BOX_PADDING_ANGSTROM)
        frozen_center = [round(v, 3) for v in correction["center"]]
        if box["center"] != frozen_center:
            raise ValueError(f"Box center {box['center']} != frozen pocket center "
                             f"{frozen_center} for {target}")
        receptor_pdb = f"evaluation/docking/prepared/{target}/receptor.pdb"
        definition = ("center = pocket-reference ligand heavy-atom centroid "
                      f"({reference['resname']} {reference['chain']}/{reference['residue_id']}); "
                      "size = reference ligand bounding box + padding on each side")
        provenance = {"receptor": correction["receptor"],
                      "pocket_reference": correction["pocket_reference"],
                      "reason": correction["reason"]}
        extra = {"input_correction": {
            key: ({**value, "source_pdb": _project_relative(value["source_pdb"])}
                  if isinstance(value, dict) else value)
            for key, value in provenance.items()}}
    cfg = {
        "target": target,
        "receptor_pdb": receptor_pdb,
        "receptor_pdbqt": f"evaluation/docking/prepared/{target}/receptor.pdbqt",
        "crystal_mol2": f"real-world_dataset/{target}/crystal.mol2",
        "box": {
            "center": box["center"],
            "size": box["size"],
            "padding_angstrom": BOX_PADDING_ANGSTROM,
            "definition": definition,
        },
        "vina": dict(VINA_PARAMS),
        **extra,
    }
    config_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = config_dir / f"{target}.json"
    cfg_path.write_text(json.dumps(cfg, indent=2) + "\n")
    return cfg_path


BACKBONE = {"N", "CA", "C", "O"}
BACKBONE_CB = BACKBONE | {"CB"}

# Covalent heavy-atom bonds within each standard side chain (CA-CB included).
SIDECHAIN_BONDS = {
    "ALA": [("CA", "CB")],
    "ARG": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"),
            ("NE", "CZ"), ("CZ", "NH1"), ("CZ", "NH2")],
    "ASN": [("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "ND2")],
    "ASP": [("CA", "CB"), ("CB", "CG"), ("CG", "OD1"), ("CG", "OD2")],
    "CYS": [("CA", "CB"), ("CB", "SG")],
    "GLN": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "NE2")],
    "GLU": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "OE1"), ("CD", "OE2")],
    "GLY": [],
    "HIS": [("CA", "CB"), ("CB", "CG"), ("CG", "ND1"), ("CG", "CD2"),
            ("ND1", "CE1"), ("CD2", "NE2"), ("CE1", "NE2")],
    "ILE": [("CA", "CB"), ("CB", "CG1"), ("CB", "CG2"), ("CG1", "CD1")],
    "LEU": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2")],
    "LYS": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE"), ("CE", "NZ")],
    "MET": [("CA", "CB"), ("CB", "CG"), ("CG", "SD"), ("SD", "CE")],
    "PHE": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
            ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ")],
    "PRO": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("N", "CD")],
    "SER": [("CA", "CB"), ("CB", "OG")],
    "THR": [("CA", "CB"), ("CB", "OG1"), ("CB", "CG2")],
    "TRP": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
            ("CD1", "NE1"), ("NE1", "CE2"), ("CD2", "CE2"), ("CD2", "CE3"),
            ("CE2", "CZ2"), ("CE3", "CZ3"), ("CZ2", "CH2"), ("CZ3", "CH2")],
    "TYR": [("CA", "CB"), ("CB", "CG"), ("CG", "CD1"), ("CG", "CD2"),
            ("CD1", "CE1"), ("CD2", "CE2"), ("CE1", "CZ"), ("CE2", "CZ"), ("CZ", "OH")],
    "VAL": [("CA", "CB"), ("CB", "CG1"), ("CB", "CG2")],
}
SIDECHAIN_ATOMS = {
    resname: {atom for bond in bonds for atom in bond if atom not in BACKBONE}
    for resname, bonds in SIDECHAIN_BONDS.items()
}
# Accepted covalent bond-length window (generous; real heavy bonds are 1.2-1.6 A).
BOND_MIN, BOND_MAX = 1.0, 2.0


def sanitize_receptor_pdb(src_pdb: Path, dst_pdb: Path, log_path: Path) -> list:
    """Resolve altlocs and truncate incomplete side chains to ALA for Meeko.

    Steps (coordinates of kept atoms never change):

    1. Alternate locations are resolved with the same policy as the PLANET RL
       receptor extraction (``planet_target_inputs.select_alternates``):
       highest mean occupancy per residue, ties -> lexicographic label.
    2. Any standard residue whose side chain is incomplete is truncated to
       backbone+CB and relabeled ALA; residues with non-standard atom names
       or missing backbone atoms are left for Meeko to report.
    3. Any complete side chain with a broken bond (outside
       [BOND_MIN, BOND_MAX] A, e.g. scrambled altloc coordinates mixing both
       conformations) is likewise truncated: broken CA-CB drops the CB as
       well and relabels GLY, otherwise backbone+CB is kept (ALA).

    Returns the truncation list (``chain:resseq RES->ALA`` or ``RES->GLY``).
    """
    lines = Path(src_pdb).read_text(encoding="utf-8").splitlines()
    atom_lines = [l for l in lines if l.startswith(("ATOM  ", "HETATM"))]
    selected = select_alternates(atom_lines)
    keep_serials = {line[6:11] for line in selected}
    lines = [l for l in lines
             if not l.startswith(("ATOM  ", "HETATM")) or l[6:11] in keep_serials]
    altloc_removed = len(atom_lines) - len(selected)

    residues: dict = {}
    for line in lines:
        if line.startswith("ATOM  "):
            residues.setdefault((line[21], line[22:27]), []).append(line)

    def _xyz(atom_line):
        return [float(atom_line[30:38]), float(atom_line[38:46]), float(atom_line[46:54])]

    truncate: dict = {}
    reasons: dict = {}
    for key, atoms in sorted(residues.items()):
        resname = atoms[0][17:20].strip()
        side = SIDECHAIN_ATOMS.get(resname)
        names = {a[12:16].strip() for a in atoms}
        if side is None or not BACKBONE <= names \
                or not names <= BACKBONE | side | {"OXT"}:
            continue  # non-standard chemistry: leave for Meeko to report
        if side - names:
            truncate[key] = "ALA"
            reasons[key] = f"incomplete side chain (missing: {' '.join(sorted(side - names))})"
            continue
        coords = {a[12:16].strip(): _xyz(a) for a in atoms if a[12:16].strip() != "OXT"}
        broken = [(a, b, math.dist(coords[a], coords[b]))
                  for a, b in SIDECHAIN_BONDS[resname]
                  if not BOND_MIN <= math.dist(coords[a], coords[b]) <= BOND_MAX]
        if broken:
            a, b, dist = broken[0]
            if (a, b) == ("CA", "CB"):
                truncate[key] = "GLY"
                reasons[key] = f"broken CA-CB bond ({dist:.2f} A); backbone kept"
            else:
                truncate[key] = "ALA"
                reasons[key] = f"broken side-chain bond {a}-{b} ({dist:.2f} A)"

    out = []
    for line in lines:
        if line.startswith("ATOM  ") and (line[21], line[22:27]) in truncate:
            target_resname = truncate[(line[21], line[22:27])]
            keep_names = BACKBONE_CB if target_resname == "ALA" else BACKBONE
            if line[12:16].strip() not in keep_names:
                continue  # drop the untypable side chain (beyond CB)
            line = line[:17] + target_resname + line[20:]
        out.append(line)
    dst_pdb.parent.mkdir(parents=True, exist_ok=True)
    dst_pdb.write_text("\n".join(out) + "\n", encoding="utf-8")

    report = [
        f"alternate locations resolved: kept {len(selected)} of {len(atom_lines)} "
        f"atom records ({altloc_removed} dropped)",
        "untypable side chains truncated (coordinates of kept atoms unchanged):",
    ]
    report += [f"{chain}:{resseq.strip()} {resname}->{truncate[(chain, resseq)]} "
               f"({reasons[(chain, resseq)]})"
               for (chain, resseq), resname in
               sorted(((k, residues[k][0][17:20].strip()) for k in truncate),
                      key=lambda item: item[0])] or ["(none)"]
    log_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    return [f"{chain}:{resseq.strip()} {resname}->{truncate[(chain, resseq)]}"
            for (chain, resseq), resname in
            sorted(((k, residues[k][0][17:20].strip()) for k in truncate),
                   key=lambda item: item[0])]


def prepare_corrected_receptor_pdb(target: str, correction: dict,
                                   prepared_dir: Path = PREPARED_DIR) -> Path:
    """Write the extracted author-chain receptor PDB for a corrected target."""
    pdb_path = prepared_dir / target / "receptor.pdb"
    pdb_path.parent.mkdir(parents=True, exist_ok=True)
    pdb_path.write_text(correction["receptor_text"], encoding="utf-8")
    return pdb_path


def prepare_sanitized_receptor(target: str, receptor_pdb: Path,
                               prepared_dir: Path = PREPARED_DIR) -> Path:
    """Write the sanitized receptor copy that Meeko reads; log relabelings."""
    sanitized = prepared_dir / target / "receptor_sanitized.pdb"
    log_path = prepared_dir / target / "receptor_sanitized.log"
    relabeled = sanitize_receptor_pdb(receptor_pdb, sanitized, log_path)
    if relabeled:
        print(f"[{target}] relabeled {len(relabeled)} truncated residue(s) to ALA"
              f" (see {log_path.relative_to(ROOT)})", flush=True)
    return sanitized


def prepare_receptor_pdbqt(receptor_pdb: Path, pdbqt_path: Path) -> None:
    """Convert a receptor PDB to PDBQT via Meeko's mk_prepare_receptor.py.

    ``--default_altloc A`` deterministically resolves alternate locations
    (present in several original receptors). ``--allow_bad_res`` is kept as a
    last-resort guard; any residue it drops is recorded in
    ``receptor.meeko.log`` and must be checked to be far from the grid box
    (truncated side chains are preserved via ALA relabeling instead, see
    ``sanitize_receptor_pdb``).
    """
    exe = shutil.which("mk_prepare_receptor.py")
    if exe is None:
        raise SystemExit(
            "mk_prepare_receptor.py not found on PATH. Activate the conda "
            "environment first (see evaluation/docking/README.md)."
        )
    pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "--read_pdb", str(receptor_pdb), "-p", str(pdbqt_path),
           "--default_altloc", "A", "--allow_bad_res"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    pdbqt_path.with_suffix(".meeko.log").write_text(
        result.stdout + result.stderr, encoding="utf-8")
    if result.returncode != 0 or not pdbqt_path.is_file():
        sys.stderr.write(result.stdout)
        sys.stderr.write(result.stderr)
        raise SystemExit(f"receptor preparation failed for {receptor_pdb}")


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--targets", nargs="*", default=None,
                        help="subset of targets (default: all discovered)")
    parser.add_argument("--box-only", action="store_true",
                        help="only (re)write the box configs, skip receptor PDBQT conversion")
    args = parser.parse_args(argv)

    targets = args.targets or discover_targets()
    if not targets:
        raise SystemExit(f"no targets found under {DATASET_DIR}")

    corrections = load_corrections()
    unknown = set(corrections) - set(discover_targets())
    if unknown:
        raise SystemExit(f"correction target(s) not in dataset: {sorted(unknown)}")
    for target in targets:
        correction = corrections.get(target)
        cfg_path = write_config(target, correction)
        print(f"[{target}] wrote {cfg_path.relative_to(ROOT)}", flush=True)
        if not args.box_only:
            if correction is not None:
                receptor_pdb = prepare_corrected_receptor_pdb(target, correction)
                print(f"[{target}] wrote {receptor_pdb.relative_to(ROOT)} (corrected)", flush=True)
            else:
                receptor_pdb = DATASET_DIR / target / "receptor_out.pdb"
            sanitized = prepare_sanitized_receptor(target, receptor_pdb)
            pdbqt_path = PREPARED_DIR / target / "receptor.pdbqt"
            prepare_receptor_pdbqt(sanitized, pdbqt_path)
            print(f"[{target}] wrote {pdbqt_path.relative_to(ROOT)}", flush=True)
    print(f"done: {len(targets)} target(s)", flush=True)


if __name__ == "__main__":
    main()
