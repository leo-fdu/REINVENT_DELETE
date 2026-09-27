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

# Standard amino-acid side-chain heavy atoms (beyond the N, CA, C, O backbone).
SIDECHAIN_ATOMS = {
    "ALA": {"CB"},
    "ARG": {"CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"},
    "ASN": {"CB", "CG", "OD1", "ND2"},
    "ASP": {"CB", "CG", "OD1", "OD2"},
    "CYS": {"CB", "SG"},
    "GLN": {"CB", "CG", "CD", "OE1", "NE2"},
    "GLU": {"CB", "CG", "CD", "OE1", "OE2"},
    "GLY": set(),
    "HIS": {"CB", "CG", "ND1", "CD2", "CE1", "NE2"},
    "ILE": {"CB", "CG1", "CG2", "CD1"},
    "LEU": {"CB", "CG", "CD1", "CD2"},
    "LYS": {"CB", "CG", "CD", "CE", "NZ"},
    "MET": {"CB", "CG", "SD", "CE"},
    "PHE": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ"},
    "PRO": {"CB", "CG", "CD"},
    "SER": {"CB", "OG"},
    "THR": {"CB", "OG1", "CG2"},
    "TRP": {"CB", "CG", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"},
    "TYR": {"CB", "CG", "CD1", "CD2", "CE1", "CE2", "CZ", "OH"},
    "VAL": {"CB", "CG1", "CG2"},
}


def sanitize_receptor_pdb(src_pdb: Path, dst_pdb: Path, log_path: Path) -> list:
    """Resolve altlocs and truncate incomplete side chains to ALA for Meeko.

    Steps (coordinates of kept atoms never change):

    1. Alternate locations are resolved with the same policy as the PLANET RL
       receptor extraction (``planet_target_inputs.select_alternates``):
       highest mean occupancy per residue, ties -> lexicographic label.
    2. Any standard residue whose side chain is incomplete is truncated to
       backbone+CB and relabeled ALA; residues with non-standard atom names
       or missing backbone atoms are left for Meeko to report.

    Returns the truncation list (``chain:resseq RES->ALA (dropped N atoms)``).
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
    truncate: dict = {}
    for key, atoms in sorted(residues.items()):
        resname = atoms[0][17:20].strip()
        side = SIDECHAIN_ATOMS.get(resname)
        names = {a[12:16].strip() for a in atoms}
        if side is None or not BACKBONE <= names \
                or not names <= BACKBONE | side | {"OXT"}:
            continue  # non-standard chemistry: leave for Meeko to report
        if side - names:
            truncate[key] = resname, len(atoms) - sum(
                a[12:16].strip() in BACKBONE_CB for a in atoms)

    out = []
    for line in lines:
        if line.startswith("ATOM  ") and (line[21], line[22:27]) in truncate:
            if line[12:16].strip() not in BACKBONE_CB:
                continue  # drop the incomplete side chain beyond CB
            line = line[:17] + "ALA" + line[20:]
        out.append(line)
    dst_pdb.parent.mkdir(parents=True, exist_ok=True)
    dst_pdb.write_text("\n".join(out) + "\n", encoding="utf-8")

    report = [
        f"alternate locations resolved: kept {len(selected)} of {len(atom_lines)} "
        f"atom records ({altloc_removed} dropped)",
        "incomplete side chains truncated to backbone+CB and relabeled ALA "
        "(coordinates of kept atoms unchanged):",
    ]
    report += [f"{chain}:{resseq.strip()} {resname}->ALA (dropped {dropped} atoms)"
               for (chain, resseq), (resname, dropped) in truncate.items()] or ["(none)"]
    log_path.write_text("\n".join(report) + "\n", encoding="utf-8")
    return [f"{chain}:{resseq.strip()} {truncate[(chain, resseq)][0]}->ALA"
            for chain, resseq in truncate]


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
