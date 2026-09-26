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

from box_utils import box_from_mol2

ROOT = Path(__file__).resolve().parents[2]
DATASET_DIR = ROOT / "real-world_dataset"
DOCKING_DIR = ROOT / "evaluation" / "docking"
CONFIG_DIR = DOCKING_DIR / "configs"
PREPARED_DIR = DOCKING_DIR / "prepared"

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


def write_config(target: str, config_dir: Path = CONFIG_DIR) -> Path:
    """Compute the grid box from the crystal ligand and freeze the config."""
    box = box_from_mol2(DATASET_DIR / target / "crystal.mol2", padding=BOX_PADDING_ANGSTROM)
    cfg = {
        "target": target,
        "receptor_pdb": f"real-world_dataset/{target}/receptor_out.pdb",
        "receptor_pdbqt": f"evaluation/docking/prepared/{target}/receptor.pdbqt",
        "crystal_mol2": f"real-world_dataset/{target}/crystal.mol2",
        "box": {
            "center": box["center"],
            "size": box["size"],
            "padding_angstrom": BOX_PADDING_ANGSTROM,
            "definition": "center = crystal ligand heavy-atom centroid; "
            "size = ligand bounding box + padding on each side",
        },
        "vina": dict(VINA_PARAMS),
    }
    config_dir.mkdir(parents=True, exist_ok=True)
    cfg_path = config_dir / f"{target}.json"
    cfg_path.write_text(json.dumps(cfg, indent=2) + "\n")
    return cfg_path


def prepare_receptor_pdbqt(receptor_pdb: Path, pdbqt_path: Path) -> None:
    """Convert a receptor PDB to PDBQT via Meeko's mk_prepare_receptor.py."""
    exe = shutil.which("mk_prepare_receptor.py")
    if exe is None:
        raise SystemExit(
            "mk_prepare_receptor.py not found on PATH. Activate the conda "
            "environment first (see evaluation/docking/README.md)."
        )
    pdbqt_path.parent.mkdir(parents=True, exist_ok=True)
    cmd = [exe, "--read_pdb", str(receptor_pdb), "-p", str(pdbqt_path)]
    result = subprocess.run(cmd, capture_output=True, text=True)
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

    for target in targets:
        cfg_path = write_config(target)
        print(f"[{target}] wrote {cfg_path.relative_to(ROOT)}", flush=True)
        if not args.box_only:
            pdbqt_path = PREPARED_DIR / target / "receptor.pdbqt"
            prepare_receptor_pdbqt(DATASET_DIR / target / "receptor_out.pdb", pdbqt_path)
            print(f"[{target}] wrote {pdbqt_path.relative_to(ROOT)}", flush=True)
    print(f"done: {len(targets)} target(s)", flush=True)


if __name__ == "__main__":
    main()
