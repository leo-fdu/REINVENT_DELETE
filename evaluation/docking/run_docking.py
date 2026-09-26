#!/usr/bin/env python3
"""Dock generated molecules into the prepared receptors with AutoDock Vina.

Input
-----
CSV file with at least the columns ``target``, ``method`` and ``smiles``
(extra columns are ignored). ``target`` must be one of the prepared targets
(``evaluation/docking/configs/<target>.json``; run ``prepare_receptors.py``
first) and ``method`` must be ``reinvent`` or ``delete``. The SMILES strings
must already be canonicalised by the upstream adapter; this module never
rewrites them in the output (internally only the largest fragment is kept
before docking, e.g. to drop counter-ions).

Output
------
``--out`` CSV with exactly five columns:

    molecule_id,method,target,smiles,vina_score

``molecule_id`` is a 1-based sequential number within each (target, method)
group, following input row order. ``vina_score`` is the best docked pose
energy in kcal/mol (lower = better); it is left empty for molecules that
failed anywhere in the pipeline, with the failure reason recorded in
``<out_stem>_errors.csv`` next to the output file. A small
``<out_stem>_summary.csv`` with per-target/method statistics is also written.

Determinism / fairness
----------------------
All targets and both methods share the same frozen settings (stored in the
per-target JSON configs): RDKit ETKDGv3 embedding with randomSeed=42 (fallback
to random coordinates), MMFF minimisation (fallback UFF), Meeko ligand
preparation, and Vina with scoring=vina, exhaustiveness=16, num_modes=1,
seed=42. Each docking worker runs Vina single-threaded (cpu=1), so the scores
do not depend on ``--n-jobs``.

Usage:
    python evaluation/docking/run_docking.py \
        --input generated_smiles.csv \
        --out evaluation/docking/results/run1/docking_scores.csv \
        --n-jobs 8
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_CONFIG_DIR = Path(__file__).resolve().parent / "configs"

REQUIRED_COLUMNS = ("target", "method", "smiles")
VALID_METHODS = ("reinvent", "delete")
EMBEDDING_SEED = 42

# Worker-process globals, initialised once per process by _init_worker().
_VINA = None
_LIGAND_PREP = None
_EXHAUSTIVENESS = None
_NUM_MODES = None


# --------------------------------------------------------------------------- #
# Input handling
# --------------------------------------------------------------------------- #
def read_input(input_csv: Path) -> list[dict]:
    """Read and validate the input CSV; assign per-(target, method) molecule ids."""
    with Path(input_csv).open(newline="") as handle:
        reader = csv.DictReader(handle)
        missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise SystemExit(
                f"input {input_csv} is missing required column(s): {', '.join(missing)} "
                f"(required: {', '.join(REQUIRED_COLUMNS)})"
            )
        rows = list(reader)
    if not rows:
        raise SystemExit(f"input {input_csv} contains no rows")

    counters: dict[tuple[str, str], int] = {}
    tasks: list[dict] = []
    for lineno, row in enumerate(rows, start=2):  # line 1 is the header
        target = (row.get("target") or "").strip()
        method = (row.get("method") or "").strip()
        smiles = (row.get("smiles") or "").strip()
        if method not in VALID_METHODS:
            raise SystemExit(f"row {lineno}: method must be one of {VALID_METHODS}, got {method!r}")
        if not target:
            raise SystemExit(f"row {lineno}: target must be non-empty")
        if not smiles:
            raise SystemExit(f"row {lineno}: smiles must be non-empty")
        key = (target, method)
        counters[key] = counters.get(key, 0) + 1
        tasks.append({"molecule_id": counters[key], "method": method,
                      "target": target, "smiles": smiles})
    return tasks


def load_config(config_dir: Path, target: str) -> dict:
    cfg_path = Path(config_dir) / f"{target}.json"
    if not cfg_path.is_file():
        raise SystemExit(
            f"missing config {cfg_path}; run evaluation/docking/prepare_receptors.py first"
        )
    return json.loads(cfg_path.read_text())


def resolve_path(maybe_relative: str) -> Path:
    path = Path(maybe_relative)
    return path if path.is_absolute() else ROOT / path


# --------------------------------------------------------------------------- #
# Docking pipeline (RDKit -> Meeko -> Vina); heavy imports stay local so the
# module remains importable (and testable) without the scientific stack.
# --------------------------------------------------------------------------- #
def _init_worker(receptor_pdbqt: str, center, size, scoring: str,
                 exhaustiveness: int, num_modes: int, seed: int) -> None:
    """Create one Vina instance with precomputed affinity maps per worker."""
    global _VINA, _LIGAND_PREP, _EXHAUSTIVENESS, _NUM_MODES
    from meeko import MoleculePreparation
    from vina import Vina

    vina = Vina(sf_name=scoring, cpu=1, seed=seed, verbosity=0)
    vina.set_receptor(receptor_pdbqt)
    # No ligand is loaded at this point, so maps are computed for all 22 Vina
    # atom types once and can be reused for every ligand of this target.
    vina.compute_vina_maps(center=list(center), box_size=list(size))
    _VINA = vina
    _LIGAND_PREP = MoleculePreparation()
    _EXHAUSTIVENESS = exhaustiveness
    _NUM_MODES = num_modes


def _smiles_to_pdbqt(smiles: str) -> str:
    """SMILES -> 3D conformer (ETKDGv3, seed=42) -> MMFF/UFF -> PDBQT string."""
    from meeko import PDBQTWriterLegacy
    from rdkit import Chem
    from rdkit.Chem import AllChem

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("RDKit could not parse SMILES")
    fragments = Chem.GetMolFrags(mol, asMols=True)
    if len(fragments) > 1:
        mol = max(fragments, key=lambda m: m.GetNumHeavyAtoms())
    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = EMBEDDING_SEED
    if AllChem.EmbedMolecule(mol, params) != 0:
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError("3D conformer embedding failed")

    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        AllChem.UFFOptimizeMolecule(mol)

    setups = _LIGAND_PREP.prepare(mol)
    written = PDBQTWriterLegacy.write_string(setups[0])
    # meeko returns (pdbqt_string, ok) on some versions, plain str on others.
    pdbqt_string = written[0] if isinstance(written, tuple) else written
    return pdbqt_string


def _dock_one(task: tuple) -> dict:
    """Dock a single molecule. Never raises: failures are returned as data."""
    molecule_id, method, target, smiles = task
    try:
        pdbqt_string = _smiles_to_pdbqt(smiles)
        _VINA.set_ligand_from_string(pdbqt_string)
        _VINA.dock(exhaustiveness=_EXHAUSTIVENESS, n_poses=_NUM_MODES)
        score = float(_VINA.energies(n_poses=_NUM_MODES)[0][0])
        if not math.isfinite(score):
            raise ValueError(f"non-finite Vina score: {score}")
        error = None
    except Exception as exc:  # failures are data, not crashes
        score = None
        error = f"{type(exc).__name__}: {exc}"
    return {"molecule_id": molecule_id, "method": method, "target": target,
            "smiles": smiles, "vina_score": score, "error": error}


def dock_target(rows: list[dict], config: dict, n_jobs: int) -> list[dict]:
    """Dock all rows of one target, optionally in parallel worker processes."""
    receptor = resolve_path(config["receptor_pdbqt"])
    if not receptor.is_file():
        raise SystemExit(
            f"missing prepared receptor {receptor}; run prepare_receptors.py first"
        )
    initargs = (
        str(receptor),
        config["box"]["center"],
        config["box"]["size"],
        config["vina"]["scoring"],
        config["vina"]["exhaustiveness"],
        config["vina"]["num_modes"],
        config["vina"]["seed"],
    )
    tasks = [(r["molecule_id"], r["method"], r["target"], r["smiles"]) for r in rows]
    if n_jobs > 1 and len(tasks) > 1:
        import multiprocessing as mp

        with mp.Pool(processes=n_jobs, initializer=_init_worker, initargs=initargs) as pool:
            return pool.map(_dock_one, tasks, chunksize=1)
    _init_worker(*initargs)
    return [_dock_one(task) for task in tasks]


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def write_scores(results: list[dict], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ordered = sorted(results, key=lambda r: (r["target"], r["method"], r["molecule_id"]))
    with out_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["molecule_id", "method", "target", "smiles", "vina_score"])
        for r in ordered:
            score = f"{r['vina_score']:.3f}" if r["vina_score"] is not None else ""
            writer.writerow([r["molecule_id"], r["method"], r["target"], r["smiles"], score])


def write_errors(results: list[dict], out_path: Path) -> Path | None:
    errors = [r for r in results if r["error"]]
    if not errors:
        return None
    err_path = out_path.with_name(f"{out_path.stem}_errors.csv")
    with err_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["target", "method", "molecule_id", "smiles", "error"])
        for r in errors:
            writer.writerow([r["target"], r["method"], r["molecule_id"], r["smiles"], r["error"]])
    return err_path


def write_summary(results: list[dict], out_path: Path) -> Path:
    summary_path = out_path.with_name(f"{out_path.stem}_summary.csv")
    groups: dict[tuple[str, str], list[dict]] = {}
    for r in results:
        groups.setdefault((r["target"], r["method"]), []).append(r)
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["target", "method", "n_molecules", "n_success", "n_failed",
                         "vina_score_best", "vina_score_mean", "vina_score_median"])
        for target, method in sorted(groups):
            rows = groups[(target, method)]
            scores = [r["vina_score"] for r in rows if r["vina_score"] is not None]
            if scores:
                stats = (f"{min(scores):.3f}", f"{statistics.mean(scores):.3f}",
                         f"{statistics.median(scores):.3f}")
            else:
                stats = ("", "", "")
            writer.writerow([target, method, len(rows), len(scores),
                             len(rows) - len(scores), *stats])
    return summary_path


# --------------------------------------------------------------------------- #
def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--input", required=True, type=Path,
                        help="CSV with columns target,method,smiles")
    parser.add_argument("--out", required=True, type=Path,
                        help="output CSV path (molecule_id,method,target,smiles,vina_score)")
    parser.add_argument("--config-dir", type=Path, default=DEFAULT_CONFIG_DIR,
                        help=f"per-target config dir (default: {DEFAULT_CONFIG_DIR})")
    parser.add_argument("--n-jobs", type=int, default=1,
                        help="worker processes per target (Vina itself stays single-threaded)")
    parser.add_argument("--limit", type=int, default=None,
                        help="dock at most N molecules per target (smoke tests)")
    args = parser.parse_args(argv)

    tasks = read_input(args.input)
    available = {p.stem for p in args.config_dir.glob("*.json")}
    unknown = sorted({t["target"] for t in tasks} - available)
    if unknown:
        raise SystemExit(f"unknown target(s) {unknown}; available configs: {sorted(available)}")

    by_target: dict[str, list[dict]] = {}
    for task in tasks:
        by_target.setdefault(task["target"], []).append(task)

    results: list[dict] = []
    for target in sorted(by_target):
        rows = by_target[target][: args.limit] if args.limit else by_target[target]
        config = load_config(args.config_dir, target)
        print(f"[{target}] docking {len(rows)} molecule(s) ...", flush=True)
        target_results = dock_target(rows, config, args.n_jobs)
        n_failed = sum(1 for r in target_results if r["error"])
        print(f"[{target}] done: {len(target_results) - n_failed} ok, {n_failed} failed", flush=True)
        results.extend(target_results)

    write_scores(results, args.out)
    err_path = write_errors(results, args.out)
    summary_path = write_summary(results, args.out)
    print(f"wrote {args.out} ({len(results)} row(s))", flush=True)
    if err_path:
        print(f"wrote {err_path} ({sum(1 for r in results if r['error'])} failure(s))", flush=True)
    print(f"wrote {summary_path}", flush=True)


if __name__ == "__main__":
    main()
