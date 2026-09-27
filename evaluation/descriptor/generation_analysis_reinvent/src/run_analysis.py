#!/usr/bin/env python3
"""Summarize REINVENT generation runs and draw batch-wise curves.

Reads : <repo>/results/reinvent-20260927/all_generated.csv.gz
        (30 tasks = 15 targets x {libinvent, linkinvent}, 100 batches x 100 molecules)
Writes: analysis/summary_table.csv
        analysis/<metric>/figure/*.png and analysis/<metric>/table/*.csv
        analysis/{summary,linkinvent,libinvent}/{figure,table}/*  (copies of the
        multi-line overview charts/tables)

Usage:
    python3 src/run_analysis.py

Metric definitions and all processing conventions are documented in README.md.
"""

from __future__ import annotations

import csv
import gzip
import shutil
import sys
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors

RDLogger.DisableLog("rdApp.*")

HERE = Path(__file__).resolve().parent
PKG_ROOT = HERE.parent  # generation_analysis_reinvent/
REPO_ROOT = PKG_ROOT.parents[2]  # repository root
DATA_PATH = REPO_ROOT / "results" / "reinvent-20260927" / "all_generated.csv.gz"
CRYSTAL_DIR = REPO_ROOT / "real-world_dataset"
ANALYSIS_ROOT = PKG_ROOT / "analysis"

# Reuse the repository's validated MOL2 loader (charge normalization, aromaticity
# fixes, stereo perceived from 3D) so lead molecules match the project convention.
sys.path.insert(0, str(REPO_ROOT / "evaluation" / "descriptor" / "similarity_test" / "src"))
from crystal_actives_similarity import load_crystal_ligand  # noqa: E402

N_STEPS = 100
MODES = ("libinvent", "linkinvent")

# Exactly the six requested per-batch metrics; order defines plot/table order.
METRICS = [
    "average_oracle_score_origin",
    "average_oracle_score_sigmoid",
    "average_atom_number",
    "average_MW",
    "repetition_count",
    "noncanonical_count",
]
METRIC_LABELS = {
    "average_oracle_score_origin": "Average oracle score (original)",
    "average_oracle_score_sigmoid": "Average oracle score (sigmoid)",
    "average_atom_number": "Average heavy atom number",
    "average_MW": "Average molecular weight",
    "repetition_count": "Repetition count",
    "noncanonical_count": "Invalid SMILES count",
}

SUMMARY_COLUMNS = [
    "oracle_score_original",
    "final_oracle_score_original",
    "best_oracle_score_original",
    "oracle_score_sigmoid",
    "final_oracle_score_sigmoid",
    "best_oracle_score_sigmoid",
    "noncanonical_rate",
    "final_noncanonical_rate",
    "repetition_rate",
    "final_repetition_rate",
    "recall_count",
]

for _style in ("seaborn-v0_8-whitegrid", "ggplot"):
    try:
        plt.style.use(_style)
        break
    except OSError:
        continue


class Batch:
    """Accumulators for one (task, batch) cell."""

    __slots__ = (
        "aff_sum", "aff_n", "sig_sum", "sig_n",
        "heavy_sum", "heavy_n", "mw_sum", "mw_n",
        "smiles", "invalid", "total",
    )

    def __init__(self) -> None:
        self.aff_sum = 0.0
        self.aff_n = 0
        self.sig_sum = 0.0
        self.sig_n = 0
        self.heavy_sum = 0.0
        self.heavy_n = 0
        self.mw_sum = 0.0
        self.mw_n = 0
        self.smiles: Counter = Counter()
        self.invalid = 0
        self.total = 0


def load_leads(targets):
    """Return {target: strict canonical SMILES of the co-crystal ligand}."""
    leads = {}
    for target in targets:
        mol = load_crystal_ligand(CRYSTAL_DIR / target / "crystal.mol2")
        leads[target] = Chem.MolToSmiles(mol)  # isomeric (stereo from 3D)
    return leads


def accumulate(rows_reader, leads):
    """Single pass over all generated molecules."""
    batches = defaultdict(Batch)  # (task, step) -> Batch
    global_smiles = defaultdict(Counter)  # task -> Counter over all batches
    recall = defaultdict(int)  # task -> lead occurrences
    task_of = {}  # task -> (target, mode)

    n_rows = 0
    for row in rows_reader:
        n_rows += 1
        if n_rows % 50000 == 0:
            print(f"  processed {n_rows} rows ...")
        target = row["target"]
        mode = row["mode"]
        task = f"{target}_{mode}"
        step = int(row["step"])
        task_of[task] = (target, mode)
        batch = batches[(task, step)]
        smiles = row["smiles"]
        batch.total += 1
        batch.smiles[smiles] += 1
        global_smiles[task][smiles] += 1

        if row["planet_status"] == "ok":
            try:
                batch.aff_sum += float(row["planet_affinity"])
                batch.aff_n += 1
                batch.sig_sum += float(row["sigmoid_reward"])
                batch.sig_n += 1
            except ValueError:
                pass

        # SmilesState (REINVENT4 sample_batch.py): 0 = INVALID, 1 = VALID,
        # 2 = DUPLICATE (valid molecule repeated within the sampling batch).
        mol = Chem.MolFromSmiles(smiles)
        if mol is not None:
            if Chem.MolToSmiles(mol) == leads[target]:
                recall[task] += 1
            if row["smiles_state"] != "0":  # structurally legal molecule
                batch.heavy_sum += mol.GetNumHeavyAtoms()
                batch.heavy_n += 1
                batch.mw_sum += Descriptors.MolWt(mol)
                batch.mw_n += 1
        if row["smiles_state"] == "0":
            batch.invalid += 1

    return batches, global_smiles, recall, task_of


def repetition_count(counter: Counter) -> int:
    """Number of duplicate occurrences: a SMILES seen n times contributes n-1."""
    return sum(n - 1 for n in counter.values() if n >= 2)


def batch_series(batches, task):
    """Return {metric: [value_batch1 .. value_batch100]} for one task."""
    series = {m: [None] * N_STEPS for m in METRICS}
    for step in range(1, N_STEPS + 1):
        batch = batches.get((task, step))
        if batch is None:
            continue
        if batch.aff_n:
            series["average_oracle_score_origin"][step - 1] = batch.aff_sum / batch.aff_n
            series["average_oracle_score_sigmoid"][step - 1] = batch.sig_sum / batch.sig_n
        if batch.heavy_n:
            series["average_atom_number"][step - 1] = batch.heavy_sum / batch.heavy_n
            series["average_MW"][step - 1] = batch.mw_sum / batch.mw_n
        series["repetition_count"][step - 1] = repetition_count(batch.smiles)
        series["noncanonical_count"][step - 1] = batch.invalid
    return series


def summary_row(task, batches, global_smiles, recall):
    """Return the 11 summary-table values for one task."""
    aff_means, sig_means = [], []
    aff_sum = sig_sum = aff_n = sig_n = 0
    total = invalid = 0
    for step in range(1, N_STEPS + 1):
        batch = batches.get((task, step))
        if batch is None:
            continue
        if batch.aff_n:
            aff_means.append(batch.aff_sum / batch.aff_n)
            sig_means.append(batch.sig_sum / batch.sig_n)
            aff_sum += batch.aff_sum
            sig_sum += batch.sig_sum
            aff_n += batch.aff_n
            sig_n += batch.sig_n
        total += batch.total
        invalid += batch.invalid

    last = batches.get((task, N_STEPS))
    final_aff = last.aff_sum / last.aff_n if last is not None and last.aff_n else None
    final_sig = last.sig_sum / last.sig_n if last is not None and last.sig_n else None
    final_total = last.total if last is not None else None
    final_invalid = last.invalid if last is not None else None
    final_rep = repetition_count(last.smiles) if last is not None else None

    return {
        "oracle_score_original": aff_sum / aff_n if aff_n else None,
        "final_oracle_score_original": final_aff,
        "best_oracle_score_original": max(aff_means) if aff_means else None,
        "oracle_score_sigmoid": sig_sum / sig_n if sig_n else None,
        "final_oracle_score_sigmoid": final_sig,
        "best_oracle_score_sigmoid": max(sig_means) if sig_means else None,
        "noncanonical_rate": invalid / total if total else None,
        "final_noncanonical_rate": final_invalid / final_total if final_total else None,
        "repetition_rate": repetition_count(global_smiles[task]) / total if total else None,
        "final_repetition_rate": final_rep / final_total if final_total else None,
        "recall_count": recall[task],
    }


def fmt(value):
    if value is None:
        return ""
    if isinstance(value, int):
        return str(value)
    return f"{value:.6f}"


def write_summary_table(rows, path):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["task"] + SUMMARY_COLUMNS)
        for name, values in rows:
            writer.writerow([name] + [fmt(values[c]) for c in SUMMARY_COLUMNS])


def write_series_table(names, series_map, metric, path):
    """Wide CSV: first column `batch`, one column per task."""
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["batch"] + names)
        for step in range(N_STEPS):
            writer.writerow(
                [step + 1]
                + [fmt(series_map[name][metric][step]) for name in names]
            )


def plot_single(xs, ys, metric, task, path):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(xs, ys, color="#1f77b4", linewidth=1.8)
    ax.set_title(f"{task} - {metric}", fontsize=13)
    ax.set_xlabel("Batch", fontsize=11)
    ax.set_ylabel(METRIC_LABELS[metric], fontsize=11)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def plot_multi(xs, series_map, metric, names, title, path, legend_ncol):
    fig, ax = plt.subplots(figsize=(12, 7))
    colors = plt.cm.turbo(np.linspace(0.05, 0.95, len(names)))
    for color, name in zip(colors, names):
        ax.plot(xs, series_map[name][metric], color=color, linewidth=1.1,
                alpha=0.85, label=name)
    ax.set_title(f"{title} - {metric}", fontsize=13)
    ax.set_xlabel("Batch", fontsize=11)
    ax.set_ylabel(METRIC_LABELS[metric], fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="center left", bbox_to_anchor=(1.02, 0.5),
              fontsize=7, ncol=legend_ncol, frameon=True)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    if not DATA_PATH.is_file():
        raise SystemExit(f"data file not found: {DATA_PATH}")

    print(f"reading {DATA_PATH} ...")
    with gzip.open(DATA_PATH, "rt") as handle:
        targets = sorted({row["target"] for row in csv.DictReader(handle)})
    task_order = [f"{t}_{m}" for t in targets for m in MODES]

    leads = load_leads(targets)
    print(f"loaded {len(leads)} lead molecules")

    print("computing metrics ...")
    with gzip.open(DATA_PATH, "rt") as handle:
        batches, global_smiles, recall, _ = accumulate(csv.DictReader(handle), leads)

    series_map = {task: batch_series(batches, task) for task in task_order}
    rows = [(task, summary_row(task, batches, global_smiles, recall)) for task in task_order]

    mode_rows = {}
    for mode in MODES:
        names = [f"{t}_{mode}" for t in targets]
        mode_rows[f"{mode}_average"] = {
            c: sum(d[c] for task, d in rows if task.endswith(f"_{mode}")) / len(names)
            for c in SUMMARY_COLUMNS
        }
    mode_rows["average"] = {
        c: sum(d[c] for _, d in rows) / len(task_order) for c in SUMMARY_COLUMNS
    }

    # ---- write outputs -------------------------------------------------
    for metric in METRICS:
        (ANALYSIS_ROOT / metric / "figure").mkdir(parents=True, exist_ok=True)
        (ANALYSIS_ROOT / metric / "table").mkdir(parents=True, exist_ok=True)
    for folder in ("summary", "linkinvent", "libinvent"):
        (ANALYSIS_ROOT / folder / "figure").mkdir(parents=True, exist_ok=True)
        (ANALYSIS_ROOT / folder / "table").mkdir(parents=True, exist_ok=True)

    write_summary_table(
        rows + [(k, mode_rows[k]) for k in ("libinvent_average", "linkinvent_average", "average")],
        ANALYSIS_ROOT / "summary_table.csv",
    )
    print(f"wrote {ANALYSIS_ROOT / 'summary_table.csv'}")

    xs = list(range(1, N_STEPS + 1))
    for metric in METRICS:
        fig_dir = ANALYSIS_ROOT / metric / "figure"
        tab_dir = ANALYSIS_ROOT / metric / "table"
        for task in task_order:
            plot_single(xs, series_map[task][metric], metric, task,
                        fig_dir / f"{task}.png")
            write_series_table([task], series_map, metric, tab_dir / f"{task}.csv")

        link_names = [f"{t}_linkinvent" for t in targets]
        lib_names = [f"{t}_libinvent" for t in targets]
        groups = [
            ("all_tasks", task_order, "All tasks", 2, "summary"),
            ("linkinvent", link_names, "LinkINVENT tasks", 1, "linkinvent"),
            ("libinvent", lib_names, "LibINVENT tasks", 1, "libinvent"),
        ]
        for stem, names, title, ncol, copy_folder in groups:
            fig_path = fig_dir / f"{stem}.png"
            tab_path = tab_dir / f"{stem}.csv"
            plot_multi(xs, series_map, metric, names, title, fig_path, ncol)
            write_series_table(names, series_map, metric, tab_path)
            shutil.copyfile(fig_path, ANALYSIS_ROOT / copy_folder / "figure" / f"{metric}.png")
            shutil.copyfile(tab_path, ANALYSIS_ROOT / copy_folder / "table" / f"{metric}.csv")
        print(f"finished metric {metric}")

    print("done")


if __name__ == "__main__":
    main()
