#!/usr/bin/env python3
"""Analyze the independent capped high=9 size-penalty rerun with the existing six-metric protocol."""

import csv
import gzip
import hashlib
import importlib.util
from importlib.metadata import version
import json
import platform
from collections import Counter, defaultdict
import shutil
from pathlib import Path


def extra_curves(analysis, rows):
    """Keep signed training rewards separate from pre-penalty sigmoid curves."""
    metrics = {
        "average_generated_heavy_atoms": ("generated_heavy_atoms", "Generated heavy atoms (scored samples)"),
        "average_base_reward": ("base_reward", "Base reward before repeat penalty (all attempts)"),
        "average_training_reward": ("training_reward", "Training reward after repeat penalty (all attempts)"),
    }
    buckets = defaultdict(list)
    for row in rows:
        buckets[(row['target']+'_'+row['mode'],int(row['step']))].append(row)
    tasks = sorted({task for task, _ in buckets})
    series = {task:{} for task in tasks}
    for metric,(column,label) in metrics.items():
        analysis.METRIC_LABELS[metric]=label
        for task in tasks:
            values=[]
            for step in range(1,11):
                selected=[float(row[column]) for row in buckets[task,step] if row[column]!='']
                values.append(sum(selected)/len(selected) if selected else None)
            series[task][metric]=values
        figdir=analysis.ANALYSIS_ROOT/metric/'figure';figdir.mkdir(parents=True,exist_ok=True)
        tabdir=analysis.ANALYSIS_ROOT/metric/'table';tabdir.mkdir(parents=True,exist_ok=True)
        for task in tasks:
            analysis.plot_single(list(range(1,11)),series[task][metric],metric,task,figdir/(task+'.png'))
            analysis.write_series_table([task],series,metric,tabdir/(task+'.csv'))
        groups=[('all_tasks',tasks,'summary'),('libinvent',[t for t in tasks if t.endswith('_libinvent')],'libinvent'),('linkinvent',[t for t in tasks if t.endswith('_linkinvent')],'linkinvent')]
        for stem,names,folder in groups:
            analysis.plot_multi(list(range(1,11)),series,metric,names,stem,figdir/(stem+'.png'),2 if stem=='all_tasks' else 1)
            analysis.write_series_table(names,series,metric,tabdir/(stem+'.csv'))
            shutil.copyfile(figdir/(stem+'.png'),analysis.ANALYSIS_ROOT/folder/'figure'/(metric+'.png'))
            shutil.copyfile(tabdir/(stem+'.csv'),analysis.ANALYSIS_ROOT/folder/'table'/(metric+'.csv'))


def main():
    package = Path(__file__).resolve().parents[1]
    repo = package.parents[3]
    data = repo / "results/reinvent-high9-cap9-size001-1000-20261025/all_generated.csv.gz"
    experiment_path = data.parent / "experiment.json"
    experiment = json.loads(experiment_path.read_text())
    assert experiment["sigmoid"] == {"low": 4.0, "high": 9.0, "k": 0.25}
    assert experiment["reward_adjustment"] == {"affinity_cap": 9.0, "generated_heavy_atom_penalty": 0.01}
    assert (experiment["batch_size"], experiment["steps"], experiment["seed"]) == (100, 10, 42)
    with gzip.open(data, "rt") as handle:
        counts = Counter((row["target"], row["mode"], int(row["step"]))
                         for row in csv.DictReader(handle))
    tasks = {(target, mode) for target, mode, _ in counts}
    expected = {(target, mode, step): 100 for target, mode in tasks for step in range(1, 11)}
    if len(tasks) != 30 or counts != expected:
        raise ValueError("Expected exactly 30 tasks x 10 steps x 100 attempts; refusing to truncate")
    source = package.parent / "generation_analysis_reinvent_1000/src/run_analysis.py"
    spec = importlib.util.spec_from_file_location("existing_generation_analysis", source)
    analysis = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(analysis)
    analysis.DATA_PATH = data
    analysis.ANALYSIS_ROOT = package / "analysis"
    analysis.main()
    with gzip.open(data, "rt") as handle:
        extra_curves(analysis, list(csv.DictReader(handle)))
    baseline = package.parent / "high9_1000_seed42_20261025/analysis/summary_table.csv"
    with baseline.open() as handle:
        old = {row["task"]: row for row in csv.DictReader(handle)}
    with (package / "analysis/summary_table.csv").open() as handle:
        current = list(csv.DictReader(handle))
    if set(old) != {row["task"] for row in current}:
        raise ValueError("Baseline and rerun task sets differ")
    with (package / "analysis/comparison_high9_baseline.csv").open("w", newline="") as handle:
        fields = ["task", "metric", "high9_baseline", "cap9_size001", "difference_size_minus_high9"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in current:
            for metric in analysis.SUMMARY_COLUMNS:
                before, after = float(old[row["task"]][metric]), float(row[metric])
                writer.writerow(dict(zip(fields, [row["task"], metric, before, after, after - before])))
    provenance = {
        "independent_rerun": True, "rows": sum(counts.values()), "tasks": len(tasks),
        "experiment": experiment,
        "input": str(data.relative_to(repo)),
        "input_sha256": hashlib.sha256(data.read_bytes()).hexdigest(),
        "shared_analysis_source": str(source.relative_to(repo)),
        "shared_analysis_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
        "baseline_summary": str(baseline.relative_to(repo)),
        "baseline_summary_sha256": hashlib.sha256(baseline.read_bytes()).hexdigest(),
        "analysis_versions": {name: version(name) for name in ("rdkit", "numpy", "matplotlib")},
        "analysis_python": platform.python_version(),
    }
    (package / "analysis/provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
