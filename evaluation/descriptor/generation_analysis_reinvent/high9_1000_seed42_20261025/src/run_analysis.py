#!/usr/bin/env python3
"""Analyze the independent high=9 rerun with the existing six-metric protocol."""

import csv
import gzip
import hashlib
import importlib.util
from importlib.metadata import version
import json
import platform
from collections import Counter
from pathlib import Path


def main():
    package = Path(__file__).resolve().parents[1]
    repo = package.parents[3]
    data = repo / "results/reinvent-high9-1000-20261025/all_generated.csv.gz"
    experiment_path = data.parent / "experiment.json"
    experiment = json.loads(experiment_path.read_text())
    assert experiment["sigmoid"] == {"low": 4.0, "high": 9.0, "k": 0.25}
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
    baseline = package.parent / "generation_analysis_reinvent_1000/analysis/summary_table.csv"
    with baseline.open() as handle:
        old = {row["task"]: row for row in csv.DictReader(handle)}
    with (package / "analysis/summary_table.csv").open() as handle:
        current = list(csv.DictReader(handle))
    if set(old) != {row["task"] for row in current}:
        raise ValueError("Baseline and rerun task sets differ")
    with (package / "analysis/comparison_high16_first1000.csv").open("w", newline="") as handle:
        fields = ["task", "metric", "high16_first1000", "high9_rerun", "difference_high9_minus_high16"]
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
