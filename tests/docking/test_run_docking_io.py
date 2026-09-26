"""Tests for the pure-stdlib parts of run_docking.py (input contract, id
numbering, output writers). RDKit/Meeko/Vina are imported lazily inside
functions, so these tests run without the scientific stack installed."""

import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "evaluation" / "docking"))

from run_docking import read_input, write_errors, write_scores, write_summary


def _write_csv(path, rows, header=("target", "method", "smiles")):
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)
    return path


def test_read_input_assigns_per_group_ids(tmp_path):
    csv_path = _write_csv(tmp_path / "in.csv", [
        ["adrb1", "reinvent", "CCO"],
        ["cdk2", "delete", "c1ccccc1"],
        ["adrb1", "reinvent", "CCN"],
        ["adrb1", "delete", "O=C=O"],
    ])
    tasks = read_input(csv_path)
    ids = {(t["target"], t["method"]): [] for t in tasks}
    for t in tasks:
        ids[(t["target"], t["method"])].append(t["molecule_id"])
    assert ids[("adrb1", "reinvent")] == [1, 2]
    assert ids[("adrb1", "delete")] == [1]
    assert ids[("cdk2", "delete")] == [1]


def test_read_input_rejects_bad_method_and_missing_columns(tmp_path):
    bad_method = _write_csv(tmp_path / "bad_method.csv", [["adrb1", "other", "CCO"]])
    with pytest.raises(SystemExit):
        read_input(bad_method)
    missing_col = _write_csv(tmp_path / "missing.csv", [["adrb1", "CCO"]], header=("target", "smiles"))
    with pytest.raises(SystemExit):
        read_input(missing_col)
    empty = _write_csv(tmp_path / "empty.csv", [])
    with pytest.raises(SystemExit):
        read_input(empty)


def _result(target, method, mid, score, error=None):
    return {"molecule_id": mid, "method": method, "target": target,
            "smiles": "CCO", "vina_score": score, "error": error}


def test_write_scores_exact_five_columns_sorted(tmp_path):
    out = tmp_path / "scores.csv"
    results = [
        _result("cdk2", "delete", 1, -5.0),
        _result("adrb1", "reinvent", 2, None, "ValueError: boom"),
        _result("adrb1", "reinvent", 1, -7.1234),
    ]
    write_scores(results, out)
    with out.open() as handle:
        rows = list(csv.reader(handle))
    assert rows[0] == ["molecule_id", "method", "target", "smiles", "vina_score"]
    assert [r[0] for r in rows[1:]] == ["1", "2", "1"]  # adrb1 first, then cdk2
    assert rows[1][4] == "-7.123"
    assert rows[2][4] == ""  # failed molecule -> empty score
    assert all(len(r) == 5 for r in rows)


def test_write_errors_only_when_failures(tmp_path):
    out = tmp_path / "scores.csv"
    assert write_errors([_result("adrb1", "reinvent", 1, -7.0)], out) is None
    err_path = write_errors([_result("adrb1", "reinvent", 1, None, "ValueError: boom")], out)
    assert err_path is not None and err_path.name == "scores_errors.csv"
    with err_path.open() as handle:
        rows = list(csv.reader(handle))
    assert rows[0] == ["target", "method", "molecule_id", "smiles", "error"]
    assert "boom" in rows[1][4]


def test_write_summary_statistics(tmp_path):
    out = tmp_path / "scores.csv"
    results = [
        _result("adrb1", "reinvent", 1, -6.0),
        _result("adrb1", "reinvent", 2, -8.0),
        _result("adrb1", "reinvent", 3, None, "ValueError: boom"),
        _result("adrb1", "delete", 1, -7.0),
    ]
    summary_path = write_summary(results, out)
    with summary_path.open() as handle:
        rows = {r[1]: r for r in csv.reader(handle) if r[0] == "adrb1"}
    # reinvent: scores -6, -8 -> best -8, mean -7, median -7, 1 failure
    assert rows["reinvent"][2:5] == ["3", "2", "1"]
    assert rows["reinvent"][5:] == ["-8.000", "-7.000", "-7.000"]
    # delete: single score
    assert rows["delete"][2:5] == ["1", "1", "0"]
    assert rows["delete"][5:] == ["-7.000", "-7.000", "-7.000"]
