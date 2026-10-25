"""Independently audit completed manual PLANET RL outputs before transfer/analysis."""

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path

from prepare_manual_planet_rl import read_experiment, sha256, write_json


def validate(root, experiment_path, report_path):
    root = Path(root)
    cfg = read_experiment(experiment_path)
    manifest = json.loads((root / "run.json").read_text())
    status = json.loads((root / "status.json").read_text())
    if manifest["experiment"] != cfg or json.loads((root / "experiment.json").read_text()) != cfg:
        raise ValueError("Run configuration does not match the requested experiment")
    for filename, digest in manifest["input_hashes"].items():
        if sha256(filename) != digest:
            raise ValueError(f"Source input changed after preparation: {filename}")
    bound_clients = {entry["client_config"] for entry in manifest["targets"]}
    for filename, digest in manifest["prepared_hashes"].items():
        # The runner intentionally binds each client's oracle ID at startup.
        if filename not in bound_clients and sha256(root / filename) != digest:
            raise ValueError(f"Prepared file changed after preparation: {filename}")
    expected = {(e["target"], t["mode"]): (e, t)
                for e in manifest["targets"] for t in e["tasks"]}
    names = {f"{target}/{mode}" for target, mode in expected}
    if (status["state"] != "complete" or set(status["tasks"]) != names
            or any(value != "complete" for value in status["tasks"].values())):
        raise ValueError("Not all requested tasks completed")
    records = {}
    with (root / "all_generated.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            key = (row["target"], row["mode"])
            if key not in expected:
                raise ValueError(f"Unexpected task: {key}")
            records.setdefault(key, []).append(row)
    if set(records) != set(expected):
        raise ValueError("Missing task outputs")
    results = {}
    low, high, k = (cfg["sigmoid"][name] for name in ("low", "high", "k"))
    adjustment = cfg.get("reward_adjustment")
    for key, (entry, task) in expected.items():
        rows = records[key]
        folder = root / task["directory"]
        with (folder / "all_generated.csv").open(newline="") as handle:
            if list(csv.DictReader(handle)) != rows:
                raise ValueError(f"Combined output differs from task output: {key}")
        counts = Counter(int(row["step"]) for row in rows)
        if counts != Counter({step: cfg["batch_size"] for step in range(1, cfg["steps"] + 1)}):
            raise ValueError(f"Wrong sampling budget: {key}")
        sample_ids = Counter((int(row["step"]), int(row["sample_index"])) for row in rows)
        if sample_ids != Counter({(step, index): 1 for step in range(1, cfg["steps"] + 1)
                                  for index in range(1, cfg["batch_size"] + 1)}):
            raise ValueError(f"Missing/duplicate sample indices: {key}")
        oracle = json.loads((root / entry["oracle_manifest"]).read_text())
        seen, states, statuses, affinities = set(), Counter(), Counter(), []
        for row in rows:
            state = int(row["smiles_state"])
            if state not in (0, 1, 2) or int(row["seed"]) != cfg["seed"]:
                raise ValueError(f"Unexpected state/seed: {key}")
            if row["input_smiles"] != task["input_smiles"]:
                raise ValueError(f"Unexpected design fragment: {key}")
            states[state] += 1
            statuses[row["planet_status"]] += 1
            duplicate = state != 0 and row["smiles"] in seen
            if bool(int(row["is_repeat"])) != duplicate:
                raise ValueError(f"Incorrect repeat annotation: {key}")
            if state != 0:
                seen.add(row["smiles"])
            reward = float(row["sigmoid_reward"])
            training = float(row["training_reward"])
            if not (math.isfinite(reward) and math.isfinite(training) and 0 <= reward <= 1
                    and (adjustment is not None or 0 <= training <= 1)):
                raise ValueError(f"Non-finite/out-of-range reward: {key}")
            base_reward = float(row["base_reward"]) if adjustment is not None else reward
            if not math.isfinite(base_reward):
                raise ValueError(f"Non-finite base reward: {key}")
            if not math.isclose(training, base_reward * (cfg["penalty_multiplier"] if duplicate else 1), abs_tol=1e-6):
                raise ValueError(f"Incorrect training reward: {key}")
            if row["planet_status"] == "ok":
                score = float(row["planet_affinity"])
                if not math.isfinite(score) or state == 0:
                    raise ValueError(f"Invalid successful oracle score: {key}")
                if row["oracle_id"] != oracle["oracle_id"] or row["planet_target_id"] != key[0]:
                    raise ValueError(f"Oracle identity mismatch: {key}")
                scoring_value = min(score, adjustment["affinity_cap"]) if adjustment else score
                z = math.log(10) * 10 * k / (high - low) * (scoring_value - (high + low) / 2)
                sigmoid = 1 / (1 + math.exp(-z)) if z >= 0 else math.exp(z) / (1 + math.exp(z))
                if not math.isclose(reward, sigmoid, abs_tol=2e-6):
                    raise ValueError(f"Sigmoid disagrees with requested parameters: {key}")
                if adjustment:
                    # Recompute from the two molecular graphs, not from reported counts.
                    from rdkit import Chem
                    product = Chem.MolFromSmiles(row["smiles"])
                    fixed = Chem.MolFromSmiles(task["input_smiles"].replace("|", "."))
                    if product is None or fixed is None:
                        raise ValueError(f"Cannot audit generated atoms: {key}")
                    added = (sum(a.GetAtomicNum() > 1 for a in product.GetAtoms())
                             - sum(a.GetAtomicNum() > 1 for a in fixed.GetAtoms()))
                    penalty = added * adjustment["generated_heavy_atom_penalty"]
                    if (added < 0 or int(row["generated_heavy_atoms"]) != added
                            or not math.isclose(float(row["capped_planet_affinity"]), scoring_value, abs_tol=1e-6)
                            or not math.isclose(float(row["size_penalty"]), penalty, abs_tol=1e-6)
                            or not math.isclose(base_reward, sigmoid - penalty, abs_tol=2e-6)):
                        raise ValueError(f"Size-adjusted reward disagrees with formula: {key}")
                affinities.append(score)
            elif row["planet_status"] not in ("not_scored", "invalid_input") or row["planet_affinity"]:
                raise ValueError(f"Unexpected failed scoring record: {key}")
            if state == 0 and training != 0:
                raise ValueError(f"Invalid molecule has nonzero reward: {key}")
            if adjustment and row["planet_status"] != "ok":
                if (reward != 0 or base_reward != 0 or training != 0
                        or row["generated_heavy_atoms"] != "" or row["capped_planet_affinity"] != ""
                        or float(row["size_penalty"]) != 0):
                    raise ValueError(f"Unscored molecule has fabricated reward terms: {key}")
        checkpoint = folder / "agent.chkpt"
        if not checkpoint.is_file() or checkpoint.stat().st_size == 0:
            raise ValueError(f"Missing final checkpoint: {key}")
        results[f"{key[0]}/{key[1]}"] = {
            "attempts": len(rows), "state_counts": dict(states), "status_counts": dict(statuses),
            "valid_unique_smiles": len(seen), "repeats": sum(int(r["is_repeat"]) for r in rows),
            "raw_mean": sum(affinities) / len(affinities) if affinities else None,
            "raw_min": min(affinities) if affinities else None,
            "raw_max": max(affinities) if affinities else None,
            "checkpoint_sha256": sha256(checkpoint), "checkpoint_bytes": checkpoint.stat().st_size,
        }
    report = {"passed": True, "tasks": len(results), "rows": sum(len(r) for r in records.values()),
              "experiment": cfg, "all_generated_sha256": sha256(root / "all_generated.csv"),
              "task_results": results}
    write_json(report_path, report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--experiment", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = validate(args.run_dir, args.experiment, args.report)
    print(f"PASS: {report['tasks']} tasks, {report['rows']} attempts; budgets, rewards, identities and checkpoints verified")


if __name__ == "__main__":
    main()
