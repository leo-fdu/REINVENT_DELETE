"""Check finite likelihoods and learned weights in our own trusted RL checkpoints."""

import argparse
import csv
import json
import math
from pathlib import Path

import torch


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    torch.set_num_threads(2)
    manifest = json.loads((args.run_dir / "run.json").read_text())
    priors, results = {}, {}
    for entry in manifest["targets"]:
        for task in entry["tasks"]:
            mode = task["mode"]
            if mode not in priors:
                path = Path(manifest["project"]) / "REINVENT4/priors" / f"{mode}_transformer_pubchem.prior"
                priors[mode] = torch.load(path, map_location="cpu", weights_only=False)["network_state"]
            folder = args.run_dir / task["directory"]
            checkpoint = torch.load(folder / "agent.chkpt", map_location="cpu", weights_only=False)
            weights = checkpoint["network_state"]
            if set(weights) != set(priors[mode]):
                raise ValueError(f"Checkpoint architecture mismatch: {task['directory']}")
            if not all(bool(torch.isfinite(tensor).all()) for tensor in weights.values()):
                raise ValueError(f"Non-finite learned parameters: {task['directory']}")
            changed = sum(not torch.equal(tensor, priors[mode][name]) for name, tensor in weights.items())
            if changed == 0:
                raise ValueError(f"Agent did not change from prior: {task['directory']}")
            rows = 0
            with (folder / "summary_1.csv").open(newline="") as handle:
                for row in csv.DictReader(handle):
                    if not all(math.isfinite(float(row[name])) for name in ("Agent", "Prior", "Target", "Score")):
                        raise ValueError(f"Non-finite training likelihood/reward: {task['directory']}")
                    rows += 1
            results[task["directory"]] = {"finite_parameters": True, "changed_tensors": changed,
                                          "total_tensors": len(weights), "finite_training_rows": rows}
    args.report.write_text(json.dumps({"passed": True, "tasks": results}, indent=2) + "\n")
    print(f"PASS: {len(results)} checkpoints have finite, updated parameters and finite training likelihoods")


if __name__ == "__main__":
    main()
