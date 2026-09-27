"""Validate every prepared target with real PLANET, without starting RL."""

import argparse
import gc
import json
from pathlib import Path
import sys

from prepare_manual_planet_rl import write_json
from run_manual_planet_rl import local_path, read_run, verify_prepared


def validate(run_dir, report_dir):
    # Imports stay lazy so --help works in a standard-library-only environment.
    import torch
    from planet_oracle.backend import PlanetPredictor
    from planet_oracle.config import server_config

    root, report_dir = Path(run_dir).resolve(), Path(report_dir).resolve()
    manifest = read_run(root)
    verify_prepared(root, manifest)
    if report_dir.is_relative_to(Path(manifest["project"])) or report_dir.exists():
        raise ValueError("Choose a new validation report directory outside the project")
    report_dir.mkdir(parents=True)
    report = {"run_dir": str(root), "total_tasks": manifest["total_tasks"],
              "started_rl": False, "targets": {}, "passed": False}
    for entry in manifest["targets"]:
        target = entry["target"]
        predictor = None
        try:
            cfg = server_config(local_path(root, entry["server_config"]))
            predictor = PlanetPredictor(cfg)
            molecule, error = predictor.prepare("CCO")
            if error:
                raise RuntimeError(error)
            affinity = predictor.predict([molecule])[0]
            if not torch.isfinite(torch.tensor(affinity)):
                raise RuntimeError("Non-finite probe affinity")
            write_json(report_dir / f"{target}-oracle.json", predictor.manifest)
            report["targets"][target] = {
                "passed": True, "pocket": predictor.pocket_info,
                "probe_smiles": "CCO", "probe_affinity": affinity,
                "oracle_id": predictor.oracle_id,
                "pocket_reference": entry.get("pocket_reference"),
                "input_correction": entry.get("input_correction"),
            }
            print(f"{target}: PASS ({predictor.pocket_info['residue_count']} residues)", flush=True)
        except Exception as error:
            # Keep all diagnostics, but never turn a failed pocket into a score.
            report["targets"][target] = {"passed": False,
                                         "error": f"{type(error).__name__}: {error}"}
            print(f"{target}: FAIL ({error})", flush=True)
        finally:
            del predictor
            gc.collect()
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
        write_json(report_dir / "report.json", report)
    report["passed"] = all(result["passed"] for result in report["targets"].values())
    write_json(report_dir / "report.json", report)
    return report["passed"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--report-dir", type=Path, required=True)
    args = parser.parse_args()
    try:
        passed = validate(args.run_dir, args.report_dir)
    except (ValueError, KeyError, OSError) as error:
        parser.exit(1, f"Validation failed: {error}\n")
    sys.exit(0 if passed else 1)


if __name__ == "__main__":
    main()
