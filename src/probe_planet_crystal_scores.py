"""Sample paired PDB/SDF entries from a local archive and probe raw PLANET scores."""

import argparse
import csv
import json
import random
import tarfile
import tempfile
import time
from pathlib import Path

import numpy as np
import torch
from rdkit import Chem

from planet_oracle.backend import PlanetPredictor, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--planet-root", type=Path, required=True)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--count", type=int, default=20)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.count <= 0:
        parser.error("count must be positive")
    args.output.mkdir(parents=True, exist_ok=False)
    torch.set_num_threads(1)
    root = args.planet_root.resolve()
    rows, manifests = [], {}
    started = time.monotonic()
    with tarfile.open(args.archive) as archive, tempfile.TemporaryDirectory() as tmp:
        members = {m.name: m for m in archive.getmembers() if m.isfile()}
        pairs = sorted(name for name in members if name.endswith("_protein.pdb")
                       and name.replace("_protein.pdb", "_ligand.sdf") in members)
        chosen = sorted(random.Random(args.seed).sample(pairs, min(args.count, len(pairs))))
        for name in chosen:
            pdb_id = Path(name).name.removesuffix("_protein.pdb")
            row = {"pdb_id": pdb_id, "smiles": "", "heavy_atoms": "", "planet_score": "",
                   "sigmoid_reward": "", "seconds": "", "error": ""}
            tick = time.monotonic()
            try:
                paths = []
                for member in (name, name.replace("_protein.pdb", "_ligand.sdf")):
                    path = Path(tmp) / Path(member).name
                    path.write_bytes(archive.extractfile(members[member]).read())
                    paths.append(path)
                ligand = Chem.SDMolSupplier(str(paths[1]))[0]
                if ligand is None:
                    raise ValueError("Cannot parse crystal ligand")
                row["smiles"] = Chem.MolToSmiles(Chem.RemoveHs(ligand))
                row["heavy_atoms"] = ligand.GetNumHeavyAtoms()
                cfg = {"planet_root": str(root), "checkpoint": str(root / "PLANET.param"),
                       "device": "cpu", "seed": args.seed, "protein_pdb": str(paths[0]),
                       "ligand_sdf": str(paths[1])}
                predictor = PlanetPredictor(cfg)
                manifests[pdb_id] = predictor.manifest
                molecule, error = predictor.prepare(row["smiles"])
                if error:
                    raise ValueError(error)
                score = predictor.predict([molecule])[0]
                if not np.isfinite(score):
                    raise ValueError("Non-finite score")
                row["planet_score"] = score
                row["sigmoid_reward"] = float(1 / (1 + 10 ** (-2.5 / 12 * (score - 10))))
            except Exception as exc:
                row["error"] = f"{type(exc).__name__}: {exc}"
            row["seconds"] = time.monotonic() - tick
            rows.append(row)
            print(json.dumps(row), flush=True)
    with (args.output / "scores.csv").open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    values = [r["planet_score"] for r in rows if not r["error"]]
    report = {"archive_sha256": sha256(args.archive), "seed": args.seed,
              "available_pairs": len(pairs), "attempted": len(rows), "success": len(values),
              "elapsed_seconds": time.monotonic() - started,
              "quantiles_min_q25_median_q75_max": np.quantile(values, [0, .25, .5, .75, 1]).tolist() if values else [],
              "mean": float(np.mean(values)) if values else None,
              "sigmoid": {"low": 4, "high": 16, "k": .25}, "manifests": manifests}
    (args.output / "report.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "manifests"}), flush=True)


if __name__ == "__main__":
    main()
