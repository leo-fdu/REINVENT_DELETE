"""Prepare portable manual-design PLANET RL runs without loading either model."""

import argparse
from collections import Counter
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shlex


ROOT = Path(__file__).resolve().parents[1]
MODES = ("libinvent", "linkinvent")


def sha256(path):
    with Path(path).open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def write_json(path, value):
    Path(path).write_text(
        json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def read_experiment(path):
    cfg = json.loads(Path(path).read_text(encoding="utf-8"))
    keys = {"schema_version", "seed", "batch_size", "steps", "sample_strategy",
            "temperature", "learning_strategy", "sigmoid", "penalty_multiplier",
            "planet_batch_size", "request_timeout", "startup_timeout"}
    if set(cfg) != keys or cfg["schema_version"] != 1:
        raise ValueError("Unexpected experiment configuration fields or schema")
    for key in ("seed", "batch_size", "steps", "planet_batch_size"):
        low = 0 if key == "seed" else 1
        if type(cfg[key]) is not int or not low <= cfg[key] < 2**32:
            raise ValueError(f"Invalid {key}")
    for key in ("temperature", "request_timeout", "startup_timeout"):
        if type(cfg[key]) not in (int, float) or not math.isfinite(cfg[key]) or cfg[key] <= 0:
            raise ValueError(f"Invalid {key}")
    penalty = cfg["penalty_multiplier"]
    if type(penalty) not in (int, float) or not math.isfinite(penalty) or not 0 <= penalty <= 1:
        raise ValueError("penalty_multiplier must be between 0 and 1")
    learning = cfg["learning_strategy"]
    if set(learning) != {"type", "sigma", "rate"} or learning["type"] != "dap":
        raise ValueError("This experiment uses DAP")
    if type(learning["sigma"]) is not int or learning["sigma"] <= 0:
        raise ValueError("sigma must be a positive integer")
    if type(learning["rate"]) not in (int, float) or not math.isfinite(learning["rate"]) or learning["rate"] <= 0:
        raise ValueError("rate must be positive and finite")
    curve = cfg["sigmoid"]
    if set(curve) != {"low", "high", "k"} or any(
        type(v) not in (int, float) or not math.isfinite(v) for v in curve.values()
    ) or curve["high"] <= curve["low"] or curve["k"] <= 0:
        raise ValueError("Invalid increasing sigmoid")
    if cfg["sample_strategy"] != "multinomial":
        raise ValueError("Fixed attempt budgets require multinomial sampling")
    return cfg


def crystal_center(path):
    """Geometric centroid of original MOL2 heavy atoms, without coordinate edits."""
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    if lines.count("@<TRIPOS>MOLECULE") != 1:
        raise ValueError(f"Expected one crystal ligand: {path}")
    molecule = lines.index("@<TRIPOS>MOLECULE")
    expected = int(lines[molecule + 2].split()[0])
    start = lines.index("@<TRIPOS>ATOM") + 1
    coordinates, ids = [], set()
    for line in lines[start:]:
        if line.startswith("@<TRIPOS>"):
            break
        if not line.strip():
            continue
        fields = line.split()
        atom_id = int(fields[0])
        xyz = [float(value) for value in fields[2:5]]
        if atom_id <= 0 or atom_id in ids or len(xyz) != 3 or not all(map(math.isfinite, xyz)):
            raise ValueError(f"Invalid MOL2 atom coordinates: {path}")
        ids.add(atom_id)
        if fields[5].split(".")[0] not in ("H", "D", "T"):
            coordinates.append(xyz)
    if len(ids) != expected or not coordinates:
        raise ValueError(f"Incomplete crystal ligand: {path}")
    return [math.fsum(xyz[axis] for xyz in coordinates) / len(coordinates) for axis in range(3)]


def interpreter(path):
    # Do not dereference a venv's Python symlink: that would lose its environment.
    path = Path(os.path.abspath(Path(path).expanduser()))
    if not path.is_file() or not os.access(path, os.X_OK):
        raise ValueError(f"Python executable does not exist or is not executable: {path}")
    return str(path)


def rl_toml(cfg, mode, project, task_dir, client, python, device):
    quote = json.dumps
    prior = project / "REINVENT4" / "priors" / f"{mode}_transformer_pubchem.prior"
    learning, curve = cfg["learning_strategy"], cfg["sigmoid"]
    args = shlex.join(["-m", "planet_oracle.client", "--config", str(client)])
    return f'''run_type = "staged_learning"
device = {quote(device)}
seed = {cfg["seed"]}
tb_logdir = {quote(str(task_dir / "tensorboard"))}
json_out_config = {quote(str(task_dir / "resolved.json"))}

[parameters]
prior_file = {quote(str(prior))}
agent_file = {quote(str(prior))}
smiles_file = {quote(str(task_dir / "input.smi"))}
summary_csv_prefix = {quote(str(task_dir / "summary"))}
use_checkpoint = false
purge_memories = false
batch_size = {cfg["batch_size"]}
sample_strategy = "multinomial"
temperature = {cfg["temperature"]}
randomize_smiles = false
isomeric_smiles = true

[learning_strategy]
type = "dap"
sigma = {learning["sigma"]}
rate = {learning["rate"]}

[diversity_filter]
type = "PenalizeSameSmiles"
penalty_multiplier = {cfg["penalty_multiplier"]}

[[stage]]
chkpt_file = {quote(str(task_dir / "agent.chkpt"))}
termination = "null"
max_score = 1.0
min_steps = 0
max_steps = {cfg["steps"]}

[stage.scoring]
type = "arithmetic_mean"
parallel = 1
use_pumas = false

[[stage.scoring.component]]
[stage.scoring.component.ExternalProcess]
[[stage.scoring.component.ExternalProcess.endpoint]]
name = "PLANET reward"
weight = 1.0
params.executable = {quote(python)}
params.args = {quote(args)}
params.property = "planet_affinity_for_scoring"
transform.type = "sigmoid"
transform.low = {curve["low"]}
transform.high = {curve["high"]}
transform.k = {curve["k"]}
'''


def prepare(project, input_manifest, experiment, run_dir, reinvent_python,
            planet_python, device="cuda:0", port=8765):
    project, run_dir = Path(project).resolve(), Path(run_dir).resolve()
    if run_dir.is_relative_to(project) or run_dir.exists():
        raise ValueError("Choose a new run directory outside the project")
    if type(port) is not int or not 1 <= port <= 65535:
        raise ValueError("Invalid port")
    if not re.fullmatch(r"cpu|cuda(?::[0-9]+)?", device):
        raise ValueError("device must be cpu, cuda, or cuda:N")
    reinvent_python, planet_python = interpreter(reinvent_python), interpreter(planet_python)
    cfg = read_experiment(experiment)
    input_manifest = Path(input_manifest).resolve()
    inputs = json.loads(input_manifest.read_text(encoding="utf-8"))
    if inputs["schema_version"] != 1:
        raise ValueError("Unexpected input manifest schema")
    records, skipped, seen = {}, [], set()
    counts = Counter()
    for task in inputs["tasks"]:
        target, mode, status = task["target"], task["mode"], task["status"]
        if not re.fullmatch(r"[a-zA-Z0-9_-]+", target) or mode not in MODES:
            raise ValueError("Invalid target or mode")
        if (target, mode) in seen or status not in ("designed", "skipped"):
            raise ValueError("Duplicate or unfinished input task")
        seen.add((target, mode))
        counts[status] += 1
        if status == "skipped":
            skipped.append(task)
            continue
        expected_source = f"real-world_dataset/{target}/crystal.mol2"
        if task["source_mol2"] != expected_source:
            raise ValueError(f"Unexpected source ligand: {target}")
        source = project / expected_source
        if sha256(source) != task["source_sha256"]:
            raise ValueError(f"Crystal ligand changed: {target}")
        if task["validation"].get("connectivity_match") is not True:
            raise ValueError(f"Unvalidated design: {target}/{mode}")
        smiles_path = (input_manifest.parent / task["smiles_file"]).resolve()
        if not smiles_path.is_relative_to(input_manifest.parent):
            raise ValueError("Input file must belong to the input manifest directory")
        smiles = task["input_smiles"]
        if not smiles or any(c.isspace() for c in smiles) or smiles_path.read_text() != smiles + "\n":
            raise ValueError(f"Input SMILES changed: {target}/{mode}")
        protein = project / "real-world_dataset" / target / "receptor_out.pdb"
        if not protein.is_file() or protein.stat().st_size == 0:
            raise ValueError(f"Missing protein: {target}")
        record = records.setdefault(target, {
            "target": target, "protein_pdb": str(protein), "source_mol2": str(source),
            "center": crystal_center(source), "tasks": [],
        })
        record["tasks"].append(task)
    actual = {s: counts[s] for s in ("designed", "skipped", "incomplete", "not_started")}
    if actual != inputs["status_counts"] or not records:
        raise ValueError("Input task counts are inconsistent or empty")
    model_paths = {mode: project / "REINVENT4/priors" / f"{mode}_transformer_pubchem.prior"
                   for mode in MODES}
    checkpoint = project / "PLANET/PLANET.param"
    source_files = [input_manifest, *model_paths.values(), checkpoint]
    source_files += [Path(record[key]) for record in records.values()
                     for key in ("protein_pdb", "source_mol2")]
    source_hashes = {str(path): sha256(path) for path in source_files}
    manifest = {
        "schema_version": 1, "project": str(project), "run_dir": str(run_dir),
        "reinvent_python": reinvent_python, "planet_python": planet_python,
        "experiment": cfg, "model_ids": inputs["model_ids"], "input_hashes": source_hashes,
        "prepared_hashes": {}, "targets": [], "skipped": skipped,
        "attempts_per_task": cfg["batch_size"] * cfg["steps"],
        "total_tasks": counts["designed"],
    }
    run_dir.mkdir(parents=True)
    write_json(run_dir / "experiment.json", cfg)
    for target, record in sorted(records.items()):
        folder = run_dir / target
        folder.mkdir()
        server = {
            "target_id": target, "planet_root": str(project / "PLANET"),
            "checkpoint": str(checkpoint), "protein_pdb": record["protein_pdb"],
            "center": record["center"], "device": device,
            "batch_size": cfg["planet_batch_size"], "port": port, "seed": cfg["seed"],
        }
        write_json(folder / "server.json", server)
        write_json(folder / "client.json", {
            "url": f"http://127.0.0.1:{port}", "target_id": target,
            "oracle_id": "UNBOUND_UNTIL_SERVER_STARTS", "timeout": cfg["request_timeout"],
        })
        entry = {key: value for key, value in record.items() if key != "tasks"}
        entry.update(server_config=f"{target}/server.json", client_config=f"{target}/client.json",
                     oracle_manifest=f"{target}/oracle-manifest.json", tasks=[])
        for task in sorted(record["tasks"], key=lambda t: MODES.index(t["mode"])):
            mode = task["mode"]
            task_dir = folder / mode
            task_dir.mkdir()
            (task_dir / "input.smi").write_text(task["input_smiles"] + "\n", encoding="utf-8")
            (task_dir / "rl.toml").write_text(
                rl_toml(cfg, mode, project, task_dir, folder / "client.json", reinvent_python, device),
                encoding="utf-8",
            )
            entry["tasks"].append({
                "mode": mode, "directory": f"{target}/{mode}",
                "validation": task["validation"], "input_smiles": task["input_smiles"],
            })
        manifest["targets"].append(entry)
    for path in sorted(run_dir.rglob("*")):
        if path.is_file():
            manifest["prepared_hashes"][str(path.relative_to(run_dir))] = sha256(path)
    write_json(run_dir / "run.json", manifest)
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project", type=Path, default=ROOT)
    parser.add_argument("--inputs", type=Path)
    parser.add_argument("--experiment", type=Path)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--reinvent-python", required=True)
    parser.add_argument("--planet-python", required=True)
    parser.add_argument("--device", default="cuda:0")
    parser.add_argument("--port", type=int, default=8765)
    args = parser.parse_args()
    project = args.project.resolve()
    try:
        manifest = prepare(
            project, args.inputs or project / "configs/manual_design/manifest.json",
            args.experiment or project / "configs/manual_planet_rl/experiment.json",
            args.run_dir, args.reinvent_python, args.planet_python, args.device, args.port,
        )
    except (ValueError, KeyError, IndexError, OSError) as error:
        parser.exit(1, f"Preparation failed: {error}\n")
    print(f"Prepared {manifest['total_tasks']} tasks, {manifest['attempts_per_task']} attempts each: {args.run_dir}")


if __name__ == "__main__":
    main()
