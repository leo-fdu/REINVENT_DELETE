"""Run prepared tasks sequentially and retain every sampling attempt."""

import argparse
from collections import Counter
import csv
import json
import math
import os
from pathlib import Path
import signal
import socket
import subprocess
import time
from urllib.error import URLError
from urllib.request import ProxyHandler, build_opener

from prepare_manual_planet_rl import sha256, write_json


FIELDS = ["target", "mode", "seed", "step", "sample_index", "input_smiles",
          "generated_smiles", "smiles", "smiles_state", "is_repeat",
          "planet_affinity", "planet_status", "planet_error", "sigmoid_reward",
          "training_reward", "planet_target_id", "oracle_id"]


def local_path(root, relative):
    path = (root / relative).resolve()
    if not path.is_relative_to(root) or path == root:
        raise ValueError(f"Path escapes run directory: {relative}")
    return path


def read_run(root):
    manifest = json.loads((root / "run.json").read_text())
    if manifest["schema_version"] != 1 or manifest["run_dir"] != str(root):
        raise ValueError("Run directory was moved or its manifest is invalid; prepare it again")
    return manifest


def verify_prepared(root, manifest):
    for filename, digest in manifest["input_hashes"].items():
        if sha256(filename) != digest:
            raise ValueError(f"Source input changed: {filename}")
    for filename, digest in manifest["prepared_hashes"].items():
        if sha256(local_path(root, filename)) != digest:
            raise ValueError(f"Prepared configuration/input changed: {filename}")


def stop_process(process):
    if process is None:
        return
    # Every child starts a new session. Its subprocesses (including scoring
    # clients) belong to this group, not to the user's shell or another run.
    try:
        os.killpg(process.pid, signal.SIGTERM)
    except ProcessLookupError:
        pass
    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()


def bind_oracle(process, root, entry, startup_timeout):
    client_path = local_path(root, entry["client_config"])
    client = json.loads(client_path.read_text())
    deadline = time.monotonic() + startup_timeout
    opener = build_opener(ProxyHandler({}))
    while time.monotonic() < deadline:
        if process.poll() is not None:
            raise RuntimeError(f"PLANET exited for {entry['target']}; see server.log")
        try:
            with opener.open(client["url"] + "/health", timeout=2) as response:
                health = json.load(response)
        except (URLError, TimeoutError):
            time.sleep(0.5)
            continue
        oracle = json.loads(local_path(root, entry["oracle_manifest"]).read_text())
        server = json.loads(local_path(root, entry["server_config"]).read_text())
        if (health.get("version") != 1 or health.get("ready") is not True
                or health.get("target_id") != entry["target"]
                or not health.get("oracle_id") or health["oracle_id"] != oracle["oracle_id"]
                or oracle["config"] != server):
            raise RuntimeError("PLANET service identity/configuration mismatch")
        client["oracle_id"] = health["oracle_id"]
        write_json(client_path, client)
        return
    raise TimeoutError(f"PLANET did not start within {startup_timeout}s: {entry['target']}")


def numeric(value, name):
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Non-finite {name}")
    return result


def collect_task(root, manifest, entry, task):
    folder = local_path(root, task["directory"])
    cfg = manifest["experiment"]
    oracle = json.loads(local_path(root, entry["oracle_manifest"]).read_text())
    counts, seen, output = Counter(), set(), []
    with (folder / "summary_1.csv").open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"step", "SMILES", "SMILES_state", "Score", "Input_SMILES", "Generated_SMILES",
                    "PLANET reward", "planet_affinity (PLANET reward)",
                    "planet_status (PLANET reward)", "planet_error (PLANET reward)",
                    "planet_target_id (PLANET reward)", "planet_oracle_id (PLANET reward)"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Expected CSV produced by the manual_planet_rl entry point")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("Malformed or truncated REINVENT CSV row")
            step, state = int(row["step"]), int(row["SMILES_state"])
            if state not in (0, 1, 2):
                raise ValueError("Unknown SMILES state")
            counts[step] += 1
            smiles = row["SMILES"]
            repeat = state != 0 and smiles in seen
            if state != 0:
                seen.add(smiles)
            status = row["planet_status (PLANET reward)"]
            affinity = row["planet_affinity (PLANET reward)"]
            error = row["planet_error (PLANET reward)"]
            target_id = row["planet_target_id (PLANET reward)"]
            oracle_id = row["planet_oracle_id (PLANET reward)"]
            if status in ("", "None", "null"):
                status, affinity, error, target_id, oracle_id = "not_scored", "", "", "", ""
            elif status == "ok":
                affinity = numeric(affinity, "PLANET affinity")
                if target_id != entry["target"] or oracle_id != oracle["oracle_id"]:
                    raise ValueError("CSV contains a score from a different PLANET oracle")
            elif status == "invalid_input":
                if affinity not in ("", "None", "null"):
                    raise ValueError("Invalid molecule has a real affinity in CSV")
                affinity = ""
            else:
                raise ValueError(f"Unknown PLANET status: {status}")
            reward = numeric(row["PLANET reward"], "sigmoid reward")
            final_reward = numeric(row["Score"], "training reward")
            if not 0 <= reward <= 1 or not 0 <= final_reward <= 1:
                raise ValueError("Reward outside [0, 1]")
            expected_reward = reward * (cfg["penalty_multiplier"] if repeat else 1)
            if not math.isclose(final_reward, expected_reward, rel_tol=1e-5, abs_tol=1e-7):
                raise ValueError("Repeat reward does not match the configured multiplier")
            if state == 0 and (status == "ok" or final_reward != 0):
                raise ValueError("Invalid generated molecule has a successful score or reward")
            output.append(dict(zip(FIELDS, [
                entry["target"], task["mode"], cfg["seed"], step, counts[step],
                row["Input_SMILES"], row["Generated_SMILES"], smiles, state, int(repeat),
                affinity, status, error, reward, final_reward, target_id, oracle_id,
            ])))
    expected_counts = Counter({step: cfg["batch_size"] for step in range(1, cfg["steps"] + 1)})
    if counts != expected_counts:
        raise ValueError(f"Incomplete attempt budget for {entry['target']}/{task['mode']}: {dict(counts)}")
    with (folder / "all_generated.csv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(output)
    return len(output)


def combine(root, entries):
    paths = [local_path(root, task["directory"]) / "all_generated.csv"
             for entry in entries for task in entry["tasks"]]
    if not all(path.is_file() for path in paths):
        raise ValueError("Some tasks have no collected output")
    with (root / "all_generated.csv").open("x", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for path in paths:
            with path.open(newline="", encoding="utf-8") as source:
                reader = csv.DictReader(source)
                if reader.fieldnames != FIELDS:
                    raise ValueError(f"Unexpected collected CSV columns: {path}")
                writer.writerows(reader)


def run(root, targets=None):
    root = Path(root).resolve()
    manifest = read_run(root)
    entries = manifest["targets"]
    if targets:
        if len(set(targets)) != len(targets) or set(targets) - {e["target"] for e in entries}:
            raise ValueError("Unknown or duplicate target selection")
        entries = [entry for entry in entries if entry["target"] in targets]
    # A claimed directory stays claimed. Report it before the hash check: a run
    # legitimately rewrites client.json (oracle_id), so verifying a used
    # directory would misreport the runner's own change as tampering.
    if (root / ".started").exists():
        raise FileExistsError("Run directory was already used; prepare a new one")
    verify_prepared(root, manifest)
    cfg = manifest["experiment"]
    project = Path(manifest["project"])
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(project / "src"), str(project / "REINVENT4"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    env["PYTHONHASHSEED"] = str(cfg["seed"])
    state = {"state": "running", "targets": [e["target"] for e in entries], "tasks": {}}
    # This marker permanently claims the run, including failed/interrupted runs.
    with (root / ".started").open("x") as handle:
        handle.write(f"{os.getpid()}\n")
    write_json(root / "status.json", state)
    server_process = rl_process = None
    try:
        for entry in entries:
            folder = root / entry["target"]
            server = json.loads(local_path(root, entry["server_config"]).read_text())
            # Do not connect to or replace an already running service on this port.
            # SO_REUSEADDR keeps TIME_WAIT sockets left by the previous target's
            # stopped server from blocking the bind; an active listener still does.
            with socket.socket() as probe:
                probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
                probe.bind(("127.0.0.1", server["port"]))
            command = [manifest["planet_python"], "-m", "planet_oracle.server",
                       "--config", str(local_path(root, entry["server_config"])),
                       "--manifest", str(local_path(root, entry["oracle_manifest"]))]
            with (folder / "server.log").open("x") as server_log:
                server_process = subprocess.Popen(
                    command, cwd=project, env=env, stdout=server_log,
                    stderr=subprocess.STDOUT, start_new_session=True,
                )
                try:
                    bind_oracle(server_process, root, entry, cfg["startup_timeout"])
                    for task in entry["tasks"]:
                        name = f"{entry['target']}/{task['mode']}"
                        task_dir = local_path(root, task["directory"])
                        state["tasks"][name] = "running"
                        write_json(root / "status.json", state)
                        print(f"Starting {name}: {manifest['attempts_per_task']} attempts", flush=True)
                        command = [manifest["reinvent_python"], "-m", "manual_planet_rl",
                                   "-l", str(task_dir / "run.log"), str(task_dir / "rl.toml")]
                        with (task_dir / "process.log").open("x") as process_log:
                            rl_process = subprocess.Popen(
                                command, cwd=project, env=env, stdout=process_log,
                                stderr=subprocess.STDOUT, start_new_session=True,
                            )
                            return_code = rl_process.wait()
                        if return_code != 0:
                            raise RuntimeError(f"RL failed ({return_code}): {name}; see run.log/process.log")
                        rl_process = None
                        if server_process.poll() is not None:
                            raise RuntimeError(f"PLANET unexpectedly stopped: {entry['target']}")
                        checkpoint = task_dir / "agent.chkpt"
                        if not checkpoint.is_file() or checkpoint.stat().st_size == 0:
                            raise RuntimeError(f"RL did not save its final agent: {name}")
                        collect_task(root, manifest, entry, task)
                        state["tasks"][name] = "complete"
                        write_json(root / "status.json", state)
                finally:
                    stop_process(rl_process)
                    rl_process = None
                    stop_process(server_process)
                    server_process = None
        combine(root, entries)
        state["state"] = "complete"
    except BaseException as error:
        state.update(state="interrupted" if isinstance(error, KeyboardInterrupt) else "failed",
                     error=f"{type(error).__name__}: {error}")
        state["tasks"] = {name: "interrupted" if value == "running" else value
                          for name, value in state["tasks"].items()}
        raise
    finally:
        stop_process(rl_process)
        stop_process(server_process)
        write_json(root / "status.json", state)
    return state


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--targets", nargs="+", help="Optional target subset; default: all prepared targets")
    args = parser.parse_args()

    def interrupted(signum, frame):
        raise KeyboardInterrupt(f"Received signal {signum}")

    previous = signal.signal(signal.SIGTERM, interrupted)
    try:
        run(args.run_dir, args.targets)
    except KeyboardInterrupt:
        parser.exit(130, "Run interrupted; partial outputs are retained.\n")
    except (ValueError, KeyError, OSError, RuntimeError) as error:
        parser.exit(1, f"Run failed: {error}\n")
    finally:
        signal.signal(signal.SIGTERM, previous)
    print(f"All selected tasks complete: {args.run_dir / 'all_generated.csv'}")


if __name__ == "__main__":
    main()
