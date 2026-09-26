"""Exercise preparation, repeat semantics, CSV accounting, and child lifecycle.

The subprocess fixture speaks the production startup/CSV protocols but never
loads PLANET or a generator; it is only available inside this test file.
"""

import csv
import json
import os
from pathlib import Path
import socket
import sys
import textwrap
import tomllib
from types import SimpleNamespace

import numpy as np
import pytest

from prepare_manual_planet_rl import ROOT, crystal_center, prepare, sha256, write_json
from run_manual_planet_rl import FIELDS, collect_task, run


FAKE_CHILD = '''
import csv, json, sys, time, tomllib
from pathlib import Path
from http.server import BaseHTTPRequestHandler, HTTPServer

module = sys.argv[2]
if module == 'planet_oracle.server':
    cfg = json.loads(Path(sys.argv[sys.argv.index('--config') + 1]).read_text())
    oid = 'fixture-' + cfg['target_id']
    manifest = {'oracle_id': oid, 'config': cfg}
    Path(sys.argv[sys.argv.index('--manifest') + 1]).write_text(json.dumps(manifest))
    class Handler(BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.end_headers()
            self.wfile.write(json.dumps({'version': 1, 'ready': True,
                'target_id': cfg['target_id'], 'oracle_id': oid}).encode())
    HTTPServer(('127.0.0.1', cfg['port']), Handler).serve_forever()
elif module == 'manual_planet_rl':
    cfg = tomllib.loads(Path(sys.argv[-1]).read_text())
    folder = Path(cfg['stage'][0]['chkpt_file']).parent
    client = json.loads((folder.parent / 'client.json').read_text())
    header = ['step','SMILES','SMILES_state','Score','Input_SMILES','Generated_SMILES',
        'PLANET reward','planet_affinity (PLANET reward)','planet_status (PLANET reward)',
        'planet_error (PLANET reward)','planet_target_id (PLANET reward)',
        'planet_oracle_id (PLANET reward)']
    seen = set()
    with (folder / 'summary_1.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        for step in range(1, cfg['stage'][0]['max_steps'] + 1):
            for index in range(cfg['parameters']['batch_size']):
                smi = 'invalid' if index == 2 else 'CCO'
                state = 0 if index == 2 else (2 if index == 1 else 1)
                reward = 0 if state == 0 else .7
                final = reward * (.5 if smi in seen else 1)
                if state != 0:
                    seen.add(smi)
                writer.writerow(dict(zip(header, [step,smi,state,final,'*C','*O',reward,
                    'None' if state == 0 else 11.0, 'None' if state == 0 else 'ok',
                    'None' if state == 0 else '', client['target_id'],client['oracle_id']])))
    (folder / 'agent.chkpt').write_bytes(b'fixture checkpoint')
else:
    raise RuntimeError(module)
'''


def build_run(tmp_path, target_names):
    """Prepare a fixture run whose targets share one port, like production."""
    project = tmp_path / "project"
    inputs = project / "configs/manual_design"
    for mode in ("libinvent", "linkinvent"):
        prior = project / "REINVENT4/priors" / f"{mode}_transformer_pubchem.prior"
        prior.parent.mkdir(parents=True, exist_ok=True)
        prior.write_bytes(b"fixture prior")
    checkpoint = project / "PLANET/PLANET.param"
    checkpoint.parent.mkdir(parents=True, exist_ok=True)
    checkpoint.write_bytes(b"fixture PLANET")
    tasks = []
    for target in target_names:
        ligand_dir = project / "real-world_dataset" / target
        ligand_dir.mkdir(parents=True)
        ligand = ligand_dir / "crystal.mol2"
        ligand.write_text("@<TRIPOS>MOLECULE\nligand\n3 0 0\nSMALL\nNO_CHARGES\n"
                          "@<TRIPOS>ATOM\n1 C 0 0 0 C.3\n2 N 3 0 0 N.3\n3 H 99 0 0 H\n")
        (ligand_dir / "receptor_out.pdb").write_text("ATOM fixture\n")
        (inputs / target).mkdir(parents=True)
        for mode, smi in (("libinvent", "*C"), ("linkinvent", "*C|*N")):
            (inputs / target / f"{mode}.smi").write_text(smi + "\n")
            tasks.append({"target": target, "mode": mode, "status": "designed",
                          "smiles_file": f"{target}/{mode}.smi", "input_smiles": smi,
                          "source_mol2": f"real-world_dataset/{target}/crystal.mol2",
                          "source_sha256": sha256(ligand),
                          "validation": {"connectivity_match": True, "stereochemistry_match": False}})
    write_json(inputs / "manifest.json", {"schema_version": 1, "model_ids": {}, "tasks": tasks,
               "status_counts": {"designed": len(tasks), "skipped": 0, "incomplete": 0,
                                 "not_started": 0}})
    cfg = json.loads((ROOT / "configs/manual_planet_rl/experiment.json").read_text())
    cfg.update(steps=2, batch_size=3, startup_timeout=10)
    experiment = project / "experiment.json"
    write_json(experiment, cfg)
    executable = tmp_path / "fixture-python"
    executable.write_text(f"#!{sys.executable}\n" + textwrap.dedent(FAKE_CHILD))
    executable.chmod(0o755)
    with socket.socket() as probe:
        probe.bind(("127.0.0.1", 0))
        port = probe.getsockname()[1]
    directory = tmp_path / "run with spaces"
    manifest = prepare(project, inputs / "manifest.json", experiment, directory,
                       executable, executable, "cpu", port)
    return directory, manifest, port


@pytest.fixture
def prepared(tmp_path):
    return build_run(tmp_path, ["example"])


def test_prepare_budget_center_no_early_stop(prepared):
    directory, manifest, _ = prepared
    assert manifest["targets"][0]["center"] == [1.5, 0, 0]
    task = directory / "example/libinvent"
    cfg = tomllib.loads((task / "rl.toml").read_text())
    assert cfg["parameters"]["prior_file"] == cfg["parameters"]["agent_file"]
    assert cfg["stage"][0]["termination"] == "null"
    assert cfg["diversity_filter"]["penalty_multiplier"] == .5
    assert cfg["stage"][0]["max_steps"] * cfg["parameters"]["batch_size"] == 6
    assert manifest["targets"][0]["tasks"][0]["validation"]["stereochemistry_match"] is False


def test_runner_retains_duplicates_invalids_and_stops_children(prepared):
    directory, _, port = prepared
    assert run(directory)["state"] == "complete"
    with (directory / "all_generated.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 12 and set(rows[0]) == set(FIELDS)
    assert sum(row["smiles_state"] == "0" for row in rows) == 4
    assert sum(row["is_repeat"] == "1" for row in rows) == 6
    assert rows[1]["training_reward"] == "0.35"
    assert rows[2]["planet_affinity"] == "" and rows[2]["planet_status"] == "not_scored"
    assert rows[3]["training_reward"] == "0.35"  # repeat across batches
    assert rows[6]["training_reward"] == "0.7"  # fresh memory for the other mode
    # No active listener may remain; TIME_WAIT leftovers must not fail this bind.
    with socket.socket() as probe:
        probe.setsockopt(socket.SOL_SOCKET, socket.SO_REUSEADDR, 1)
        probe.bind(("127.0.0.1", port))
    with pytest.raises(FileExistsError):
        run(directory)


def test_sequential_targets_share_one_port(tmp_path):
    # Each stopped server leaves TIME_WAIT sockets on the shared port; the next
    # target's probe and server must still bind without waiting for them to age.
    directory, manifest, _ = build_run(tmp_path, ["alpha", "beta"])
    assert [entry["target"] for entry in manifest["targets"]] == ["alpha", "beta"]
    assert run(directory)["state"] == "complete"
    with (directory / "all_generated.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 24
    assert {row["target"] for row in rows} == {"alpha", "beta"}
    for target in ("alpha", "beta"):
        oracle = json.loads((directory / target / "oracle-manifest.json").read_text())
        assert oracle["oracle_id"] == f"fixture-{target}"


def test_input_tampering_rejected_before_process_start(prepared):
    directory, _, _ = prepared
    (directory / "example/libinvent/input.smi").write_text("*N\n")
    with pytest.raises(ValueError, match="changed"):
        run(directory)
    assert not (directory / ".started").exists()


def test_incomplete_budget_is_not_exported(prepared):
    directory, manifest, _ = prepared
    run(directory)
    entry = manifest["targets"][0]
    task = entry["tasks"][0]
    folder = directory / task["directory"]
    lines = (folder / "summary_1.csv").read_text().splitlines()
    (folder / "summary_1.csv").write_text("\n".join(lines[:-1]) + "\n")
    with pytest.raises(ValueError, match="Incomplete attempt budget"):
        collect_task(directory, manifest, entry, task)


def test_same_batch_repeats_receive_half_reward_and_enter_update():
    from manual_planet_rl import LibinventRepeatLearning
    from reinvent.models.model_factory.sample_batch import SmilesState
    from reinvent.runmodes.RL.memories.penalize_same_smiles import PenalizeSameSmiles

    learning = object.__new__(LibinventRepeatLearning)
    learning.sampled = SimpleNamespace(smilies=["CCO", "CCO", "invalid"],
        states=np.array([SmilesState.VALID, SmilesState.DUPLICATE, SmilesState.INVALID]))
    learning.invalid_mask = np.array([True, True, False])
    learning.duplicate_mask = np.array([True, False, True])

    def scorer(smiles, invalid, duplicates, **kwargs):
        return SimpleNamespace(total_scores=.8 * (invalid & duplicates), smilies=smiles)

    learning.scoring_function = scorer
    results = learning.score()
    df = PenalizeSameSmiles(bucket_size=25, minscore=.4, minsimilarity=.4,
                           penalty_multiplier=.5, rdkit_smiles_flags={"isomericSmiles": True})
    df.update_score(results.total_scores, results.smilies, learning.invalid_mask)
    np.testing.assert_allclose(results.total_scores, [.8, .4, 0])
    adapter = SimpleNamespace(likelihood_smiles=lambda sampled: SimpleNamespace(likelihood="nll"))
    learning._state = SimpleNamespace(agent=adapter)
    learning.prior, learning.inception = adapter, None
    learning.reward_nlls = lambda *args: args[4]
    np.testing.assert_array_equal(learning._update_common(results, []), [0, 1])


def test_stable_reporter_keeps_metadata_columns_when_all_invalid(caplog):
    from manual_planet_rl import StableCSVReporter
    from reinvent.models.model_factory.sample_batch import SmilesState

    result = SimpleNamespace(fetch_scores=lambda *args, **kwargs: [[0]],
                             fetch_metadata=lambda *args: {})
    component = SimpleNamespace(component_names=["PLANET reward"], component_result=result,
                                transformed_scores=[[0]])
    data = SimpleNamespace(score_results=SimpleNamespace(completed_components=[component],
        smilies=["invalid"], total_scores=[0]), agent_nll=[1], prior_nll=[1], augmented_nll=[1],
        sampled=SimpleNamespace(states=[SmilesState.INVALID], items1=["*C"], items2=["invalid"]),
        step=1, mean_score=0, fraction_valid_smiles=0, fraction_duplicate_smiles=0)
    with caplog.at_level("INFO", logger="csv"):
        StableCSVReporter().submit(data)
    records = [record for record in caplog.records if record.name == "csv"]
    assert "planet_affinity (PLANET reward)" in records[0].msg
    assert len(records[0].msg) == len(records[1].msg)
