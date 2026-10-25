"""Exercise signed size rewards through scoring, repeats, DAP and exported audits."""

import csv
import json
import logging
from types import SimpleNamespace

import numpy as np
import pytest
import torch

from manual_planet_rl import LibinventRepeatLearning, LinkinventRepeatLearning, StableCSVReporter
from planet_size_reward import generated_heavy_atoms, reward_terms
from prepare_manual_planet_rl import ROOT, read_experiment, write_json
from reinvent.models.model_factory.sample_batch import SmilesState
from reinvent.runmodes.RL.memories.penalize_same_smiles import PenalizeSameSmiles
from reinvent.runmodes.RL.reward import RLReward
from reinvent.utils.logmon import CsvFormatter
from run_manual_planet_rl import collect_task, combine
from validate_manual_planet_results import validate
from test_workflow import build_run


CONFIG = ROOT / "configs/manual_planet_rl/experiment_high9_cap9_size001_1000.json"


@pytest.mark.parametrize("product,fixed,expected", [
    ("CCO", "*C", 2), ("CCN", "*C|*N", 1),
    ("CCO", "[H]C(*)", 2), ("c1ccccc1N", "*c1ccccc1", 1),
])
def test_only_generated_real_heavy_atoms(product, fixed, expected):
    assert generated_heavy_atoms(product, fixed) == expected


def test_capping_and_negative_reward():
    cfg = read_experiment(CONFIG)
    at9 = reward_terms(9, "CCO", "*C", cfg)
    assert reward_terms(99, "CCO", "*C", cfg) == at9
    assert at9["sigmoid_reward"] == pytest.approx(1 / (1 + 10 ** -1.25))
    assert at9["base_reward"] == pytest.approx(at9["sigmoid_reward"] - .02)
    assert reward_terms(0, "CCO", "*C", cfg)["base_reward"] < 0
    with pytest.raises(ValueError):
        generated_heavy_atoms("C", "*CCC")


@pytest.mark.parametrize("value", [-.01, float("nan"), True])
def test_reject_bad_penalty(tmp_path, value):
    cfg = json.loads(CONFIG.read_text())
    cfg["reward_adjustment"]["generated_heavy_atom_penalty"] = value
    p = tmp_path / "bad.json"
    p.write_text(json.dumps(cfg))
    with pytest.raises(ValueError):
        read_experiment(p)


def learning_fixture(mode, cfg):
    cls = LibinventRepeatLearning if mode == "libinvent" else LinkinventRepeatLearning
    learning = object.__new__(cls)
    learning.reward_experiment = cfg
    product, fixed = ("CCO", "*C") if mode == "libinvent" else ("CCN", "*C|*N")
    learning.sampled = SimpleNamespace(smilies=[product, product, "invalid"],
        states=np.array([SmilesState.VALID, SmilesState.DUPLICATE, SmilesState.INVALID]),
        items1=[fixed] * 3, items2=["*CO", "*CO", "invalid"])
    learning.invalid_mask = np.array([True, True, False])
    learning.duplicate_mask = np.array([True, False, True])

    def scorer(smiles, invalid, duplicates, *args, **kwargs):
        raw = [learning.affinity, learning.affinity, None]
        meta = {"planet_affinity": raw, "planet_status": ["ok", "ok", None],
                "planet_error": ["", "", None], "planet_target_id": ["example"] * 3,
                "planet_oracle_id": ["fixture-example"] * 3}
        component = SimpleNamespace(component_names=["PLANET reward"],
            transformed_scores=[np.array([.99, .99, 0.], dtype=np.float32)],
            component_result=SimpleNamespace(fetch_metadata=lambda *a: meta,
                                             fetch_scores=lambda *a, **kw: [raw]))
        return SimpleNamespace(smilies=smiles, total_scores=np.array([.99, .99, 0.], dtype=np.float32),
                               completed_components=[component])

    learning.scoring_function = scorer
    return learning


def test_both_modes_report_collect_and_independently_validate_signed_rewards(tmp_path):
    root, manifest, _ = build_run(tmp_path, ["example"])
    cfg = read_experiment(CONFIG)
    cfg.update(batch_size=3, steps=2)
    manifest["experiment"] = cfg
    write_json(root / "experiment.json", cfg)
    # This synthetic run tests the audit contract; no model is launched.
    from prepare_manual_planet_rl import sha256
    manifest["prepared_hashes"]["experiment.json"] = sha256(root / "experiment.json")
    write_json(root / "run.json", manifest)
    entry = manifest["targets"][0]
    write_json(root / entry["oracle_manifest"], {"oracle_id": "fixture-example"})
    for task in entry["tasks"]:
        learning = learning_fixture(task["mode"], cfg)
        df = PenalizeSameSmiles(bucket_size=25, minscore=.4, minsimilarity=.4,
                               penalty_multiplier=.5, rdkit_smiles_flags={})
        folder = root / task["directory"]
        logger = logging.getLogger("csv")
        handler = logging.FileHandler(folder / "summary_1.csv")
        handler.setFormatter(CsvFormatter())
        old_level = logger.level
        logger.setLevel(logging.INFO)
        logger.addHandler(handler)
        reporter = StableCSVReporter()
        try:
            for step, affinity in enumerate([12., 0.], 1):
                learning.affinity = affinity
                results = learning.score()
                base = results.size_reward_details[0]["base_reward"]
                df.update_score(results.total_scores, results.smilies, learning.invalid_mask)
                expected = [base, .5 * base, 0] if step == 1 else [.5 * base, .5 * base, 0]
                np.testing.assert_allclose(results.total_scores, expected, atol=1e-7)
                assert results.completed_components[0].transformed_scores[0][0] >= 0
                # Actual production DAP optimizer receives signed rewards unchanged.
                agent_nll = torch.nn.Parameter(torch.tensor([2., 2., 2.]))
                optimizer = torch.optim.SGD([agent_nll], lr=.001)
                reward = RLReward(optimizer, sigma=128)
                _, _, augmented, loss = reward([], results.total_scores, agent_nll,
                                               torch.tensor([2., 2., 2.]), np.array([0, 1]), None, None)
                assert torch.isfinite(loss)
                np.testing.assert_allclose(augmented.detach().numpy(), -2 + 128 * np.array(expected), atol=1e-5)
                if step == 2:
                    assert agent_nll[0].item() > 2  # Negative reward lowers sequence probability.
                reporter.submit(SimpleNamespace(score_results=results, sampled=learning.sampled,
                    agent_nll=[2.] * 3, prior_nll=[2.] * 3, augmented_nll=augmented.detach().tolist(),
                    step=step, mean_score=float(np.mean(results.total_scores)),
                    fraction_valid_smiles=2/3, fraction_duplicate_smiles=1/3))
        finally:
            logger.removeHandler(handler)
            logger.setLevel(old_level)
            handler.close()
        (folder / "agent.chkpt").write_bytes(b"fixture checkpoint")
        assert collect_task(root, manifest, entry, task) == 6
    combine(root, manifest["targets"])
    write_json(root / "status.json", {"state": "complete", "tasks": {
        "example/libinvent": "complete", "example/linkinvent": "complete"}})
    validate(root, root / "experiment.json", root / "audit.json")
    assert json.loads((root / "audit.json").read_text())["passed"]
    # A wrong generated-atom count must not pass even if all CSV copies agree.
    path = root / "example/libinvent/all_generated.csv"
    rows = list(csv.DictReader(path.open()))
    rows[0]["generated_heavy_atoms"] = "100"
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
    (root / "all_generated.csv").unlink()
    combine(root, manifest["targets"])
    with pytest.raises(ValueError, match="Size-adjusted reward"):
        validate(root, root / "experiment.json", root / "bad-audit.json")
