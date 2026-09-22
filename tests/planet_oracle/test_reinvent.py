import csv
import json
import logging
from pathlib import Path
import shlex
import sys
from types import SimpleNamespace

import numpy as np
import pytest
from rdkit import Chem
import torch

from reinvent.models.model_factory.sample_batch import SampleBatch, SmilesState
from reinvent.models.transformer.core.vocabulary import SMILESTokenizer
from reinvent.runmodes.RL.validation import RLConfig
from reinvent.scoring import Scorer
from reinvent.utils import config_parse
from reinvent.utils.helpers import get_tokens_from_vocabulary
from reinvent.validation import ReinventConfig
from reinvent_plugins.components.comp_external_process import ExternalProcess, Parameters

PROJECT = Path(__file__).resolve().parents[2]
CONFIGS = PROJECT / "configs/planet_oracle"


def scoring_config(mode, tmp_path, cfg):
    config = config_parse.read_config(CONFIGS / f"{mode}_rl.toml", "toml")
    ReinventConfig(**config)
    rl = RLConfig(**{k: config[k] for k in ("parameters", "stage", "learning_strategy")})
    assert rl.parameters.batch_size == 16
    assert rl.parameters.randomize_smiles and rl.parameters.isomeric_smiles
    assert config["seed"] == 42 and len(rl.stage) == 1
    assert rl.stage[0].min_steps == rl.stage[0].max_steps == 5
    assert rl.learning_strategy.sigma == 128 and rl.learning_strategy.rate == 0.0001
    assert config["parameters"]["prior_file"] == config["parameters"]["agent_file"]
    client = tmp_path / "client.json"
    client.write_text(json.dumps(cfg))
    scoring = config["stage"][0]["scoring"]
    endpoint = scoring["component"][0]["ExternalProcess"]["endpoint"][0]
    assert len(scoring["component"]) == 1 and endpoint["weight"] == 1
    assert "transform" not in endpoint
    endpoint["params"].update(executable=sys.executable, args=f"-m planet_oracle.client --config {shlex.quote(str(client))}")
    return scoring, endpoint["params"]


def test_real_external_process(running_oracle, tmp_path):
    _, cfg = running_oracle
    _, params = scoring_config("libinvent", tmp_path, cfg)
    component = ExternalProcess(Parameters(**{key: [value] for key, value in params.items()}))
    results = component(["CCO", "invalid", "CCO"])
    assert results.scores[0].tolist() == [0.5, 0.0, 0.5]
    assert results.metadata["planet_affinity"] == [7., None, 7.]
    assert results.metadata["planet_status"] == ["ok", "invalid_input", "ok"]
    cfg["oracle_id"] = "wrong"
    (tmp_path / "client.json").write_text(json.dumps(cfg))
    with pytest.raises(ValueError, match="failed with exit"):
        component(["CCO"])


@pytest.mark.parametrize("mode", ["libinvent", "linkinvent"])
def test_templates_fragments_full_molecule_scoring_and_csv(mode, running_oracle, tmp_path):
    from reinvent.runmodes.RL.libinvent import LibinventLearning
    from reinvent.runmodes.RL.linkinvent import LinkinventLearning
    from reinvent.runmodes.samplers.libinvent import LibinventSampler
    from reinvent.runmodes.samplers.linkinvent import LinkinventSampler
    from reinvent.runmodes.RL.reports.csv_summmary import write_summary
    from reinvent.utils.logmon import CsvFormatter

    oracle, cfg = running_oracle
    scoring, _ = scoring_config(mode, tmp_path, cfg)
    inputs = (CONFIGS / f"{mode}.smi").read_text().splitlines()
    for part in inputs[0].split("|"):
        mol = Chem.MolFromSmiles(part)
        dummy = [a for a in mol.GetAtoms() if a.GetAtomicNum() == 0]
        assert len(dummy) == 1 and dummy[0].GetDegree() == 1
    prior = PROJECT / "REINVENT4/priors" / f"{mode}_transformer_pubchem.prior"
    saved = torch.load(prior, map_location="cpu", weights_only=False)
    allowed = get_tokens_from_vocabulary(saved["vocabulary"])
    assert all(token in allowed[0] for token in SMILESTokenizer().tokenize(inputs[0]))
    assert config_parse.read_smiles_csv_file(str(CONFIGS / f"{mode}.smi"), 0, allowed) == inputs
    del saved

    sampler_cls, learning_cls = ((LibinventSampler, LibinventLearning) if mode == "libinvent"
                                 else (LinkinventSampler, LinkinventLearning))
    generated = ["*C", "*CC", "*C"] if mode == "libinvent" else ["*C*", "*CC*", "*C*"]
    sampled = SampleBatch(inputs * 3, generated, torch.zeros(3))
    sampler = object.__new__(sampler_cls)
    molecules = sampler._join_fragments(sampled)
    assert all(mol is not None for mol in molecules)
    sampled.smilies = [Chem.MolToSmiles(mol, isomericSmiles=True) for mol in molecules]
    sampled.states = np.array([SmilesState.VALID] * 3)
    learning = object.__new__(learning_cls)
    learning.sampled = sampled
    learning.invalid_mask = np.ones(3, dtype=bool)
    learning.duplicate_mask = np.ones(3, dtype=bool)
    learning.scoring_function = Scorer(scoring)
    results = learning.score()
    received = [smi for batch in oracle.predictor.batches for smi in batch]
    assert len(received) == 3
    assert all("*" not in smi and "|" not in smi for smi in received)
    assert sorted(Chem.MolFromSmiles(s).GetNumHeavyAtoms() for s in received) == sorted(m.GetNumHeavyAtoms() for m in molecules)
    assert np.isfinite(results.total_scores).all()
    n_batches = len(oracle.predictor.batches)
    learning.score()  # REINVENT's cache still avoids repeated HTTP inference.
    assert len(oracle.predictor.batches) == n_batches

    data = SimpleNamespace(agent_nll=[0.] * 3, prior_nll=[0.] * 3, augmented_nll=[0.] * 3,
                           score_results=results, sampled=sampled, scaffolds=[], step=1,
                           model_type=mode.capitalize())
    csv_path = tmp_path / f"{mode}.csv"
    logger = logging.getLogger("csv")
    handler = logging.FileHandler(csv_path)
    handler.setFormatter(CsvFormatter())
    level = logger.level
    logger.setLevel(logging.INFO)
    logger.addHandler(handler)
    try:
        write_summary(data, write_header=True)
    finally:
        logger.removeHandler(handler)
        logger.setLevel(level)
        handler.close()
    with csv_path.open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 3
    for row in rows:
        assert row["planet_status (PLANET reward)"] == "ok"
        assert row["planet_oracle_id (PLANET reward)"] == cfg["oracle_id"]
        assert row["planet_target_id (PLANET reward)"] == cfg["target_id"]
        assert float(row["planet_affinity (PLANET reward)"]) > 0


def test_scoring_invalid_zero_and_metadata(running_oracle, tmp_path):
    _, cfg = running_oracle
    scoring, _ = scoring_config("libinvent", tmp_path, cfg)
    scorer = Scorer(scoring)
    smiles = ["CCO", "C.C"]  # syntactically valid RDKit input rejected by oracle chemistry policy
    results = scorer(smiles, np.ones(2, dtype=bool), np.ones(2, dtype=bool))
    assert results.total_scores.tolist() == [0.5, 0.0]
    metadata = results.completed_components[0].component_result.fetch_metadata(smiles)
    # REINVENT stringifies metadata before CSV serialization (the client retains JSON null).
    assert list(metadata["planet_affinity"]) == ["7.0", "None"]
    assert list(metadata["planet_status"]) == ["ok", "invalid_input"]
