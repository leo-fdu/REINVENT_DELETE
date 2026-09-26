"""REINVENT entry point with the same repeat penalty within and across batches.

The upstream sampler still marks duplicates for reporting. Only the two
Transformer learning classes used by this experiment are replaced, in this
process; the upstream source and all other generation modes are untouched.
"""

import logging

import numpy as np

from reinvent.models.model_factory.sample_batch import SmilesState
from reinvent.runmodes import RL
from reinvent.runmodes.RL.libinvent import LibinventLearning
from reinvent.runmodes.RL.linkinvent import LinkinventLearning
from reinvent.runmodes.RL.reports.csv_summmary import RLCSVReporter


class StableCSVReporter:
    """Keep all generated rows and stable columns even if a batch is all invalid."""

    def __init__(self):
        self.header_written = False

    def submit(self, data):
        components = data.score_results.completed_components
        if len(components) != 1 or components[0].component_names != ["PLANET reward"]:
            raise ValueError("This entry point requires exactly one PLANET reward component")
        component = components[0]
        raw = component.component_result.fetch_scores(data.score_results.smilies, transpose=True)[0]
        metadata = component.component_result.fetch_metadata(data.score_results.smilies)
        keys = ("planet_affinity", "planet_status", "planet_error", "planet_target_id", "planet_oracle_id")
        columns = {key: metadata.get(key, [None] * len(raw)) for key in keys}
        header = ["Agent", "Prior", "Target", "Score", "SMILES", "SMILES_state",
                  "Input_SMILES", "Generated_SMILES", "PLANET reward", "PLANET reward (raw)"]
        header += [f"{key} (PLANET reward)" for key in keys] + ["step"]
        csv_logger = logging.getLogger("csv")
        if not self.header_written:
            csv_logger.info(header)
            self.header_written = True
        for index, smiles in enumerate(data.score_results.smilies):
            row = [data.agent_nll[index], data.prior_nll[index], data.augmented_nll[index],
                   data.score_results.total_scores[index], smiles, data.sampled.states[index].value,
                   data.sampled.items1[index], data.sampled.items2[index],
                   component.transformed_scores[0][index], raw[index]]
            row += [columns[key][index] for key in keys] + [data.step]
            csv_logger.info(row)
        logging.getLogger("reinvent").info(
            "Step %s: reward %.4f, valid %.1f%%, batch duplicates %.1f%%",
            data.step, data.mean_score, 100 * data.fraction_valid_smiles,
            100 * data.fraction_duplicate_smiles,
        )


class RepeatRewards:
    def _setup_reporters(self, tb_logdir):
        super()._setup_reporters(tb_logdir)
        self.reporters = [StableCSVReporter() if isinstance(reporter, RLCSVReporter) else reporter
                          for reporter in self.reporters]

    def score(self):
        # Keep invalid masking, but let the diversity filter handle every repeat.
        self.duplicate_mask = np.ones(len(self.sampled.states), dtype=bool)
        return super().score()

    def _update_common(self, results, orig_smilies):
        agent_nlls = self._state.agent.likelihood_smiles(self.sampled).likelihood
        prior_nlls = self.prior.likelihood_smiles(self.sampled).likelihood
        valid_indices = np.flatnonzero(self.sampled.states != SmilesState.INVALID)
        return self.reward_nlls(
            orig_smilies,
            results.total_scores,
            agent_nlls,
            prior_nlls,
            valid_indices,
            self.inception,
            self._state.agent,
        )


class LibinventRepeatLearning(RepeatRewards, LibinventLearning):
    pass


class LinkinventRepeatLearning(RepeatRewards, LinkinventLearning):
    pass


def main():
    RL.LibinventTransformerLearning = LibinventRepeatLearning
    RL.LinkinventTransformerLearning = LinkinventRepeatLearning
    from reinvent.Reinvent import main_script

    main_script()


if __name__ == "__main__":
    main()
