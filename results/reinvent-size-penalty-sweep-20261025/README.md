# Generated-fragment size penalty sweep

Two independent reruns: generated_heavy_atom_penalty=0.015 and 0.02.
Each has 30 tasks × 10 RL steps × batch 100 = 30,000 sampling attempts;
seed 42, sigmoid low=4/high=9/k=0.25, raw affinity capped at 9 only for reward,
DAP sigma=128 and learning rate=0.0001. Duplicate multiplier 0.5 applies to the
whole signed reward. No other scientific setting changes from the 0.01 run.

The groups run sequentially. Each group uses four target workers, each running
its LibINVENT and LinkINVENT tasks in sequence. Thus at most four RL tasks and
four PLANET services are active. Each target has a distinct port; OMP/MKL/OpenBLAS
threads are 2 per process. PLANET CUDA caching is disabled, as in the 0.01 run.

Cloud root: `/home/ubuntu/experiments/size_penalty_sweep_20261025`.
The independent `size0015/run` and `size002/run` trees retain raw results and
checkpoints. `run_cloud.sh` records the commands and validates every result and
trained checkpoint before marking a group accepted.

Analysis wrappers are under `evaluation/descriptor/generation_analysis_reinvent/`
in `high9_cap9_size0015_1000_seed42_20261025` and
`high9_cap9_size002_1000_seed42_20261025`. They reproduce all nine metrics and
compare with the 0.01 run. This is a same-seed comparison, not a multi-seed test.

## Completed results

Both groups passed the independent budget/reward/oracle/checkpoint audit and
finite/updated-parameter audit. Each has exactly 30,000 sampling attempts.
Four-worker training wall times were 552.21 s (0.015) and 517.22 s (0.02).
The cloud driver exited 0. The baseline 0.01 is an earlier independent run.

The table uses equal-weight means across 30 tasks at step 10. Repeats are counted
across each task's full 1,000-attempt history, using the training duplicate flag.

| Penalty | Final whole heavy atoms | Final generated heavy atoms | Final MW | Final raw PLANET | Full-run repeats | Invalid |
|---|---:|---:|---:|---:|---:|---:|
| 0.01 | 41.21 | 22.34 | 583.70 | 9.249 | 13.18% | 0.93% |
| 0.015 | 38.55 | 19.69 | 545.63 | 8.894 | 14.29% | 0.82% |
| 0.02 | 34.86 | 15.99 | 492.39 | 8.478 | 19.45% | 0.77% |

Increasing the penalty reduced size growth. With 0.02 the average size curve
levels off near step 4, but the raw PLANET score is lower and repetition higher.
All three groups share exactly the same first 3,000 SMILES and raw PLANET scores.
These are same-seed observations; they do not establish multi-seed robustness.
Signed negative rewards are intentional: 300 / 548 / 1,013 attempts for
0.01 / 0.015 / 0.02 respectively. The duplicate multiplier still applies to the
whole signed reward, as requested.

Each new analysis contains all nine metric groups, 324 verified PNGs, 33 summary
rows, and per-task/per-mode tables. `comparison_curves.png` plots task-equal means;
`comparison.csv` and `comparison.json` retain the numeric comparison.
`source-checksums.json` records 397 cloud source/config/weight files, all identical
to their local counterparts. `environment.json` and `timing.json` preserve runtime
provenance. Reproduce the analysis by running each experiment's `src/run_analysis.py`,
then this folder's `summarize.py` with the local `reinvent4` Python environment.

Full raw mirrors, including model checkpoints, are under
`results/server-evidence/size_penalty_sweep_20261025/size0015/run` and `size002/run`.
They are covered by the existing raw-evidence ignore rule. Compact result tables
and audit reports are versioned in `results/reinvent-high9-cap9-size0015-1000-20261025`
and `results/reinvent-high9-cap9-size002-1000-20261025`.
`verify_transfer.py` checks every raw file against the server's finalized manifest.

Validation: 97 existing tests and two new scheduler tests passed; the full GPU
runs additionally validated the four-worker execution on all 60 real tasks.
