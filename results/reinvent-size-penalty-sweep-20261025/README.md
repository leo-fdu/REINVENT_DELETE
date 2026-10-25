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
