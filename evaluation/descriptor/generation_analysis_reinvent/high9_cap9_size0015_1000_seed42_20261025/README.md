# Capped PLANET size penalty .015

30 tasks × 10 RL steps × batch 100; seed 42; low=4, high=9, k=0.25;
raw PLANET capped at 9 only for reward. Duplicate multiplier 0.5 applies after
the signed generated-heavy-atom penalty.

Run `src/run_analysis.py` with the local reinvent4 environment to recreate all
nine analysis groups. `analysis/comparison_size001_baseline.csv` compares with
penalty 0.01; `analysis/provenance.json` records inputs and shared analysis hashes.
`analysis/artifact-validation.json` records the complete image/table checks.

The cross-penalty comparison and acceptance details are in
`results/reinvent-size-penalty-sweep-20261025/README.md` at the repository root.
