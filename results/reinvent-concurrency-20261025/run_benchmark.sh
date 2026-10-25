#!/usr/bin/env bash
set -euo pipefail
PROJECT=/home/ubuntu/experiments/high9_cap9_size001_1000_seed42_20261025/project
BENCH_ROOT=/home/ubuntu/benchmarks/reinvent_concurrency_20261025
PY=/home/ubuntu/experiments/high9_1000_seed42_20261025/project/REINVENT4/.venv/bin/python
export PYTHONPATH="$PROJECT/src:$PROJECT/REINVENT4"
export PYTHONHASHSEED=42
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cd "$PROJECT"
"$PY" src/benchmark_manual_planet_concurrency.py \
  --project "$PROJECT" \
  --experiment "$PROJECT/configs/manual_planet_rl/experiment_high9_cap9_size001_1000.json" \
  --python "$PY" \
  --planet-python /home/ubuntu/experiments/high9_cap9_size001_1000_seed42_20261025/planet_python_nocache.sh \
  --output "$BENCH_ROOT/sweep_t2" --levels 1 2 3 4 --cpu-threads 2 \
  --targets cp3a4 egfr jak2 aa2ar
