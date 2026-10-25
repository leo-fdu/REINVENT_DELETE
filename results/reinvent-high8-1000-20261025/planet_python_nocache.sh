#!/usr/bin/env bash
set -euo pipefail
export PYTORCH_NO_CUDA_MEMORY_CACHING=1
exec /home/ubuntu/experiments/high9_1000_seed42_20261025/project/REINVENT4/.venv/bin/python "$@"
