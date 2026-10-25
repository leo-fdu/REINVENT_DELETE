#!/usr/bin/env bash
set -euo pipefail
BASE=/home/ubuntu/experiments/size_penalty_sweep_20261025
PROJECT="$BASE/project"
PY=/home/ubuntu/experiments/high9_1000_seed42_20261025/project/REINVENT4/.venv/bin/python
export PYTHONPATH="$PROJECT/src:$PROJECT/REINVENT4" PYTHONHASHSEED=42
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 MKL_NUM_THREADS=2
cd "$PROJECT"
for TAG in 0015 002; do
  RUN="$BASE/size${TAG}/run"
  CFG="$PROJECT/configs/manual_planet_rl/experiment_high9_cap9_size${TAG}_1000.json"
  "$PY" src/prepare_manual_planet_rl.py --project "$PROJECT" --experiment "$CFG" --run-dir "$RUN" --reinvent-python "$PY" --planet-python "$BASE/planet_python_nocache.sh" --device cuda:0 --separate-ports --port 8900
  "$PY" src/validate_manual_planet_inputs.py --run-dir "$RUN" --report-dir "$BASE/size${TAG}/pocket-validation"
  "$PY" src/run_manual_planet_rl.py --run-dir "$RUN" --workers 4
  "$PY" src/validate_manual_planet_results.py --run-dir "$RUN" --experiment "$CFG" --report "$RUN/acceptance-report.json"
  "$PY" src/validate_manual_planet_training.py --run-dir "$RUN" --report "$RUN/training-numerics.json"
  "$PY" - "$RUN" <<'PY'
import gzip,hashlib,json,shutil,sys
from pathlib import Path
root=Path(sys.argv[1])
with (root/'all_generated.csv').open('rb') as source,gzip.open(root/'all_generated.csv.gz','wb') as dest:shutil.copyfileobj(source,dest)
checks={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.rglob('*')) if p.is_file() and p.name!='checksums.json'}
(root/'checksums.json').write_text(json.dumps(checks,indent=2)+'\n')
print('ACCEPTED',root,flush=True)
PY
done
