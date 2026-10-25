#!/usr/bin/env bash
set -euo pipefail
EXPERIMENT_ROOT=/home/ubuntu/experiments/high8_1000_seed42_20261025
PROJECT=/home/ubuntu/experiments/high9_1000_seed42_20261025/project
export EXPERIMENT_ROOT
PY="$PROJECT/REINVENT4/.venv/bin/python"
export PYTHONPATH="$PROJECT/src:$PROJECT/REINVENT4"
export PYTHONHASHSEED=42
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
cd "$PROJECT"
"$PROJECT/../bootstrap/bin/uv" pip freeze --python "$PY" > "$EXPERIMENT_ROOT/environment-freeze.txt"
"$PY" - <<'PY'
import json,torch,rdkit,numpy,scipy,pandas,sys,os
from pathlib import Path
assert torch.cuda.is_available()
x=torch.randn(128,128,device='cuda'); y=(x@x).sum(); assert torch.isfinite(y)
info={'python':sys.version,'torch':torch.__version__,'rdkit':rdkit.__version__,'numpy':numpy.__version__,'scipy':scipy.__version__,'pandas':pandas.__version__,'cuda_runtime':torch.version.cuda,'device':torch.cuda.get_device_name(0),'omp_threads':2}
(Path(os.environ['EXPERIMENT_ROOT'])/'environment.json').write_text(json.dumps(info,indent=2)+'\n')
print(info,flush=True)

PY
"$PY" src/prepare_manual_planet_rl.py --project "$PROJECT" --experiment "$PROJECT/configs/manual_planet_rl/experiment_high8_1000.json" --run-dir "$EXPERIMENT_ROOT/run" --reinvent-python "$PY" --planet-python "$PY" --device cuda:0
"$PY" src/validate_manual_planet_inputs.py --run-dir "$EXPERIMENT_ROOT/run" --report-dir "$EXPERIMENT_ROOT/pocket-validation"
"$PY" src/run_manual_planet_rl.py --run-dir "$EXPERIMENT_ROOT/run"
"$PY" src/validate_manual_planet_results.py --run-dir "$EXPERIMENT_ROOT/run" --experiment "$PROJECT/configs/manual_planet_rl/experiment_high8_1000.json" --report "$EXPERIMENT_ROOT/run/acceptance-report.json"
"$PY" src/validate_manual_planet_training.py --run-dir "$EXPERIMENT_ROOT/run" --report "$EXPERIMENT_ROOT/run/training-numerics.json"
"$PY" - <<'PY'
import gzip,hashlib,json,shutil,os
from pathlib import Path
base=Path(os.environ['EXPERIMENT_ROOT'])
root=base/'run'
with (root/'all_generated.csv').open('rb') as src, gzip.open(root/'all_generated.csv.gz','wb') as dst:
 shutil.copyfileobj(src,dst)
shutil.copyfile(base/'environment.json',root/'environment.json')
shutil.copyfile(base/'environment-freeze.txt',root/'environment-freeze.txt')
shutil.copyfile(Path.cwd() / '../bootstrap/runtime-requirements.txt',root/'runtime-requirements.txt')
shutil.copyfile(base/'run_high8_cloud.sh',root/'run_high8_cloud.sh')
shutil.copyfile('src/validate_manual_planet_results.py',root/'validate_manual_planet_results.py')
shutil.copyfile(base/'pocket-validation/report.json',root/'pocket-validation-report.json')
shutil.copyfile(base/'source-transfer-manifest.json',root/'source-transfer-manifest.json')
checks={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.rglob('*')) if p.is_file() and p.name!='checksums.json'}
(root/'checksums.json').write_text(json.dumps(checks,indent=2)+'\n')
print('FORMAL RUN ACCEPTED; transfer manifest written',flush=True)
PY
