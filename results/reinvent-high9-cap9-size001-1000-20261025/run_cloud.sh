#!/usr/bin/env bash
set -euo pipefail
EXPERIMENT_ROOT=/home/ubuntu/experiments/high9_cap9_size001_1000_seed42_20261025
PROJECT="$EXPERIMENT_ROOT/project"
export EXPERIMENT_ROOT
PY=/home/ubuntu/experiments/high9_1000_seed42_20261025/project/REINVENT4/.venv/bin/python
export PYTHONPATH="$PROJECT/src:$PROJECT/REINVENT4"
export PYTHONHASHSEED=42
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
cd "$PROJECT"
/home/ubuntu/experiments/high9_1000_seed42_20261025/bootstrap/bin/uv pip freeze --python "$PY" > "$EXPERIMENT_ROOT/environment-freeze.txt"
"$PY" - <<'PY'
import json,torch,rdkit,numpy,scipy,pandas,sys,os
from pathlib import Path
assert torch.cuda.is_available()
x=torch.randn(128,128,device='cuda'); y=(x@x).sum(); assert torch.isfinite(y)
info={'python':sys.version,'torch':torch.__version__,'rdkit':rdkit.__version__,'numpy':numpy.__version__,'scipy':scipy.__version__,'pandas':pandas.__version__,'cuda_runtime':torch.version.cuda,'device':torch.cuda.get_device_name(0),'omp_threads':2,'planet_environment':{'PYTORCH_NO_CUDA_MEMORY_CACHING':'1'},'reinvent_cuda_allocator':'default'}
(Path(os.environ['EXPERIMENT_ROOT'])/'environment.json').write_text(json.dumps(info,indent=2)+'\n')
print(info,flush=True)

PY
"$PY" src/prepare_manual_planet_rl.py --project "$PROJECT" --experiment "$PROJECT/configs/manual_planet_rl/experiment_high9_cap9_size001_1000.json" --run-dir "$EXPERIMENT_ROOT/run" --reinvent-python "$PY" --planet-python "$EXPERIMENT_ROOT/planet_python_nocache.sh" --device cuda:0
"$PY" src/validate_manual_planet_inputs.py --run-dir "$EXPERIMENT_ROOT/run" --report-dir "$EXPERIMENT_ROOT/pocket-validation"
"$PY" - <<'PYSMOKE'
import os,json
from pathlib import Path
base=Path(os.environ['EXPERIMENT_ROOT'])
cfg=json.loads((base/'run/experiment.json').read_text());cfg.update(batch_size=10,steps=2)
(base/'smoke-experiment.json').write_text(json.dumps(cfg,indent=2)+'\n')
PYSMOKE
"$PY" src/prepare_manual_planet_rl.py --project "$PROJECT" --experiment "$EXPERIMENT_ROOT/smoke-experiment.json" --run-dir "$EXPERIMENT_ROOT/smoke-run" --reinvent-python "$PY" --planet-python "$EXPERIMENT_ROOT/planet_python_nocache.sh" --device cuda:0
"$PY" src/run_manual_planet_rl.py --run-dir "$EXPERIMENT_ROOT/smoke-run" --targets aa2ar cp3a4
"$PY" - <<'PYSMOKE'
import os,csv,math,json
from pathlib import Path
from collections import Counter
base=Path(os.environ['EXPERIMENT_ROOT'])
rows=list(csv.DictReader((base/'smoke-run/all_generated.csv').open()))
assert len(rows)==80 and len(Counter((r['target'],r['mode']) for r in rows))==4
for r in rows:
 if r['planet_status']=='ok':
  capped=min(float(r['planet_affinity']),9)
  sig=1/(1+10**(-.5*(capped-6.5)))
  reward=sig-.01*int(r['generated_heavy_atoms'])
  expected=reward*(.5 if r['is_repeat']=='1' else 1)
  assert abs(float(r['training_reward'])-expected)<2e-6
 else:assert float(r['training_reward'])==0
(base/'smoke-validation.json').write_text(json.dumps({'passed':True,'tasks':4,'attempts':80})+'\n')
print('SMOKE PASSED: actual models, four tasks, 80 attempts, signed size reward verified',flush=True)
PYSMOKE
"$PY" src/run_manual_planet_rl.py --run-dir "$EXPERIMENT_ROOT/run"
"$PY" src/validate_manual_planet_results.py --run-dir "$EXPERIMENT_ROOT/run" --experiment "$PROJECT/configs/manual_planet_rl/experiment_high9_cap9_size001_1000.json" --report "$EXPERIMENT_ROOT/run/acceptance-report.json"
"$PY" src/validate_manual_planet_training.py --run-dir "$EXPERIMENT_ROOT/run" --report "$EXPERIMENT_ROOT/run/training-numerics.json"
"$PY" - <<'PY'
import gzip,hashlib,json,shutil,os
from pathlib import Path
base=Path(os.environ['EXPERIMENT_ROOT'])
root=base/'run'
with (root/'all_generated.csv').open('rb') as src, gzip.open(root/'all_generated.csv.gz','wb') as dst:
 shutil.copyfileobj(src,dst)
shutil.copyfile(base/'smoke-validation.json',root/'smoke-validation.json')
shutil.copyfile(base/'environment.json',root/'environment.json')
shutil.copyfile(base/'environment-freeze.txt',root/'environment-freeze.txt')
shutil.copyfile('/home/ubuntu/experiments/high9_1000_seed42_20261025/bootstrap/runtime-requirements.txt',root/'runtime-requirements.txt')
shutil.copyfile(base/'run_cloud.sh',root/'run_cloud.sh')
shutil.copyfile(base/'planet_python_nocache.sh',root/'planet_python_nocache.sh')
shutil.copyfile('src/validate_manual_planet_results.py',root/'validate_manual_planet_results.py')
shutil.copyfile(base/'pocket-validation/report.json',root/'pocket-validation-report.json')
shutil.copyfile(base/'source-transfer-manifest.json',root/'source-transfer-manifest.json')
checks={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.rglob('*')) if p.is_file() and p.name!='checksums.json'}
(root/'checksums.json').write_text(json.dumps(checks,indent=2)+'\n')
print('FORMAL RUN ACCEPTED; transfer manifest written',flush=True)
PY
