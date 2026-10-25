#!/usr/bin/env bash
set -euo pipefail
EXPERIMENT_ROOT=/home/ubuntu/experiments/high9_1000_seed42_20261025
PROJECT="$EXPERIMENT_ROOT/project"
PY="$PROJECT/REINVENT4/.venv/bin/python"
export PYTHONPATH="$PROJECT/src:$PROJECT/REINVENT4"
export PYTHONHASHSEED=42
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
cd "$PROJECT"
if test -e "$EXPERIMENT_ROOT/run/.started"; then echo "Formal run already used"; exit 2; fi
"$EXPERIMENT_ROOT/bootstrap/bin/uv" pip freeze --python "$PY" > "$EXPERIMENT_ROOT/environment-freeze.txt"
"$PY" - <<'PY'
import json,torch,rdkit,numpy,scipy,pandas,sys
from pathlib import Path
assert torch.cuda.is_available()
x=torch.randn(128,128,device='cuda'); y=(x@x).sum(); assert torch.isfinite(y)
info={'python':sys.version,'torch':torch.__version__,'rdkit':rdkit.__version__,'numpy':numpy.__version__,'scipy':scipy.__version__,'pandas':pandas.__version__,'cuda_runtime':torch.version.cuda,'device':torch.cuda.get_device_name(0),'omp_threads':2}
Path('../environment.json').write_text(json.dumps(info,indent=2)+'\n')
print(info,flush=True)
cfg=json.loads(Path('configs/manual_planet_rl/experiment_high9_1000.json').read_text())
cfg.update(batch_size=10,steps=2)
Path('../smoke-experiment.json').write_text(json.dumps(cfg,indent=2)+'\n')
PY
"$PY" src/prepare_manual_planet_rl.py --project "$PROJECT" --experiment "$EXPERIMENT_ROOT/smoke-experiment.json" --run-dir "$EXPERIMENT_ROOT/smoke-run-v2" --reinvent-python "$PY" --planet-python "$PY" --device cuda:0
"$PY" src/run_manual_planet_rl.py --run-dir "$EXPERIMENT_ROOT/smoke-run-v2" --targets aa2ar jak2 drd3
"$PY" - <<'PY'
import csv,json
from collections import Counter
from pathlib import Path
p=Path('../smoke-run-v2'); s=json.loads((p/'status.json').read_text())
assert s['state']=='complete' and len(s['tasks'])==6
rows=list(csv.DictReader((p/'all_generated.csv').open()))
assert len(rows)==120
assert set(Counter((r['target'],r['mode']) for r in rows).values())=={20}
print('SMOKE PASSED: six tasks, 120 attempts',flush=True)
PY
"$PY" src/run_manual_planet_rl.py --run-dir "$EXPERIMENT_ROOT/run"
"$PY" src/validate_manual_planet_results.py --run-dir "$EXPERIMENT_ROOT/run" --experiment "$PROJECT/configs/manual_planet_rl/experiment_high9_1000.json" --report "$EXPERIMENT_ROOT/run/acceptance-report.json"
"$PY" - <<'PY'
import gzip,hashlib,json,shutil
from pathlib import Path
root=Path('../run')
with (root/'all_generated.csv').open('rb') as src, gzip.open(root/'all_generated.csv.gz','wb') as dst:
 shutil.copyfileobj(src,dst)
shutil.copyfile('../environment.json',root/'environment.json')
shutil.copyfile('../environment-freeze.txt',root/'environment-freeze.txt')
shutil.copyfile('../bootstrap/runtime-requirements.txt',root/'runtime-requirements.txt')
shutil.copyfile('../run_high9_cloud.sh',root/'run_high9_cloud.sh')
shutil.copyfile('../resume_high9_cloud.sh',root/'resume_high9_cloud.sh')
shutil.copyfile('src/validate_manual_planet_results.py',root/'validate_manual_planet_results.py')
shutil.copyfile('../pocket-validation/report.json',root/'pocket-validation-report.json')
shutil.copyfile('transfer-manifest.json',root/'source-transfer-manifest.json')
checks={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sorted(root.rglob('*')) if p.is_file() and p.name!='checksums.json'}
(root/'checksums.json').write_text(json.dumps(checks,indent=2)+'\n')
print('FORMAL RUN ACCEPTED; transfer manifest written',flush=True)
PY
