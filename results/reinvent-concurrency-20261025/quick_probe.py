import os,time,json,subprocess,csv
from pathlib import Path
from prepare_manual_planet_rl import prepare,write_json
from run_manual_planet_rl import stop_process
from benchmark_manual_planet_concurrency import gpu_sample,cpu_ticks
project=Path('/home/ubuntu/experiments/high9_cap9_size001_1000_seed42_20261025/project')
root=Path('/home/ubuntu/benchmarks/reinvent_concurrency_20261025/quick');root.mkdir()
py='/home/ubuntu/experiments/high9_1000_seed42_20261025/project/REINVENT4/.venv/bin/python'
planet='/home/ubuntu/experiments/high9_cap9_size001_1000_seed42_20261025/planet_python_nocache.sh'
cfg=project/'configs/manual_planet_rl/experiment_high9_cap9_size001_1000.json'
env=dict(os.environ,PYTHONPATH=f'{project}/src:{project}/REINVENT4',PYTHONHASHSEED='42',OMP_NUM_THREADS='2',MKL_NUM_THREADS='2',OPENBLAS_NUM_THREADS='2')
cohort=[('cp3a4','linkinvent'),('egfr','linkinvent'),('jak2','libinvent'),('aa2ar','libinvent')]
results=[]
for c in [2,3,4]:
 folder=root/f'c{c}';folder.mkdir();jobs=[]
 for i,(target,mode) in enumerate(cohort[:c]):
  run=folder/f'{target}_{mode}';m=prepare(project,project/'configs/manual_design/manifest.json',cfg,run,py,planet,'cuda:0',19760+i)
  e=next(e for e in m['targets'] if e['target']==target);e['tasks']=[t for t in e['tasks'] if t['mode']==mode];m['targets']=[e];m['total_tasks']=1;write_json(run/'run.json',m);jobs.append(run)
 active=[];handles=[];samples=[];start=time.monotonic();prev=cpu_ticks()
 print('START quick concurrency',c,flush=True)
 try:
  for run in jobs:
   h=(run/'quick-worker.log').open('w');handles.append(h)
   active.append(subprocess.Popen([py,str(project/'src/run_manual_planet_rl.py'),'--run-dir',str(run)],cwd=project,env=env,stdout=h,stderr=subprocess.STDOUT,start_new_session=True))
  while time.monotonic()-start<35:
   s=gpu_sample();ticks=cpu_ticks();total=ticks[0]-prev[0];idle=ticks[1]-prev[1];prev=ticks
   s.update(elapsed_seconds=time.monotonic()-start,cpu_busy_pct=100*(1-idle/total) if total else 0,alive=sum(p.poll() is None for p in active));samples.append(s)
   if any(p.poll() not in (None,0) for p in active):break
   time.sleep(1)
  codes=[p.poll() for p in active]
 finally:
  for p in active:stop_process(p)
  for h in handles:h.close()
 warm=[s for s in samples if s['elapsed_seconds']>=10]
 out={'concurrency':c,'cpu_threads_per_process':2,'warmup_seconds':10,'observed_seconds':time.monotonic()-start,'measured_samples':len(warm),'gpu_mean':sum(s['gpu_util_pct'] for s in warm)/len(warm),'gpu_peak':max(s['gpu_util_pct'] for s in warm),'gpu_ge90_fraction':sum(s['gpu_util_pct']>=90 for s in warm)/len(warm),'vram_peak_mib':max(s['gpu_memory_mib'] for s in samples),'cpu_mean':sum(s['cpu_busy_pct'] for s in warm)/len(warm),'pre_stop_returncodes':codes,'note':'Short load probe; incomplete RL tasks stopped intentionally after measurement, not production results.'}
 with (folder/'samples.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(samples[0]));w.writeheader();w.writerows(samples)
 results.append(out);write_json(root/'summary.json',results);print(json.dumps(out),flush=True)
 if any(code not in (None,0) for code in codes):break
