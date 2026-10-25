"""Benchmark isolated RL workers on a fixed cohort without changing training code.

Run on Linux with nvidia-smi. Each worker owns one PLANET service and one RL
task. Preparation, result validation and comparisons are outside measured time.
"""

import argparse
import csv
import json
import os
from pathlib import Path
import statistics
import subprocess
import time

from prepare_manual_planet_rl import prepare, read_experiment, sha256, write_json
from run_manual_planet_rl import stop_process
from validate_manual_planet_results import validate


def cpu_ticks():
    values = [int(v) for v in Path('/proc/stat').read_text().splitlines()[0].split()[1:9]]
    return sum(values), values[3] + values[4]


def gpu_sample():
    result = subprocess.run([
        'nvidia-smi', '--query-gpu=utilization.gpu,utilization.memory,memory.used,power.draw',
        '--format=csv,noheader,nounits', '--id=0',
    ], capture_output=True, text=True, check=True, timeout=10)
    return dict(zip(['gpu_util_pct', 'memory_util_pct', 'gpu_memory_mib', 'power_w'],
                    [float(v.strip()) for v in result.stdout.strip().split(',')]))


def compare_outputs(jobs, baseline):
    result = {'rows': 0, 'smiles_mismatches': 0, 'max_raw_score_difference': 0.0,
              'max_training_reward_difference': 0.0}
    for job in jobs:
        target, mode = job['target'], job['mode']
        with (job['root'] / 'all_generated.csv').open() as f:
            current = list(csv.DictReader(f))
        with (baseline / f'{target}_{mode}' / 'all_generated.csv').open() as f:
            old = list(csv.DictReader(f))
        if len(current) != len(old):
            raise ValueError('Benchmark cohorts have different row counts')
        for a, b in zip(current, old):
            result['rows'] += 1
            if (a['step'], a['sample_index']) != (b['step'], b['sample_index']):
                raise ValueError('Benchmark sample indices differ')
            if a['smiles'] != b['smiles']:
                result['smiles_mismatches'] += 1
                continue
            for field, key in [('planet_affinity', 'max_raw_score_difference'),
                               ('training_reward', 'max_training_reward_difference')]:
                if a[field] and b[field]:
                    result[key] = max(result[key], abs(float(a[field]) - float(b[field])))
    return result


def run_round(args, concurrency, baseline):
    folder = args.output / f'c{concurrency}_t{args.cpu_threads}'
    folder.mkdir()
    jobs = []
    # All levels process exactly the same jobs in the same queue order.
    for index, (target, mode) in enumerate((t, m) for t in args.targets
                                          for m in ('libinvent', 'linkinvent')):
        root = folder / f'{target}_{mode}'
        manifest = prepare(args.project, args.project / 'configs/manual_design/manifest.json',
                           args.experiment, root, args.python, args.planet_python,
                           'cuda:0', args.port_base + index)
        entry = next(e for e in manifest['targets'] if e['target'] == target)
        entry['tasks'] = [t for t in entry['tasks'] if t['mode'] == mode]
        # run.json is the selection manifest; all prepared input hashes stay intact.
        manifest['targets'] = [entry]
        manifest['total_tasks'] = 1
        manifest['benchmark_selection'] = {'target': target, 'mode': mode}
        write_json(root / 'run.json', manifest)
        jobs.append({'root': root, 'target': target, 'mode': mode})

    env = os.environ.copy()
    env.update(PYTHONHASHSEED='42', OMP_NUM_THREADS=str(args.cpu_threads),
               OPENBLAS_NUM_THREADS=str(args.cpu_threads), MKL_NUM_THREADS=str(args.cpu_threads),
               PYTHONPATH=f'{args.project}/src:{args.project}/REINVENT4')
    active, finished, samples, handles = [], [], [], []
    pending = list(jobs)
    failed = False
    previous_ticks = cpu_ticks()
    started = time.monotonic()
    start_epoch = time.time()
    try:
        while pending or active:
            for item in list(active):
                code = item['process'].poll()
                if code is not None:
                    item['returncode'] = code
                    item['wall_seconds'] = time.monotonic() - item['started']
                    finished.append(item)
                    active.remove(item)
                    failed |= code != 0
            if failed:
                break
            while pending and len(active) < concurrency:
                job = pending.pop(0)
                handle = (job['root'] / 'benchmark-worker.log').open('w')
                handles.append(handle)
                process = subprocess.Popen(
                    [args.python, str(args.project / 'src/run_manual_planet_rl.py'),
                     '--run-dir', str(job['root'])], env=env, cwd=args.project,
                    stdout=handle, stderr=subprocess.STDOUT, start_new_session=True)
                active.append({**job, 'process': process, 'started': time.monotonic()})
            if not active:
                break
            sample = gpu_sample()
            ticks = cpu_ticks()
            total, idle = ticks[0] - previous_ticks[0], ticks[1] - previous_ticks[1]
            previous_ticks = ticks
            mem = dict(line.split(':', 1) for line in Path('/proc/meminfo').read_text().splitlines())
            sample.update(elapsed_seconds=time.monotonic() - started, active_workers=len(active),
                          cpu_busy_pct=100 * (1 - idle / total) if total else 0,
                          available_ram_mib=int(mem['MemAvailable'].split()[0]) / 1024)
            samples.append(sample)
            time.sleep(1)
    finally:
        elapsed = time.monotonic() - started
        for item in active:
            stop_process(item['process'])
        for handle in handles:
            handle.close()

    with (folder / 'samples.csv').open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(samples[0]) if samples else [])
        writer.writeheader()
        writer.writerows(samples)
    cfg = read_experiment(args.experiment)
    count = len(jobs) * cfg['batch_size'] * cfg['steps']
    saturated = [s for s in samples if s['active_workers'] == concurrency]
    result = {
        'concurrency': concurrency, 'cpu_threads_per_process': args.cpu_threads,
        'requested_tasks': len(jobs),
        'completed_tasks': sum(j['returncode'] == 0 for j in finished),
        'interrupted_tasks': len(active), 'failed': failed,
        'attempts': count if not failed else None, 'wall_seconds': elapsed,
        'start_epoch': start_epoch, 'end_epoch': start_epoch + elapsed,
        'attempts_per_second': count / elapsed if not failed else None,
        'gpu_util_mean_pct': statistics.mean(s['gpu_util_pct'] for s in samples),
        'gpu_util_median_pct': statistics.median(s['gpu_util_pct'] for s in samples),
        'gpu_util_ge90_fraction': sum(s['gpu_util_pct'] >= 90 for s in samples) / len(samples),
        'full_worker_count_gpu_util_mean_pct': statistics.mean(s['gpu_util_pct'] for s in saturated),
        'gpu_memory_peak_mib': max(s['gpu_memory_mib'] for s in samples),
        'cpu_busy_mean_pct': statistics.mean(s['cpu_busy_pct'] for s in samples),
        'power_mean_w': statistics.mean(s['power_w'] for s in samples),
        'available_ram_min_mib': min(s['available_ram_mib'] for s in samples),
        'sampling_interval_seconds': 1, 'samples': len(samples),
        'workers': [{'target': j['target'], 'mode': j['mode'], 'returncode': j['returncode'],
                     'wall_seconds': j['wall_seconds']} for j in finished],
        'notes': 'One independent PLANET service per worker; startup and task transitions included. Preparation and audits excluded. GPU-memory peak is sampled, not allocator peak.',
    }
    write_json(folder / 'summary.json', result)
    if not failed:
        for job in jobs:
            validate(job['root'], args.experiment, job['root'] / 'acceptance-report.json')
        result['all_reward_audits_passed'] = True
        if baseline:
            result['comparison_to_c1'] = compare_outputs(jobs, baseline)
        write_json(folder / 'summary.json', result)
    print(json.dumps(result), flush=True)
    return result, folder


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project', type=Path, required=True)
    parser.add_argument('--experiment', type=Path, required=True)
    parser.add_argument('--python', required=True)
    parser.add_argument('--planet-python', required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--levels', nargs='+', type=int, default=[1, 2, 3, 4])
    parser.add_argument('--targets', nargs='+', default=['cp3a4', 'egfr', 'jak2', 'aa2ar'])
    parser.add_argument('--cpu-threads', type=int, default=2)
    parser.add_argument('--port-base', type=int, default=18760)
    parser.add_argument('--baseline', type=Path)
    args = parser.parse_args()
    if (args.output.exists() or args.cpu_threads < 1 or len(set(args.targets)) != len(args.targets)
            or len(set(args.levels)) != len(args.levels)
            or any(n < 1 or n > 2 * len(args.targets) for n in args.levels)):
        raise ValueError('Use a new output directory and valid unique levels/targets')
    args.project = args.project.resolve()
    args.output = args.output.resolve()
    args.output.mkdir(parents=True)
    write_json(args.output / 'benchmark-config.json', {
        'project': str(args.project), 'experiment': str(args.experiment),
        'experiment_sha256': sha256(args.experiment), 'levels': args.levels,
        'targets': args.targets, 'cpu_threads_per_process': args.cpu_threads,
        'cohort_tasks': 2 * len(args.targets), 'script_sha256': sha256(__file__),
    })
    baseline = args.baseline
    results = []
    for level in args.levels:
        print(f'START concurrency={level}, cpu_threads={args.cpu_threads}', flush=True)
        result, folder = run_round(args, level, baseline)
        results.append(result)
        write_json(args.output / 'summary.json', results)
        if result['failed']:
            print('STOP: worker failed; higher concurrency not attempted', flush=True)
            break
        if level == 1:
            baseline = folder


if __name__ == '__main__':
    main()
