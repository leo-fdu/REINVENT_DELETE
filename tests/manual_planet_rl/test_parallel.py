"""Exercise the bounded scheduler without neural models."""
import csv
import json
import time

import pytest
import run_manual_planet_rl as runner


def fake_target(root, targets, _worker):
    assert _worker
    target = targets[0]
    (root / f'{target}.start').write_text(str(time.monotonic()))
    time.sleep(.1)
    with (root / target / 'all_generated.csv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=runner.FIELDS)
        writer.writeheader()
    (root / f'{target}.end').write_text(str(time.monotonic()))
    (root / f'status-{target}.json').write_text(json.dumps({'tasks': {target: 'complete'}}))


def prepared(tmp_path, monkeypatch, same_port=False):
    entries=[]
    for i in range(6):
        name=f't{i}'
        (tmp_path/name).mkdir()
        (tmp_path/name/'server.json').write_text(json.dumps({'port': 8900 if same_port else 8900+i}))
        entries.append({'target':name,'server_config':f'{name}/server.json','tasks':[{'directory':name}]})
    manifest={'schema_version':1,'run_dir':str(tmp_path),'targets':entries,'input_hashes':{},'prepared_hashes':{}}
    (tmp_path/'run.json').write_text(json.dumps(manifest))
    monkeypatch.setattr(runner,'run',fake_target)


def test_parallel_four_workers_collects_every_task(tmp_path,monkeypatch):
    prepared(tmp_path,monkeypatch)
    result=runner.run_parallel(tmp_path,4)
    assert result['state']=='complete' and len(result['tasks'])==6
    events=[]
    for i in range(6):
        events.extend([(float((tmp_path/f't{i}.start').read_text()),1),(float((tmp_path/f't{i}.end').read_text()),-1)])
    active=peak=0
    for _,delta in sorted(events):
        active+=delta
        peak=max(peak,active)
    assert peak==4 and active==0
    assert (tmp_path/'all_generated.csv').is_file()
    with pytest.raises(FileExistsError):runner.run_parallel(tmp_path,4)


def test_parallel_rejects_shared_ports_before_claim(tmp_path,monkeypatch):
    prepared(tmp_path,monkeypatch,same_port=True)
    with pytest.raises(ValueError,match='separate-ports'):runner.run_parallel(tmp_path,4)
    assert not (tmp_path/'.started').exists()
