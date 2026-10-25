"""Verify the complete local raw mirror against the finalized server manifests."""
import hashlib
import json
from pathlib import Path

repo=Path(__file__).resolve().parents[2]
report={}
for tag in ['0015','002']:
    root=repo/f'results/server-evidence/size_penalty_sweep_20261025/size{tag}/run'
    hashes=json.loads((root/'checksums.json').read_text())
    errors=[]
    for name,expected in hashes.items():
        path=root/name
        if not path.is_file():errors.append(f'missing: {name}');continue
        digest=hashlib.sha256()
        with path.open('rb') as handle:
            for chunk in iter(lambda:handle.read(1024*1024),b''):digest.update(chunk)
        if digest.hexdigest()!=expected:errors.append(f'hash mismatch: {name}')
    report[tag]={'passed':not errors,'files_verified':len(hashes),'bytes_verified':sum((root/p).stat().st_size for p in hashes if (root/p).is_file()),'errors':errors}
    print(tag,report[tag],flush=True)
Path(__file__).with_name('transfer-verification.json').write_text(json.dumps(report,indent=2)+'\n')
assert all(r['passed'] for r in report.values())
