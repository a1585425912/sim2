"""Merge the four disjoint formal seed groups; never rerun or drop samples."""
from pathlib import Path
import json,csv,shutil
HERE=Path(__file__).resolve().parent


def main():
    out=HERE/'outputs';chunks=out/'chunks';rows=[];manifests=[]
    for seed in (20261028,20261033,20261038,20261043):
        folder=chunks/str(seed);manifest=json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
        if manifest['status']!='completed' or manifest['seeds']!=list(range(seed,seed+5)):raise ValueError('Incomplete or wrong seed group')
        with (folder/'per_seed.csv').open(encoding='utf-8-sig') as f:part=list(csv.DictReader(f))
        if len(part)!=8000:raise ValueError('Formal seed group must have 8000 records')
        rows.extend(part);manifests.append(manifest)
    expected=set(range(20261028,20261048))
    if set(int(r['seed']) for r in rows)!=expected or len(rows)!=32000:raise ValueError('Missing/duplicated full batch')
    signatures=[m['source_sha256'] for m in manifests]
    if not all(s==signatures[0] for s in signatures):raise ValueError('Workers ran different source versions')
    with (out/'per_seed.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    for p in (chunks/'20261028').glob('first_seed_*.npz'):shutil.copy2(p,out/p.name)
    combined=manifests[0].copy();combined.update(trials=20,seeds=sorted(expected),elapsed_s=max(m['elapsed_s'] for m in manifests),worker_elapsed_s=[m['elapsed_s'] for m in manifests],workers=4,baseline_checks=sum([m['baseline_checks'] for m in manifests],[]),records=len(rows))
    combined['implementation']='sequential.py: persistent joint Gaussian arrival factor and stochastic target cloning; all past measurements retained in frozen linearized factor, each current raw observation counted once'
    combined['boundary']='Centralized augmented iterative MAP with persistent Gaussian arrival and target cloning, not a finite-round distributed GP or RI-FLS. Measured-CPU timing diagnostic supplements matched fixed delay; power allocation is a separate total-power-only variant.'
    (out/'manifest.json').write_text(json.dumps(combined,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='completed',records=len(rows),seeds=20,workers=4),indent=2))


if __name__=='__main__':main()
