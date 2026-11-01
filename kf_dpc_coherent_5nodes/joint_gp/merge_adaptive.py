"""Combine complete Gaussian and T-Trend studies, preserving original CSVs."""
from pathlib import Path
import json,csv,time,shutil
import numpy as np
HERE=Path(__file__).resolve().parent


def main():
    out=HERE/'outputs';extra=HERE/'adaptive_outputs';base=json.loads((out/'manifest.json').read_text(encoding='utf-8'))
    with (out/'per_seed.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
    if len(rows)!=32000:raise ValueError('Complete Gaussian study required')
    manuscripts=[];adaptive_rows=[]
    if (extra/'schedule.json').exists():
        schedule=json.loads((extra/'schedule.json').read_text(encoding='utf-8'))
        folders=[(extra/'chunks'/str(p['seed']),p['seed'],p['count']) for p in schedule['parts']]
    else:folders=[(extra,20261028,20)] if (extra/'manifest.json').exists() else [(extra/'chunks'/str(seed),seed,10) for seed in (20261028,20261038)]
    for folder,seed,count in folders:
        m=json.loads((folder/'manifest.json').read_text(encoding='utf-8'))
        if m['status'] not in ('completed','checkpointed_complete') or m['seeds']!=list(range(seed,seed+count)):raise ValueError('Incomplete adaptive seed group')
        with (folder/'per_seed.csv').open(encoding='utf-8-sig') as f:part=list(csv.DictReader(f))
        if len(part)!=640*count:raise ValueError('Adaptive group has wrong record count')
        manuscripts.append(m);adaptive_rows.extend(part)
    if any(m['source_sha256']!=manuscripts[0]['source_sha256'] for m in manuscripts):raise ValueError('Adaptive source mismatch')
    if len(adaptive_rows)!=12800 or sorted(s for m in manuscripts for s in m['seeds'])!=list(range(20261028,20261048)):raise ValueError('Adaptive study must cover exactly twenty disjoint seeds')
    for p,h in manuscripts[0]['source_sha256'].items():
        if p in base['source_sha256'] and base['source_sha256'][p]!=h:raise ValueError('Shared Gaussian/adaptive source changed')
    with (extra/'per_seed.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(adaptive_rows[0]));w.writeheader();w.writerows(adaptive_rows)
    for hs in (.5,1.5):
        path=out/f'first_seed_hs{hs:g}.npz';backup=out/f'first_seed_gaussian_hs{hs:g}.npz'
        if not backup.exists():shutil.copy2(path,backup)
        with np.load(backup) as old:data={k:old[k] for k in old.files}
        with np.load(folders[0][0]/f'first_seed_hs{hs:g}.npz') as new:
            for key in new.files:
                if key in data and not np.array_equal(data[key],new[key]):raise ValueError(f'Mismatched shared truth: {key}')
                data[key]=new[key]
        np.savez_compressed(path,**data)
    allrows=rows+adaptive_rows
    with (out/'per_seed_all.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(allrows)
    manifest=base.copy();manifest['methods']=list(base['methods'])+manuscripts[0]['methods'];manifest['records']=len(allrows)
    manifest['adaptive_parameters']=manuscripts[0]['parameters'];manifest['adaptive_worker_elapsed_s']=[m['elapsed_s'] for m in manuscripts]
    manifest['source_sha256'].update(manuscripts[0]['source_sha256'])
    logs=list((out/'chunks').glob('*/full_run.log'));manifest['overall_elapsed_s']=time.time()-min(p.stat().st_ctime for p in logs)
    manifest['boundary']+='; T-Trend extension uses Student-t acceleration impulses plus background GP, Gaussian matched-unconditional-covariance ablation, no robust measurement rejection.'
    (out/'manifest_all.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='completed',records=len(allrows),methods=7,seeds=20),indent=2))


if __name__=='__main__':main()
