"""Finish two preassigned seed prefixes before redistributing remaining work.

This is scheduling only: complete seeds, all scenes/frames/methods retained.
Original process logs and any overlapping progress are preserved explicitly.
Run only after each of the two existing processes has completed two seeds.
"""
from pathlib import Path
import csv,json,hashlib,shutil
from dataclasses import asdict
HERE=Path(__file__).resolve().parent;RAW=HERE.parent/'raw_sensor_maritime'


def main():
    root=HERE/'adaptive_outputs';parts=[]
    source_files=[HERE/p for p in ('adaptive.py','run_adaptive.py','sequential.py','solver.py','navigation.py','run_joint.py')]+[RAW/p for p in ('physics.py','estimation.py','target.py','run_full.py')]
    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
    parameters=dict(background_jerk_psd=.001,jump_scale_std=.02,degrees_of_freedom=3.,iterations=3,starts=[1.,.001])
    for seed in (20261028,20261038):
        folder=root/'chunks'/str(seed);path=folder/'per_seed.csv'
        with path.open(encoding='utf-8-sig') as f:allrows=list(csv.DictReader(f))
        rows=[r for r in allrows if int(r['seed']) in (seed,seed+1)]
        if len(rows)!=1280:raise ValueError('Both complete preassigned seeds required')
        shutil.copy2(path,folder/'original_progress_before_resharding.csv')
        with path.open('w',encoding='utf-8-sig',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
        progress=json.loads((folder/'progress.json').read_text(encoding='utf-8'))
        m=dict(status='checkpointed_complete',trials=2,seeds=[seed,seed+1],elapsed_s=progress['elapsed_s'],parameters=parameters,source_sha256=hashes,methods=['JointGP_Gaussian_background','JointGP_TTrend'],boundary='Explicit scheduling checkpoint of two complete seeds from originally ten-seed worker; all remaining seeds assigned disjointly before results, no scene/frame selection')
        (folder/'manifest.json').write_text(json.dumps(m,indent=2),encoding='utf-8');parts.append(dict(seed=seed,count=2))
    parts.extend([dict(seed=s,count=n) for s,n in ((20261030,3),(20261033,3),(20261036,2),(20261040,3),(20261043,3),(20261046,2))])
    (root/'schedule.json').write_text(json.dumps(dict(parts=parts,reason='Redistribute remaining independent seeds after Gaussian workers release memory; estimator and parameters unchanged'),indent=2),encoding='utf-8')
    print('Two complete prefixes checkpointed; remaining sixteen seeds assigned to six disjoint formal workers.')


if __name__=='__main__':main()
