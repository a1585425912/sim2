"""Apply all twenty seeds of corrected runtime-CI evaluation transparently."""
from pathlib import Path
import csv,json,hashlib,shutil
HERE=Path(__file__).resolve().parent;OUT=HERE/'outputs'


def read(path):
    with path.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))


def key(r):return (int(r['seed']),float(r['hs']),r['trajectory'],r['method'],float(r['frequency_mhz']),r['phase'],r['timing'])


def main():
    corrections=[]
    for seed in (20261028,20261033,20261038,20261043):
        part=read(OUT/'corrections'/f'rows_{seed}.csv')
        if len(part)!=800:raise ValueError('Each correction group requires every complete scene')
        corrections.extend(part)
    if len({key(r) for r in corrections})!=3200:raise ValueError('Wrong correction coverage')
    path=OUT/'per_seed_all.csv';backup=OUT/'corrections/per_seed_all_before_runtime_ci_fix.csv'
    if not backup.exists():shutil.copy2(path,backup)
    rows=read(backup);replacement={key(r):r for r in corrections};changed=0
    for i,r in enumerate(rows):
        if key(r) in replacement:rows[i]=replacement[key(r)];changed+=1
    if changed!=3200 or len(rows)!=44800:raise ValueError('Unexpected result correction scope')
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    manifest=OUT/'manifest_all.json';original=OUT/'corrections/manifest_all_before_runtime_ci_fix.json'
    if not original.exists():shutil.copy2(manifest,original)
    m=json.loads(original.read_text(encoding='utf-8'));m['original_source_sha256']=m['source_sha256'].copy()
    for p in (HERE/'run_joint.py',HERE/'correct_baseline.py',HERE/'apply_correction.py'):
        m['source_sha256'][str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    m['corrections']=[dict(id='runtime_ci_forecast',records=3200,seeds=20,scope='IMM_CI_linear, measured_cpu only: future_pos_rmse, path_rmse and power/gain metrics',reason='Original runtime evaluation mistakenly used local forecasts; recomputed complete raw-sensor IMM branch forecasts and information-space CI at recorded original CPU horizons',preserved='original child results, original source and merged CSV/manifest snapshots retained; fixed-delay and all GP trajectory estimates unchanged')]
    manifest.write_text(json.dumps(m,indent=2),encoding='utf-8')
    print(json.dumps(dict(status='corrected',records=changed,total=len(rows)),indent=2))


if __name__=='__main__':main()
