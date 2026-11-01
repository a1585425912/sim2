"""Full matched Gaussian-background / T-Trend process comparison, 20 seeds."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(name,'1')
import sys,json,time,hashlib,argparse,gc
from pathlib import Path
from dataclasses import asdict
HERE=Path(__file__).resolve().parent;RAW=HERE.parent/'raw_sensor_maritime';ROOT=HERE.parents[1]
for p in (ROOT,HERE.parent/'imm_ci',HERE.parent/'imm_dpc',RAW,HERE):sys.path.insert(0,str(p))
import numpy as np
from physics import Config,SeaTruth,sensor_record,rotate,spherical
from navigation import navigation_factors
from adaptive import fit,AdaptiveParameters
from sequential import predict
from target import wrap
import run_joint as shared

METHODS=('JointGP_Gaussian_background','JointGP_TTrend')


def run_case(cfg,seed,out,keep):
    tic=time.perf_counter();scene=SeaTruth(cfg,seed);packets,evaluation=sensor_record(cfg,scene,seed+100000)
    nav=navigation_factors(cfg,packets);nav_seconds=time.perf_counter()-tic;t=packets['frame_time'];h=cfg.horizons()
    std=np.array([cfg.range_std,*([np.deg2rad(cfg.angle_std_deg)]*2)]);noise=np.random.default_rng(seed+200000).normal(size=(len(t),cfg.nodes,3))*std
    masks=dict(all=np.ones(len(t),bool),startup=t<5,after_startup=t>=5,acceleration=(t>=15)&(t<20),turn_entry=(t>=35)&(t<40),turn_steady=(t>=40)&(t<55),braking=(t>=65)&(t<70),hover_end=(t>=75)&(t<90))
    saved=dict(time_s=t) if keep else None;rows=[]
    for kind in shared.KINDS:
        p,v=shared.reference.trajectory(t,kind,cfg)
        radar=spherical(rotate(evaluation['frame']['R'].swapaxes(-1,-2),p[:,None]-evaluation['frame']['p']))+noise;radar[...,1:]=wrap(radar[...,1:])
        for method,student in zip(METHODS,(False,True)):
            hist=fit(cfg,t,radar,nav,packets,student=student);seconds=hist['seconds']+nav_seconds/len(t)
            for timing in ('fixed_delay','measured_cpu'):
                hh=np.broadcast_to(h,(len(t),len(h))) if timing=='fixed_delay' else cfg.rounds*cfg.round_s+np.maximum(cfg.processing_s,seconds)[:,None]+np.linspace(0,cfg.slot_s,cfg.slot_samples)
                future,_=shared.reference.trajectory(t[:,None]+hh,kind,cfg);true_d=np.linalg.norm(future[:,:,None]-scene.at(t[:,None]+hh)['p'],axis=-1)
                forecasts=[predict(hist,hh[:,j],cfg,acceleration=True,covariance=False) for j in range(len(h))]
                predictions=np.stack([a[0] for a in forecasts],1);cmd=np.stack([a[2] for a in forecasts],1)
                shared.score(cfg,seed,kind,t,p,v,future,true_d,method,hist['state'],hist['covariance'],predictions,cmd,None,seconds,masks,saved,rows,timing=timing,horizons=hh)
            if keep:
                saved[f'{kind}__truth_p']=p;saved[f'{kind}__truth_v']=v;saved[f'{kind}__{method}__scale']=hist['precision_scale'];saved[f'{kind}__{method}__cost']=hist['cost']
            print(f'  {kind} {method} seed={seed} Hs={cfg.hs:g}; solver={hist["seconds"].mean()*1000:.2f}ms elapsed={time.perf_counter()-tic:.1f}s',flush=True)
            del hist,forecasts;gc.collect()
    if keep:np.savez_compressed(out/f'first_seed_hs{cfg.hs:g}.npz',**saved)
    return rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=20);ap.add_argument('--seed-start',type=int,default=20261028);ap.add_argument('--output',type=Path,default=HERE/'adaptive_outputs');args=ap.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True);rows=[];start=time.perf_counter()
    for i in range(args.trials):
        seed=args.seed_start+i
        for hs in (.5,1.5):
            rows+=run_case(Config(hs=hs),seed,out,seed==20261028)
            shared.reference.csv_write(out/'per_seed.csv',rows)
            (out/'progress.json').write_text(json.dumps(dict(completed_cases=2*i+int(hs==1.5)+1,total_cases=2*args.trials,elapsed_s=time.perf_counter()-start),indent=2),encoding='utf-8')
        print(f'SEED {i+1}/{args.trials} complete {time.perf_counter()-start:.1f}s',flush=True)
    files=[HERE/p for p in ('adaptive.py','run_adaptive.py','sequential.py','solver.py','navigation.py','run_joint.py')]+[RAW/p for p in ('physics.py','estimation.py','target.py','run_full.py')]
    manifest=dict(status='completed',trials=args.trials,seeds=list(range(args.seed_start,args.seed_start+args.trials)),elapsed_s=time.perf_counter()-start,parameters=asdict(AdaptiveParameters()),source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},methods=METHODS,boundary='Gaussian background GP plus Student-t acceleration impulse, persistent joint arrival; Gaussian-matched-covariance ablation; raw Gaussian sensors unchanged; NAV linearized once then radar MAP; conditional scale covariance, not automatic calibration')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')


if __name__=='__main__':main()
