"""Complete paired FMCW experiment: 20 seeds, two seas, five full 90s paths."""
import bootstrap
import sys,time,json,csv,hashlib,argparse,gc
from pathlib import Path
from dataclasses import asdict,replace
import numpy as np
from physics import Config,SeaTruth,sensor_record,rotate,spherical
from navigation import navigation_factors
from solver import Parameters
import run_joint as original
from target import wrap
from waveform import Waveform,buoy_bursts,observe
from joint_filter import fit,predict

HERE=Path(__file__).resolve().parent
KINDS=original.reference.KINDS
TRAJECTORY=original.reference.trajectory
METHODS={
    'GP_noD_q1':(1.,False,False),
    'GP_Doppler_q1':(1.,True,False),
    'GP_noD_q001':(.01,False,False),
    'GP_Doppler_q001':(.01,True,False),
    'GP_Doppler_ignore_buoy':(1.,True,True),
}


def write(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def run_case(cfg,seed,out,keep):
    tic=time.perf_counter();w=Waveform();scene=SeaTruth(cfg,seed)
    packets,evaluation=sensor_record(cfg,scene,seed+100000)
    ntic=time.perf_counter();nav=navigation_factors(cfg,packets);nav_s=time.perf_counter()-ntic
    t=packets['frame_time'];truth_nav=evaluation['frame']
    burst_p,burst_v=buoy_bursts(scene,t,w)
    radar_std=np.array([cfg.range_std,*([np.deg2rad(cfg.angle_std_deg)]*2)])
    # Exactly the previous radar draws for the same seed; supplementary FMCW
    # channel has independent thermal noise and supplies Doppler only.
    noise=np.random.default_rng(seed+200000).normal(size=(len(t),cfg.nodes,3))*radar_std
    masks=dict(all=np.ones(len(t),bool),startup=t<5,after_startup=t>=5,
               acceleration=(t>=15)&(t<20),turn_entry=(t>=35)&(t<40),
               turn_steady=(t>=40)&(t<55),braking=(t>=65)&(t<70),hover_end=(t>=75)&(t<90))
    saved=dict(time_s=t) if keep else None;rows=[];sensor_rows=[]
    for ik,kind in enumerate(KINDS):
        p,v=TRAJECTORY(t,kind,cfg)
        radar=spherical(rotate(truth_nav['R'].swapaxes(-1,-2),p[:,None]-truth_nav['p']))+noise
        radar[...,1:]=wrap(radar[...,1:])
        doppler=observe(cfg,t,kind,TRAJECTORY,burst_p,burst_v,seed+400000+ik*10000,w,keep)
        for phase,mask in masks.items():
            e=doppler['error'][mask]
            sensor_rows.append(dict(seed=seed,hs=cfg.hs,trajectory=kind,phase=phase,
                velocity_rmse=float(np.sqrt(np.mean(e**2))),velocity_bias=float(e.mean()),
                absolute_error_p95=float(np.quantile(np.abs(e),.95)),
                nominal95_coverage=float(np.mean(np.abs(e)<=1.96*w.velocity_std)),
                adc_snr_mean_db=float(doppler['snr_db'][mask].mean()),
                processing_mean_ms=float(doppler['seconds'][mask].mean()*1000)))
        if keep:
            saved[f'{kind}__truth_p']=p;saved[f'{kind}__truth_v']=v
            saved[f'{kind}__doppler_measured']=doppler['velocity']
            saved[f'{kind}__doppler_truth']=doppler['truth_velocity']
            saved[f'{kind}__doppler_snr_db']=doppler['snr_db']
            for name,value in doppler['sample'].items():saved[f'{kind}__{name}']=value
        print(f'SIGNAL seed={seed} Hs={cfg.hs:g} {kind}: vRMSE {np.sqrt(np.mean(doppler["error"]**2)):.4f} m/s, wall {time.perf_counter()-tic:.1f}s',flush=True)
        for method,(q,useD,ignore) in METHODS.items():
            # Two retained target clones: all older Gaussian information and
            # cross-covariance remain in the persistent arrival distribution.
            par=Parameters(window=2,jerk_psd=q)
            history=fit(cfg,t,radar,nav,par=par,packets=packets,
                        doppler=doppler if useD else None,ignore_buoy=ignore)
            # Paired no-D controls collect/process the identical supplementary
            # waveform too; only the likelihood changes. No timing advantage.
            seconds=history['seconds']+nav_s/len(t)+doppler['seconds']
            for timing in ('fixed_delay','measured_cpu'):
                compute=np.full(len(t),cfg.processing_s) if timing=='fixed_delay' else np.maximum(cfg.processing_s,seconds)
                hh=w.half_span_s+cfg.rounds*cfg.round_s+compute[:,None]+np.linspace(0,cfg.slot_s,cfg.slot_samples)
                truth_future,_=TRAJECTORY(t[:,None]+hh,kind,cfg)
                distances=np.linalg.norm(truth_future[:,:,None]-scene.at(t[:,None]+hh)['p'],axis=-1)
                forecasts=[predict(history,hh[:,j],cfg,covariance=False) for j in range(cfg.slot_samples)]
                predictions=np.stack([a[0] for a in forecasts],1);cmd=np.stack([a[2] for a in forecasts],1)
                original.score(cfg,seed,kind,t,p,v,truth_future,distances,method,
                    history['state'],history['covariance'],predictions,cmd,None,
                    seconds,masks,saved,rows,timing=timing,horizons=hh)
                if keep and timing=='fixed_delay':
                    saved[f'{kind}__{method}__compute_s']=seconds
                    saved[f'{kind}__{method}__future_error']=predictions-truth_future
            print(f'FIT seed={seed} Hs={cfg.hs:g} {kind} {method}: {history["seconds"].mean()*1000:.2f}ms/frame wall {time.perf_counter()-tic:.1f}s',flush=True)
            del history,forecasts;gc.collect()
    if keep:np.savez_compressed(out/f'first_seed_hs{cfg.hs:g}.npz',**saved)
    return rows,sensor_rows


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=20)
    ap.add_argument('--seed-start',type=int,default=20261028)
    ap.add_argument('--output',type=Path,default=HERE/'outputs');args=ap.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True);tic=time.perf_counter();rows=[];sensors=[]
    # Hash inputs at start as well as end, detect edits during a formal run.
    sources=[HERE/f for f in ('bootstrap.py','waveform.py','joint_filter.py','run_full.py')]+[bootstrap.PARENT/'joint_gp'/f for f in ('solver.py','navigation.py','run_joint.py','sequential.py')]+[bootstrap.PARENT/'raw_sensor_maritime'/f for f in ('physics.py','estimation.py','target.py','run_full.py')]
    source_hash={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    for i in range(args.trials):
        seed=args.seed_start+i
        for hs in (.5,1.5):
            r,s=run_case(Config(hs=hs),seed,out,seed==20261028);rows+=r;sensors+=s
            write(out/'per_seed.csv',rows);write(out/'sensor_per_seed.csv',sensors)
            progress=dict(status='running',completed_cases=2*i+1+int(hs==1.5),total_cases=2*args.trials,last_seed=seed,last_hs=hs,elapsed_s=time.perf_counter()-tic)
            (out/'progress.json').write_text(json.dumps(progress,indent=2),encoding='utf-8')
        print(f'SEED {i+1}/{args.trials} complete',flush=True)
    unchanged=all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in source_hash.items())
    if not unchanged:raise RuntimeError('Source changed during formal execution')
    manifest=dict(status='completed',seeds=list(range(args.seed_start,args.seed_start+args.trials)),trials=args.trials,rows=len(rows),sensor_rows=len(sensors),elapsed_s=time.perf_counter()-tic,
        methods=METHODS,waveform=Waveform().metadata(),configs=[asdict(Config(hs=h)) for h in (.5,1.5)],source_sha256=source_hash,
        scope='Centralized joint Gaussian, point body supplementary FMCW; same range/angles and nav as old experiment; all frames, uniform 20 nodes, no truth gating/initialization; CPI centered with causal availability included; raw ADC processing measured, waveform synthesis is simulator time.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')


if __name__=='__main__':main()
