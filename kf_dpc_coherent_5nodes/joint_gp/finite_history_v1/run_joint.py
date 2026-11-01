"""Full matched research experiment. No smoke, gating or truth tuning."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(name,'1')
import sys,time,json,csv,hashlib,argparse,gc
from pathlib import Path
from dataclasses import asdict
HERE=Path(__file__).resolve().parent;RAW=HERE.parent/'raw_sensor_maritime';ROOT=HERE.parents[1]
for path in (ROOT,HERE.parent/'imm_ci',HERE.parent/'imm_dpc',RAW,HERE):sys.path.insert(0,str(path))
import numpy as np
from physics import Config,SeaTruth,sensor_record,rotate,spherical
from navigation import navigation_factors
from solver import fit,predict,Parameters as GPParameters
from target import fit_radar,wrap
from imm_filter import Parameters as IMMParameters,predict as imm_predict
from fusion import forecast,fuse,kernel
import run_full as reference

METHODS=('IMM_CI_linear','GP_independent_nav','JointGP_linear','JointGP_accel','JointGP_power')
KINDS=reference.KINDS


def power(cmd,true_d,frequency,total_power,amplitude_weights=None):
    lam=299792458./(frequency*1e6);N=true_d.shape[-1]
    channel=lam/(4*np.pi*true_d)
    if amplitude_weights is None:amplitude_weights=np.full((len(cmd),N),1/np.sqrt(N))
    amplitude=np.sqrt(total_power)*amplitude_weights[:,None]*channel
    field=amplitude*np.exp(1j*2*np.pi/lam*(cmd-true_d))
    rx=np.abs(field.sum(-1))**2
    inc=np.sum(amplitude**2,-1)
    uniform_ideal=total_power/N*channel.sum(-1)**2
    own_ideal=amplitude.sum(-1)**2
    optimal_ideal=total_power*np.sum(channel**2,-1)
    return dict(rx=rx,inc=inc,ideal=uniform_ideal,own=own_ideal,optimal=optimal_ideal)


def covariance_beamformer(cmd,C,frequency):
    """One power allocation per frame; total power fixed, no node gating.

    Principal eigenvector of time-averaged conditional channel second moment.
    All phase commands remain predicted geometry. No actual channel/truth used.
    """
    k=2*np.pi*frequency*1e6/299792458.
    diagonal=np.diagonal(C,axis1=-2,axis2=-1)
    pairvar=diagonal[..., :,None]+diagonal[...,None,:]-2*C
    # Tiny negative values can arise only from floating point cancellation.
    pairvar=np.maximum(pairvar,0.)
    mean_channel=(299792458./(frequency*1e6))/(4*np.pi*cmd)
    second=mean_channel[..., :,None]*mean_channel[...,None,:]*np.exp(-.5*k*k*pairvar)
    second=second.mean(axis=1)
    _,vectors=np.linalg.eigh(second)
    a=np.abs(vectors[...,-1]);return a/np.linalg.norm(a,axis=-1,keepdims=True)


def score(cfg,seed,kind,t,truth_p,truth_v,pfuture,true_d,method,state,P,predictions,cmd,cov_distance,seconds,masks,saved,rows,timing='fixed_delay',horizons=None):
    if state.ndim==2:state=np.broadcast_to(state[:,None],(len(t),cfg.nodes,state.shape[-1]))
    if P.ndim==3:P=np.broadcast_to(P[:,None],(len(t),cfg.nodes,P.shape[-1],P.shape[-1]))
    ep=state[...,:3]-truth_p[:,None];ev=state[...,3:6]-truth_v[:,None]
    err=np.concatenate([ep,ev],-1)
    if predictions.ndim==3:predictions=np.broadcast_to(predictions[:,:,None,:3],(len(t),cfg.slot_samples,cfg.nodes,3))
    future_error=predictions[...,:3]-pfuture[:,:,None]
    path=cmd-true_d;centered_path=path-path.mean(-1,keepdims=True)
    error_rate=np.diff(ep,axis=0)/cfg.dt
    velocity_error_rate=np.diff(ev,axis=0)/cfg.dt
    consistency=np.diff(state[...,:3],axis=0)/cfg.dt-.5*(state[1:,...,3:6]+state[:-1,...,3:6])
    if saved is not None and timing=='fixed_delay':
        saved[f'{kind}__{method}__state']=state[...,:6].copy()
        saved[f'{kind}__{method}__sigma']=np.sqrt(np.diagonal(P[...,:6,:6],axis1=-2,axis2=-1))
    for f in cfg.frequencies_mhz:
        weights=covariance_beamformer(cmd,cov_distance,f) if method=='JointGP_power' else None
        pw=power(cmd,true_d,f,cfg.total_power,weights)
        if saved is not None and timing=='fixed_delay':
            saved[f'{kind}__{method}__eta{int(f)}']=pw['rx'].sum(1)/pw['ideal'].sum(1)
            if weights is not None:saved[f'{kind}__{method}__weights{int(f)}']=weights
        for phase,mask in masks.items():
            nees,coverage=reference.metrics(err[mask],P[mask,:,:6,:6])
            neesp,coveragep=reference.metrics(ep[mask],P[mask,:,:3,:3])
            pairs=mask[1:]&mask[:-1]
            rms=lambda a:float(np.sqrt(np.mean(np.sum(a*a,axis=-1))))
            gain=float(pw['rx'][mask].sum()/pw['inc'][mask].sum())
            row=dict(seed=seed,hs=cfg.hs,trajectory=kind,method=method,frequency_mhz=f,phase=phase,timing=timing,
                     pos_rmse=rms(ep[mask]),vel_rmse=rms(ev[mask]),future_pos_rmse=rms(future_error[mask]),
                     path_rmse=float(np.sqrt(np.mean(centered_path[mask]**2))),position_error_rate=rms(error_rate[pairs]),velocity_error_rate=rms(velocity_error_rate[pairs]),kinematic_defect=rms(consistency[pairs]),
                     eta=float(pw['rx'][mask].sum()/pw['ideal'][mask].sum()),eta_own=float(pw['rx'][mask].sum()/pw['own'][mask].sum()),eta_optimal_bound=float(pw['rx'][mask].sum()/pw['optimal'][mask].sum()),gain_linear=gain,gain_db=float(10*np.log10(gain)),rx_mean_w=float(pw['rx'][mask].mean()),
                     nees6=nees,coverage95_6=coverage,nees_position=neesp,coverage95_position=coveragep,
                     compute_mean_ms=float(np.mean(seconds[mask])*1000),compute_p95_ms=float(np.quantile(seconds[mask],.95)*1000),
                     effective_nodes=float(cfg.nodes if weights is None else np.mean(1/np.sum(weights[mask]**4,axis=-1))),
                     tx_end_mean_ms=float(np.mean(horizons[mask,-1])*1000) if horizons is not None else float(cfg.horizons()[-1]*1000))
            for axis,name in enumerate('XYZ'):
                row[f'p{name}_rmse']=float(np.sqrt(np.mean(ep[mask,:,axis]**2)));row[f'v{name}_rmse']=float(np.sqrt(np.mean(ev[mask,:,axis]**2)))
            rows.append(row)


def run_case(cfg,seed,out,keep):
    tic=time.perf_counter();scene=SeaTruth(cfg,seed);packets,evaluation=sensor_record(cfg,scene,seed+100000)
    nav=navigation_factors(cfg,packets);nav_seconds=time.perf_counter()-tic
    t=packets['frame_time'];horizons=cfg.horizons();truth_nav=evaluation['frame']
    true_buoy=scene.at(t[:,None]+horizons)['p']
    radar_std=np.array([cfg.range_std,*([np.deg2rad(cfg.angle_std_deg)]*2)])
    rng=np.random.default_rng(seed+200000);noise=rng.normal(size=(len(t),cfg.nodes,3))*radar_std
    C=kernel(reference.previous.ring_weights(cfg.nodes),cfg.rounds)
    masks=dict(all=np.ones(len(t),bool),startup=t<5,after_startup=t>=5,acceleration=(t>=15)&(t<20),turn_entry=(t>=35)&(t<40),turn_steady=(t>=40)&(t<55),braking=(t>=65)&(t<70),hover_end=(t>=75)&(t<90))
    saved=dict(time_s=t) if keep else None;rows=[];diagnostics=[]
    for kind in KINDS:
        p,v=reference.trajectory(t,kind,cfg);pfuture,_=reference.trajectory(t[:,None]+horizons,kind,cfg)
        radar=spherical(rotate(truth_nav['R'].swapaxes(-1,-2),p[:,None]-truth_nav['p']))+noise;radar[...,1:]=wrap(radar[...,1:])
        true_d=np.linalg.norm(pfuture[:,:,None]-true_buoy,axis=-1)
        if keep:saved[f'{kind}__truth_p']=p;saved[f'{kind}__truth_v']=v
        baseline_start=time.perf_counter()
        imm=fit_radar(t,radar,nav,np.diag(radar_std**2),par=IMMParameters())
        forecasts=[forecast(imm,float(h)) for h in (0.,*horizons)]
        states,Ps=fuse(np.stack([a[0] for a in forecasts],1),np.stack([a[1] for a in forecasts],1),C,'ci')
        baseline_seconds=np.full(len(t),(time.perf_counter()-baseline_start+nav_seconds)/len(t))
        buoy=nav['p'][:,None]+horizons[None,:,None,None]*nav['v'][:,None]
        cmd=np.linalg.norm(states[:,1:,:,:3]-buoy,axis=-1)
        score(cfg,seed,kind,t,p,v,pfuture,true_d,'IMM_CI_linear',states[:,0],Ps[:,0],states[:,1:,:,:3],cmd,None,baseline_seconds,masks,saved,rows)
        runtime_h=cfg.rounds*cfg.round_s+np.maximum(cfg.processing_s,baseline_seconds)[:,None]+np.linspace(0,cfg.slot_s,cfg.slot_samples)
        runtime_predictions=np.empty((len(t),cfg.slot_samples,cfg.nodes,3))
        for k in range(len(t)):
            local={name:imm[name][k] for name in ('branches','probabilities')};local['enabled']=imm['enabled']
            for j in range(cfg.slot_samples):runtime_predictions[k,j]=imm_predict(local,float(runtime_h[k,j]))[:,:3]
        runtime_buoy=nav['p'][:,None]+runtime_h[:,:,None,None]*nav['v'][:,None]
        runtime_cmd=np.linalg.norm(runtime_predictions-runtime_buoy,axis=-1)
        runtime_truth,_=reference.trajectory(t[:,None]+runtime_h,kind,cfg)
        runtime_d=np.linalg.norm(runtime_truth[:,:,None]-scene.at(t[:,None]+runtime_h)['p'],axis=-1)
        score(cfg,seed,kind,t,p,v,runtime_truth,runtime_d,'IMM_CI_linear',states[:,0],Ps[:,0],runtime_predictions,runtime_cmd,None,baseline_seconds,masks,saved,rows,timing='measured_cpu',horizons=runtime_h)
        if keep:
            with np.load(RAW/'outputs'/f'first_seed_hs{cfg.hs:g}.npz') as old:
                diagnostics.append(dict(seed=seed,hs=cfg.hs,trajectory=kind,baseline_max_state_difference=float(np.max(np.abs(states[:,0]-old[f'{kind}__IMM_CI_linear__state'])))))
        del imm,states,Ps,forecasts
        for joint in (False,True):
            label='JOINT' if joint else 'GP independent'
            history=fit(cfg,t,radar,nav,joint=joint)
            seconds=history['seconds']+nav_seconds/len(t)
            for timing in ('fixed_delay','measured_cpu'):
                hh=np.broadcast_to(horizons,(len(t),len(horizons))) if timing=='fixed_delay' else cfg.rounds*cfg.round_s+np.maximum(cfg.processing_s,seconds)[:,None]+np.linspace(0,cfg.slot_s,cfg.slot_samples)
                if timing=='fixed_delay':truth_future,distances=pfuture,true_d
                else:
                    truth_future,_=reference.trajectory(t[:,None]+hh,kind,cfg)
                    distances=np.linalg.norm(truth_future[:,:,None]-scene.at(t[:,None]+hh)['p'],axis=-1)
                forecasts=[predict(history,hh[:,j],cfg,acceleration=joint,covariance=joint) for j in range(len(horizons))]
                predictions=np.stack([a[0] for a in forecasts],axis=1)
                cmd=np.stack([a[2] for a in forecasts],axis=1);distance_cov=np.stack([a[3] for a in forecasts],axis=1) if joint else None
                if joint:
                    linear=[predict(history,hh[:,j],cfg,acceleration=False,covariance=False)[2] for j in range(len(horizons))]
                    linear_cmd=np.stack(linear,axis=1)
                    score(cfg,seed,kind,t,p,v,truth_future,distances,'JointGP_linear',history['state'],history['covariance'],predictions,linear_cmd,None,seconds,masks,saved,rows,timing=timing,horizons=hh)
                    for method in ('JointGP_accel','JointGP_power'):
                        score(cfg,seed,kind,t,p,v,truth_future,distances,method,history['state'],history['covariance'],predictions,cmd,distance_cov,seconds,masks,saved,rows,timing=timing,horizons=hh)
                else:score(cfg,seed,kind,t,p,v,truth_future,distances,'GP_independent_nav',history['state'],history['covariance'],predictions,cmd,None,seconds,masks,saved,rows,timing=timing,horizons=hh)
                if keep and joint and timing=='fixed_delay':saved[f'{kind}__distance_covariance']=distance_cov
            if keep:
                saved[f'{kind}__{label}__cost']=history['cost'];saved[f'{kind}__{label}__compute_s']=seconds
                if joint:
                    saved[f'{kind}__joint_buoy_p']=np.array([n['p'] for n in history['nav']])
            print(f'  {kind} {label} seed={seed} Hs={cfg.hs:g}; solver mean {history["seconds"].mean()*1000:.2f}ms; elapsed {time.perf_counter()-tic:.1f}s',flush=True)
            del history,forecasts;gc.collect()
    if keep:np.savez_compressed(out/f'first_seed_hs{cfg.hs:g}.npz',**saved)
    return rows,diagnostics


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=20);ap.add_argument('--seed-start',type=int,default=20261028);ap.add_argument('--output',type=Path,default=HERE/'outputs');args=ap.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True);start=time.perf_counter();rows=[];diagnostics=[]
    for i in range(args.trials):
        seed=args.seed_start+i
        for hs in (.5,1.5):
            r,d=run_case(Config(hs=hs),seed,out,seed==20261028);rows+=r;diagnostics+=d
            reference.csv_write(out/'per_seed.csv',rows)
            (out/'progress.json').write_text(json.dumps(dict(completed_cases=2*i+int(hs==1.5)+1,total_cases=2*args.trials,elapsed_s=time.perf_counter()-start,last_seed=seed,last_hs=hs),indent=2),encoding='utf-8')
        print(f'SEED {i+1}/{args.trials} complete {time.perf_counter()-start:.1f}s',flush=True)
    files=list(HERE.glob('*.py'))+[RAW/p for p in ('physics.py','estimation.py','target.py','run_full.py')]+[HERE.parent/'imm_dpc/imm_filter.py',HERE.parent/'imm_ci/fusion.py']
    manifest=dict(status='completed',trials=args.trials,seeds=list(range(args.seed_start,args.seed_start+args.trials)),configs=[asdict(Config(hs=h)) for h in (.5,1.5)],gp=asdict(GPParameters()),elapsed_s=time.perf_counter()-start,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},baseline_checks=diagnostics,python=sys.version,numpy=np.__version__,methods=METHODS,
                  boundary='Centralized joint finite-history MAP, not finite-round distributed GP; same raw sensors/seeds, no post-hoc smoothing, no gating. Communication/CPU timing assumed, measured compute reported. Power variant maximizes Gaussian-posterior expected RX power with variable per-node power and fixed total, not guaranteed actual global optimum.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')


if __name__=='__main__':main()
