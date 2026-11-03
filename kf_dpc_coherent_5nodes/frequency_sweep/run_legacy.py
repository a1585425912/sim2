"""Replay comparable legacy models fully; assess tracking before low-f gain."""
import sys,importlib.util,time,json,csv,hashlib,argparse,gc
from pathlib import Path
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frequency_core',HERE/'run_full.py')
core=importlib.util.module_from_spec(spec);spec.loader.exec_module(core)
import numpy as np
from physics import Config,SeaTruth,sensor_record,rotate,spherical
from navigation import navigation_factors
from target import fit_radar,wrap
from fusion import forecast,fuse,kernel
from imm_filter import Parameters as IMMParameters
from estimation import WavePredictor
from sequential import fit as independent_fit,predict
from adaptive import fit as adaptive_fit
from solver import Parameters
import run_joint as original

METHODS=('CV_CI_linear','CV_CI_wave','IMM_DPC_linear','IMM_CI_linear','IMM_CI_wave',
         'GP_independent_nav','JointGP_Gaussian_background','JointGP_TTrend')

def masks(t):
    return dict(all=np.ones(len(t),bool),startup=t<5,after_startup=t>=5,
        acceleration=(t>=15)&(t<20),turn_entry=(t>=35)&(t<40),turn_steady=(t>=40)&(t<55),
        braking=(t>=65)&(t<70),hover_end=(t>=75)&(t<90))

def tracking(cfg,seed,kind,method,t,p,v,state,P):
    if state.ndim==2:state=state[:,None]
    if P.ndim==3:P=P[:,None]
    ep=state[...,:3]-p[:,None];ev=state[...,3:6]-v[:,None];err=np.concatenate([ep,ev],-1)
    derivative=np.diff(ep,axis=0)/cfg.dt;dv=np.diff(ev,axis=0)/cfg.dt
    defect=np.diff(state[...,:3],axis=0)/cfg.dt-.5*(state[1:,...,3:6]+state[:-1,...,3:6])
    rms=lambda x:float(np.sqrt(np.mean(np.sum(x*x,axis=-1))))
    out=[]
    for phase,mask in masks(t).items():
        pairs=mask[1:]&mask[:-1]
        n6,c6=original.reference.metrics(err[mask],P[mask,...,:6,:6])
        n3,c3=original.reference.metrics(ep[mask],P[mask,...,:3,:3])
        row=dict(seed=seed,hs=cfg.hs,trajectory=kind,method=method,phase=phase,
            pos_rmse=rms(ep[mask]),vel_rmse=rms(ev[mask]),position_error_rate=rms(derivative[pairs]),
            velocity_error_rate=rms(dv[pairs]),kinematic_defect=rms(defect[pairs]),
            nees6=n6,coverage95_6=c6,nees_position=n3,coverage95_position=c3)
        for j,a in enumerate('XYZ'):
            row['p'+a+'_rmse']=float(np.sqrt(np.mean(ep[mask,...,j]**2)))
            row['v'+a+'_rmse']=float(np.sqrt(np.mean(ev[mask,...,j]**2)))
        out.append(row)
    return out

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=20);ap.add_argument('--seed-start',type=int,default=20261028)
    ap.add_argument('--output',type=Path,default=HERE/'outputs/legacy');args=ap.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True);cache=out/'geometry';cache.mkdir(exist_ok=True)
    parent=HERE.parent;raw=parent/'raw_sensor_maritime';joint=parent/'joint_gp'
    sources=[Path(__file__),HERE/'run_full.py',parent/'fmcw_doppler/bootstrap.py',parent/'fmcw_doppler/waveform.py',
        *[raw/f for f in ('physics.py','estimation.py','target.py','run_full.py')],
        *[joint/f for f in ('solver.py','navigation.py','adaptive.py','sequential.py','run_joint.py')],
        parent/'imm_ci/fusion.py',parent/'imm_dpc/imm_filter.py',parent/'imm_dpc/run_experiment.py',parent/'gm_dpc_unified/experiment.py']
    source_hash={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    old={}
    for path in (raw/'outputs/per_seed.csv',joint/'outputs/per_seed_all.csv'):
        source_hash[str(path)]=hashlib.sha256(path.read_bytes()).hexdigest()
        with path.open(encoding='utf-8-sig') as f:
            for r in csv.DictReader(f):
                if r['phase']=='all' and float(r['frequency_mhz'])==300 and r.get('timing','fixed_delay')=='fixed_delay':
                    old[(int(r['seed']),float(r['hs']),r['trajectory'],r['method'])]=r
    begin=time.perf_counter();rows=0;trackrows=[];reproductions=[]
    with (out/'per_seed.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.writer(f);writer.writerow(core.FIELDS)
        for trial in range(args.trials):
            seed=args.seed_start+trial
            for hs in (.5,1.5):
                cfg=Config(hs=hs);scene=SeaTruth(cfg,seed);packets,evaluation=sensor_record(cfg,scene,seed+100000)
                nav=navigation_factors(cfg,packets);t=packets['frame_time'];wave=WavePredictor(cfg);wave_history=wave.fit(nav)
                h_common=core.Waveform().half_span_s+cfg.horizons();h_native=cfg.horizons()
                std=np.array([cfg.range_std,*([np.deg2rad(cfg.angle_std_deg)]*2)])
                noise=np.random.default_rng(seed+200000).normal(size=(len(t),cfg.nodes,3))*std
                C=kernel(original.reference.previous.ring_weights(cfg.nodes),cfg.rounds)
                for kind in original.reference.KINDS:
                    p,v=original.reference.trajectory(t,kind,cfg)
                    radar=spherical(rotate(evaluation['frame']['R'].swapaxes(-1,-2),p[:,None]-evaluation['frame']['p']))+noise
                    radar[...,1:]=wrap(radar[...,1:])
                    cv=fit_radar(t,radar,nav,np.diag(std**2),enabled=(1,),par=IMMParameters())
                    imm=fit_radar(t,radar,nav,np.diag(std**2),par=IMMParameters())
                    for method in METHODS:
                        hist=None
                        if method.startswith(('CV_','IMM_')):
                            base=cv if method.startswith('CV_') else imm;fusion='dpc' if 'DPC' in method else 'ci'
                            state,P=fuse(*forecast(base,0.),C,fusion)
                            def commands(h):
                                forecasts=[forecast(base,float(hh)) for hh in h]
                                means,Ps=fuse(np.stack([a[0] for a in forecasts],1),np.stack([a[1] for a in forecasts],1),C,fusion)
                                buoy=np.stack([wave.predict(wave_history,float(hh)) for hh in h],1) if method.endswith('wave') else nav['p'][:,None]+h[None,:,None,None]*nav['v'][:,None]
                                return np.linalg.norm(means[...,:3]-buoy,axis=-1)
                        else:
                            if method=='GP_independent_nav':hist=independent_fit(cfg,t,radar,nav,par=Parameters(window=2),joint=False,packets=packets)
                            else:hist=adaptive_fit(cfg,t,radar,nav,packets,student=method=='JointGP_TTrend')
                            state,P=hist['state'],hist['covariance']
                            def commands(h):
                                return np.stack([predict(hist,float(hh),cfg,acceleration=method!='GP_independent_nav',covariance=False)[2] for hh in h],1)
                        tr=tracking(cfg,seed,kind,method,t,p,v,state,P);trackrows+=tr
                        future,_=original.reference.trajectory(t[:,None]+h_common,kind,cfg)
                        distances=np.linalg.norm(future[:,:,None]-scene.at(t[:,None]+h_common)['p'],axis=-1)
                        cmd=commands(h_common);error=cmd-distances
                        core.evaluate(seed,cfg,kind,method,'fixed_delay',t,error,distances,writer)
                        rows+=len(core.FREQUENCIES)*2
                        np.savez_compressed(cache/f'{seed}_hs{hs:g}_{kind}_{method}.npz',time_s=t,state=state,covariance=P,
                            fixed_delay__error=error,fixed_delay__distance=distances,fixed_delay__horizons=h_common)
                        old_method='CV_CI_wave' if method=='CV_CI_linear' else method
                        previous=old[(seed,hs,kind,old_method)];now=tr[0]
                        native_eta_difference=''
                        if method!='CV_CI_linear':
                            future_native,_=original.reference.trajectory(t[:,None]+h_native,kind,cfg)
                            dn=np.linalg.norm(future_native[:,:,None]-scene.at(t[:,None]+h_native)['p'],axis=-1)
                            power=original.power(commands(h_native),dn,300,cfg.total_power)
                            native_eta_difference=float(power['rx'].sum()/power['ideal'].sum())-float(previous['eta'])
                        reproductions.append(dict(seed=seed,hs=hs,trajectory=kind,method=method,
                            pos_rmse_difference=now['pos_rmse']-float(previous['pos_rmse']),
                            vel_rmse_difference=now['vel_rmse']-float(previous['vel_rmse']),native_eta_difference=native_eta_difference))
                        print(f'DONE seed={seed} Hs={hs:g} {kind} {method} elapsed={time.perf_counter()-begin:.1f}s',flush=True)
                        del hist;gc.collect()
                    del cv,imm;f.flush()
                for name,items in (('tracking_per_seed.csv',trackrows),('reproduction.csv',reproductions)):
                    with (out/name).open('w',newline='',encoding='utf-8-sig') as g:
                        w=csv.DictWriter(g,fieldnames=list(items[0]));w.writeheader();w.writerows(items)
                (out/'progress.json').write_text(json.dumps(dict(completed_cases=2*trial+1+int(hs==1.5),total_cases=2*args.trials,elapsed_s=time.perf_counter()-begin),indent=2),encoding='utf-8')
    if not all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in source_hash.items()):raise RuntimeError('Source changed during replay')
    (out/'manifest.json').write_text(json.dumps(dict(status='completed',seeds=list(range(args.seed_start,args.seed_start+args.trials)),rows=rows,methods=METHODS,
        source_sha256=source_hash,elapsed_s=time.perf_counter()-begin,
        scope='Legacy full raw-sensor models, common183.395ms horizon; CV linear forecast is added controlled ablation; native177ms used only for reproduction, no frequency-dependent trajectory changes'),indent=2),encoding='utf-8')

if __name__=='__main__':main()
