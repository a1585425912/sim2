"""Full paired raw-sensor experiment: 20 seeds x two seas x five UAV paths."""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(name,'1')
import sys,json,time,hashlib,csv,argparse
from pathlib import Path
from dataclasses import asdict,replace
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT));sys.path.insert(0,str(HERE.parent/'imm_ci'));sys.path.insert(0,str(HERE.parent/'imm_dpc'))
import numpy as np
from scipy.spatial.transform import Rotation
from scipy.stats import chi2
import run_experiment as previous
from fusion import forecast,fuse,kernel
from imm_filter import Parameters
from physics import Config,SeaTruth,sensor_record,rotate,spherical
from estimation import navigate,WavePredictor
from target import fit_radar,wrap

KINDS=('hover','straight','circle','mixed','climb')
METHODS=('CV_CI_wave','IMM_DPC_linear','IMM_CI_linear','IMM_CI_wave')


def trajectory(t,kind,cfg):
    if kind!='climb':return previous.trajectory(t,kind,cfg)
    p,v=previous.base.uav_trajectory(t,'straight')
    positive=lambda a:np.maximum(a,0.)
    p[...,2]+=.5*(positive(t-20)**2-positive(t-22)**2-positive(t-40)**2+positive(t-42)**2)
    v[...,2]=np.clip(t-20,0,2)-np.clip(t-40,0,2)
    return p,v


def csv_write(path,rows):
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def metrics(e,P):
    nees=np.sum(e*np.linalg.solve(P,e[...,None])[...,0],axis=-1)
    return float(nees.mean()),float(np.mean(nees<=chi2.ppf(.95,e.shape[-1])))


def run_case(cfg,seed,out,keep):
    tic=time.perf_counter();scene=SeaTruth(cfg,seed)
    packets,eval_nav=sensor_record(cfg,scene,seed+100000)
    nav=navigate(cfg,packets);t=packets['frame_time'];h=cfg.horizons()
    wave_model=WavePredictor(cfg);wave_history=wave_model.fit(nav)
    buoy_linear=nav['p'][:,None]+h[None,:,None,None]*nav['v'][:,None]
    buoy_wave=np.stack([wave_model.predict(wave_history,float(a)) for a in h],axis=1)
    true_buoy=scene.at(t[:,None]+h)['p']
    truth_nav=eval_nav['frame'];ep=nav['p']-truth_nav['p'];ev=nav['v']-truth_nav['v']
    attitude=Rotation.from_matrix((truth_nav['R'].swapaxes(-1,-2)@nav['R']).reshape(-1,3,3)).as_rotvec().reshape(ep.shape)
    nav_error=np.concatenate([ep,ev],-1)
    nav_nees,nav_coverage=metrics(nav_error,nav['P'][...,:6,:6])
    rmse=lambda error:float(np.sqrt(np.mean(np.sum(error**2,axis=-1))))
    navrow=dict(seed=seed,hs=cfg.hs,pos_rmse=rmse(ep),vel_rmse=rmse(ev),
                attitude_rmse_deg=float(np.rad2deg(rmse(attitude))),nees6=nav_nees,coverage95_6=nav_coverage,
                linear_future_rmse=rmse(buoy_linear-true_buoy),wave_future_rmse=rmse(buoy_wave-true_buoy),
                accel_bias_rmse=rmse(nav['ba']-eval_nav['bias_accel']),gyro_bias_rmse=rmse(nav['bg']-eval_nav['bias_gyro']))
    if keep:
        np.savez_compressed(out/f'sensors_hs{cfg.hs:g}.npz',**packets)
        np.savez_compressed(out/f'navigation_hs{cfg.hs:g}.npz',**nav,wave_history=wave_history,buoy_linear=buoy_linear,buoy_wave=buoy_wave)
        np.savez_compressed(out/f'truth_navigation_eval_only_hs{cfg.hs:g}.npz',time=t,**truth_nav,future_p=true_buoy,**{k:v for k,v in eval_nav.items() if k!='frame'})
    print(f'NAV seed={seed} Hs={cfg.hs:g} p={navrow["pos_rmse"]:.3f}m v={navrow["vel_rmse"]:.3f}m/s attitude={navrow["attitude_rmse_deg"]:.3f}deg ({time.perf_counter()-tic:.1f}s)',flush=True)
    C=kernel(previous.ring_weights(cfg.nodes),cfg.rounds)
    radar_std=np.array([cfg.range_std,*([np.deg2rad(cfg.angle_std_deg)]*2)])
    rng=np.random.default_rng(seed+200000);noise=rng.normal(size=(len(t),cfg.nodes,3))*radar_std
    phases=dict(all=np.ones(len(t),bool),startup=t<5,after_startup=t>=5,acceleration=(t>=15)&(t<20),turn_entry=(t>=35)&(t<40),turn_steady=(t>=40)&(t<55),braking=(t>=65)&(t<70),hover_end=(t>=75)&(t<90))
    rows=[];saved=dict(time_s=t);par=Parameters()
    for kind in KINDS:
        p,v=trajectory(t,kind,cfg);future,_=trajectory(t[:,None]+h,kind,cfg)
        radar=spherical(rotate(truth_nav['R'].swapaxes(-1,-2),p[:,None]-truth_nav['p']))+noise
        radar[...,1:]=wrap(radar[...,1:])
        true_d=np.linalg.norm(future[:,:,None]-true_buoy,axis=-1)
        if keep:
            np.savez_compressed(out/f'radar_{kind}_hs{cfg.hs:g}.npz',time=t,measurement=radar,covariance=np.diag(radar_std**2))
            saved[f'{kind}__truth_p']=p;saved[f'{kind}__truth_v']=v
        for target_method,enabled in (('CV',(1,)),('IMM',(0,1,2))):
            hist=fit_radar(t,radar,nav,np.diag(radar_std**2),enabled,par)
            forecasts=[forecast(hist,float(a),par) for a in (0.,*h)]
            means=np.stack([f[0] for f in forecasts],1);covs=np.stack([f[1] for f in forecasts],1)
            selected=('CV_CI_wave',) if target_method=='CV' else METHODS[1:]
            fused={}
            for method in selected:
                fusion='dpc' if 'DPC' in method else 'ci'
                if fusion not in fused:fused[fusion]=fuse(means,covs,C,fusion)
                x,P=fused[fusion];current=x[:,0];currentP=P[:,0]
                buoy=buoy_wave if method.endswith('wave') else buoy_linear
                cmd=np.linalg.norm(x[:,1:,:,:3]-buoy,axis=-1)[...,None]
                e=current-np.concatenate([p,v],-1)[:,None]
                future_error=x[:,1:,:,:3]-future[:,:,None]
                path=cmd[...,0]-true_d;path-=path.mean(-1,keepdims=True)
                if keep:
                    saved[f'{kind}__{method}__state']=current;saved[f'{kind}__{method}__sigma']=np.sqrt(np.diagonal(currentP,axis1=-2,axis2=-1))
                    saved[f'{kind}__{method}__probability']=hist['probabilities']
                for frequency in cfg.frequencies_mhz:
                    pw=previous.base.powers(cmd,true_d[...,None],frequency,cfg.total_power)
                    eta=pw['rx'].sum(1)/pw['ideal'].sum(1)
                    if keep:saved[f'{kind}__{method}__eta{int(frequency)}']=eta
                    for phase,mask in phases.items():
                        nees,coverage=metrics(e[mask],currentP[mask])
                        neesp,coveragep=metrics(e[mask,:,:3],currentP[mask,:,:3,:3])
                        row=dict(seed=seed,hs=cfg.hs,trajectory=kind,method=method,frequency_mhz=frequency,phase=phase,eta=float(pw['rx'][mask].sum()/pw['ideal'][mask].sum()),rx_mean_w=float(pw['rx'][mask].mean()),ideal_mean_w=float(pw['ideal'][mask].mean()),pos_rmse=float(np.sqrt(np.mean(np.sum(e[mask,:,:3]**2,-1)))),vel_rmse=float(np.sqrt(np.mean(np.sum(e[mask,:,3:]**2,-1)))),future_pos_rmse=float(np.sqrt(np.mean(np.sum(future_error[mask]**2,-1)))),path_rmse=float(np.sqrt(np.mean(path[mask]**2))),nees6=nees,coverage95_6=coverage,nees_position=neesp,coverage95_position=coveragep)
                        for axis,name in enumerate('XYZ'):
                            row[f'p{name}_rmse']=float(np.sqrt(np.mean(e[mask,:,axis]**2)));row[f'v{name}_rmse']=float(np.sqrt(np.mean(e[mask,:,axis+3]**2)))
                        rows.append(row)
        print(f'  TARGET {kind} seed={seed} Hs={cfg.hs:g} done ({time.perf_counter()-tic:.1f}s)',flush=True)
    if keep:np.savez_compressed(out/f'first_seed_hs{cfg.hs:g}.npz',**saved)
    return rows,navrow


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--trials',type=int,default=20);parser.add_argument('--output',type=Path,default=HERE/'outputs')
    args=parser.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=True)
    configs=[Config(hs=hs) for hs in (.5,1.5)];rows=[];navrows=[];start=time.perf_counter()
    for i in range(args.trials):
        seed=20261028+i
        for cfg in configs:
            r,n=run_case(cfg,seed,out,i==0);rows.extend(r);navrows.append(n)
            csv_write(out/'per_seed.csv',rows);csv_write(out/'navigation_per_seed.csv',navrows)
            (out/'progress.json').write_text(json.dumps(dict(completed_cases=len(navrows),total_cases=args.trials*len(configs),elapsed_s=time.perf_counter()-start,last_seed=seed,last_hs=cfg.hs),indent=2),encoding='utf-8')
        print(f'SEED {i+1}/{args.trials} complete; total {time.perf_counter()-start:.1f}s',flush=True)
    sources=list(HERE.glob('*.py'))+[HERE.parent/'imm_dpc/imm_filter.py',HERE.parent/'imm_ci/fusion.py',HERE.parent/'imm_dpc/run_experiment.py',HERE.parent/'gm_dpc_unified/experiment.py']
    manifest=dict(status='completed',trials=args.trials,seeds=list(range(20261028,20261028+args.trials)),configs=[asdict(c) for c in configs],elapsed_s=time.perf_counter()-start,python=sys.version,numpy=np.__version__,source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},definition='Raw pseudorange+IMU+mag navigation ESKF; raw radar IMM; fixed 25-round DPC/CI; zero levers; no truth initialization, no supplied slow/fast truth; no gating; fixed network total TX power',limitations=['Hierarchical navigation/radar fusion, not full joint MAP; temporal pose/target correlation approximated','Secondary wave filter reuses correlated navigation outputs; its covariance is not used as calibrated joint confidence','Finite spectral surrogate RAO, not measured hydrodynamics','GNSS corrected code 0.3m noise with drifting receiver clock; no carrier ambiguity/multipath/ionosphere model','Ideal RF clock/carrier synchronization; co-located sensor and single isotropic TX per buoy','Identical assumed processing and network delay, not real-time hardware timing'])
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')


if __name__=='__main__':main()
