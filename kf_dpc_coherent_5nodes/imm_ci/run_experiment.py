"""IMM local posterior shared by Local/DPC/CI; only fusion rule differs."""
import sys, argparse, json, time, hashlib
from pathlib import Path
from dataclasses import asdict, replace
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'imm_dpc'))
import run_experiment as imm_run
from imm_filter import fit, Parameters
from fusion import forecast, fuse, kernel
base=imm_run.base
KINDS=('hover','straight','circle','mixed','climb')
METHODS=('local','dpc','ci')


def trajectory(t, kind, cfg):
    if kind!='climb':return imm_run.trajectory(t,kind,cfg)
    p,v=base.uav_trajectory(t,'straight')
    positive=lambda a:np.maximum(a,0)
    p[...,2]+=.5*(positive(t-20)**2-positive(t-22)**2-positive(t-40)**2+positive(t-42)**2)
    v[...,2]=np.clip(t-20,0,2)-np.clip(t-40,0,2)
    return p,v


def run_seed(cfg, par, keep=False):
    start=time.perf_counter();t=np.arange(round(cfg.duration/cfg.dt)+1)*cfg.dt
    nt=len(t);L=cfg.main_rounds
    horizons=cfg.processing_s+.005*L+np.linspace(0,cfg.slot_s,cfg.slot_samples)
    if horizons[-1]>cfg.frame_budget_s:raise ValueError('Frame budget exceeded')
    ss=np.random.SeedSequence(cfg.seed).spawn(3)
    b,v=base.gm_buoy_truth(cfg,np.random.default_rng(ss[0]),cfg.duration+cfg.frame_budget_s)
    frame_idx=np.arange(nt)*round(cfg.dt/cfg.truth_grid_s)
    nav=base.navigation(cfg,b[frame_idx],v[frame_idx],np.random.default_rng(ss[1]),'nominal')
    std=np.array([cfg.range_std_m,*([np.deg2rad(cfg.angle_std_deg)]*2)])
    noise=np.random.default_rng(ss[2]).standard_normal((nt,cfg.nodes,3))*std
    C=kernel(imm_run.ring_weights(cfg.nodes),L);levers=base.levers_for(cfg)
    idx=frame_idx[:,None]+np.rint(horizons[None]/cfg.truth_grid_s).astype(int)
    elem=b[idx][...,None,:]+levers
    phases=dict(all=np.ones(nt,bool),after_startup=t>=5,acceleration=(t>=15)&(t<20),
                turn_entry=(t>=35)&(t<40),turn_steady=(t>=40)&(t<55),
                braking=(t>=65)&(t<70),hover_end=(t>=75)&(t<90))
    rows=[];saved={'time_s':t}
    for kind in KINDS:
        p,vtrue=trajectory(t,kind,cfg)
        pfuture,_=trajectory(t[:,None]+horizons,kind,cfg)
        true_d=np.linalg.norm(pfuture[:,:,None,None,:]-elem,axis=-1)
        sph=base.spherical(p[:,None]-b[frame_idx])+noise
        sph[...,1]=base.wrap_angle(sph[...,1])
        obs,R=base.position_observations(sph,nav,np.diag(std**2))
        hist=fit(t,obs,R,par,retain_branch_covariances=True)
        pairs=[forecast(hist,h,par) for h in (0.,*horizons)]
        means=np.stack([pair[0] for pair in pairs],axis=1)
        covs=np.stack([pair[1] for pair in pairs],axis=1)
        del pairs
        for method in METHODS:
            tic=time.perf_counter()
            states,Ps=fuse(means,covs,C,method)
            fusion_s=time.perf_counter()-tic
            current=states[:,0];currentP=Ps[:,0]
            predictions=[states[:,j] for j in range(1,len(horizons)+1)]
            cmd=imm_run.exact_commands(predictions,nav,levers,horizons,cfg)
            ep=current[...,:3]-p[:,None];ev=current[...,3:]-vtrue[:,None]
            err=np.concatenate([ep,ev],axis=-1)
            nees=np.einsum('tni,tni->tn',err,np.linalg.solve(currentP,err[...,None])[...,0])
            future=states[:,1:,:,:3]-pfuture[:,:,None]
            path=cmd-true_d;path-=path.mean(axis=(-1,-2),keepdims=True)
            if keep:
                saved[f'{kind}__{method}__state']=current.copy()
                saved[f'{kind}__{method}__sigma_z']=np.sqrt(currentP[...,2,2])
                if method=='dpc' and kind in imm_run.KINDS:
                    ref=HERE.parent/'imm_dpc/outputs/first_seed.npz'
                    with np.load(ref) as old:
                        np.testing.assert_allclose(current,old[f'{kind}__IMM__state'],rtol=1e-11,atol=1e-10)
            for frequency in cfg.frequencies_mhz:
                pw=base.powers(cmd,true_d,frequency,cfg.total_power_w)
                eta=pw['rx'].sum(1)/pw['ideal'].sum(1)
                if keep:
                    saved[f'{kind}__{method}__eta{int(frequency)}']=eta
                    if method=='dpc' and kind in imm_run.KINDS:
                        with np.load(HERE.parent/'imm_dpc/outputs/first_seed.npz') as old:
                            np.testing.assert_allclose(eta,old[f'{kind}__IMM__eta{int(frequency)}'],rtol=1e-10,atol=1e-10)
                for phase,mk in phases.items():
                    row=dict(seed=cfg.seed,trajectory=kind,method=method,frequency_mhz=frequency,phase=phase,
                             eta=float(pw['rx'][mk].sum()/pw['ideal'][mk].sum()),
                             rx_sum=float(pw['rx'][mk].sum()),ideal_sum=float(pw['ideal'][mk].sum()),
                             pos_rmse=float(np.sqrt(np.mean(np.sum(ep[mk]**2,axis=-1)))),
                             vel_rmse=float(np.sqrt(np.mean(np.sum(ev[mk]**2,axis=-1)))),
                             future_pos_rmse=float(np.sqrt(np.mean(np.sum(future[mk]**2,axis=-1)))),
                             path_rmse=float(np.sqrt(np.mean(path[mk]**2))),
                             z_bias=float(ep[mk,:,2].mean()),nees6=float(nees[mk].mean()),
                             fusion_seconds=fusion_s,
                             disagree_pos=float(np.sqrt(np.mean(np.sum((current[mk,:,:3]-current[mk,:,:3].mean(1,keepdims=True))**2,axis=-1)))))
                    for j,axis in enumerate('XYZ'):
                        row[f'p{axis}_rmse']=float(np.sqrt(np.mean(ep[mk,:,j]**2)))
                        row[f'v{axis}_rmse']=float(np.sqrt(np.mean(ev[mk,:,j]**2)))
                    rows.append(row)
            del states,Ps
        if keep:
            saved[f'{kind}__truth_p']=p;saved[f'{kind}__truth_v']=vtrue
        print(f'  {kind} done ({time.perf_counter()-start:.1f}s)',flush=True)
    return rows,saved


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=20)
    ap.add_argument('--output',type=Path,default=HERE/'outputs');args=ap.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    cfg=base.Config();par=Parameters();rows=[];start=time.perf_counter()
    for i in range(args.trials):
        r,s=run_seed(replace(cfg,seed=cfg.seed+i),par,keep=i==0);rows+=r
        if i==0:np.savez_compressed(args.output/'first_seed.npz',**s)
        base.write_csv(args.output/'per_seed.csv',rows)
        print(f'SEED {i+1}/{args.trials} saved, {time.perf_counter()-start:.1f}s',flush=True)
    files=[HERE/'fusion.py',HERE/'run_experiment.py',HERE.parent/'imm_dpc/imm_filter.py',HERE.parent/'imm_dpc/run_experiment.py',HERE.parent/'gm_dpc_unified/experiment.py']
    metadata=dict(trials=args.trials,seeds=[cfg.seed,cfg.seed+args.trials-1],config=asdict(cfg),imm=asdict(par),
                  elapsed_s=time.perf_counter()-start,python=sys.version,numpy=np.__version__,
                  source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files},
                  definition='CI: finite 25-round convex information consensus, fixed Metropolis weights; no logdet optimization; no N multiplier',
                  feedback=False,prediction='frozen-mode branch mean/covariance, then CI each slot',
                  timing='Identical assumed 12ms processing + 25*5ms communication; covariance payload and hardware processing unvalidated',
                  baseline_check='First-seed DPC states and both-frequency power curves match prior IMM outputs with atol=1e-10')
    (args.output/'config.json').write_text(json.dumps(metadata,indent=2),encoding='utf-8')

if __name__=='__main__':main()
