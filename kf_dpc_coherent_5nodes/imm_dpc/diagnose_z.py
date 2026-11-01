"""Diagnostic ablations only; baseline estimator/files remain unchanged."""
from pathlib import Path
import sys, csv, json
import numpy as np
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import run_experiment as run
import imm_filter as imm
base=run.base


def main():
    rows=[]; original=imm.process_covariance
    for seed in range(20261028,20261031):
        cfg=base.Config(seed=seed); t=np.arange(451)*cfg.dt; mk=t>=5
        ss=np.random.SeedSequence(seed).spawn(3)
        b,v=base.gm_buoy_truth(cfg,np.random.default_rng(ss[0]),90.2)
        idx=np.arange(len(t))*200; b,v=b[idx],v[idx]
        std=np.array([.8,np.deg2rad(.35),np.deg2rad(.35)])
        noise=np.random.default_rng(ss[2]).standard_normal((len(t),20,3))*std
        nav=base.navigation(cfg,b,v,np.random.default_rng(ss[1]),'nominal')
        ideal=base.navigation(cfg,b,v,np.random.default_rng(ss[1]),'ideal')
        for kind in ('hover','straight','circle','mixed','climb'):
            if kind=='climb':
                p,vt=base.uav_trajectory(t,'straight')
                positive=lambda a:np.maximum(a,0)
                p[:,2]+=.5*(positive(t-20)**2-positive(t-22)**2-positive(t-40)**2+positive(t-42)**2)
                vt[:,2]=np.clip(t-20,0,2)-np.clip(t-40,0,2)
            else:p,vt=run.trajectory(t,kind,cfg)
            sph=base.spherical(p[:,None]-b)+noise
            obs,R=base.position_observations(sph,nav,np.diag(std**2))
            diag=np.zeros_like(R)
            diag[...,range(3),range(3)]=np.diagonal(R,axis1=-2,axis2=-1)
            for variant in ('original','qz_x0.01','diagonal_R_diagnostic','ideal_navigation','marginal_Z_KF_q0.001'):
                imm.process_covariance=original
                zz,rr=obs,R
                if variant=='diagonal_R_diagnostic':rr=diag
                if variant=='ideal_navigation':zz,rr=base.position_observations(sph,ideal,np.diag(std**2))
                if variant=='qz_x0.01':
                    def low(h,model,par):
                        Q=original(h,model,par)
                        for a,c in ((2,2),(2,5),(5,2),(5,5)):Q[a,c]*=.01
                        return Q
                    imm.process_covariance=low
                try:
                    if variant=='marginal_Z_KF_q0.001':
                        local=base.local_kf(t,obs,diag,.001,20.)
                    else:local=imm.fit(t,zz,rr)['state'][...,:6]
                finally:imm.process_covariance=original
                state=run.dpc(local,run.ring_weights(20),25)
                ep=state[...,:3]-p[:,None]; ev=state[...,3:]-vt[:,None]
                row=dict(seed=seed,trajectory=kind,variant=variant,
                         z_bias=float(ep[mk,:,2].mean()),z_rmse=float(np.sqrt(np.mean(ep[mk,:,2]**2))),
                         vz_rmse=float(np.sqrt(np.mean(ev[mk,:,2]**2))),
                         xy_rmse=float(np.sqrt(np.mean(np.sum(ep[mk,:,:2]**2,-1)))),
                         raw_z_rmse=float(np.sqrt(np.mean((obs[mk,:,2]-p[mk,None,2])**2))),
                         median_z_sigma=float(np.median(np.sqrt(R[mk,:,2,2]))),
                         mean_abs_xz_correlation=float(np.mean(np.abs(R[mk,:,0,2])/np.sqrt(R[mk,:,0,0]*R[mk,:,2,2]))))
                rows.append(row)
            print(seed,kind,'done',flush=True)
    out=HERE/'outputs/z_diagnostics';out.mkdir(exist_ok=True)
    with (out/'ablations.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    for kind in ('hover','straight','circle','mixed','climb'):
        for variant in ('original','qz_x0.01','diagonal_R_diagnostic','ideal_navigation','marginal_Z_KF_q0.001'):
            rs=[r for r in rows if r['trajectory']==kind and r['variant']==variant]
            print(kind,variant,'Z_RMSE',round(np.sqrt(np.mean([r['z_rmse']**2 for r in rs])),4),
                  'bias',round(np.mean([r['z_bias'] for r in rs]),4),flush=True)

if __name__=='__main__':main()
