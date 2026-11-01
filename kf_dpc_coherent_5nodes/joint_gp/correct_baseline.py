"""Full-seed correction of runtime baseline evaluation, no GP refit/tuning.

The recorded CPU delays are retained. The original code mistakenly evaluated
unfused local IMM forecasts for this one timing variant. Apply the same W^25
CI information fusion as the fixed-delay baseline to all future states.

IMPORTANT: CI fusion of means requires covariances, not arithmetic DPC. This
script forecasts each branch covariance at each recorded scalar horizon and
uses the original CI rule, exactly matching the fixed-delay path.
"""
import os
for name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ.setdefault(name,'1')
import argparse,json,csv,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;RAW=HERE.parent/'raw_sensor_maritime';ROOT=HERE.parents[1]
for p in (ROOT,HERE.parent/'imm_ci',HERE.parent/'imm_dpc',RAW,HERE):sys.path.insert(0,str(p))
import numpy as np
from physics import Config,SeaTruth,sensor_record,rotate,spherical
from navigation import navigation_factors
from target import fit_radar,wrap
from imm_filter import Parameters
from fusion import forecast,fuse,kernel
import run_joint as shared


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--seed-start',type=int,required=True);ap.add_argument('--trials',type=int,default=5);args=ap.parse_args()
    out=HERE/'outputs'
    with (out/'per_seed_all.csv').open(encoding='utf-8-sig') as f:original=list(csv.DictReader(f))
    source=[r for r in original if args.seed_start<=int(r['seed'])<args.seed_start+args.trials and r['method']=='IMM_CI_linear' and r['timing']=='measured_cpu']
    updated=[]
    for seed in range(args.seed_start,args.seed_start+args.trials):
        for hs in (.5,1.5):
            cfg=Config(hs=hs);scene=SeaTruth(cfg,seed);packets,evaluation=sensor_record(cfg,scene,seed+100000);nav=navigation_factors(cfg,packets)
            t=packets['frame_time'];C=kernel(shared.reference.previous.ring_weights(cfg.nodes),cfg.rounds)
            std=np.array([cfg.range_std,*([np.deg2rad(cfg.angle_std_deg)]*2)]);noise=np.random.default_rng(seed+200000).normal(size=(len(t),cfg.nodes,3))*std
            masks=dict(all=np.ones(len(t),bool),startup=t<5,after_startup=t>=5,acceleration=(t>=15)&(t<20),turn_entry=(t>=35)&(t<40),turn_steady=(t>=40)&(t<55),braking=(t>=65)&(t<70),hover_end=(t>=75)&(t<90))
            for kind in shared.KINDS:
                rows=[r for r in source if int(r['seed'])==seed and float(r['hs'])==hs and r['trajectory']==kind]
                if len(rows)!=16:raise ValueError('Incomplete original baseline records')
                seconds=float(rows[0]['compute_mean_ms'])/1000
                h=cfg.rounds*cfg.round_s+max(cfg.processing_s,seconds)+np.linspace(0,cfg.slot_s,cfg.slot_samples)
                p,_=shared.reference.trajectory(t,kind,cfg);pfuture,_=shared.reference.trajectory(t[:,None]+h,kind,cfg)
                radar=spherical(rotate(evaluation['frame']['R'].swapaxes(-1,-2),p[:,None]-evaluation['frame']['p']))+noise;radar[...,1:]=wrap(radar[...,1:])
                history=fit_radar(t,radar,nav,np.diag(std**2),par=Parameters())
                pairs=[forecast(history,float(hh)) for hh in h]
                means=np.stack([a[0] for a in pairs],axis=1);covs=np.stack([a[1] for a in pairs],axis=1)
                states,_=fuse(means,covs,C,'ci')
                predicted=states[...,:3];buoy=nav['p'][:,None]+h[None,:,None,None]*nav['v'][:,None]
                cmd=np.linalg.norm(predicted-buoy,axis=-1)
                true_d=np.linalg.norm(pfuture[:,:,None]-scene.at(t[:,None]+h)['p'],axis=-1)
                error=predicted-pfuture[:,:,None];path=cmd-true_d;path-=path.mean(-1,keepdims=True)
                for frequency in cfg.frequencies_mhz:
                    power=shared.power(cmd,true_d,frequency,cfg.total_power)
                    for old in rows:
                        if float(old['frequency_mhz'])!=frequency:continue
                        mask=masks[old['phase']];r=old.copy();gain=float(power['rx'][mask].sum()/power['inc'][mask].sum())
                        r.update(eta=float(power['rx'][mask].sum()/power['ideal'][mask].sum()),eta_own=float(power['rx'][mask].sum()/power['own'][mask].sum()),eta_optimal_bound=float(power['rx'][mask].sum()/power['optimal'][mask].sum()),gain_linear=gain,gain_db=float(10*np.log10(gain)),rx_mean_w=float(power['rx'][mask].mean()),future_pos_rmse=float(np.sqrt(np.mean(np.sum(error[mask]**2,axis=-1)))),path_rmse=float(np.sqrt(np.mean(path[mask]**2))))
                        updated.append(r)
            print(f'Corrected all five runtime-baseline tracks seed={seed} Hs={hs:g}',flush=True)
    path=out/'corrections'/f'rows_{args.seed_start}.csv';shared.reference.csv_write(path,updated)
    print(f'Saved {len(updated)} corrected records to {path}',flush=True)


if __name__=='__main__':main()
