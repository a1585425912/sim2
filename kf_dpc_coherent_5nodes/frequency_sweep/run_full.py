"""Full paired frequency sensitivity study; communication f only is varied."""
import os,sys
from pathlib import Path
HERE=Path(__file__).resolve().parent;BASE=HERE.parent/'fmcw_doppler'
sys.path.insert(0,str(BASE))
import bootstrap
import time,json,csv,hashlib,argparse,gc
import numpy as np
from physics import Config,SeaTruth,sensor_record,rotate,spherical
from navigation import navigation_factors
from solver import Parameters
from target import wrap
from waveform import Waveform,buoy_bursts,observe
from joint_filter import fit,predict
import run_joint as original

# Log grid resolves LF/MF/HF extrapolation; integer MHz grid covers 1..300.
# Formula in THEORY defines all intermediate positive frequencies, not merely
# these numerical samples. 300MHz is retained as the reference endpoint.
FREQUENCIES=np.r_[np.logspace(-3,0,61),np.arange(2.,301.)]
METHODS=('GP_noD_q1','GP_Doppler_q1')
FIELDS=('seed','hs','trajectory','method','timing','phase','frequency_mhz',
        'eta_energy','eta_frame_mean','eta_frame_p05','gain_db','rx_ideal_model_w')


def evaluate(seed,cfg,kind,method,timing,t,error,distance,writer):
    amplitude=1/distance
    ideal=cfg.total_power/cfg.nodes*amplitude.sum(-1)**2
    incoherent=cfg.total_power/cfg.nodes*(amplitude**2).sum(-1)
    masks={'all':np.ones(len(t),bool),'late75_90':(t>=75)&(t<90)}
    for start in range(0,len(FREQUENCIES),12):
        fs=FREQUENCIES[start:start+12]
        phase=(2*np.pi/299792458.)*fs[:,None,None,None]*1e6*error[None]
        field=(amplitude[None]*np.exp(1j*phase)).sum(-1)
        rx=cfg.total_power/cfg.nodes*np.abs(field)**2
        frame_eta=rx.sum(-1)/ideal.sum(-1)[None]
        for phase_name,mask in masks.items():
            for i,f in enumerate(fs):
                energy=float(rx[i,mask].sum()/ideal[mask].sum())
                gain=float(10*np.log10(rx[i,mask].sum()/incoherent[mask].sum()))
                # Below30MHz only normalized phase sensitivity is interpreted;
                # extrapolated Friis absolute power is deliberately not output.
                power=float(rx[i,mask].mean()*(299792458./(4*np.pi*f*1e6))**2) if f>=30 else ''
                writer.writerow((seed,cfg.hs,kind,method,timing,phase_name,float(f),energy,
                    float(frame_eta[i,mask].mean()),float(np.quantile(frame_eta[i,mask],.05)),gain,power))


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--trials',type=int,default=20)
    ap.add_argument('--seed-start',type=int,default=20261028);ap.add_argument('--output',type=Path,default=HERE/'outputs');args=ap.parse_args()
    out=args.output;out.mkdir(parents=True,exist_ok=True);cache=out/'geometry';cache.mkdir(exist_ok=True)
    source_files=[Path(__file__),BASE/'bootstrap.py',BASE/'waveform.py',BASE/'joint_filter.py',
                  HERE.parent/'joint_gp/navigation.py',HERE.parent/'joint_gp/solver.py',HERE.parent/'joint_gp/run_joint.py',
                  HERE.parent/'raw_sensor_maritime/physics.py',HERE.parent/'raw_sensor_maritime/estimation.py',
                  HERE.parent/'raw_sensor_maritime/target.py',HERE.parent/'raw_sensor_maritime/run_full.py',
                  HERE.parent/'imm_dpc/run_experiment.py',HERE.parent/'gm_dpc_unified/experiment.py']
    sources={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
    old_path=BASE/'outputs/per_seed_all.csv'
    with old_path.open(encoding='utf-8-sig') as f:
        old={(int(r['seed']),float(r['hs']),r['trajectory'],r['method']):r for r in csv.DictReader(f)
             if r['phase']=='all' and r['timing']=='fixed_delay' and float(r['frequency_mhz'])==300 and r['method'] in METHODS}
    sources[str(old_path)]=hashlib.sha256(old_path.read_bytes()).hexdigest()
    begin=time.perf_counter();reproduction=[];row_count=0
    with (out/'per_seed.csv').open('w',newline='',encoding='utf-8-sig') as f:
        writer=csv.writer(f);writer.writerow(FIELDS)
        for trial in range(args.trials):
            seed=args.seed_start+trial
            for hs in (.5,1.5):
                cfg=Config(hs=hs);wave=Waveform();scene=SeaTruth(cfg,seed)
                packets,evaluation=sensor_record(cfg,scene,seed+100000);nav=navigation_factors(cfg,packets)
                t=packets['frame_time'];buoy_p,buoy_v=buoy_bursts(scene,t,wave)
                std=np.array([cfg.range_std,*([np.deg2rad(cfg.angle_std_deg)]*2)])
                noise=np.random.default_rng(seed+200000).normal(size=(len(t),cfg.nodes,3))*std
                for ik,kind in enumerate(original.reference.KINDS):
                    p,v=original.reference.trajectory(t,kind,cfg)
                    radar=spherical(rotate(evaluation['frame']['R'].swapaxes(-1,-2),p[:,None]-evaluation['frame']['p']))+noise
                    radar[...,1:]=wrap(radar[...,1:])
                    doppler=observe(cfg,t,kind,original.reference.trajectory,buoy_p,buoy_v,seed+400000+ik*10000,wave)
                    for method in METHODS:
                        history=fit(cfg,t,radar,nav,par=Parameters(window=2,jerk_psd=1.),packets=packets,
                                    doppler=doppler if method=='GP_Doppler_q1' else None)
                        previous=old[(seed,hs,kind,method)]
                        rms=lambda e:float(np.sqrt(np.mean(np.sum(e**2,axis=-1))))
                        reproduction.append(dict(seed=seed,hs=hs,trajectory=kind,method=method,
                            position_rmse_difference=rms(history['state'][:,:3]-p)-float(previous['pos_rmse']),
                            velocity_rmse_difference=rms(history['state'][:,3:6]-v)-float(previous['vel_rmse'])))
                        saved={'time_s':t,'state':history['state'],'covariance':history['covariance']}
                        for timing,compute in (('fixed_delay',cfg.processing_s),('recorded_mean_cpu',float(previous['compute_mean_ms'])/1000)):
                            # Freeze old recorded track-mean processing time for
                            # every frequency. This is not the old per-frame CPU
                            # schedule; fixed_delay reproduces the previous study.
                            h=wave.half_span_s+cfg.rounds*cfg.round_s+max(cfg.processing_s,compute)+np.linspace(0,cfg.slot_s,cfg.slot_samples)
                            future,_=original.reference.trajectory(t[:,None]+h,kind,cfg)
                            distance=np.linalg.norm(future[:,:,None]-scene.at(t[:,None]+h)['p'],axis=-1)
                            command=np.stack([predict(history,hj,cfg,covariance=False)[2] for hj in h],axis=1)
                            error=command-distance
                            saved[timing+'__error']=error;saved[timing+'__distance']=distance;saved[timing+'__horizons']=h
                            evaluate(seed,cfg,kind,method,timing,t,error,distance,writer)
                            row_count+=len(FREQUENCIES)*2
                        np.savez_compressed(cache/f'{seed}_hs{hs:g}_{kind}_{method}.npz',**saved)
                        print(f'DONE seed={seed} Hs={hs:g} {kind} {method} elapsed={time.perf_counter()-begin:.1f}s',flush=True)
                        del history,saved;gc.collect()
                    f.flush()
                progress=dict(completed_cases=2*trial+1+int(hs==1.5),total_cases=2*args.trials,elapsed_s=time.perf_counter()-begin)
                (out/'progress.json').write_text(json.dumps(progress,indent=2),encoding='utf-8')
    if not all(hashlib.sha256(Path(p).read_bytes()).hexdigest()==s for p,s in sources.items()):raise RuntimeError('Source changed during execution')
    with (out/'reproduction.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(reproduction[0]));w.writeheader();w.writerows(reproduction)
    manifest=dict(status='completed',seeds=list(range(args.seed_start,args.seed_start+args.trials)),rows=row_count,
        frequencies_mhz=FREQUENCIES.tolist(),source_sha256=sources,elapsed_s=time.perf_counter()-begin,
        timing_note='recorded_mean_cpu uses original completed study mean per track, constant across frames/frequencies, not original varying per-frame latency',
        scope='Communication frequency only; ideal phase/channel sensitivity, not hardware/sea propagation feasibility or spectrum authorization. Below30MHz no absolute received power reported.')
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')

if __name__=='__main__':main()
