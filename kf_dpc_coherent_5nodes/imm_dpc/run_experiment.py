"""Paired CV/CT/IMM coherent-transmission experiment; baseline files untouched."""
import sys, time, json, argparse, csv, hashlib
from pathlib import Path
from dataclasses import replace, asdict
import numpy as np
from scipy.integrate import cumulative_trapezoid
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT))
from kf_dpc_coherent_5nodes.gm_dpc_unified import experiment as base
from kf_dpc_coherent_5nodes.pipeline import dpc, ring_weights
from imm_filter import fit, predict, Parameters

KINDS = ('hover', 'straight', 'circle', 'mixed')
METHODS = ('CV_original', 'CV_exact', 'CT', 'IMM_SC', 'IMM')


def trajectory(t, kind, cfg):
    if kind != 'mixed':
        return base.uav_trajectory(t, kind)
    # Continuous position and velocity. Hidden switch times only define truth.
    grid = np.arange(round((cfg.duration+cfg.frame_budget_s)/.001)+1)*.001
    speed = 3*np.clip(grid-15, 0, 5)-3*np.clip(grid-65, 0, 5)
    heading = .1*np.clip(grid-35, 0, 20)
    v = np.c_[speed*np.cos(heading), speed*np.sin(heading), np.zeros(len(grid))]
    p = np.array([0., -200., 120.])+cumulative_trapezoid(v, grid, axis=0, initial=0)
    shape = np.asarray(t).shape+(3,)
    return tuple(np.stack([np.interp(np.asarray(t).ravel(), grid, a[:, j]) for j in range(3)], -1).reshape(shape) for a in (p, v))


def exact_commands(predictions, nav, levers, horizons, cfg):
    out = []
    for h, pu in zip(horizons, predictions):
        b, _ = base.buoy_prediction(nav, h, cfg, 'gm')
        elem = b[:, :, None]+np.einsum('tnij,mj->tnmi', nav['C'], levers)
        out.append(np.linalg.norm(pu[..., None, :3]-elem, axis=-1))
    return np.stack(out, axis=1)


def run_seed(cfg, par, keep=False):
    start = time.perf_counter()
    t = np.arange(round(cfg.duration/cfg.dt)+1)*cfg.dt
    nt = len(t); L = cfg.main_rounds
    tau = cfg.processing_s+.005*L
    slot = np.linspace(0, cfg.slot_s, cfg.slot_samples); horizons = tau+slot
    if horizons[-1] > cfg.dt: raise ValueError('Frame budget exceeded')
    ss = np.random.SeedSequence(cfg.seed).spawn(3)
    b, v = base.gm_buoy_truth(cfg, np.random.default_rng(ss[0]), cfg.duration+cfg.frame_budget_s)
    idx = np.arange(nt)*round(cfg.dt/cfg.truth_grid_s)
    nav = base.navigation(cfg, b[idx], v[idx], np.random.default_rng(ss[1]), 'nominal')
    std = np.array([cfg.range_std_m, *([np.deg2rad(cfg.angle_std_deg)]*2)])
    noise = np.random.default_rng(ss[2]).standard_normal((nt, cfg.nodes, 3))*std
    W = ring_weights(cfg.nodes); levers = base.levers_for(cfg)
    future_idx = idx[:, None]+np.rint(horizons[None]/cfg.truth_grid_s).astype(int)
    elem = b[future_idx][..., None, :]+levers
    phases = {'all': np.ones(nt, bool), 'after_startup': t >= cfg.warmup}
    if cfg.duration >= 80:
        phases.update(hover_start=(t>=5)&(t<15), acceleration=(t>=15)&(t<20),
                      straight_segment=(t>=25)&(t<35), turn_entry=(t>=35)&(t<40),
                      turn_steady=(t>=40)&(t<55), braking=(t>=65)&(t<70),
                      hover_end=(t>=75)&(t<90))
    rows, saved = [], {'time_s': t}
    for kind in KINDS:
        pt, vt = trajectory(t, kind, cfg)
        ptfuture, _ = trajectory(t[:, None]+horizons, kind, cfg)
        true_d = np.linalg.norm(ptfuture[:, :, None, None, :]-elem, axis=-1)
        sph = base.spherical(pt[:, None]-b[idx])+noise
        sph[..., 1] = base.wrap_angle(sph[..., 1])
        obs, R = base.position_observations(sph, nav, np.diag(std**2))
        tic = time.perf_counter()
        cv = base.local_kf(t, obs, R, cfg.acceleration_psd, cfg.initial_velocity_std)
        cv_seconds = time.perf_counter()-tic
        histories = {}; filter_seconds = {}
        for method, enabled in (('CT', (2,)), ('IMM_SC', (0,1)), ('IMM', (0,1,2))):
            tic = time.perf_counter()
            histories[method] = fit(t, obs, R, par, enabled)
            filter_seconds[method] = time.perf_counter()-tic
        for method in METHODS:
            if method.startswith('CV'):
                current = dpc(cv, W, L)
                preds = [current.copy() for h in horizons]
                for h, pred in zip(horizons, preds): pred[..., :3] += h*pred[..., 3:]
                elapsed = cv_seconds
            else:
                hist = histories[method]
                current = dpc(hist['state'][..., :6], W, L)
                preds = [dpc(predict(hist, h, par), W, L) for h in horizons]
                elapsed = filter_seconds[method]
            cmd = (base.commands(current, nav, levers, tau, slot, cfg)
                   if method == 'CV_original' else exact_commands(preds, nav, levers, horizons, cfg))
            ep = current[..., :3]-pt[:, None]; ev = current[..., 3:]-vt[:, None]
            ef = np.stack([p[..., :3] for p in preds], 1)-ptfuture[:, :, None]
            path = cmd-true_d; path -= path.mean(axis=(-1,-2), keepdims=True)
            for f in cfg.frequencies_mhz:
                pw = base.powers(cmd, true_d, f, cfg.total_power_w)
                if keep:
                    saved[f'{kind}__{method}__eta{int(f)}'] = pw['rx'].sum(1)/pw['ideal'].sum(1)
                for phase, mask in phases.items():
                    if not mask.any(): continue
                    row = dict(seed=cfg.seed, trajectory=kind, method=method, phase=phase,
                               frequency_mhz=f, eta=float(pw['rx'][mask].sum()/pw['ideal'][mask].sum()),
                               rx_sum=float(pw['rx'][mask].sum()), ideal_sum=float(pw['ideal'][mask].sum()),
                               gain_db=float(10*np.log10(pw['rx'][mask].sum()/pw['elem'][mask].sum())),
                               pos_rmse=float(np.sqrt(np.mean(np.sum(ep[mask]**2, -1)))),
                               vel_rmse=float(np.sqrt(np.mean(np.sum(ev[mask]**2, -1)))),
                               future_pos_rmse=float(np.sqrt(np.mean(np.sum(ef[mask]**2, -1)))),
                               path_rmse=float(np.sqrt(np.mean(path[mask]**2))),
                               filter_seconds=elapsed, frames=nt)
                    if method == 'IMM':
                        means = histories[method]['probabilities'][mask].mean(axis=(0,1))
                        row.update(dict(zip(('prob_S','prob_CV','prob_CT'),map(float,means))))
                    else: row.update(prob_S='',prob_CV='',prob_CT='')
                    rows.append(row)
            if keep:
                saved[f'{kind}__{method}__state'] = current
        if keep:
            saved[f'{kind}__truth_p'] = pt; saved[f'{kind}__truth_v'] = vt
            saved[f'{kind}__probabilities'] = histories['IMM']['probabilities']
            saved[f'{kind}__omega'] = histories['IMM']['branches'][..., 2, 6]
        print(f'  {kind} complete ({time.perf_counter()-start:.1f}s)', flush=True)
    return rows, saved


def main():
    ap=argparse.ArgumentParser(); ap.add_argument('--trials',type=int,default=20)
    ap.add_argument('--seed',type=int,default=20261028)
    ap.add_argument('--output',type=Path,default=HERE/'outputs')
    args=ap.parse_args(); args.output.mkdir(parents=True,exist_ok=True)
    cfg=base.Config(seed=args.seed); par=Parameters(); rows=[]; start=time.perf_counter()
    for i in range(args.trials):
        a,s=run_seed(replace(cfg,seed=cfg.seed+i),par,keep=i==0); rows+=a
        if i==0: np.savez_compressed(args.output/'first_seed.npz',**s)
        base.write_csv(args.output/'per_seed.csv',rows)
        print(f'SEED {i+1}/{args.trials} saved; elapsed {time.perf_counter()-start:.1f}s',flush=True)
    hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in (HERE/'imm_filter.py',HERE/'run_experiment.py',ROOT/'kf_dpc_coherent_5nodes/gm_dpc_unified/experiment.py')}
    (args.output/'config.json').write_text(json.dumps(dict(config=asdict(cfg),imm=asdict(par),
        trials=args.trials,seeds=[cfg.seed,cfg.seed+args.trials-1],elapsed_s=time.perf_counter()-start,
        source_sha256=hashes,parameter_policy='Exploratory fixed parameters; not calibrated or held-out validated',
        prediction='frozen mode over <=177ms; discrete additive CT process noise',
        timing='Assumed identical 12ms processing +25x5ms; hardware and larger payload unvalidated'),indent=2),encoding='utf-8')

if __name__=='__main__': main()
