"""GM修订版文档的统一实验：一阶GM漂移 → 位置观测KF → DPC → 几何相位补偿。

python kf_dpc_coherent_5nodes/gm_dpc_unified/experiment.py --trials 20

只有本文件的仿真/评价部分读取真值；KF、DPC和发射命令只用带噪雷达、导航输出和邻居状态。
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for path in (HERE.parent, ROOT):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))
import argparse, csv, json, time
from dataclasses import dataclass, asdict, replace
import numpy as np
from scipy.signal import lfilter
from pipeline import ring_weights, dpc
from position_only_velocity_demo.model import spherical, cartesian, wrap_angle, cv_matrices

C = 299792458.
TRAJECTORIES = ('hover', 'straight', 'circle')
NAV_MODES = ('ideal', 'nominal')
# 原网络每轮5ms：12+5L+40<=200ms → L<=29；另评价假设每轮1ms的更快网络。
TIMINGS = {'5ms': .005, '1ms': .001}
ROUNDS = (0, 1, 3, 5, 10, 15, 20, 25, 29, 50, 100, 140)
VARIANTS = ('buoy_cv', 'buoy_frozen', 'uav_hold')


@dataclass(frozen=True)
class Config:
    seed: int = 20261028
    nodes: int = 20
    dt: float = .2
    duration: float = 90.
    warmup: float = 5.
    deployment_radius_m: float = 100.
    buoy_height_m: float = 1.5
    # 一阶GM（dfpc_experiment correlated_drift 分支参数）
    mean_current_mps: tuple = (.30, 0.)
    gm_tau_s: float = 20.
    gm_sigma_mps: float = .05
    truth_grid_s: float = .001
    # 雷达与UAV滤波
    range_std_m: float = .8
    angle_std_deg: float = .35
    acceleration_psd: float = .5
    initial_velocity_std: float = 20.
    # 导航输出（nominal：独立白误差；取 buoy_cbf_simple NavConfig 的独立分量）
    nav_position_std_m: float = .02
    nav_velocity_std_mps: float = .02
    nav_roll_pitch_std_deg: float = .10
    nav_yaw_std_deg: float = .30
    # 通信、发射
    main_rounds: int = 25
    processing_s: float = .012
    slot_s: float = .04
    slot_samples: int = 9
    frame_budget_s: float = .2
    subarray: int = 2
    subarray_reference_mhz: float = 300.
    frequencies_mhz: tuple = (300., 1000.)
    total_power_w: float = 1.


def uav_trajectory(t, kind):
    t = np.asarray(t, dtype=float)
    p = np.zeros(t.shape + (3,)); v = np.zeros_like(p)
    if kind == 'hover':
        p[...] = [0., 0., 20.]
    elif kind == 'straight':
        p[..., 0] = 15*t; p[..., 1] = -200.; p[..., 2] = 120.; v[..., 0] = 15.
    elif kind == 'circle':
        w = .075
        p[..., 0] = 200*np.sin(w*t); p[..., 1] = -200*np.cos(w*t); p[..., 2] = 120.
        v[..., 0] = 15*np.cos(w*t); v[..., 1] = 15*np.sin(w*t)
    else:
        raise ValueError(kind)
    return p, v


def gm_buoy_truth(cfg, rng, horizon_s):
    """式(5)(6)：GM速度精确离散 + 位置梯形积分，在1ms细网格上生成。

    初始速度偏差为0；垂向保持 buoy_height_m；姿态固定为单位阵、角速度为0。
    """
    h = cfg.truth_grid_s
    n = int(round(horizon_s/h)) + 1
    rho = np.exp(-h/cfg.gm_tau_s)
    eps = rng.standard_normal((n-1, cfg.nodes, 2))
    drive = np.concatenate([np.zeros((1, cfg.nodes, 2)), cfg.gm_sigma_mps*np.sqrt(1-rho**2)*eps])
    u = lfilter([1.], [1., -rho], drive, axis=0)
    vbar = np.array(cfg.mean_current_mps)
    az = 2*np.pi*np.arange(cfg.nodes)/cfg.nodes
    b0 = np.c_[cfg.deployment_radius_m*np.cos(az), cfg.deployment_radius_m*np.sin(az)]
    steps = h*(vbar + .5*(u[1:] + u[:-1]))
    xy = b0[None] + np.concatenate([np.zeros((1, cfg.nodes, 2)), np.cumsum(steps, axis=0)])
    pos = np.concatenate([xy, np.full((n, cfg.nodes, 1), cfg.buoy_height_m)], axis=-1)
    vel = np.concatenate([vbar + u, np.zeros((n, cfg.nodes, 1))], axis=-1)
    return pos, vel


def so3_exp(w):
    """Rodrigues，w (...,3) → (...,3,3)。"""
    th = np.linalg.norm(w, axis=-1)[..., None, None]
    K = skew(w)
    small = th < 1e-12
    a = np.where(small, 1., np.sin(th)/np.where(small, 1., th))
    b = np.where(small, .5, (1-np.cos(th))/np.where(small, 1., th**2))
    return np.eye(3) + a*K + b*(K@K)


def skew(w):
    w = np.asarray(w)
    z = np.zeros(w.shape[:-1])
    return np.stack([np.stack([z, -w[..., 2], w[..., 1]], -1),
                     np.stack([w[..., 2], z, -w[..., 0]], -1),
                     np.stack([-w[..., 1], w[..., 0], z], -1)], -2)


def cartesian_jacobian(z):
    r, a, e = np.moveaxis(z, -1, 0)
    ca, sa, ce, se = np.cos(a), np.sin(a), np.cos(e), np.sin(e)
    return np.stack([np.stack([ce*ca, -r*ce*sa, -r*se*ca], -1),
                     np.stack([ce*sa, r*ce*ca, -r*se*sa], -1),
                     np.stack([se, np.zeros_like(r), r*ce], -1)], -2)


def navigation(cfg, b, v, rng, mode):
    """式(8)(9)：导航输出 = 真值 + 独立白误差；ideal 模式误差为零（同一组随机数）。"""
    nt, nn, _ = b.shape
    e_b = rng.standard_normal((nt, nn, 3)); e_v = rng.standard_normal((nt, nn, 3))
    e_a = rng.standard_normal((nt, nn, 3))
    att_std = np.deg2rad([cfg.nav_roll_pitch_std_deg, cfg.nav_roll_pitch_std_deg, cfg.nav_yaw_std_deg])
    scale = 0. if mode == 'ideal' else 1.
    p_std = scale*cfg.nav_position_std_m
    return dict(b=b + p_std*e_b, v=v + scale*cfg.nav_velocity_std_mps*e_v,
                C=so3_exp(scale*att_std*e_a),
                P_b=p_std**2*np.eye(3), P_att=np.diag((scale*att_std)**2))


def position_observations(z, nav, R_sph):
    """式(12)(13)：世界系位置观测与一阶协方差（含导航位置/姿态边际项）。"""
    g = cartesian(z)
    Ch = nav['C']
    p = nav['b'] + np.einsum('tnij,tnj->tni', Ch, g)
    A = Ch @ cartesian_jacobian(z)
    B = -Ch @ skew(g)
    R = A @ R_sph @ np.swapaxes(A, -1, -2) + nav['P_b'] + B @ nav['P_att'] @ np.swapaxes(B, -1, -2)
    return p, R


def local_kf(t, z, R, q, sigma_v0):
    """式(14)–(22)：各节点独立的六维CV线性KF（对节点批量计算，与 position_filter 相同）。"""
    nt, nn, _ = z.shape
    F, Q = cv_matrices(t[1]-t[0], q)
    x = np.concatenate([z[0], np.zeros((nn, 3))], axis=-1)
    P = np.zeros((nn, 6, 6)); P[:, :3, :3] = R[0]; P[:, 3:, 3:] = sigma_v0**2*np.eye(3)
    H = np.c_[np.eye(3), np.zeros((3, 3))]
    out = np.empty((nt, nn, 6)); out[0] = x
    for k in range(1, nt):
        x = x @ F.T; P = F @ P @ F.T + Q
        S = P[:, :3, :3] + R[k]
        K = np.swapaxes(np.linalg.solve(S, np.swapaxes(P[:, :, :3], 1, 2)), 1, 2)
        x = x + np.einsum('nij,nj->ni', K, z[k] - x[:, :3])
        IKH = np.eye(6) - K @ H
        P = IKH @ P @ np.swapaxes(IKH, 1, 2) + K @ R[k] @ np.swapaxes(K, 1, 2)
        P = .5*(P + np.swapaxes(P, 1, 2))
        out[k] = x
    return out


def buoy_prediction(nav, tau, cfg, mode):
    """式(28)(29)GM条件均值；buoy_cv=常速度外推；buoy_frozen=保持当前位置、速度0。"""
    b, v = nav['b'], nav['v']
    if mode == 'gm':
        vbar = np.r_[cfg.mean_current_mps, 0.]
        u = v - vbar; u[..., 2] = 0.
        decay = np.exp(-tau/cfg.gm_tau_s)
        b_tx = b + vbar*tau + cfg.gm_tau_s*(1-decay)*u
        v_tx = vbar + decay*u
        b_tx[..., 2] = b[..., 2]
    elif mode == 'cv':
        b_tx = b + tau*v; v_tx = v.copy()
    elif mode == 'frozen':
        b_tx = b.copy(); v_tx = np.zeros_like(v)
    else:
        raise ValueError(mode)
    return b_tx, v_tx


def commands(states, nav, levers, tau, slot, cfg, buoy_mode='gm', use_uav_velocity=True):
    """式(27)(31)–(33)：只用DPC状态与导航输出。固定姿态、ω=0，故阵元速度=浮标速度。返回 (T,S,N,M)。"""
    vu = states[..., 3:] if use_uav_velocity else np.zeros_like(states[..., 3:])
    pu = states[..., :3] + tau*vu
    b_tx, v_tx = buoy_prediction(nav, tau, cfg, buoy_mode)
    elem = b_tx[:, :, None] + np.einsum('tnij,mj->tnmi', nav['C'], levers)
    q = pu[:, :, None] - elem
    d = np.linalg.norm(q, axis=-1)
    rate = np.einsum('tnmi,tnmi->tnm', q/d[..., None], vu[:, :, None] - v_tx[:, :, None])
    return d[:, None] + slot[None, :, None, None]*rate[:, None]


def powers(cmd, true_d, freq_mhz, total_power):
    """式(34)–(37)，同步残差为0。"""
    n, m = true_d.shape[-2:]
    lam = C/(freq_mhz*1e6)
    amp = np.sqrt(total_power/(n*m))*lam/(4*np.pi*true_d)
    field = amp*np.exp(1j*2*np.pi/lam*(cmd - true_d))
    node = field.sum(-1)
    return dict(rx=np.abs(node.sum(-1))**2, node=np.sum(np.abs(node)**2, -1),
                elem=np.sum(amp**2, (-1, -2)), ideal=amp.sum((-1, -2))**2)


def levers_for(cfg):
    d = C/(cfg.subarray_reference_mhz*1e6)/2
    s = (np.arange(cfg.subarray) - (cfg.subarray-1)/2)*d
    X, Y = np.meshgrid(s, s, indexing='ij')
    return np.c_[X.ravel(), Y.ravel(), np.zeros(X.size)]


def run_seed(cfg, keep=False):
    t = np.arange(round(cfg.duration/cfg.dt) + 1)*cfg.dt
    nt = len(t)
    phases = {'all': np.ones(nt, bool), 'startup': t < cfg.warmup, 'after_startup': t >= cfg.warmup}
    ss = np.random.SeedSequence(cfg.seed).spawn(3)
    horizon = cfg.duration + cfg.frame_budget_s
    b_fine, v_fine = gm_buoy_truth(cfg, np.random.default_rng(ss[0]), horizon)
    step = int(round(cfg.dt/cfg.truth_grid_s))
    frame_idx = np.arange(nt)*step
    b_k, v_k = b_fine[frame_idx], v_fine[frame_idx]
    nav_rng_state = np.random.default_rng(ss[1]).bit_generator.state
    std = np.array([cfg.range_std_m, np.deg2rad(cfg.angle_std_deg), np.deg2rad(cfg.angle_std_deg)])
    radar_noise = np.random.default_rng(ss[2]).standard_normal((nt, cfg.nodes, 3))*std
    R_sph = np.diag(std**2)
    W = ring_weights(cfg.nodes)
    levers = levers_for(cfg)
    slot = np.linspace(0, cfg.slot_s, cfg.slot_samples)
    slot_ms = np.rint(slot/cfg.truth_grid_s).astype(int)
    track, power_rows, saved = [], [], {}

    def add_power(kind, mode, timing, L, variant, cmd, true_d):
        err = cmd - true_d
        err = err - err.mean(axis=(-1, -2), keepdims=True)
        for f in cfg.frequencies_mhz:
            pw = powers(cmd, true_d, f, cfg.total_power_w)
            for ph, mk in phases.items():
                power_rows.append(dict(seed=cfg.seed, trajectory=kind, nav_mode=mode, timing=timing, rounds=L,
                    variant=variant, frequency_mhz=f, phase=ph, samples=int(mk.sum()*len(slot)),
                    **{k+'_sum': float(v[mk].sum()) for k, v in pw.items()},
                    diff_path_sq_sum=float(np.sum(err[mk]**2)), diff_path_count=int(err[mk].size)))

    for mode in NAV_MODES:
        rng = np.random.default_rng(); rng.bit_generator.state = nav_rng_state   # 两种导航模式用同一组创新
        nav = navigation(cfg, b_k, v_k, rng, mode)
        for kind in TRAJECTORIES:
            p_true, v_true = uav_trajectory(t, kind)
            # 真姿态为单位阵：雷达系相对位置 q = p_u - b
            z = spherical(p_true[:, None] - b_k) + radar_noise
            z[..., 1] = wrap_angle(z[..., 1])
            obs, R = position_observations(z, nav, R_sph)
            local = local_kf(t, obs, R, cfg.acceleration_psd, cfg.initial_velocity_std)
            truth = np.concatenate([p_true, v_true], -1)[:, None]
            fused = {L: dpc(local, W, L) for L in ROUNDS}
            for L, x in fused.items():
                dev = x - x.mean(axis=1, keepdims=True)
                e = x - truth
                for ph, mk in phases.items():
                    track.append(dict(seed=cfg.seed, trajectory=kind, nav_mode=mode, rounds=L, phase=ph,
                        points=int(mk.sum()*cfg.nodes),
                        **{f'{q}_{a}_sq_sum': float(np.sum(e[mk][..., o+j]**2))
                           for q, o in (('pos', 0), ('vel', 3)) for j, a in enumerate('XYZ')},
                        dis_pos_sq_sum=float(np.sum(dev[mk][..., :3]**2)),
                        dis_vel_sq_sum=float(np.sum(dev[mk][..., 3:]**2))))
            cache = {}
            for timing, per_round in TIMINGS.items():
                for L in ROUNDS:
                    tau = cfg.processing_s + per_round*L
                    if tau + cfg.slot_s > cfg.frame_budget_s + 1e-12:
                        continue
                    tau_ms = int(round(tau/cfg.truth_grid_s))
                    if tau_ms not in cache:
                        idx = frame_idx[:, None] + tau_ms + slot_ms[None]
                        pu_tx, _ = uav_trajectory(t[:, None] + tau + slot[None], kind)
                        elem = b_fine[idx][:, :, :, None] + levers[None, None, None]
                        cache[tau_ms] = np.linalg.norm(pu_tx[:, :, None, None] - elem, axis=-1)
                    true_d = cache[tau_ms]
                    add_power(kind, mode, timing, L, 'main', commands(fused[L], nav, levers, tau, slot, cfg), true_d)
                    if timing == '5ms' and L == cfg.main_rounds:
                        x = fused[L]
                        add_power(kind, mode, timing, L, 'buoy_cv', commands(x, nav, levers, tau, slot, cfg, 'cv'), true_d)
                        add_power(kind, mode, timing, L, 'buoy_frozen', commands(x, nav, levers, tau, slot, cfg, 'frozen'), true_d)
                        add_power(kind, mode, timing, L, 'uav_hold', commands(x, nav, levers, tau, slot, cfg, 'gm', False), true_d)
                        add_power(kind, mode, timing, L, 'local_kf_same_delay', commands(fused[0], nav, levers, tau, slot, cfg), true_d)
                        if keep and mode == 'nominal':
                            pw = powers(commands(x, nav, levers, tau, slot, cfg), true_d, 300., cfg.total_power_w)
                            pw0 = powers(commands(fused[0], nav, levers, tau, slot, cfg), true_d, 300., cfg.total_power_w)
                            saved[f'{kind}__eta300_dpc'] = pw['rx'].sum(1)/pw['ideal'].sum(1)
                            saved[f'{kind}__eta300_local'] = pw0['rx'].sum(1)/pw0['ideal'].sum(1)
            if keep and mode == 'nominal':
                saved[f'{kind}__truth_p'] = p_true; saved[f'{kind}__truth_v'] = v_true
                saved[f'{kind}__local'] = local; saved[f'{kind}__dpc25'] = fused[cfg.main_rounds]
    if keep:
        saved['time_s'] = t; saved['buoy_truth_frames'] = b_k
    return track, power_rows, saved


def pool(rows, keys, fn):
    groups = {}
    for r in rows:
        groups.setdefault(tuple(r[k] for k in keys), []).append(r)
    return [dict(zip(keys, g), trials=len(v), **fn(v)) for g, v in groups.items()]


def track_stats(v):
    n = sum(r['points'] for r in v); s = lambda k: sum(r[k] for r in v)
    out = {f'{q}_{a}_rmse': float(np.sqrt(s(f'{q}_{a}_sq_sum')/n)) for q in ('pos', 'vel') for a in 'XYZ'}
    out['pos_3d_rmse'] = float(np.sqrt(sum(s(f'pos_{a}_sq_sum') for a in 'XYZ')/n))
    out['vel_3d_rmse'] = float(np.sqrt(sum(s(f'vel_{a}_sq_sum') for a in 'XYZ')/n))
    out['disagree_pos_rms'] = float(np.sqrt(s('dis_pos_sq_sum')/n))
    out['disagree_vel_rms'] = float(np.sqrt(s('dis_vel_sq_sum')/n))
    return out


def power_stats(v):
    s = lambda k: sum(r[k] for r in v)
    rx, node, elem, ideal = s('rx_sum'), s('node_sum'), s('elem_sum'), s('ideal_sum')
    return dict(eta=rx/ideal, G_sub=node/elem, G_inter=rx/node, G_total=rx/elem,
                G_total_db=float(10*np.log10(rx/elem)), G_ideal=ideal/elem,
                rx_power_w=rx/s('samples'), diff_path_rmse_m=float(np.sqrt(s('diff_path_sq_sum')/s('diff_path_count'))))


def write_csv(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument('--trials', type=int, default=20); args = ap.parse_args()
    cfg = Config(); out = HERE/'outputs'; out.mkdir(parents=True, exist_ok=True)
    start = time.perf_counter(); track, power, saved = [], [], None
    for i in range(args.trials):
        a, b, s = run_seed(replace(cfg, seed=cfg.seed + i), keep=i == 0)
        track += a; power += b; saved = saved or s
        print(f'seed {i+1}/{args.trials} done ({time.perf_counter()-start:.1f}s)', flush=True)
    tracking = pool(track, ('trajectory', 'nav_mode', 'rounds', 'phase'), track_stats)
    transmit = pool(power, ('trajectory', 'nav_mode', 'timing', 'rounds', 'variant', 'frequency_mhz', 'phase'), power_stats)
    write_csv(out/'tracking_summary.csv', tracking); write_csv(out/'transmit_summary.csv', transmit)
    write_csv(out/'tracking_per_seed.csv', track); write_csv(out/'transmit_per_seed.csv', power)
    np.savez_compressed(out/'first_seed.npz', **saved)
    (out/'config.json').write_text(json.dumps(dict(config=asdict(cfg), trials=args.trials,
        seeds=[cfg.seed, cfg.seed + args.trials - 1], elapsed_s=time.perf_counter()-start), indent=2), encoding='utf-8')
    print('saved', out)


if __name__ == '__main__':
    main()
