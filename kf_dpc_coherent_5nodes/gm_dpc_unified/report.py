"""由 outputs/ 中已保存的统计与首种子数据画图。python kf_dpc_coherent_5nodes/gm_dpc_unified/report.py"""
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = Path(__file__).resolve().parent/'outputs'
TRAJ = {'hover': '悬停', 'straight': '直线匀速', 'circle': '匀速圆周'}
COLOR = {'hover': '#2a78d6', 'straight': '#eb6834', 'circle': '#1baf7a'}
INK, MUTED, GRID = '#1f1f1e', '#6b6a64', '#e4e3df'
plt.rcParams.update({'font.sans-serif': ['Microsoft YaHei', 'SimHei', 'DejaVu Sans'],
                     'axes.unicode_minus': False, 'font.size': 10, 'axes.spines.top': False,
                     'axes.spines.right': False, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
                     'xtick.color': MUTED, 'ytick.color': MUTED, 'axes.grid': True,
                     'grid.color': GRID, 'grid.linewidth': .6, 'lines.linewidth': 2})


def read(name):
    return list(csv.DictReader(open(OUT/name, encoding='utf-8-sig')))


def save(fig, name):
    for ext in ('png', 'svg'):
        fig.savefig(OUT/f'{name}.{ext}', dpi=180, bbox_inches='tight')
    plt.close(fig)


def rounds_scan(tx):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.4), sharey=False, layout='constrained')
    for ax, f in zip(axes, ('300.0', '1000.0')):
        for kind in TRAJ:
            for timing, ls in (('5ms', '-'), ('1ms', '--')):
                pts = sorted((int(r['rounds']), float(r['eta'])) for r in tx
                             if r['trajectory'] == kind and r['timing'] == timing and r['frequency_mhz'] == f
                             and r['variant'] == 'main' and r['nav_mode'] == 'nominal' and r['phase'] == 'after_startup')
                L, eta = np.array(pts).T
                ax.plot(L, eta, ls, color=COLOR[kind], marker='o', ms=4 if timing == '1ms' else 5,
                        label=f'{TRAJ[kind]}，{timing}/轮')
                if timing == '1ms':
                    ax.annotate(TRAJ[kind], (L[-1], eta[-1]), xytext=(4, 0), textcoords='offset points',
                                va='center', color=INK, fontsize=9)
        ax.axvline(25, color=MUTED, lw=1, ls=':')
        ax.text(25, 0.02, ' L=25', color=MUTED, fontsize=9, va='bottom')
        ax.set_xscale('symlog', linthresh=10)
        ax.set_xticks([0, 1, 3, 5, 10, 25, 50, 140]); ax.set_xticklabels(['0', '1', '3', '5', '10', '25', '50', '140'])
        ax.set_xlabel('每帧DPC轮数 L'); ax.set_ylabel('相干效率 η'); ax.set_ylim(0, 1 if f == '300.0' else .5)
        ax.set_title(f'{float(f):.0f} MHz（实线：5ms/轮，L≤29；虚线：假设1ms/轮）', fontsize=10, color=INK)
    axes[0].legend(fontsize=8, ncol=2, frameon=False, loc='upper left')
    fig.suptitle('共识轮数对相干效率的影响（20浮标、2×2子阵、nominal导航、t≥5s、20种子）', color=INK)
    save(fig, 'rounds_scan')


def velocity(d):
    t = d['time_s']
    fig, axes = plt.subplots(3, 3, figsize=(14, 8.5), sharex=True, layout='constrained')
    for row, kind in enumerate(TRAJ):
        tv, loc, fus = d[f'{kind}__truth_v'], d[f'{kind}__local'], d[f'{kind}__dpc25']
        for j, a in enumerate('XYZ'):
            ax = axes[row, j]
            ax.plot(t, loc[:, 0, 3+j], color='#86b6ef', lw=1, label='节点1本地KF（L=0）')
            ax.plot(t, fus[:, 0, 3+j], color='#eb6834', lw=1.4, label='节点1 KF+DPC（L=25）')
            ax.plot(t, tv[:, j], color=INK, lw=1.4, ls='--', label='真值（仅评价）')
            ax.axvspan(0, 5, color=GRID, alpha=.6, lw=0)
            lo, hi = tv[t >= 5, j].min(), tv[t >= 5, j].max()
            pad = max(1.5, .25*(hi-lo))
            ax.set_ylim(lo-pad, hi+pad)
            if row == 0: ax.set_title(f'v{a}', color=INK)
            if j == 0: ax.set_ylabel(f'{TRAJ[kind]}\n速度 m/s')
            if row == 2: ax.set_xlabel('时间 s')
    axes[0, 0].legend(fontsize=8, frameon=False, loc='upper right')
    fig.suptitle('首个种子的三轴速度估计（灰底为启动段0–5s，纵轴截掉了首帧零初速的大误差）', color=INK)
    save(fig, 'velocity_xyz')


def eta_time(d):
    t = d['time_s']; k = 11
    smooth = lambda y: np.convolve(y, np.ones(k), 'same')/np.convolve(np.ones_like(y), np.ones(k), 'same')
    fig, axes = plt.subplots(3, 1, figsize=(11, 7.5), sharex=True, sharey=True, layout='constrained')
    for ax, kind in zip(axes, TRAJ):
        ax.plot(t, smooth(d[f'{kind}__eta300_local']), color='#86b6ef', lw=1.4, label='本地KF，同样137ms延迟')
        ax.plot(t, smooth(d[f'{kind}__eta300_dpc']), color='#eb6834', lw=1.6, label='KF+DPC（L=25）')
        ax.set_ylim(0, 1); ax.set_ylabel(f'{TRAJ[kind]}\nη')
    axes[0].legend(fontsize=8, frameon=False, ncol=2, loc='upper right')
    axes[-1].set_xlabel('时间 s')
    fig.suptitle('首个种子300MHz相干效率随时间变化（11帧滑动平均仅用于显示）', color=INK)
    save(fig, 'eta_time')


def trajectories(d):
    fig, ax = plt.subplots(figsize=(8.5, 6.5), layout='constrained')
    b = d['buoy_truth_frames']
    for i in range(b.shape[1]):
        ax.plot(b[:, i, 0], b[:, i, 1], color=MUTED, lw=1, label='浮标GM漂移真值（90s）' if i == 0 else None)
        ax.plot(b[0, i, 0], b[0, i, 1], 'o', color=MUTED, ms=3)
    for kind in TRAJ:
        p = d[f'{kind}__truth_p']; e = d[f'{kind}__dpc25'][:, 0, :3]
        ax.plot(p[:, 0], p[:, 1], color=COLOR[kind], lw=2, label=f'{TRAJ[kind]}真值')
        ax.plot(e[:, 0], e[:, 1], color=COLOR[kind], lw=.8, ls=':', alpha=.9)
    ax.plot(0, 0, '*', color=COLOR['hover'], ms=12)
    ax.set_xlim(-260, 420); ax.set_ylim(-260, 260); ax.set_aspect('equal')
    ax.set_xlabel('X m'); ax.set_ylabel('Y m')
    ax.legend(fontsize=8, frameon=False, loc='upper right')
    ax.set_title('XY平面几何（直线轨迹截至X=420m，实际飞到1350m；点线为节点1的DPC估计）', fontsize=10, color=INK)
    save(fig, 'geometry_xy')


if __name__ == '__main__':
    tx = read('transmit_summary.csv'); d = np.load(OUT/'first_seed.npz')
    rounds_scan(tx); velocity(d); eta_time(d); trajectories(d)
    print('figures saved to', OUT)
