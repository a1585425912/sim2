"""Plot saved IMM+DPC coordinates/velocities and errors; no simulation rerun."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

HERE = Path(__file__).resolve().parent
OUT = HERE/'outputs'
KINDS = {'hover':'悬停', 'straight':'直线', 'circle':'圆周', 'mixed':'混合运动'}
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],
                     'axes.unicode_minus':False, 'font.size':10,
                     'axes.spines.top':False, 'axes.spines.right':False})


def main():
    data = np.load(OUT/'first_seed.npz')
    t = data['time_s']; stable = t >= 5
    rows = []
    for kind, name in KINDS.items():
        state = data[f'{kind}__IMM__state']
        truth = np.concatenate([data[f'{kind}__truth_p'], data[f'{kind}__truth_v']], axis=-1)
        assert state.shape == (len(t), 20, 6) and truth.shape == (len(t), 6)
        error = state-truth[:, None]
        assert np.isfinite(error).all()
        fig, axes = plt.subplots(6, 2, figsize=(14, 15), sharex='col', layout='constrained')
        for j in range(6):
            axis = 'XYZ'[j % 3]
            unit = 'm' if j < 3 else 'm/s'
            label = ('位置 ' if j < 3 else '速度 v')+axis
            rmse = float(np.sqrt(np.mean(error[stable, :, j]**2)))
            rows.append(dict(trajectory=kind, component=label, unit=unit,
                             all_nodes_rmse_after5=rmse,
                             node1_rmse_after5=float(np.sqrt(np.mean(error[stable, 0, j]**2))),
                             all_nodes_rmse_full=float(np.sqrt(np.mean(error[..., j]**2)))))
            left, right = axes[j]
            left.plot(t, state[:, 0, j], color='#c4504f', lw=1.25, label='IMM+DPC：节点1')
            left.plot(t, truth[:, j], color='#222222', lw=1.5, ls='--', label='真值')
            left.axvspan(0, 5, color='#cccccc', alpha=.25)
            lo, hi = np.quantile(error[stable, :, j], [.1, .9], axis=1)
            right.fill_between(t[stable], lo, hi, color='#4686b8', alpha=.2,
                               label='20节点误差的10%–90%分位范围')
            right.plot(t[stable], error[stable, 0, j], color='#c4504f', lw=1.1, label='节点1：估计−真值')
            right.axhline(0, color='#222222', ls='--', lw=1)
            right.text(.98, .94, f'20节点RMSE = {rmse:.3f} {unit}', transform=right.transAxes,
                       ha='right', va='top', fontsize=9,
                       bbox=dict(facecolor='white', edgecolor='none', alpha=.85))
            left.set_ylabel(f'{label} ({unit})'); right.set_ylabel(f'{label}误差 ({unit})')
            for ax in (left, right):
                ax.grid(alpha=.2)
                fmt = ScalarFormatter(useOffset=False); fmt.set_scientific(False)
                ax.yaxis.set_major_formatter(fmt)
                if kind == 'mixed':
                    for boundary in (15, 20, 35, 55, 65, 70):
                        ax.axvline(boundary, color='gray', lw=.5, ls=':')
        axes[0, 0].set_title('完整0–90s：三轴位置与速度对照（灰底为启动0–5s）')
        axes[0, 1].set_title('5–90s：误差细节（分位范围不是置信区间）')
        axes[0, 0].legend(fontsize=8, loc='lower right')
        axes[0, 1].legend(fontsize=8, loc='lower right')
        for ax in axes[-1]: ax.set_xlabel('时间 (s)')
        axes[0, 0].set_xlim(0, 90); axes[0, 1].set_xlim(5, 90)
        fig.suptitle(f'{name}：IMM+DPC三轴坐标/速度与真值的差距\n'
                     '首种子20261028；20浮标、25轮DPC；节点1曲线与全节点误差统计', fontsize=14)
        fig.savefig(OUT/f'xyz_{kind}.png', dpi=160)
        plt.close(fig)
    with (OUT/'xyz_rmse_first_seed.csv').open('w', encoding='utf-8-sig', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0])); writer.writeheader(); writer.writerows(rows)
    print('Generated four XYZ figures and first-seed per-axis RMSE table.')
    for r in rows:
        if r['trajectory']=='mixed': print(r['component'], round(r['all_nodes_rmse_after5'], 4), r['unit'])


if __name__=='__main__': main()
