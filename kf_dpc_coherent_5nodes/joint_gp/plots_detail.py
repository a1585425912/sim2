"""Readable full-time errors without the huge initial velocity prior band."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':10})


def main():
    out=HERE/'outputs';d=np.load(out/'first_seed_hs1.5.npz');t=d['time_s'];kind='mixed'
    methods={'IMM_CI_linear':('旧IMM+CI','#888888'),'JointGP_accel':('高斯联合GP','#bf4859'),'JointGP_TTrend':('T-Trend联合估计','#228a76')}
    truth=np.c_[d[f'{kind}__truth_p'],d[f'{kind}__truth_v']]
    fig,axes=plt.subplots(3,2,figsize=(14,9),sharex=True,layout='constrained')
    for axis,name in enumerate('XYZ'):
        for method,(label,color) in methods.items():
            value=d[f'{kind}__{method}__state'][:,0]
            axes[axis,0].plot(t,value[:,axis]-truth[:,axis],color=color,label=label,lw=1)
            axes[axis,1].plot(t,value[:,axis+3]-truth[:,axis+3],color=color,label=label,lw=1)
        for col,unit in ((0,'m'),(1,'m/s')):
            axes[axis,col].axhline(0,color='black',ls='--',lw=.8);axes[axis,col].set_ylabel(f'{name}误差 ({unit})');axes[axis,col].grid(alpha=.2)
    axes[0,0].set_title('位置估计减真值');axes[0,1].set_title('速度估计减真值');axes[0,0].legend(ncol=3,fontsize=8)
    for ax in axes[-1]:ax.set_xlabel('时间 (s)')
    fig.suptitle('混合机动：首种子、节点1、全部0–90s；无平滑、无裁剪、无样本筛选')
    fig.savefig(out/'mixed_errors_readable.png',dpi=170);plt.close(fig)


if __name__=='__main__':main()
