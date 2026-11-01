"""Full-run CI comparison. No gating, smoothing, outlier rejection or rerun."""
from pathlib import Path
from collections import defaultdict
import csv, json, argparse
import numpy as np
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter, MaxNLocator

HERE=Path(__file__).resolve().parent
LABELS={'hover':'悬停','straight':'直线','circle':'圆周','mixed':'混合运动','climb':'直线爬升'}
METHODS={'local':('本地IMM','#aaaaaa'),'dpc':('IMM+DPC','#3786ba'),'ci':('IMM+CI','#c25054')}
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],
                     'axes.unicode_minus':False,'font.size':10,'axes.spines.top':False,'axes.spines.right':False})


def mean_ci(values):
    a=np.asarray(values)
    return float(a.mean()),float(student_t.ppf(.975,len(a)-1)*a.std(ddof=1)/np.sqrt(len(a))) if len(a)>1 else float('nan')


def plots(out,summary):
    d=np.load(out/'first_seed.npz');t=d['time_s']
    def result(kind,method,frequency=300.,phase='all'):
        return next(r for r in summary if (r['trajectory'],r['method'],r['frequency_mhz'],r['phase'])==(kind,method,frequency,phase))
    kinds=list(LABELS);x=np.arange(len(kinds))
    fig,axes=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    for ax,metric,title in zip(axes.ravel(),('pos_rmse','vel_rmse','pZ_rmse','eta_mean'),
                               ('三维位置 RMSE (m)','三维速度 RMSE (m/s)','高度Z RMSE (m)','300MHz 相干效率 η')):
        for j,(method,(label,color)) in enumerate(METHODS.items()):
            vals=[result(k,method)[metric] for k in kinds]
            ax.bar(x+(j-1)*.25,vals,width=.23,color=color,label=label)
        ax.set_xticks(x,[LABELS[k] for k in kinds]);ax.set_title(title);ax.grid(axis='y',alpha=.2)
        if metric=='eta_mean':ax.set_ylim(0,1)
    axes[0,0].legend(fontsize=9)
    fig.suptitle('CI完整配对实验：20个种子、0–90s全部帧、相同时延（无筛选）')
    fig.savefig(out/'ci_comparison.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(5,2,figsize=(14,13),sharex=True,layout='constrained')
    for row,(kind,label) in enumerate(LABELS.items()):
        truth=d[f'{kind}__truth_p'][:,2]
        axes[row,0].plot(t,truth,'k--',label='真值',lw=1.4)
        for method in ('dpc','ci'):
            name,color=METHODS[method];z=d[f'{kind}__{method}__state'][:,0,2]
            axes[row,0].plot(t,z,color=color,label=name,lw=1)
            axes[row,1].plot(t,z-truth,color=color,label=name,lw=1)
        axes[row,1].axhline(0,color='black',ls='--',lw=1)
        axes[row,0].set_ylabel(label+'\nZ (m)');axes[row,1].set_ylabel('Z误差 (m)')
        for ax in axes[row]:
            ax.grid(alpha=.2);ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    axes[0,0].set_title('高度坐标与真值');axes[0,1].set_title('估计高度减真值')
    axes[0,0].legend(fontsize=9)
    for ax in axes[-1]:ax.set_xlabel('时间 (s)')
    fig.suptitle('首种子20261028，节点1：Z轴完整时间曲线，无平滑')
    fig.savefig(out/'z_comparison.png',dpi=155);plt.close(fig)
    kind='mixed';truth=np.concatenate([d[f'{kind}__truth_p'],d[f'{kind}__truth_v']],axis=-1)
    fig,axes=plt.subplots(6,2,figsize=(14,15),sharex=True,layout='constrained')
    for j in range(6):
        axis='XYZ'[j%3];unit='m' if j<3 else 'm/s';name=('位置 ' if j<3 else '速度 v')+axis
        axes[j,0].plot(t,truth[:,j],'k--',lw=1.3,label='真值')
        for method in ('dpc','ci'):
            label,color=METHODS[method];value=d[f'{kind}__{method}__state'][:,0,j]
            axes[j,0].plot(t,value,color=color,lw=1,label=label)
            axes[j,1].plot(t,value-truth[:,j],color=color,lw=1,label=label)
        axes[j,1].axhline(0,color='black',ls='--',lw=.8)
        axes[j,0].set_ylabel(f'{name} ({unit})');axes[j,1].set_ylabel(f'误差 ({unit})')
        for ax in axes[j]:ax.grid(alpha=.2);ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    axes[0,0].legend();axes[0,0].set_title('真值与估计');axes[0,1].set_title('估计减真值')
    for ax in axes[-1]:ax.set_xlabel('时间 (s)')
    fig.suptitle('混合运动：相同IMM，DPC与CI对照；首种子，节点1，0–90s全部帧')
    fig.savefig(out/'xyz_mixed_comparison.png',dpi=155);plt.close(fig)
    fig=plt.figure(figsize=(16,10))
    for i,(kind,label) in enumerate(LABELS.items()):
        ax=fig.add_subplot(2,3,i+1,projection='3d')
        for method in ('dpc','ci'):
            name,color=METHODS[method];p=d[f'{kind}__{method}__state'][:,0,:3]
            ax.plot(*p.T,color=color,lw=1,label=name)
        ax.plot(*d[f'{kind}__truth_p'].T,'k--',lw=1.8,label='真值')
        ax.set_title(label);ax.set_xlabel('X (m)');ax.set_ylabel('Y (m)');ax.set_zlabel('Z (m)',labelpad=7)
        ax.set_box_aspect((1,1,.65));ax.view_init(23,-58)
        for axis in (ax.xaxis,ax.yaxis,ax.zaxis):
            axis.set_major_locator(MaxNLocator(4));axis.get_major_formatter().set_useOffset(False)
    fig.legend(handles=[Line2D([0],[0],color='black',ls='--',label='真值')]+[Line2D([0],[0],color=METHODS[m][1],label=METHODS[m][0]) for m in ('dpc','ci')],
               loc='lower center',ncol=3)
    fig.suptitle('3D轨迹对照：首种子、节点1、全部时间\n各轴独立缩放，非等比例几何')
    fig.subplots_adjust(left=.02,right=.95,bottom=.08,top=.9,hspace=.23,wspace=.24)
    fig.savefig(out/'fits_3d.png',dpi=170);plt.close(fig)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'outputs');args=ap.parse_args();out=args.output
    with (out/'per_seed.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
    config=json.loads((out/'config.json').read_text(encoding='utf-8'))
    assert len({r['seed'] for r in rows})==config['trials']
    groups=defaultdict(list)
    for r in rows:groups[r['trajectory'],r['method'],r['phase'],r['frequency_mhz']].append(r)
    summary=[]
    rms_keys=('pos_rmse','vel_rmse','future_pos_rmse','path_rmse','pX_rmse','pY_rmse','pZ_rmse','vX_rmse','vY_rmse','vZ_rmse','disagree_pos')
    for (kind,method,phase,frequency),rs in groups.items():
        mean,half=mean_ci([float(r['eta']) for r in rs])
        s=dict(trajectory=kind,method=method,phase=phase,frequency_mhz=float(frequency),trials=len(rs),eta_mean=mean,eta_ci95_half=half,
               eta_power_pooled=sum(float(r['rx_sum']) for r in rs)/sum(float(r['ideal_sum']) for r in rs))
        for key in rms_keys:s[key]=float(np.sqrt(np.mean([float(r[key])**2 for r in rs])))
        for key in ('z_bias','nees6','fusion_seconds'):s[key]=float(np.mean([float(r[key]) for r in rs]))
        summary.append(s)
    with (out/'summary.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    plots(out,summary)
    def result(kind,method,f=300.,phase='all'):
        return next(r for r in summary if (r['trajectory'],r['method'],r['frequency_mhz'],r['phase'])==(kind,method,f,phase))
    text=['# IMM：DPC与CI融合完整实验','', '## Material Passport',
          '- 阶段：已实现并运行的配对融合实验。',
          f"- {config['trials']}个种子，范围{config['seeds']}；5条轨迹；全部节点、全部帧保留。",
          '- 主表统计完整0–90s，不设置误差/功率达标门槛，不剔除起动、机动或异常样本。分段表仅用于解释。',
          '- 同一份IMM后验，位置/速度及模型参数不变；只有融合规则不同。',
          '', '## 定义与公平性',
          '- CI=Covariance Intersection：固定Metropolis权重的信息域凸共识，不是独立信息求和。',
          '- 每轮Y_i=sum_j w_ij P_j^-1，y_i=sum_j w_ij P_j^-1 x_j；输出P_i=Y_i^-1，x_i=P_i y_i。',
          '- 20节点环网、25轮、137ms起始延迟、40ms时隙。local也按同样发射时刻评价。',
          '- W^25是有限轮邻居交换的等效计算，不代替为全网集中均值。不回灌IMM。',
          '- 各分支预测均值及协方差，再保留模型间散布做矩匹配；先对omega边际化，再形成六维信息。',
          '- 未做CI权重logdet/trace优化。所有协方差及节点权重来自观测与滤波，不使用真值。',
          '', '## 全时段相干效率','',
          '|轨迹|DPC 300MHz|CI 300MHz|CI−DPC及95%配对区间|DPC 1GHz|CI 1GHz|',
          '|---|---:|---:|---:|---:|---:|']
    for kind in LABELS:
        a=sorted(groups[kind,'dpc','all','300.0'],key=lambda r:int(r['seed']))
        b=sorted(groups[kind,'ci','all','300.0'],key=lambda r:int(r['seed']))
        assert [r['seed'] for r in a]==[r['seed'] for r in b]
        delta,half=mean_ci([float(y['eta'])-float(x['eta']) for x,y in zip(a,b)])
        text.append(f"|{LABELS[kind]}|{result(kind,'dpc')['eta_mean']:.3f}|{result(kind,'ci')['eta_mean']:.3f}|{delta:+.3f} ± {half:.3f}|{result(kind,'dpc',1000.)['eta_mean']:.3f}|{result(kind,'ci',1000.)['eta_mean']:.3f}|")
    text+=['','区间反映固定轨迹下随机种子变化，未作多重比较校正，不是任意飞行任务的保证。',
           '', '## 全时段跟踪与差分路径误差','',
           '|轨迹|方法|三维位置m|三维速度m/s|高度Z m|垂直速度m/s|差分路径m|',
           '|---|---|---:|---:|---:|---:|---:|']
    for kind in LABELS:
        for method in METHODS:
            r=result(kind,method)
            text.append(f"|{LABELS[kind]}|{method}|{r['pos_rmse']:.3f}|{r['vel_rmse']:.3f}|{r['pZ_rmse']:.3f}|{r['vZ_rmse']:.3f}|{r['path_rmse']:.3f}|")
    text+=['','## 混合运动分段（所有该阶段帧，不按效果筛选）','',
           '|阶段|DPC位置m|CI位置m|DPC 300MHz η|CI 300MHz η|','|---|---:|---:|---:|---:|']
    for phase in ('acceleration','turn_entry','turn_steady','braking','hover_end'):
        a=result('mixed','dpc',phase=phase);b=result('mixed','ci',phase=phase)
        text.append(f"|{phase}|{a['pos_rmse']:.3f}|{b['pos_rmse']:.3f}|{a['eta_mean']:.3f}|{b['eta_mean']:.3f}|")
    text+=['','## 混合运动协方差诊断：六维NEES的均值','',
           '|阶段|DPC|CI|','|---|---:|---:|']
    for phase in ('all','acceleration','turn_steady','braking'):
        text.append(f"|{phase}|{result('mixed','dpc',phase=phase)['nees6']:.3f}|{result('mixed','ci',phase=phase)['nees6']:.3f}|")
    text+=['','NEES来自误差与输出协方差的比较，只用于解释，不决定保留哪些数据。机动段数值很大时，说明误差相对于估计不确定度偏大；CI不能保证修复这种输入模型失配。全程平均改善可能伴随机动局部退化，应同时读取上表。']
    text+=['','## 去掉启动5s的辅助对照（不是主结果）','',
           '|轨迹|DPC位置m|CI位置m|DPC Z m|CI Z m|DPC 300MHz η|CI 300MHz η|',
           '|---|---:|---:|---:|---:|---:|---:|']
    for kind in LABELS:
        a=result(kind,'dpc',phase='after_startup');b=result(kind,'ci',phase='after_startup')
        text.append(f"|{LABELS[kind]}|{a['pos_rmse']:.3f}|{b['pos_rmse']:.3f}|{a['pZ_rmse']:.3f}|{b['pZ_rmse']:.3f}|{a['eta_mean']:.3f}|{b['eta_mean']:.3f}|")
    text+=['','## 边界与可解释性',
           '- 本次混合运动的起飞加速和刹停阶段，CI相干效率低于DPC；刹停阶段位置误差也更大。全程平均改善不能掩盖这些局部退化。',
           '- 起飞加速阶段CI三维位置RMSE虽下降，相干效率却下降，说明绝对位置精度与跨节点差分路径精度并不等价。',
           '- CI的保守性以输入协方差可信为前提；IMM/EKF与坐标转换是近似，不宣称严格一致性。NEES仅保存为诊断量，无门控。',
           '- 同一后验的不确定度存在方向差异时，信息融合可以优于坐标等权平均；输入协方差相同时，CI均值退化为DPC均值。',
           '- CI不能凭空补充未观测的速度，也不自动消除所有节点共有的机动滞后、偏置或错误运动模型。',
           '- 固定相同时延隔离算法影响。CI需要交换协方差，消息量及硬件处理开销未实测；实际带宽受限时可能改变结果。',
           '- 参数未重新调优；模型窗内冻结、固定姿态、独立白导航误差、理想载波同步及自由空间假设与旧实验一致。',
           '- 只增加旧IMM的可选分支协方差返回值，没有改变滤波递推。首种子DPC状态与两频率功率曲线已与旧数据核对。',
           '', '## 来源',
           '- Julier & Uhlmann (1997): https://doi.org/10.1109/ACC.1997.609105 （CI原始出处）',
           '- MathWorks CI: https://www.mathworks.com/help/fusion/ref/fusecovint.html （官方实现为优化权重，区别于本固定权重实验）',
           '', '![完整对照](outputs/ci_comparison.png)', '', '![Z轴](outputs/z_comparison.png)', '', '![3D轨迹](outputs/fits_3d.png)']
    (HERE/'RESULTS.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print('Saved summary.csv, full report, four comparison figures.')

if __name__=='__main__':main()
