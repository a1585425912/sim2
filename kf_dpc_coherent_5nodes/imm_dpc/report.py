"""Recompute summaries and plots from saved paired trials; no simulation."""
from pathlib import Path
import argparse, csv, json
from collections import defaultdict
import numpy as np
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
LABELS={'hover':'悬停','straight':'直线','circle':'圆周','mixed':'混合运动'}
COLORS={'S':'#25836c','CV':'#4782c2','CT':'#d47e29'}
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],
                     'axes.unicode_minus':False,'font.size':10,'axes.spines.top':False,
                     'axes.spines.right':False})


def interval(x):
    x=np.asarray(x); half=(student_t.ppf(.975,len(x)-1)*x.std(ddof=1)/np.sqrt(len(x))) if len(x)>1 else np.nan
    return float(x.mean()),float(half)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=HERE/'outputs');args=ap.parse_args()
    out=args.output
    with (out/'per_seed.csv').open(encoding='utf-8-sig') as f: rows=list(csv.DictReader(f))
    config=json.loads((out/'config.json').read_text(encoding='utf-8'))
    groups=defaultdict(list)
    for r in rows:groups[(r['trajectory'],r['method'],r['phase'],r['frequency_mhz'])].append(r)
    summary=[]
    for (kind,method,phase,f),rs in groups.items():
        mean,ci=interval([float(r['eta']) for r in rs])
        row=dict(trajectory=kind,method=method,phase=phase,frequency_mhz=f,trials=len(rs),eta_mean=mean,eta_ci95_half=ci,
                 eta_power_pooled=sum(float(r['rx_sum']) for r in rs)/sum(float(r['ideal_sum']) for r in rs))
        for metric in ('pos_rmse','vel_rmse','future_pos_rmse','path_rmse'):
            row[metric]=float(np.sqrt(np.mean([float(r[metric])**2 for r in rs])))
        if method=='IMM':
            for key in ('prob_S','prob_CV','prob_CT'):row[key]=float(np.mean([float(r[key]) for r in rs]))
        else:row.update(prob_S='',prob_CV='',prob_CT='')
        summary.append(row)
    with (out/'summary.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    def result(kind,method,f=300.,phase='after_startup'):
        return next(r for r in summary if (r['trajectory'],r['method'],float(r['frequency_mhz']),r['phase'])==(kind,method,f,phase))
    d=np.load(out/'first_seed.npz');time=d['time_s']
    fig,axes=plt.subplots(4,3,figsize=(15,12),sharex=True,layout='constrained')
    for i,kind in enumerate(LABELS):
        prob=d[f'{kind}__probabilities'].mean(axis=1)
        for j,name in enumerate(('S','CV','CT')):axes[i,0].plot(time,prob[:,j],label=name,color=COLORS[name])
        axes[i,0].set_ylim(0,1);axes[i,0].set_ylabel(LABELS[kind]+'\n模型概率（节点均值）')
        truth=np.linalg.norm(d[f'{kind}__truth_v'],axis=-1)
        axes[i,1].plot(time,truth,color='black',ls='--',label='真实速度')
        for method,color in (('CV_original','#4782c2'),('IMM','#b34b53')):
            axes[i,1].plot(time,np.linalg.norm(d[f'{kind}__{method}__state'][:,0,3:6],axis=-1),color=color,label=method,lw=1)
            axes[i,2].plot(time,d[f'{kind}__{method}__eta300'],color=color,label=method,lw=.9)
        axes[i,1].set_ylabel('节点1速度模长 m/s');axes[i,2].set_ylim(0,1);axes[i,2].set_ylabel('300 MHz 相干效率')
        for ax in axes[i]:
            ax.grid(alpha=.2)
            if kind=='mixed':
                for boundary in (15,20,35,55,65,70):ax.axvline(boundary,color='gray',ls=':',lw=.6)
    for ax in axes[0]:ax.legend(fontsize=8,loc='upper right')
    for ax in axes[-1]:ax.set_xlabel('时间 s')
    fig.suptitle('IMM如何工作：首个种子，概率与原始相干效率曲线（无平滑）')
    fig.savefig(out/'imm_explained.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,4.6),layout='constrained')
    kinds=list(LABELS);x=np.arange(4)
    for ax,f in zip(axes,(300.,1000.)):
        for j,method in enumerate(('CV_original','CT','IMM_SC','IMM')):
            vals=[result(k,method,f) for k in kinds]
            ax.bar(x+(j-1.5)*.19,[r['eta_mean'] for r in vals],width=.18,
                   yerr=[r['eta_ci95_half'] for r in vals],capsize=2,label=method)
        ax.set_xticks(x,[LABELS[k] for k in kinds]);ax.set_ylim(0,1);ax.set_title(f'{f:.0f} MHz');ax.set_ylabel('相干效率：种子均值与95% t区间')
    axes[0].legend(fontsize=8);fig.savefig(out/'comparison.png',dpi=170);plt.close(fig)
    text=['# IMM首版实验报告','',
          '## Material Passport','- 阶段：首版探索性实现与配对实验。',
          '- 输入：固定的GM浮标、nominal导航、位置观测；无多普勒、真值标签或未来测量。',
          f"- 输出：{config['trials']}个配对种子，原始CSV、首种子NPZ和图。",
          '- 判定：已运行；参数未经独立校准，不能声称泛化或全局最优。','',
          '## 相干效率（t≥5s，先按每种子求η，再平均）','',
          '|轨迹|CV 300MHz|IMM 300MHz|配对提升及95%区间|CV 1GHz|IMM 1GHz|',
          '|---|---:|---:|---:|---:|---:|']
    for kind in kinds:
        aa=sorted(groups[kind,'CV_original','after_startup','300.0'],key=lambda r:int(r['seed']))
        bb=sorted(groups[kind,'IMM','after_startup','300.0'],key=lambda r:int(r['seed']))
        assert [a['seed'] for a in aa]==[b['seed'] for b in bb]
        delta,ci=interval([float(b['eta'])-float(a['eta']) for a,b in zip(aa,bb)])
        text.append(f"|{LABELS[kind]}|{result(kind,'CV_original')['eta_mean']:.3f}|{result(kind,'IMM')['eta_mean']:.3f}|{delta:+.3f} ± {ci:.3f}|{result(kind,'CV_original',1000.)['eta_mean']:.3f}|{result(kind,'IMM',1000.)['eta_mean']:.3f}|")
    text+=['','区间只反映固定轨迹与噪声模型下的种子变化；未作多重比较校正，不代表所有飞行任务的置信范围。',
           '','## 位置、速度和差分路径误差（300MHz对应行，跨种子合并平方误差）','',
           '|轨迹|方法|位置RMSE m|速度RMSE m/s|未来位置RMSE m|差分路径RMSE m|',
           '|---|---|---:|---:|---:|---:|']
    for kind in kinds:
        for method in ('CV_original','CT','IMM_SC','IMM'):
            r=result(kind,method)
            text.append(f"|{LABELS[kind]}|{method}|{r['pos_rmse']:.3f}|{r['vel_rmse']:.3f}|{r['future_pos_rmse']:.3f}|{r['path_rmse']:.3f}|")
    text+=['','## 模型权重（各节点、帧与种子平均）','','|轨迹|近悬停S|CV|CT|','|---|---:|---:|---:|']
    for kind in kinds:
        r=result(kind,'IMM');text.append(f"|{LABELS[kind]}|{r['prob_S']:.3f}|{r['prob_CV']:.3f}|{r['prob_CT']:.3f}|")
    text+=['','## 混合轨迹分段：300MHz η','','|阶段|CV|IMM|','|---|---:|---:|']
    for phase in ('hover_start','acceleration','straight_segment','turn_entry','turn_steady','braking','hover_end'):
        text.append(f"|{phase}|{result('mixed','CV_original',phase=phase)['eta_mean']:.3f}|{result('mixed','IMM',phase=phase)['eta_mean']:.3f}|")
    text+=['','## 如何读结果与限制','',
           '- S把近零速度作为软先验；CV覆盖直线；CT从位置历史估计未知转弯率。所有分支共用七维状态，协方差混合含模型间分歧。',
           '- 每帧按预测创新及其协方差更新概率。CT的零转弯率极限也能描述直线，因此概率不是分类准确率。',
           '- CT过程噪声强度比CV更小，所以单CT/IMM的改善同时可能包含模型和噪声参数收益，不能全部归于“识别了转弯”。需后续同q消融与独立校准。',
           '- IMM_SC对照揭示仅加入悬停分支是否足够；CT单模型揭示纯圆周可能更适合单模型，IMM不是每种情况下都最好。',
           '- CV_exact只改变逐时隙几何计算；与CV_original的差异隔离预测接口收益。完整IMM先逐分支预测、概率加权，再做25轮DPC；不回灌本地滤波。',
           '- 预测窗≤177ms内冻结模式，是近似；模式切换仍有检测延迟，不能预知起飞、刹停或转弯。混合轨迹速度连续，切向加速度有阶段变化。',
           '- 相同12ms处理与5ms/轮属于链路假设。9时隙预测提高消息负载，额外运算、真实硬件时延与带宽尚未闭环验证。filter_seconds仅为本机20节点批处理的滤波耗时。',
           '- 载波频相已同步、导航独立白误差、固定姿态、自由空间；不涵盖多径及复杂三维机动。',
           '- 参数固定后直接做首版探索：没有独立校准集和未见轨迹测试，不据此宣称最终性能。原始CV文件与结果未覆盖。',
           '', '## 理论来源',
           '- Blom & Bar-Shalom (1988), https://doi.org/10.1109/9.1299 （IMM原始出处；本次未取得IEEE正文）',
           '- Yuan et al. (2014), https://doi.org/10.1155/2014/649276 （未知转弯率CT模型；不是复现其EM算法）',
           '- Huang et al. (2017), https://isif.org/files/isif/2024-01/JAIF_Vol12_1_complete_170828_7.pdf （IMM机动跟踪应用，传感器不同）',
           '', '![IMM解释](outputs/imm_explained.png)','', '![对照](outputs/comparison.png)']
    (HERE/'RESULTS.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    print('Saved summary.csv, figures and RESULTS.md')

if __name__=='__main__':main()
