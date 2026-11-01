"""Read complete formal results, paired summaries and scientific figures."""
from pathlib import Path
from collections import defaultdict
import csv,json,hashlib
import numpy as np
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import ScalarFormatter,MaxNLocator
HERE=Path(__file__).resolve().parent
METHODS={'IMM_CI_linear':('旧IMM+CI','#929292'),'GP_independent_nav':('GP+独立导航近似','#d7a243'),'JointGP_linear':('联合GP+直线浮标预测','#668fc2'),'JointGP_accel':('联合GP+惯性预测','#bf4859'),'JointGP_power':('联合GP+期望功率优化','#44886c')}
KINDS={'hover':'悬停','straight':'直线','circle':'圆周','mixed':'混合','climb':'爬升'}
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':9,'axes.spines.top':False,'axes.spines.right':False})


def read(path):
    with path.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))


def write(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def mean_interval(values):
    a=np.asarray(values,float)
    return float(a.mean()),float(student_t.ppf(.975,len(a)-1)*a.std(ddof=1)/np.sqrt(len(a)))


def figures(out,summary):
    def find(hs,kind,method,frequency=300.,timing='fixed_delay',phase='all'):
        return next(r for r in summary if (r['hs'],r['trajectory'],r['method'],r['frequency_mhz'],r['timing'],r['phase'])==(hs,kind,method,frequency,timing,phase))
    fig,axes=plt.subplots(2,3,figsize=(16,8),layout='constrained');xx=np.arange(5)
    selected=('IMM_CI_linear','GP_independent_nav','JointGP_accel','JointGP_power')
    for row,hs in enumerate((.5,1.5)):
        for col,(metric,title) in enumerate((('pos_rmse','位置RMSE (m)'),('position_error_rate','位置误差变化率RMS (m/s)'),('eta','300MHz相干效率，统一功率参考'))):
            ax=axes[row,col]
            for j,method in enumerate(selected):
                label,color=METHODS[method];ax.bar(xx+(j-1.5)*.2,[find(hs,k,method)[metric] for k in KINDS],width=.19,label=label,color=color)
            ax.set_xticks(xx,list(KINDS.values()));ax.set_title(f'Hs={hs:g}m · {title}');ax.grid(axis='y',alpha=.2)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('完整20种子：0–90s全部帧，固定时延对照；GP直接估计轨迹，无输出端平滑')
    fig.savefig(out/'overview.png',dpi=165);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(14,8),layout='constrained')
    for row,hs in enumerate((.5,1.5)):
        for col,freq in enumerate((300.,1000.)):
            ax=axes[row,col]
            for j,(method,timing,style) in enumerate((('IMM_CI_linear','fixed_delay',''),('JointGP_accel','fixed_delay',''),('JointGP_accel','measured_cpu','//'),('JointGP_power','measured_cpu','//'))):
                label,color=METHODS[method];label+=' · '+('固定延迟' if timing=='fixed_delay' else '计入求解耗时')
                ax.bar(xx+(j-1.5)*.2,[find(hs,k,method,freq,timing)['eta'] for k in KINDS],width=.19,label=label,color=color,hatch=style,edgecolor='white')
            ax.set_title(f'Hs={hs:g}m · {freq:g}MHz');ax.set_xticks(xx,list(KINDS.values()));ax.grid(axis='y',alpha=.2)
    axes[0,0].legend(fontsize=7)
    fig.suptitle('相同固定总功率1W下的实际接收功率/均匀功率理想参考；斜线含测得求解耗时')
    fig.savefig(out/'coherent_timing.png',dpi=170);plt.close(fig)
    d=np.load(out/'first_seed_hs1.5.npz');t=d['time_s'];kind='mixed'
    truth=np.concatenate([d[f'{kind}__truth_p'],d[f'{kind}__truth_v']],-1)
    fig,axes=plt.subplots(6,2,figsize=(14,15),sharex=True,layout='constrained')
    for j in range(6):
        axes[j,0].plot(t,truth[:,j],'k--',label='真值',lw=1.3)
        for method in ('IMM_CI_linear','GP_independent_nav','JointGP_accel'):
            label,color=METHODS[method];value=d[f'{kind}__{method}__state'][:,0,j]
            axes[j,0].plot(t,value,color=color,label=label,lw=.9);axes[j,1].plot(t,value-truth[:,j],color=color,lw=.9)
        sig=d[f'{kind}__JointGP_accel__sigma'][:,0,j]
        axes[j,1].fill_between(t,-1.96*sig,1.96*sig,color=METHODS['JointGP_accel'][1],alpha=.12,label='联合GP名义±1.96σ')
        axes[j,1].axhline(0,color='black',ls='--',lw=.7)
        name=('位置' if j<3 else '速度')+'XYZ'[j%3];unit='m' if j<3 else 'm/s'
        axes[j,0].set_ylabel(f'{name} ({unit})');axes[j,1].set_ylabel(f'误差 ({unit})')
        for ax in axes[j]:ax.grid(alpha=.2);ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    axes[0,0].legend(fontsize=8);axes[0,1].legend(fontsize=8)
    for ax in axes[-1]:ax.set_xlabel('时间 (s)')
    fig.suptitle('混合运动，Hs=1.5m，首种子节点1：当前时刻的因果输出，未事后平滑')
    fig.savefig(out/'xyz_mixed.png',dpi=150);plt.close(fig)
    fig=plt.figure(figsize=(20,10))
    chosen=('IMM_CI_linear','GP_independent_nav','JointGP_accel')
    for row,hs in enumerate((.5,1.5)):
        data=np.load(out/f'first_seed_hs{hs:g}.npz')
        for col,(kind,label) in enumerate(KINDS.items()):
            ax=fig.add_subplot(2,5,row*5+col+1,projection='3d')
            for method in chosen:ax.plot(*data[f'{kind}__{method}__state'][:,0,:3].T,color=METHODS[method][1],lw=.9)
            truthp=data[f'{kind}__truth_p'];ax.plot(*truthp.T,'k--',lw=1.3)
            if kind=='hover':ax.scatter(*truthp[0],color='black',s=24,marker='x')
            ax.set_title(f'Hs={hs:g}m · {label}');ax.set_xlabel('X (m)');ax.set_ylabel('Y (m)');ax.set_zlabel('Z (m)');ax.view_init(23,-58);ax.set_box_aspect((1,1,.8))
            for axis in (ax.xaxis,ax.yaxis,ax.zaxis):axis.set_major_formatter(ScalarFormatter(useOffset=False));axis.set_major_locator(MaxNLocator(3))
    handles=[Line2D([0],[0],color='black',ls='--',label='真值')]+[Line2D([0],[0],color=METHODS[m][1],label=METHODS[m][0]) for m in chosen]
    fig.legend(handles=handles,loc='lower center',ncol=4);fig.suptitle('首种子、完整时间的3D轨迹；各轴独立缩放，非等比例几何')
    fig.subplots_adjust(left=.02,right=.96,bottom=.09,top=.9,wspace=.3,hspace=.25)
    fig.savefig(out/'trajectories_3d.png',dpi=165);plt.close(fig)
    fig,axes=plt.subplots(2,5,figsize=(19,7),sharex=True,layout='constrained')
    for row,hs in enumerate((.5,1.5)):
        data=np.load(out/f'first_seed_hs{hs:g}.npz');t=data['time_s']
        for col,(kind,name) in enumerate(KINDS.items()):
            ax=axes[row,col]
            for method in ('IMM_CI_linear','JointGP_accel','JointGP_power'):
                label,color=METHODS[method];ax.plot(t,data[f'{kind}__{method}__eta300'],color=color,lw=.8,label=label)
            ax.set_title(f'Hs={hs:g}m · {name}');ax.grid(alpha=.2)
            if row:ax.set_xlabel('时间 (s)')
            if col==0:ax.set_ylabel('300MHz统一参考相干效率')
    axes[0,0].legend(fontsize=7);fig.suptitle('首种子：完整固定时延发射曲线，未筛选或平滑；功率优化可能集中功率')
    fig.savefig(out/'coherent_curves.png',dpi=160);plt.close(fig)


def main():
    out=HERE/'outputs';manifest=json.loads((out/'manifest.json').read_text(encoding='utf-8'));rows=read(out/'per_seed.csv')
    if manifest['trials']!=20 or len(rows)!=32000:raise ValueError('Report requires entire 20-seed formal experiment')
    dimensions=('hs','trajectory','method','frequency_mhz','phase','timing');groups=defaultdict(list)
    numeric=[k for k in rows[0] if k not in (*dimensions,'seed')]
    if not np.isfinite(np.array([[r[k] for k in numeric] for r in rows],float)).all():raise ValueError('Nonfinite formal results; no samples removed')
    for r in rows:groups[float(r['hs']),r['trajectory'],r['method'],float(r['frequency_mhz']),r['phase'],r['timing']].append(r)
    summary=[]
    for key,rr in groups.items():
        if len(rr)!=20 or {int(r['seed']) for r in rr}!=set(manifest['seeds']):raise ValueError('Missing or duplicate paired seed')
        item=dict(zip(dimensions,key));item['seeds']=20
        for metric in numeric:item[metric],item[metric+'_ci95_half']=mean_interval([r[metric] for r in rr])
        summary.append(item)
    write(out/'summary.csv',summary)
    pairs=[]
    for hs in (.5,1.5):
        for kind in KINDS:
            for freq in (300.,1000.):
                for timing in ('fixed_delay','measured_cpu'):
                    for phase in ('all','after_startup','acceleration','braking'):
                        for left,right in (('IMM_CI_linear','JointGP_accel'),('GP_independent_nav','JointGP_accel'),('JointGP_linear','JointGP_accel'),('JointGP_accel','JointGP_power')):
                            a={int(r['seed']):r for r in groups[hs,kind,left,freq,phase,timing]};b={int(r['seed']):r for r in groups[hs,kind,right,freq,phase,timing]}
                            item=dict(hs=hs,trajectory=kind,frequency_mhz=freq,timing=timing,phase=phase,left=left,right=right)
                            for metric in ('pos_rmse','future_pos_rmse','position_error_rate','velocity_error_rate','eta','gain_db','rx_mean_w','coverage95_position'):
                                item[metric+'_difference'],item[metric+'_ci95_half']=mean_interval([float(b[s][metric])-float(a[s][metric]) for s in sorted(a)])
                            pairs.append(item)
    write(out/'paired_comparisons.csv',pairs)
    figures(out,summary)
    lines=['# 文献驱动的联合GP轨迹估计：完整实验结果','','## Material Passport','- 20种子、两海况、五轨迹、五方法、两载频、两种时延评价，正式32000行种子级分段指标。','- 全时段0–90s为主结果；无输出后处理平滑、门控或样本剔除。','- 同样原始伪距/IMU/磁场/雷达、同样物理真值与随机数；杆臂为零、20单阵元浮标、固定总功率1W。','- 新方法为集中有限历史联合MAP，不是已经完成的有限轮分布式GP。','- 完整公式、正式文献、计数和适用边界见RESEARCH.md。','', '## 全时段固定时延对照','', '|Hs m|轨迹|方法|位置RMSE m|速度RMSE m/s|位置误差变化率 m/s|未来位置RMSE m|η300MHz|η1GHz|位置95%名义区域覆盖率|','|---|---|---|---:|---:|---:|---:|---:|---:|---:|']
    def find(hs,kind,m,f=300.,phase='all',timing='fixed_delay'):
        return next(r for r in summary if (r['hs'],r['trajectory'],r['method'],r['frequency_mhz'],r['phase'],r['timing'])==(hs,kind,m,f,phase,timing))
    for hs in (.5,1.5):
        for kind in KINDS:
            for m in METHODS:
                r=find(hs,kind,m);r1=find(hs,kind,m,1000.)
                lines.append(f'|{hs:g}|{KINDS[kind]}|{METHODS[m][0]}|{r["pos_rmse"]:.4f}|{r["vel_rmse"]:.4f}|{r["position_error_rate"]:.4f}|{r["future_pos_rmse"]:.4f}|{r["eta"]:.4f}|{r1["eta"]:.4f}|{r["coverage95_position"]:.3f}|')
    lines+=['','## 求解耗时与计入耗时后的相干结果','', '|Hs m|轨迹|方法|计算均值ms|计算P95 ms|时隙末端平均ms|η300MHz|G300MHz dB|η1GHz|G1GHz dB|有效节点数300MHz|','|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for hs in (.5,1.5):
        for kind in KINDS:
            for m in ('IMM_CI_linear','JointGP_accel','JointGP_power'):
                r=find(hs,kind,m,timing='measured_cpu');r1=find(hs,kind,m,1000.,timing='measured_cpu')
                lines.append(f'|{hs:g}|{KINDS[kind]}|{METHODS[m][0]}|{r["compute_mean_ms"]:.2f}|{r["compute_p95_ms"]:.2f}|{r["tx_end_mean_ms"]:.2f}|{r["eta"]:.4f}|{r["gain_db"]:.3f}|{r1["eta"]:.4f}|{r1["gain_db"]:.3f}|{r["effective_nodes"]:.2f}|')
    lines+=['','计算时间为当前电脑并行负载下的求解加导航平均耗时，未计全预测、功率优化、通信编解码和真实带宽开销；不是实时硬件保证。若时隙末端超过200ms，则不满足原先帧内结束整个时隙的约束，未剔除这些样本。','', '## 混合机动的局部置信度','', '|Hs m|阶段|方法|位置RMSE m|位置95%名义区域覆盖率|六维NEES|','|---|---|---|---:|---:|---:|']
    for hs in (.5,1.5):
        for phase in ('all','acceleration','turn_steady','braking'):
            for m in ('IMM_CI_linear','GP_independent_nav','JointGP_accel'):
                r=find(hs,'mixed',m,phase=phase)
                lines.append(f'|{hs:g}|{phase}|{METHODS[m][0]}|{r["pos_rmse"]:.4f}|{r["coverage95_position"]:.3f}|{r["nees6"]:.3f}|')
    maxdiff=max(d['baseline_max_state_difference'] for d in manifest['baseline_checks'])
    lines+=['','## 解释和复现边界',f'- 基线首种子两海况五轨迹与原完整数据的最大状态差：{maxdiff:.3g}。','- RMSE等列为先按种子计算、再跨种子平均；配对CSV给20种子差值的Student-t 95%区间，未作多重比较校正。时序/节点不作为独立显著性样本。','- 位置误差变化率针对估计减真值的变化，不是单纯奖励小jerk；完整启动和机动均保留。','- η主列统一使用均匀功率理想参考；功率重分配可能使该η超过1。另存eta_optimal_bound为任意分配最优理想上界归一化，eta_own只作自身分配同相效率。','- 主特征向量仅在条件高斯/预测幅度/固定总功率下最大化预测期望接收功率；不保证实际全局最优、不保证所有浮标保持高相干，必须同时看G及有效节点数。','- 波浪分解不再由真实分量提供，也没有再次对导航后验进行独立噪声假设的波浪滤波；惯性预测只使用当前传感器和联合导航修正。','- 导航高斯压缩、一阶误差坐标、三次Gauss–Newton、弱始端先验与有限2s窗口均为近似，经验置信度以数据为准。窗口外旧雷达信息被舍弃，上轮解只用于初始化，不当作重复观测先验。','- 集中观测也可能带来收益，GP_independent_nav用于分离这一贡献；不能把所有改善归因于联合相关性或WNOJ先验。','- 未作独立留出参数标定，不宣称任意海况普适；物理和传感器简化与上一版相同，RF载波同步仍理想。','', '![总体对照](outputs/overview.png)','![计算时延与相干](outputs/coherent_timing.png)','![三轴坐标速度](outputs/xyz_mixed.png)','![三维轨迹](outputs/trajectories_3d.png)','![相干曲线](outputs/coherent_curves.png)']
    (HERE/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    check=dict(status='complete',records=len(rows),seeds=20,source_matches={p:hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in manifest['source_sha256'].items()},baseline_max_state_difference=maxdiff,figures=[p.name for p in out.glob('*.png')])
    (out/'report_audit.json').write_text(json.dumps(check,indent=2),encoding='utf-8');print(json.dumps(check,indent=2))


if __name__=='__main__':main()
