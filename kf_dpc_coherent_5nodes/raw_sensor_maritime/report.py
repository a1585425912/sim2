"""Describe the entire full-run dataset; no rerunning or sample selection."""
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
METHODS={'CV_CI_wave':('CV+CI+波浪预测','#999999'),'IMM_DPC_linear':('IMM+DPC+直线预测','#6b91bb'),'IMM_CI_linear':('IMM+CI+直线预测','#dc9a44'),'IMM_CI_wave':('IMM+CI+波浪预测','#bd4454')}
KINDS={'hover':'悬停','straight':'直线','circle':'圆周','mixed':'混合运动','climb':'爬升'}
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':9,'axes.spines.top':False,'axes.spines.right':False})


def read(path):
    with path.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))


def write(path,rows):
    with path.open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def mean_interval(values):
    a=np.asarray(values,float);half=student_t.ppf(.975,len(a)-1)*a.std(ddof=1)/np.sqrt(len(a))
    return float(a.mean()),float(half)


def add_power_reference(rows,manifest,out):
    """Reconstruct evaluation-only channel reference; never rerun estimation."""
    from physics import Config,SeaTruth
    from run_full import trajectory
    groups=defaultdict(list)
    for r in rows:groups[int(r['seed']),float(r['hs'])].append(r)
    for (seed,hs),rr in groups.items():
        cfg=Config(**next(c for c in manifest['configs'] if c['hs']==hs))
        t=np.arange(round(cfg.duration/cfg.dt)+1)*cfg.dt;h=cfg.horizons()
        buoy=SeaTruth(cfg,seed).at(t[:,None]+h)['p']
        phases=dict(all=np.ones(len(t),bool),startup=t<5,after_startup=t>=5,acceleration=(t>=15)&(t<20),turn_entry=(t>=35)&(t<40),turn_steady=(t>=40)&(t<55),braking=(t>=65)&(t<70),hover_end=(t>=75)&(t<90))
        reference={}
        for kind in KINDS:
            future,_=trajectory(t[:,None]+h,kind,cfg)
            d=np.linalg.norm(future[:,:,None]-buoy,axis=-1)
            for f in cfg.frequencies_mhz:
                lam=299792458./(f*1e6);a=np.sqrt(cfg.total_power/cfg.nodes)*lam/(4*np.pi*d)
                inc=np.sum(a*a,axis=-1)
                for phase,mask in phases.items():reference[kind,float(f),phase]=float(inc[mask].mean())
        for r in rr:
            inc=reference[r['trajectory'],float(r['frequency_mhz']),r['phase']]
            r['incoherent_mean_w']=inc;r['gain_linear']=float(r['rx_mean_w'])/inc;r['gain_db']=10*np.log10(r['gain_linear'])
    write(out/'per_seed_with_gain.csv',rows)


def plots(out,summary):
    def get(hs,kind,method,f=300,phase='all'):
        return next(r for r in summary if (r['hs'],r['trajectory'],r['method'],r['frequency_mhz'],r['phase'])==(hs,kind,method,f,phase))
    fig,axes=plt.subplots(2,3,figsize=(16,8),layout='constrained');x=np.arange(5)
    for row,hs in enumerate((.5,1.5)):
        for col,(metric,f,title) in enumerate((('pos_rmse',300,'位置RMSE (m)'),('eta',300,'300MHz 相干效率'),('eta',1000,'1GHz 相干效率'))):
            ax=axes[row,col]
            for j,(method,(name,color)) in enumerate(METHODS.items()):ax.bar(x+(j-1.5)*.2,[get(hs,k,method,f)[metric] for k in KINDS],width=.19,color=color,label=name)
            ax.set_xticks(x,list(KINDS.values()));ax.set_title(f'Hs={hs:g}m · {title}');ax.grid(axis='y',alpha=.2)
            if metric=='eta':ax.set_ylim(0,1)
    axes[0,0].legend(fontsize=8)
    fig.suptitle('原始传感器驱动完整实验：20种子，20浮标，五条UAV轨迹，全部0–90s')
    fig.savefig(out/'overview.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,5,figsize=(19,7),sharex=True,sharey=True,layout='constrained')
    for row,hs in enumerate((.5,1.5)):
        d=np.load(out/f'first_seed_hs{hs:g}.npz');t=d['time_s']
        for col,(kind,name) in enumerate(KINDS.items()):
            ax=axes[row,col]
            for method,(label,color) in METHODS.items():ax.plot(t,d[f'{kind}__{method}__eta300'],color=color,lw=.8,label=label)
            ax.set_title(f'Hs={hs:g}m · {name}');ax.set_ylim(0,1);ax.grid(alpha=.2)
            if col==0:ax.set_ylabel('300MHz 相干效率 η')
            if row==1:ax.set_xlabel('时间 (s)')
    axes[0,0].legend(fontsize=7)
    fig.suptitle('首种子20261028：实际发射时隙相干效率，无平滑、无样本筛选')
    fig.savefig(out/'coherent_curves.png',dpi=160);plt.close(fig)
    fig=plt.figure(figsize=(20,10))
    handles=[Line2D([0],[0],color='black',ls='--',label='真值')]
    for method in ('CV_CI_wave','IMM_DPC_linear','IMM_CI_wave'):
        name,color=METHODS[method];handles.append(Line2D([0],[0],color=color,label=name))
    for row,hs in enumerate((.5,1.5)):
        d=np.load(out/f'first_seed_hs{hs:g}.npz')
        for col,(kind,label) in enumerate(KINDS.items()):
            ax=fig.add_subplot(2,5,row*5+col+1,projection='3d')
            for method in ('CV_CI_wave','IMM_DPC_linear','IMM_CI_wave'):
                p=d[f'{kind}__{method}__state'][:,0,:3];ax.plot(*p.T,color=METHODS[method][1],lw=.8)
            ax.plot(*d[f'{kind}__truth_p'].T,'k--',lw=1.4)
            ax.set_title(f'Hs={hs:g}m · {label}');ax.set_xlabel('X (m)');ax.set_ylabel('Y (m)');ax.set_zlabel('Z (m)')
            ax.set_box_aspect((1,1,.8));ax.view_init(23,-58)
            for axis in (ax.xaxis,ax.yaxis,ax.zaxis):axis.set_major_locator(MaxNLocator(3));axis.set_major_formatter(ScalarFormatter(useOffset=False))
    fig.legend(handles=handles,loc='lower center',ncol=4)
    fig.suptitle('首种子、节点1、全部时间的3D轨迹；各轴独立缩放，非等比例几何')
    fig.subplots_adjust(left=.02,right=.96,bottom=.09,top=.9,wspace=.3,hspace=.25)
    fig.savefig(out/'trajectories_3d.png',dpi=160);plt.close(fig)
    d=np.load(out/'first_seed_hs1.5.npz');t=d['time_s'];kind='mixed'
    truth=np.concatenate([d[f'{kind}__truth_p'],d[f'{kind}__truth_v']],axis=-1)
    fig,axes=plt.subplots(6,2,figsize=(14,15),sharex=True,layout='constrained')
    for j in range(6):
        axes[j,0].plot(t,truth[:,j],'k--',label='真值',lw=1.4)
        for method in ('CV_CI_wave','IMM_DPC_linear','IMM_CI_wave'):
            name,color=METHODS[method];value=d[f'{kind}__{method}__state'][:,0,j]
            axes[j,0].plot(t,value,color=color,label=name,lw=.9)
            axes[j,1].plot(t,value-truth[:,j],color=color,lw=.9)
        sigma=d[f'{kind}__IMM_CI_wave__sigma'][:,0,j]
        axes[j,1].fill_between(t,-1.96*sigma,1.96*sigma,color=METHODS['IMM_CI_wave'][1],alpha=.14,label='IMM+CI名义±1.96σ')
        axes[j,1].axhline(0,color='black',ls='--',lw=.7)
        unit='m' if j<3 else 'm/s';name=('位置' if j<3 else '速度')+'XYZ'[j%3]
        axes[j,0].set_ylabel(f'{name} ({unit})');axes[j,1].set_ylabel(f'误差 ({unit})')
        for ax in axes[j]:ax.grid(alpha=.2);ax.yaxis.set_major_formatter(ScalarFormatter(useOffset=False))
    axes[0,0].legend(fontsize=8);axes[0,1].legend(fontsize=8)
    axes[0,0].set_title('真值与估计');axes[0,1].set_title('估计减真值；阴影为名义边际区间，并非已校准保证')
    for ax in axes[-1]:ax.set_xlabel('时间 (s)')
    fig.suptitle('Hs=1.5m，混合运动，首种子节点1，完整三轴坐标与速度')
    fig.savefig(out/'xyz_mixed.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(15,8),sharex=True,layout='constrained')
    for row,hs in enumerate((.5,1.5)):
        nav=np.load(out/f'navigation_hs{hs:g}.npz');tr=np.load(out/f'truth_navigation_eval_only_hs{hs:g}.npz');t=tr['time']
        ep=nav['p'][:,0]-tr['p'][:,0];ev=nav['v'][:,0]-tr['v'][:,0]
        for j,axis in enumerate('XYZ'):
            axes[row,0].plot(t,ep[:,j],label=axis,lw=.9);axes[row,1].plot(t,ev[:,j],label=axis,lw=.9)
        lin=np.sqrt(np.sum((nav['buoy_linear'][:,:,0]-tr['future_p'][:,:,0])**2,-1)).mean(1)
        wav=np.sqrt(np.sum((nav['buoy_wave'][:,:,0]-tr['future_p'][:,:,0])**2,-1)).mean(1)
        axes[row,2].plot(t,lin,label='直线预测',lw=.9);axes[row,2].plot(t,wav,label='波浪预测',lw=.9)
        for col,title in enumerate(('位置误差 (m)','速度误差 (m/s)','未来浮标位置误差 (m)')):
            axes[row,col].set_title(f'Hs={hs:g}m · {title}');axes[row,col].grid(alpha=.2);axes[row,col].legend()
            if row==1:axes[row,col].set_xlabel('时间 (s)')
    fig.suptitle('浮标1：从伪距、IMU、磁场观测自行解算；未来误差在发射时隙内平均')
    fig.savefig(out/'buoy_navigation.png',dpi=165);plt.close(fig)


def main():
    out=HERE/'outputs';manifest=json.loads((out/'manifest.json').read_text(encoding='utf-8'))
    rows=read(out/'per_seed.csv');navrows=read(out/'navigation_per_seed.csv')
    expected=manifest['trials']*2*5*4*2*8
    if len(rows)!=expected or len(navrows)!=manifest['trials']*2:raise ValueError('Incomplete formal batch; report requires every case')
    add_power_reference(rows,manifest,out)
    groups=defaultdict(list)
    for r in rows:groups[float(r['hs']),r['trajectory'],r['method'],float(r['frequency_mhz']),r['phase']].append(r)
    numerical=[k for k in rows[0] if k not in ('seed','hs','trajectory','method','frequency_mhz','phase')]
    if not np.isfinite(np.array([[r[k] for k in numerical] for r in rows],dtype=float)).all():
        raise ValueError('Nonfinite formal results; no sample is silently removed')
    summary=[]
    for (hs,kind,method,f,phase),rr in groups.items():
        if {int(r['seed']) for r in rr}!=set(manifest['seeds']) or len(rr)!=manifest['trials']:
            raise ValueError('Missing or duplicated seed in formal comparison')
        item=dict(hs=hs,trajectory=kind,method=method,frequency_mhz=f,phase=phase,seeds=len(rr))
        for metric in numerical:
            item[metric],item[metric+'_ci95_half']=mean_interval([r[metric] for r in rr])
        summary.append(item)
    write(out/'summary.csv',summary)
    navsummary=[]
    for hs in (.5,1.5):
        rr=[r for r in navrows if float(r['hs'])==hs];item=dict(hs=hs,seeds=len(rr))
        for metric in rr[0]:
            if metric not in ('seed','hs'):item[metric],item[metric+'_ci95_half']=mean_interval([r[metric] for r in rr])
        navsummary.append(item)
    write(out/'navigation_summary.csv',navsummary)
    comparisons=[]
    for hs in (.5,1.5):
        for kind in KINDS:
            for frequency in (300.,1000.):
                for phase in ('all','acceleration','braking','turn_entry'):
                    for left,right,label in (('IMM_CI_linear','IMM_CI_wave','wave_minus_linear'),('IMM_DPC_linear','IMM_CI_linear','CI_minus_DPC'),('CV_CI_wave','IMM_CI_wave','IMM_minus_CV')):
                        a={int(r['seed']):r for r in groups[hs,kind,left,frequency,phase]};b={int(r['seed']):r for r in groups[hs,kind,right,frequency,phase]}
                        item=dict(hs=hs,trajectory=kind,frequency_mhz=frequency,phase=phase,comparison=label)
                        for metric in ('eta','pos_rmse','future_pos_rmse'):
                            item[metric+'_difference'],item[metric+'_ci95_half']=mean_interval([float(b[s][metric])-float(a[s][metric]) for s in sorted(a)])
                        comparisons.append(item)
    write(out/'paired_comparisons.csv',comparisons)
    plots(out,summary)
    lines=['# 原始传感器驱动海上浮标—UAV完整实验','', '## Material Passport','- 阶段：完整仿真实现及已运行结果。',f'- {manifest["trials"]}种子×2海况×5轨迹×4方法×2载频；{len(rows)}行种子级分段记录。','- 全时段0–90s为主结果；不设置门控、不剔除启动与机动、不平滑曲线。','- 每个浮标为一个零杆臂、共位的各向同性等效发射阵元；固定全网总功率1W。','- 算法读取8颗卫星伪距、50Hz加速度/角速度、5Hz磁场及雷达r/方位/俯仰。无真值初始化，无真实慢快分量输入。','- 导航为16维误差状态滤波；雷达直接采用原始球坐标更新；预测到137–177ms后的实际发射时隙。','', '## 浮标自身导航及未来位置预测','', '|Hs m|位置RMSE m|速度RMSE m/s|姿态RMSE °|直线未来RMSE m|波浪未来RMSE m|六维95%名义区域覆盖率|','|---|---:|---:|---:|---:|---:|---:|']
    for r in navsummary:lines.append(f'|{r["hs"]:g}|{r["pos_rmse"]:.4f}|{r["vel_rmse"]:.4f}|{r["attitude_rmse_deg"]:.3f}|{r["linear_future_rmse"]:.4f}|{r["wave_future_rmse"]:.4f}|{r["coverage95_6"]:.3f}|')
    lines+=['','覆盖率仅用于描述不确定度，不作样本门控或独立同分布显著性检验。','', '## UAV全时段跟踪与实际相干效率','', '|Hs m|轨迹|方法|位置RMSE m|速度RMSE m/s|未来位置RMSE m|η 300MHz|η 1GHz|G 300MHz dB|位置95%名义区域覆盖率|','|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
    for hs in (.5,1.5):
        for kind in KINDS:
            for method in METHODS:
                r=next(r for r in summary if (r['hs'],r['trajectory'],r['method'],r['frequency_mhz'],r['phase'])==(hs,kind,method,300.,'all'))
                r1=next(s for s in summary if (s['hs'],s['trajectory'],s['method'],s['frequency_mhz'],s['phase'])==(hs,kind,method,1000.,'all'))
                lines.append(f'|{hs:g}|{KINDS[kind]}|{METHODS[method][0]}|{r["pos_rmse"]:.4f}|{r["vel_rmse"]:.4f}|{r["future_pos_rmse"]:.4f}|{r["eta"]:.4f}|{r1["eta"]:.4f}|{r["gain_db"]:.3f}|{r["coverage95_position"]:.3f}|')
    lines+=['','## 混合运动：波浪预测相对直线预测的配对变化','', '|Hs m|阶段|Δη 300MHz|95%种子配对区间半宽|','|---|---|---:|---:|']
    for r in comparisons:
        if r['trajectory']=='mixed' and r['frequency_mhz']==300 and r['comparison']=='wave_minus_linear':lines.append(f'|{r["hs"]:g}|{r["phase"]}|{r["eta_difference"]:+.5f}|{r["eta_ci95_half"]:.5f}|')
    lines+=['','区间是固定场景下20个独立种子差值的Student-t区间，未做多重比较校正；不能推广为任意海况和航迹的保证。RMSE列为种子级RMSE的算术平均。G=P_rx/Σa_i²，为相同总功率及真实距离下相对非相干功率和的增益；dB列为每种子增益dB的平均。非相干参考由保存的种子/物理模型重建，不重跑估计算法，不改原始per_seed.csv。','', '## 实现和解释边界']
    lines += ['- '+s for s in manifest['limitations']]
    lines+=['- 补齐的是可运行的传感器驱动链路，不宣称已完成全局最优轨迹或最大期望相干功率的求解。','- 波浪预测相对直线预测可能无收益或退化；导航误差、短预测时间和模型失配均可能主导。具体因果需独立分析，不能只凭曲线归因。','- 新实验为20个共位单阵元浮标；旧实验80阵元，且导航假设不同，禁止将旧/新功率数值当作单变量算法比较。','', '## 公式、文献与复现','详见本目录 THEORY.md、README.md；原始数据、每种子指标、配对区间和源文件哈希在outputs。','', '![总体比较](outputs/overview.png)','![三轴坐标速度](outputs/xyz_mixed.png)','![三维轨迹](outputs/trajectories_3d.png)','![浮标导航](outputs/buoy_navigation.png)','![相干效率曲线](outputs/coherent_curves.png)']
    (HERE/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    audit=dict(status='complete',records=len(rows),expected_records=expected,nav_records=len(navrows),source_matches={p:hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in manifest['source_sha256'].items()},figures=[p.name for p in out.glob('*.png')])
    (out/'report_audit.json').write_text(json.dumps(audit,indent=2),encoding='utf-8')
    print(json.dumps(audit,indent=2))


if __name__=='__main__':main()
