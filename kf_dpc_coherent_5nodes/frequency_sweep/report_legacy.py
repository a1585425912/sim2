"""Model validity first, then matched low-frequency comparisons."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
from scipy.stats import chi2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;OUT=HERE/'outputs'
NAMES={'CV_CI_linear':'CV+CI直线预测','CV_CI_wave':'CV+CI波浪预测',
 'IMM_DPC_linear':'IMM+DPC','IMM_CI_linear':'IMM+CI直线预测','IMM_CI_wave':'IMM+CI波浪预测',
 'GP_independent_nav':'GP独立导航近似','GP_noD_q1':'联合GP','JointGP_Gaussian_background':'GP小背景高斯',
 'JointGP_TTrend':'T-Trend','GP_Doppler_q1':'联合GP+多普勒'}
KINDS={'hover':'悬停','straight':'直线','circle':'圆周','mixed':'混合','climb':'爬升'}
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':9})

def summarize(data,keys,metrics):
    group=data.groupby(keys,observed=True)
    if (group.size()!=20).any() or (group.seed.nunique()!=20).any():raise ValueError('Incomplete seed pairing')
    m=group[metrics].mean();h=group[metrics].std()*student_t.ppf(.975,19)/np.sqrt(20)
    h.columns=[x+'_ci95_half' for x in h.columns]
    return m.join(h).reset_index()

def main():
    chunks=[OUT/'legacy_chunks'/str(s) for s in (20261028,20261033,20261038,20261043)]
    manifests=[json.loads((p/'manifest.json').read_text()) for p in chunks]
    seeds=sorted(s for m in manifests for s in m['seeds'])
    if seeds!=list(range(20261028,20261048)) or any(m['status']!='completed' for m in manifests):raise ValueError('Not complete')
    dtype={k:'category' for k in ('trajectory','method','timing','phase')}
    legacy=pd.concat([pd.read_csv(p/'per_seed.csv',dtype=dtype) for p in chunks],ignore_index=True)
    if len(legacy)!=1152000:raise ValueError('Legacy frequency count mismatch')
    core=pd.read_csv(OUT/'per_seed_all.csv',dtype=dtype);core=core[core.timing=='fixed_delay']
    data=pd.concat([legacy,core],ignore_index=True)
    if len(data)!=1440000:raise ValueError('Combined count mismatch')
    data.to_csv(OUT/'models_frequency_per_seed.csv',index=False,encoding='utf-8-sig')
    metrics=['eta_energy','eta_frame_mean','eta_frame_p05','gain_db','rx_ideal_model_w']
    summary=summarize(data,['hs','trajectory','method','timing','phase','frequency_mhz'],metrics)
    summary.to_csv(OUT/'models_frequency_summary.csv',index=False,encoding='utf-8-sig')
    tracking=pd.concat([pd.read_csv(p/'tracking_per_seed.csv') for p in chunks],ignore_index=True)
    # Main Doppler and Gaussian GP tracking estimates are identical across
    # communication f; use original completed records, never invent new ones.
    previous=pd.read_csv(HERE.parent/'fmcw_doppler/outputs/per_seed_all.csv')
    previous=previous[(previous.frequency_mhz==300)&(previous.timing=='fixed_delay')&previous.method.isin(['GP_noD_q1','GP_Doppler_q1'])]
    columns=list(tracking.columns);tracking=pd.concat([tracking,previous[columns]],ignore_index=True)
    if len(tracking)!=16000:raise ValueError('Tracking count mismatch')
    track_metrics=[x for x in columns if x not in ('seed','hs','trajectory','method','phase')]
    ts=summarize(tracking,['hs','trajectory','method','phase'],track_metrics)
    tracking.to_csv(OUT/'models_tracking_per_seed.csv',index=False,encoding='utf-8-sig')
    ts.to_csv(OUT/'models_tracking_summary.csv',index=False,encoding='utf-8-sig')
    # Sharpness complements coverage: an arbitrarily broad region can cover
    # truth while conveying little useful information. Use every stored epoch.
    uncertainty=[];geometry=[p for c in chunks for p in (c/'geometry').glob('*.npz')]
    geometry += [p for s in (20261028,20261033,20261038,20261043) for p in (OUT/'chunks'/str(s)/'geometry').glob('*.npz')]
    if len(geometry)!=2000:raise ValueError('Incomplete covariance caches')
    volume_factor=4*np.pi/3*chi2.ppf(.95,3)**1.5
    for path in geometry:
        seed,hs,kind,*parts=path.stem.split('_');method='_'.join(parts)
        with np.load(path) as d:
            P=d['covariance'];times=d['time_s']
            if P.ndim==3:P=P[:,None]
            Pp=P[...,:3,:3];det=np.linalg.det(Pp);trace=np.trace(Pp,axis1=-2,axis2=-1)
            eigenvalues=np.linalg.eigvalsh(Pp)
            masks=dict(all=np.ones(len(times),bool),startup=times<5,after_startup=times>=5,
                acceleration=(times>=15)&(times<20),turn_entry=(times>=35)&(times<40),
                turn_steady=(times>=40)&(times<55),braking=(times>=65)&(times<70),hover_end=(times>=75)&(times<90))
            for phase,mask in masks.items():
                bad=int(np.sum(det[mask]<0));bad_eigen=int(np.sum(np.any(eigenvalues[mask]<0,axis=-1)))
                volume=float(volume_factor*np.mean(np.sqrt(det[mask]))) if not bad and not bad_eigen else float('nan')
                uncertainty.append(dict(seed=int(seed),hs=float(hs[2:]),trajectory=kind,method=method,phase=phase,
                    position_sigma_rss=float(np.sqrt(trace[mask].mean())),position95_volume_mean_m3=volume,negative_determinants=bad,
                    negative_eigenvalue_epochs=bad_eigen,minimum_position_eigenvalue=float(eigenvalues[mask].min())))
    uncertainty=pd.DataFrame(uncertainty)
    uncertainty.to_csv(OUT/'models_uncertainty_per_seed.csv',index=False,encoding='utf-8-sig')
    us=summarize(uncertainty,['hs','trajectory','method','phase'],['position_sigma_rss','position95_volume_mean_m3','negative_determinants','negative_eigenvalue_epochs','minimum_position_eigenvalue'])
    us.to_csv(OUT/'models_uncertainty_summary.csv',index=False,encoding='utf-8-sig')
    rep=pd.concat([pd.read_csv(p/'reproduction.csv') for p in chunks],ignore_index=True)
    rep.to_csv(OUT/'legacy_reproduction.csv',index=False,encoding='utf-8-sig')
    shown=['CV_CI_linear','IMM_DPC_linear','IMM_CI_linear','GP_noD_q1','JointGP_TTrend','GP_Doppler_q1']
    colors=plt.cm.tab10(np.linspace(0,.9,len(shown)))
    fig,axes=plt.subplots(2,5,figsize=(20,8),sharex=True,sharey=True,layout='constrained')
    for i,hs in enumerate((.5,1.5)):
        for j,(kind,name) in enumerate(KINDS.items()):
            ax=axes[i,j]
            for method,color in zip(shown,colors):
                r=summary[(summary.hs==hs)&(summary.trajectory==kind)&(summary.method==method)&(summary.phase=='all')].sort_values('frequency_mhz')
                ax.plot(r.frequency_mhz,r.eta_energy,color=color,label=NAMES[method])
            ax.axvspan(0,30,color='gray',alpha=.1);ax.set_xlim(0,300);ax.set_ylim(0,1.02);ax.grid(alpha=.2);ax.set_title(f'Hs={hs:g}m · {name}')
            if i:ax.set_xlabel('通信载频 (MHz)')
            if j==0:ax.set_ylabel('能量加权效率η')
    axes[0,0].legend(fontsize=7);fig.suptitle('先核对模型有效性，再读频率收益：共同183.395ms时隙末端；灰区仅相位外推')
    fig.savefig(OUT/'models_frequency_curves.png',dpi=165);plt.close(fig)
    methods=list(NAMES);x=np.arange(len(methods))
    fig,axes=plt.subplots(2,3,figsize=(19,10),layout='constrained')
    for ax,metric,label in zip(axes.ravel()[:4],['pos_rmse','vel_rmse','position_error_rate','velocity_error_rate'],['位置RMSE (m)','速度RMSE (m/s)','位置误差变化率 (m/s)','速度误差变化率 (m/s²)']):
        r=ts[(ts.hs==1.5)&(ts.trajectory=='mixed')&(ts.phase=='all')].set_index('method').loc[methods]
        ax.bar(x,r[metric]);ax.set_xticks(x,[NAMES[m] for m in methods],rotation=65,ha='right',fontsize=7);ax.set_title(label);ax.grid(axis='y',alpha=.2)
    for ax,metric,label in zip(axes.ravel()[4:],['coverage95_position','nees6'],['位置名义95%区域覆盖率','六维NEES（对数轴）']):
        for k,phase in enumerate(('all','acceleration','braking')):
            r=ts[(ts.hs==1.5)&(ts.trajectory=='mixed')&(ts.phase==phase)].set_index('method').loc[methods]
            ax.bar(x+(k-1)*.24,r[metric],width=.24,label={'all':'全程','acceleration':'加速','braking':'刹停'}[phase])
        ax.axhline(.95 if metric.startswith('coverage') else 6,color='black',ls='--',lw=.8)
        if metric=='nees6':ax.set_yscale('log')
        ax.set_xticks(x,[NAMES[m] for m in methods],rotation=65,ha='right',fontsize=7);ax.set_title(label);ax.legend(fontsize=8);ax.grid(axis='y',alpha=.2)
    fig.suptitle('模型有效性证据 · Hs=1.5m混合运动 · 完整20种子；参考线不是测量门控或新通过门限')
    fig.savefig(OUT/'model_validity_comparison.png',dpi=160);plt.close(fig)
    selected=summary[summary.frequency_mhz.isin([30,50,75,100,150,200,250,300])]
    selected.to_csv(OUT/'models_selected_frequencies.csv',index=False,encoding='utf-8-sig')
    lines=['# 旧模型有效性与低频效果：完整共同场景对照','',
      '20种子、两海况、五轨迹、10个策略，360频率、两时段，公共时隙末端183.395ms，共1440000行。新增重放12800条轨迹分段指标，加原GP/多普勒3200条，共16000条。',
      '','## 先看模型：Hs=1.5m混合运动完整0–90s','',
      '|策略|位置RMSE m|速度RMSE m/s|位置误差变化率 m/s|位置覆盖率|NEES6|',
      '|---|---:|---:|---:|---:|---:|']
    for method in methods:
        r=ts[(ts.hs==1.5)&(ts.trajectory=='mixed')&(ts.phase=='all')&(ts.method==method)].iloc[0]
        lines.append(f'|{NAMES[method]}|{r.pos_rmse:.4f}|{r.vel_rmse:.4f}|{r.position_error_rate:.4f}|{r.coverage95_position:.4f}|{r.nees6:.3f}|')
    lines+=['','## 再看频率：Hs=1.5m混合运动，全程能量加权η','',
      '|策略|30MHz|50MHz|100MHz|150MHz|200MHz|300MHz|','|---|---:|---:|---:|---:|---:|---:|']
    for method in methods:
        r=summary[(summary.hs==1.5)&(summary.trajectory=='mixed')&(summary.method==method)&(summary.phase=='all')].set_index('frequency_mhz')
        lines.append('|'+NAMES[method]+'|'+'|'.join(f'{r.loc[f,"eta_energy"]:.5f}' for f in (30,50,100,150,200,300))+'|')
    lines+=['','## 解释约束','',
      '- CV_CI_linear为新增受控预测对照，与CV_CI_wave共享同一CV轨迹；IMM的直线/波浪预测也共享同一估计，不视为新的独立模型或样本。',
      '- 旧原始导航/角距信息相同；多普勒模型额外获得77GHz观测，不把它的全部收益归因于换优化器。GP集中与有限轮DPC/CI信息结构不同，独立导航GP为相关性消融。',
      '- 当前10策略频率比较统一183.395ms末端，旧模型额外等待仅用于控制时延；不是证明它们实际处理均为12ms，也不作为实时部署结果。',
      '- 原始177ms时序仅用于300MHz复现；CV_CI_linear是新预测消融，没有旧相干参考。旧阵列/简化导航历史版本不混入公平排名。',
      '- 轨迹误差、波动、机动与协方差可信度不随通信f改变；低频η趋近1不能证明估计器有效。原始q001等参数诊断保留于MODEL_VALIDITY.md和原完整成果，不伪称本轮将所有废弃参数版本重新评估。',
      '- models_uncertainty_summary.csv另报sqrt(mean trace(Ppp))及位置95%椭球平均体积，区分覆盖率与区域宽度；体积为4π/3×χ²(3,.95)^(3/2)×sqrt(det(Ppp))。覆盖率高本身不证明区域精确或信息充分。',
      '- 极低频为固定幅度相位模型外推，真实天线、互耦、海面传播及带宽/数据率仍未建模。全部结果保留，没有按效果筛模型或删帧。',
      f'- 完整旧模型重放的位置RMSE最大差{rep.pos_rmse_difference.abs().max():.3g}m、速度{rep.vel_rmse_difference.abs().max():.3g}m/s、原时序300MHz η最大差{rep.native_eta_difference.abs().max():.3g}。']
    (HERE/'LEGACY_RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    hashes=manifests[0]['source_sha256']
    if any(m['source_sha256']!=hashes for m in manifests):raise ValueError('Source mismatch')
    match={p:hashlib.sha256(Path(p).read_bytes()).hexdigest()==h for p,h in hashes.items()}
    (OUT/'legacy_manifest_all.json').write_text(json.dumps(dict(status='completed',seeds=seeds,legacy_frequency_rows=len(legacy),combined_frequency_rows=len(data),tracking_rows=len(tracking),
       source_sha256=hashes,source_hash_match=match,worker_elapsed_s=[m['elapsed_s'] for m in manifests],report_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2),encoding='utf-8')
    print(json.dumps(dict(rows=len(data),tracking_rows=len(tracking),all_sources_match=all(match.values())),indent=2))

if __name__=='__main__':main()
