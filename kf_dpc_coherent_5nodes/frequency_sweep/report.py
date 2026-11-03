"""Full 20-seed frequency curves, paired uncertainty, phase-only boundaries."""
from pathlib import Path
import json,hashlib
import numpy as np
import pandas as pd
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent;OUT=HERE/'outputs'
KINDS={'hover':'悬停','straight':'直线','circle':'圆周','mixed':'混合','climb':'爬升'}
METHODS={'GP_noD_q1':('无多普勒','#999999'),'GP_Doppler_q1':('联合多普勒','#bd4d54')}
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':9})

def main():
    chunks=[OUT/'chunks'/str(x) for x in (20261028,20261033,20261038,20261043)]
    manifests=[json.loads((p/'manifest.json').read_text()) for p in chunks]
    seeds=sorted(s for m in manifests for s in m['seeds'])
    if seeds!=list(range(20261028,20261048)) or any(m['status']!='completed' for m in manifests):raise ValueError('Incomplete trials')
    categories={x:'category' for x in ('trajectory','method','timing','phase')}
    data=pd.concat([pd.read_csv(p/'per_seed.csv',dtype=categories) for p in chunks],ignore_index=True)
    if len(data)!=576000:raise ValueError('Full record count mismatch')
    keys=['hs','trajectory','method','timing','phase','frequency_mhz']
    metrics=['eta_energy','eta_frame_mean','eta_frame_p05','gain_db','rx_ideal_model_w']
    grouped=data.groupby(keys,observed=True)
    counts=grouped.seed.nunique()
    if (counts!=20).any() or (grouped.size()!=20).any():raise ValueError('Seed pairing mismatch')
    summary=grouped[metrics].mean()
    half=grouped[metrics].std()*student_t.ppf(.975,19)/np.sqrt(20)
    half.columns=[c+'_ci95_half' for c in half.columns]
    summary=summary.join(half).reset_index();summary.to_csv(OUT/'summary.csv',index=False,encoding='utf-8-sig')
    data.to_csv(OUT/'per_seed_all.csv',index=False,encoding='utf-8-sig')
    # Current position/velocity estimates are frequency invariant. Compare the
    # complete replay, not an arbitrary short check, against the old 20 trials.
    reproduction=pd.concat([pd.read_csv(p/'reproduction.csv') for p in chunks],ignore_index=True)
    reproduction.to_csv(OUT/'reproduction.csv',index=False,encoding='utf-8-sig')
    old=pd.read_csv(HERE.parent/'fmcw_doppler/outputs/per_seed_all.csv')
    old=old[(old.frequency_mhz==300)&(old.timing=='fixed_delay')&(old.phase=='all')&old.method.isin(METHODS)]
    new=data[(data.frequency_mhz==300)&(data.timing=='fixed_delay')&(data.phase=='all')]
    reference=new.merge(old,on=['seed','hs','trajectory','method'],suffixes=('_new','_old'))
    reference['eta_difference']=reference.eta_energy-reference.eta
    reference[['seed','hs','trajectory','method','eta_difference']].to_csv(OUT/'reference_300mhz.csv',index=False,encoding='utf-8-sig')
    # Paired differences of frequency response, using independent seed as unit.
    paired_keys=['seed','hs','trajectory','timing','phase','frequency_mhz']
    a=data[data.method=='GP_noD_q1'].set_index(paired_keys)
    b=data[data.method=='GP_Doppler_q1'].set_index(paired_keys)
    delta=(b[['eta_energy','eta_frame_mean','eta_frame_p05','gain_db']]-a[['eta_energy','eta_frame_mean','eta_frame_p05','gain_db']]).reset_index()
    dg=delta.groupby(paired_keys[1:],observed=True)
    diffs=dg[['eta_energy','eta_frame_mean','eta_frame_p05','gain_db']].mean()
    h=dg[['eta_energy','eta_frame_mean','eta_frame_p05','gain_db']].std()*student_t.ppf(.975,19)/np.sqrt(20)
    h.columns=[c+'_ci95_half' for c in h.columns];diffs.join(h).reset_index().to_csv(OUT/'paired_differences.csv',index=False,encoding='utf-8-sig')
    reference_keys=['seed','hs','trajectory','method','timing','phase']
    compared=['eta_energy','eta_frame_mean','eta_frame_p05','gain_db']
    ref300=data[data.frequency_mhz==300][reference_keys+compared]
    joined=data.merge(ref300,on=reference_keys,suffixes=('','_300'))
    changes=joined[keys+['seed']].copy()
    for metric in compared:changes[metric+'_difference']=joined[metric]-joined[metric+'_300']
    change_metrics=[m+'_difference' for m in compared]
    cg=changes.groupby(keys,observed=True)
    cm=cg[change_metrics].mean();ch=cg[change_metrics].std()*student_t.ppf(.975,19)/np.sqrt(20)
    ch.columns=[c+'_ci95_half' for c in ch.columns]
    cm.join(ch).reset_index().to_csv(OUT/'changes_from_300mhz.csv',index=False,encoding='utf-8-sig')
    del joined,changes
    def curve(hs,kind,method,phase='all',timing='recorded_mean_cpu'):
        return summary[(summary.hs==hs)&(summary.trajectory==kind)&(summary.method==method)&(summary.phase==phase)&(summary.timing==timing)].sort_values('frequency_mhz')
    fig,axes=plt.subplots(2,5,figsize=(20,8),sharex=True,sharey=True,layout='constrained')
    for i,hs in enumerate((.5,1.5)):
        for j,(kind,label) in enumerate(KINDS.items()):
            ax=axes[i,j]
            for method,(name,color) in METHODS.items():
                for phase,style in (('all','-'),('late75_90','--')):
                    r=curve(hs,kind,method,phase);ax.plot(r.frequency_mhz,r.eta_energy,style,color=color,label=name+(' 全程' if phase=='all' else ' 75–90s'))
            ax.axvspan(0,30,color='gray',alpha=.1);ax.set_xlim(0,300);ax.set_ylim(0,1.02);ax.grid(alpha=.2)
            ax.set_title(f'Hs={hs:g}m · {label}')
            if i:ax.set_xlabel('通信载频 (MHz)')
            if j==0:ax.set_ylabel('能量加权相干效率η')
    axes[0,0].legend(fontsize=7);fig.suptitle('完整20种子：同一轨迹误差、同一历史均值延迟；灰区仅相位模型外推')
    fig.savefig(OUT/'frequency_curves.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(15,5),layout='constrained')
    for ax,hs in zip(axes,(.5,1.5)):
        for kind,label in KINDS.items():
            r=curve(hs,kind,'GP_Doppler_q1');ax.semilogx(r.frequency_mhz,r.eta_frame_p05,label=label)
        ax.axvspan(.001,30,color='gray',alpha=.1);ax.set_xlim(.001,300);ax.set_ylim(0,1.02);ax.grid(alpha=.2);ax.set_xlabel('通信载频 (MHz，对数轴)');ax.set_ylabel('逐帧效率5%分位的种子均值');ax.set_title(f'Hs={hs:g}m');ax.legend()
    fig.suptitle('1kHz–300MHz理想相位敏感性；低频趋近1不代表真实天线/传播可用')
    fig.savefig(OUT/'phase_extrapolation.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(14,8),sharex=True,layout='constrained')
    for i,hs in enumerate((.5,1.5)):
        for kind,label in KINDS.items():
            r=curve(hs,kind,'GP_Doppler_q1');r=r[r.frequency_mhz>=30]
            axes[i,0].plot(r.frequency_mhz,r.gain_db,label=label)
            axes[i,1].plot(r.frequency_mhz,10*np.log10(1000*r.rx_ideal_model_w),label=label)
        axes[i,0].set_ylabel(f'Hs={hs:g}m · 相对非相干增益 (dB)')
        axes[i,1].set_ylabel(f'Hs={hs:g}m · 理想信道参考接收功率 (dBm)')
        for ax in axes[i]:ax.set_xlim(30,300);ax.grid(alpha=.2);ax.legend(fontsize=8)
    for ax in axes[-1]:ax.set_xlabel('通信载频 (MHz)')
    fig.suptitle('20节点总通信功率1W、固定单位增益参考；绝对功率不代表真实天线/海上传播性能')
    fig.savefig(OUT/'power_reference.png',dpi=165);plt.close(fig)
    selected=summary[(summary.frequency_mhz.isin([10,30,50,75,100,150,200,250,300]))&(summary.hs==1.5)&(summary.method=='GP_Doppler_q1')&(summary.timing=='recorded_mean_cpu')]
    selected.to_csv(OUT/'selected_frequencies.csv',index=False,encoding='utf-8-sig')
    table=selected[selected.phase=='all'].pivot(index='frequency_mhz',columns='trajectory',values='eta_energy')
    lines=['# 300MHz以下通信载频：完整配对扫描结果','',
        '20种子、两海况、五轨迹、两种估计器，360个数值频率点、两种时延、全程及75–90s，共576000行。全部0–90s帧保留。',
        '','连续频率函数由THEORY.md给出；网格不是穷尽实数频率，更不是可用/许可频段清单。77GHz测速链路不变，通信频率不改变轨迹、速度误差或覆盖率。',
        '','## Hs=1.5m联合多普勒：历史均值处理延迟，能量加权效率','',
        '|MHz|悬停|直线|圆周|混合|爬升|','|---:|---:|---:|---:|---:|---:|']
    for f,r in table.iterrows():lines.append('|'+f'{f:g}'+'|'+'|'.join(f'{r[k]:.5f}' for k in KINDS)+'|')
    lines+=['','## 解释限制','',
        '- 1kHz–1MHz为61点对数网格，2–300MHz为整数MHz网格，共360点；300MHz为原参考端点。可以直接使用geometry缓存评估任何中间正频率。',
        '- 30MHz以下仅报告理想相位/相对信道敏感性，不给绝对Friis功率；30MHz分界是展示约定，不是证明此频率以上硬件/海上传播都有效。',
        '- 30–300MHz绝对接收功率仍是假设单位增益、自由空间、固定可辐射功率下的参考，没有固定尺寸天线效率/匹配、海面多径或噪声模型，不能由此选出真实最佳频段。',
        '- recorded_mean_cpu固定使用上一轮同种子同场景的处理时间均值，逐帧不变且跨频率一致；它不是原measured_cpu逐帧时序的精确重放。fixed_delay用于复现原300MHz参考。',
        '- 通信带宽、调制、编码、噪声与数据率没有随载频重建；固定125ms通信成本用于隔离相位效应，不证明极低频能以相同速率传输原数据。',
        '- η能量加权会偏重近距离；同时给逐帧均值、逐帧5%分位及远段，不以平均相干效率代替轨迹可靠性。',
        '- 降频不消除原远段轨迹异常、置信度失配和200ms时限问题，也不证明全局最优。20种子配对区间为探索性描述。',
        f'- 全量复现位置RMSE最大差：{reproduction.position_rmse_difference.abs().max():.3g}m，速度：{reproduction.velocity_rmse_difference.abs().max():.3g}m/s；固定时延300MHz η最大差：{reference.eta_difference.abs().max():.3g}。']
    (HERE/'RESULTS.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    hashes=manifests[0]['source_sha256']
    if any(m['source_sha256']!=hashes for m in manifests):raise ValueError('Worker sources differ')
    match={p:hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in hashes.items()}
    (OUT/'manifest_all.json').write_text(json.dumps(dict(status='completed',rows=len(data),seeds=seeds,
        frequencies_mhz=manifests[0]['frequencies_mhz'],source_sha256=hashes,source_hash_match=match,
        worker_elapsed_s=[m['elapsed_s'] for m in manifests],report_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2),encoding='utf-8')
    print(json.dumps(dict(rows=len(data),seeds=len(seeds),all_sources_match=all(match.values())),indent=2))

if __name__=='__main__':main()
