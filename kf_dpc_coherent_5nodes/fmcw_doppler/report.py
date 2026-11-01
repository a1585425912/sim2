"""Summarize only complete paired trials; scientific plots, no curve smoothing."""
import csv,json,hashlib,shutil
from pathlib import Path
from collections import defaultdict
import numpy as np
from scipy.stats import t as student_t
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE=Path(__file__).resolve().parent;OUT=HERE/'outputs'
METHODS={
    'GP_noD_q1':('无多普勒 q=1','#929292'),
    'GP_Doppler_q1':('联合多普勒 q=1','#c45052'),
    'GP_noD_q001':('无多普勒 q=0.01','#b99749'),
    'GP_Doppler_q001':('联合多普勒 q=0.01','#268e79'),
    'GP_Doppler_ignore_buoy':('错误忽略浮标速度','#7564ac'),
}
KINDS={'hover':'悬停','straight':'直线','circle':'圆周','mixed':'混合','climb':'爬升'}
plt.rcParams.update({'font.sans-serif':['Microsoft YaHei','SimHei','DejaVu Sans'],'axes.unicode_minus':False,'font.size':9,'axes.spines.top':False,'axes.spines.right':False})

def read(path):
    with path.open(encoding='utf-8-sig') as f:return list(csv.DictReader(f))

def write(path,rows):
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def stats(values):
    a=np.asarray(values,float)
    return float(a.mean()),float(student_t.ppf(.975,len(a)-1)*a.std(ddof=1)/np.sqrt(len(a)))

def main():
    chunks=[OUT/'chunks'/str(s) for s in (20261028,20261033,20261038,20261043)]
    rows=[];sensors=[];manifests=[]
    for chunk in chunks:
        m=json.loads((chunk/'manifest.json').read_text(encoding='utf-8'))
        if m['status']!='completed':raise ValueError(f'Incomplete {chunk}')
        manifests.append(m);rows+=read(chunk/'per_seed.csv');sensors+=read(chunk/'sensor_per_seed.csv')
    seeds=sorted(s for m in manifests for s in m['seeds'])
    if seeds!=list(range(20261028,20261048)):raise ValueError('Seeds incomplete/duplicated')
    if len(rows)!=32000 or len(sensors)!=1600:raise ValueError('Full experiment count mismatch')
    groups=defaultdict(list)
    keys=('hs','trajectory','method','frequency_mhz','timing','phase')
    metrics=[k for k in rows[0] if k not in (*keys,'seed')]
    for r in rows:
        for name in metrics:
            if not np.isfinite(float(r[name])):raise ValueError(f'Non-finite {name}')
        groups[tuple(r[k] for k in keys)].append(r)
    summary=[]
    for key,group in groups.items():
        if sorted(int(r['seed']) for r in group)!=seeds:raise ValueError('Unpaired group')
        out=dict(zip(keys,key));out['seeds']=len(group)
        for name in metrics:out[name],out[name+'_ci95_half']=stats([r[name] for r in group])
        summary.append(out)
    sensor_groups=defaultdict(list)
    for r in sensors:sensor_groups[tuple(r[k] for k in ('hs','trajectory','phase'))].append(r)
    sensor_summary=[]
    for key,group in sensor_groups.items():
        out=dict(zip(('hs','trajectory','phase'),key));out['seeds']=len(group)
        for name in sensors[0]:
            if name not in ('hs','trajectory','phase','seed'):out[name],out[name+'_ci95_half']=stats([r[name] for r in group])
        sensor_summary.append(out)
    paired=[]
    pairs=(('GP_noD_q1','GP_Doppler_q1'),('GP_noD_q001','GP_Doppler_q001'),
           ('GP_Doppler_q1','GP_Doppler_q001'),('GP_Doppler_q1','GP_Doppler_ignore_buoy'))
    for key,group in groups.items():
        if key[2]!='GP_noD_q1':continue
        for left,right in pairs:
            kl=list(key);kr=list(key);kl[2]=left;kr[2]=right
            a={int(r['seed']):r for r in groups[tuple(kl)]};b={int(r['seed']):r for r in groups[tuple(kr)]}
            out=dict(zip(keys,key));out.pop('method');out.update(left=left,right=right)
            for name in ('pos_rmse','vel_rmse','pZ_rmse','vZ_rmse','future_pos_rmse','path_rmse','position_error_rate','velocity_error_rate','eta','gain_db','coverage95_position','nees6'):
                out[name+'_difference'],out[name+'_ci95_half']=stats([float(b[s][name])-float(a[s][name]) for s in seeds])
            paired.append(out)
    OUT.mkdir(exist_ok=True)
    write(OUT/'per_seed_all.csv',rows);write(OUT/'summary.csv',summary)
    write(OUT/'sensor_summary.csv',sensor_summary);write(OUT/'paired_comparisons.csv',paired)
    for hs in (.5,1.5):shutil.copy2(chunks[0]/f'first_seed_hs{hs:g}.npz',OUT/f'first_seed_hs{hs:g}.npz')
    def get(hs,kind,method,freq=300,timing='fixed_delay',phase='all'):
        return next(r for r in summary if float(r['hs'])==hs and r['trajectory']==kind and r['method']==method and float(r['frequency_mhz'])==freq and r['timing']==timing and r['phase']==phase)
    fig,axes=plt.subplots(2,3,figsize=(17,8),layout='constrained');xx=np.arange(5)
    for i,hs in enumerate((.5,1.5)):
        for j,(metric,title) in enumerate((('pos_rmse','位置RMSE (m)'),('vel_rmse','速度RMSE (m/s)'),('position_error_rate','位置误差变化率 (m/s)'))):
            ax=axes[i,j]
            for n,(method,(label,color)) in enumerate(METHODS.items()):ax.bar(xx+(n-2)*.16,[get(hs,k,method)[metric] for k in KINDS],width=.16,color=color,label=label)
            ax.set_xticks(xx,list(KINDS.values()));ax.set_title(f'Hs={hs:g}m · {title}');ax.grid(axis='y',alpha=.2)
    axes[0,0].legend(fontsize=7);fig.suptitle('20种子、完整0–90s：相同传感器、均匀20节点，无输出端平滑')
    fig.savefig(OUT/'overview.png',dpi=165);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(14,8),layout='constrained')
    for i,hs in enumerate((.5,1.5)):
        for j,freq in enumerate((300,1000)):
            ax=axes[i,j]
            for n,(method,(label,color)) in enumerate(METHODS.items()):ax.bar(xx+(n-2)*.16,[get(hs,k,method,freq,'measured_cpu')['eta'] for k in KINDS],width=.16,color=color,label=label)
            ax.set_xticks(xx,list(KINDS.values()));ax.set_title(f'Hs={hs:g}m · {freq}MHz');ax.set_ylabel('相干效率η');ax.grid(axis='y',alpha=.2)
    axes[0,0].legend(fontsize=7);fig.suptitle('含半CPI、实测求解和ADC处理耗时；通信仍为125ms假设')
    fig.savefig(OUT/'coherence.png',dpi=165);plt.close(fig)
    chosen=('GP_noD_q1','GP_Doppler_q1','GP_Doppler_q001','GP_Doppler_ignore_buoy')
    d=np.load(OUT/'first_seed_hs1.5.npz');t=d['time_s'];truth=np.c_[d['mixed__truth_p'],d['mixed__truth_v']]
    fig,axes=plt.subplots(3,2,figsize=(16,10),sharex=True,layout='constrained')
    for axis in range(3):
        for j in range(2):
            ax=axes[axis,j];col=axis+3*j
            for method in chosen:
                label,color=METHODS[method];value=d[f'mixed__{method}__state'][:,0,col]
                ax.plot(t,value-truth[:,col],color=color,label=label,lw=.85)
            ax.axhline(0,color='black',lw=.7,ls='--');ax.grid(alpha=.2)
            ax.set_ylabel('XYZ'[axis]+('位置误差 (m)' if j==0 else '速度误差 (m/s)'))
    axes[0,0].legend(fontsize=8)
    for ax in axes[-1]:ax.set_xlabel('时间 (s)')
    fig.suptitle('首种子混合运动 · 全时段因果输出误差，无截断或后平滑')
    fig.savefig(OUT/'mixed_errors.png',dpi=170);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(16,8),sharex=True,layout='constrained')
    for i,kind in enumerate(('straight','climb')):
        truth=np.c_[d[f'{kind}__truth_p'],d[f'{kind}__truth_v']]
        for method in chosen:
            label,color=METHODS[method];error=d[f'{kind}__{method}__state'][:,0,:6]-truth
            axes[i,0].plot(t,np.linalg.norm(error[:,:3],axis=-1),color=color,label=label,lw=.85)
            axes[i,1].plot(t,np.linalg.norm(error[:,3:6],axis=-1),color=color,label=label,lw=.85)
        axes[i,0].set_ylabel(KINDS[kind]+'位置误差范数 (m)');axes[i,1].set_ylabel(KINDS[kind]+'速度误差范数 (m/s)')
        for ax in axes[i]:ax.grid(alpha=.2)
    axes[0,0].legend(fontsize=8)
    for ax in axes[-1]:ax.set_xlabel('时间 (s)')
    fig.suptitle('完整远距离轨迹：保留测速异常及轨迹跳变，无删帧或截断')
    fig.savefig(OUT/'far_errors.png',dpi=170);plt.close(fig)
    fig=plt.figure(figsize=(20,7))
    for j,(kind,label) in enumerate(KINDS.items()):
        ax=fig.add_subplot(1,5,j+1,projection='3d');p=d[f'{kind}__truth_p']
        ax.plot(*p.T,'k--',lw=1.4)
        if kind=='hover':ax.scatter(*p[0],marker='x',color='black')
        for method in chosen:
            ax.plot(*d[f'{kind}__{method}__state'][:,0,:3].T,color=METHODS[method][1],lw=.8)
        ax.set_title(label);ax.set_xlabel('X (m)');ax.set_ylabel('Y (m)');ax.set_zlabel('Z (m)');ax.view_init(23,-58)
    handles=[Line2D([],[],color='black',ls='--',label='真值')]+[Line2D([],[],color=METHODS[m][1],label=METHODS[m][0]) for m in chosen]
    fig.legend(handles=handles,loc='lower center',ncol=5);fig.suptitle('Hs=1.5m 首种子完整3D轨迹；各轴独立缩放，非等比例几何')
    fig.subplots_adjust(left=.02,right=.97,bottom=.13,top=.9,wspace=.3);fig.savefig(OUT/'trajectories_3d.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(15,7),sharex=True,layout='constrained')
    for kind,ls in (('hover','-'),('mixed','--')):
        axes[0].plot(t,d[f'{kind}__doppler_measured'][:,0],ls,label=KINDS[kind]+'测得',lw=.9)
        axes[0].plot(t,d[f'{kind}__doppler_truth'][:,0],color='black',ls=ls,label=KINDS[kind]+'真值',lw=.8)
        axes[1].plot(t,d[f'{kind}__doppler_measured'][:,0]-d[f'{kind}__doppler_truth'][:,0],ls,label=KINDS[kind],lw=.9)
    axes[0].set_ylabel('相对径向速度 (m/s)');axes[1].set_ylabel('测速误差 (m/s)');axes[1].set_xlabel('CPI中心时刻 (s)')
    for ax in axes:ax.grid(alpha=.2);ax.legend()
    fig.suptitle('FMCW复数ADC全图选峰输出；悬停时浮标运动仍产生多普勒')
    fig.savefig(OUT/'doppler_measurements.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(14,5),layout='constrained')
    wavelength=manifests[0]['waveform']['cpi_s']*2*manifests[0]['waveform']['velocity_fft_resolution']
    vmax=manifests[0]['waveform']['velocity_unambiguous'];rperiod=manifests[0]['waveform']['range_unambiguous_m']
    for ax,kind in zip(axes,('hover','mixed')):
        power=d[f'{kind}__frame100_range_doppler']
        im=ax.imshow(10*np.log10(np.maximum(power,1e-30)),origin='lower',aspect='auto',
                     extent=(0,rperiod,-vmax,vmax),cmap='viridis')
        ax.set_title(KINDS[kind]+' · t=20s · 浮标1');ax.set_xlabel('快频FFT距离映射 (m)');ax.set_ylabel('多普勒速度 (m/s)')
        fig.colorbar(im,ax=ax,label='频谱功率 (dB，任意参考)')
    fig.suptitle('实际含噪复数ADC的完整距离—多普勒图；未用真值限定选峰区域')
    fig.savefig(OUT/'range_doppler.png',dpi=160);plt.close(fig)
    text=['# FMCW距离率完整配对实验结果','',
          '20种子×2海况×5轨迹×5方法×2载频×2时延×8时段=32000行；原始测速1600行。全部0–90s，无门控、删帧或真值选峰。',
          '','## Hs=1.5m混合运动，全部时间，公共固定时延含半CPI','',
          '|方法|位置RMSE m|速度RMSE m/s|Z速度RMSE m/s|位置误差变化率 m/s|位置95%覆盖率|六维NEES|未来位置RMSE m|差分距离误差 m|η300MHz|',
          '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for method,(label,_) in METHODS.items():
        r=get(1.5,'mixed',method)
        text.append('|'+label+'|'+'|'.join(f'{r[m]:.4f}' for m in ('pos_rmse','vel_rmse','vZ_rmse','position_error_rate','coverage95_position','nees6','future_pos_rmse','path_rmse','eta'))+'|')
    text+=['','## 两海况五轨迹全部主结果','',
           '|Hs m|轨迹|方法|位置RMSE m|速度RMSE m/s|误差变化率 m/s|覆盖率|含CPU η300|含CPU η1GHz|时隙末端ms|',
           '|---|---|---|---:|---:|---:|---:|---:|---:|---:|']
    for hs in (.5,1.5):
        for kind,label in KINDS.items():
            for method,(name,_) in METHODS.items():
                r=get(hs,kind,method);a=get(hs,kind,method,300,'measured_cpu');b=get(hs,kind,method,1000,'measured_cpu')
                text.append(f'|{hs}|{label}|{name}|{r["pos_rmse"]:.4f}|{r["vel_rmse"]:.4f}|{r["position_error_rate"]:.4f}|{r["coverage95_position"]:.4f}|{a["eta"]:.4f}|{b["eta"]:.4f}|{a["tx_end_mean_ms"]:.2f}|')
    text+=['','## Hs=1.5m直线/爬升远段75–90s（统计字段hover_end仅为沿用的时段名）','',
           '|轨迹|方法|位置RMSE m|速度RMSE m/s|位置覆盖率|六维NEES|含CPU η300|含CPU η1GHz|',
           '|---|---|---:|---:|---:|---:|---:|---:|']
    for kind in ('straight','climb'):
        for method,(name,_) in METHODS.items():
            a=get(1.5,kind,method,300,'measured_cpu','hover_end');b=get(1.5,kind,method,1000,'measured_cpu','hover_end')
            text.append(f'|{KINDS[kind]}|{name}|{a["pos_rmse"]:.4f}|{a["vel_rmse"]:.4f}|{a["coverage95_position"]:.4f}|{a["nees6"]:.4f}|{a["eta"]:.4f}|{b["eta"]:.4f}|')
    text+=['','## 解释边界','',
           '- 所有方法均处理同一补充FMCW通道，无D版只不使用速度似然。q=1与q=0.01为事先固定的结构对照，不是寻优或真值选参。',
           '- 旧距离/角度传感器保持一致，补充FMCW仅供速度；不是同一ADC同时获得全部角距观测的整机仿真。',
           '- 77GHz为ADC起点瞬时频率；50μs完整扫频含10μs前置保护、末40μs采样，有效采样带宽5MHz、TX总带宽6.25MHz，256扫频、64快时间样本。点机体、单目标全图最大峰；无旋翼、海杂波、多径或互扰。',
           '- 时延计半CPI与实测ADC频谱处理及求解，未计嵌入式最坏耗时、完整编码、真实网络共享负载。',
           '- 当前状态协方差是局部Gaussian近似；测速std固定0.15m/s，不能把FFT分辨率或条件协方差当严格置信保证。',
           '- 250m处ADC SNR=5dB是固定实验假设，未通过实际雷达功率/天线/RCS/噪声系数闭合链路预算；1W约束是通信功率，不证明77GHz感知硬件能达到该距离。',
           '- 图为首种子说明，统计为20种子；32000行不是32000个独立重复。配对95%均值区间见paired_comparisons.csv。']
    text+=['- η全程为ΣP_rx/ΣP_ideal，即按理想接收功率加权的效率，不是每帧η的算术平均；远距离帧功率较低，不能用全程η掩盖远段轨迹跳变。']
    (HERE/'RESULTS.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
    # Current-state reproduction uses the full completed batch, not a short
    # test. Forecasts intentionally differ because half-CPI is now included.
    old=read(HERE.parent/'joint_gp'/'outputs'/'summary.csv')
    reproduction=[]
    for hs in (.5,1.5):
        for kind in KINDS:
            current=get(hs,kind,'GP_noD_q1')
            prior=next(r for r in old if float(r['hs'])==hs and r['trajectory']==kind and r['method']=='JointGP_accel' and float(r['frequency_mhz'])==300 and r['timing']=='fixed_delay' and r['phase']=='all')
            reproduction.append(dict(hs=hs,trajectory=kind,position_rmse_difference=current['pos_rmse']-float(prior['pos_rmse']),velocity_rmse_difference=current['vel_rmse']-float(prior['vel_rmse'])))
    write(OUT/'historical_state_reproduction.csv',reproduction)
    # Physical-domain audit of ALL trials, not a measurement gate. The emitted
    # linear sweep spans the 50us repetition: first 10us unrecorded guard,
    # last 40us ADC. The dechirp expression is valid only for same-sweep RX.
    # A global analytic speed bound proves the full intra-CPI range domain
    # without rerunning FFT or selecting/removing any observation.
    import bootstrap
    from physics import Config,SeaTruth
    import run_joint as original
    wave=manifests[0]['waveform'];guard=wave['repetition_s']-wave['chirp_s']
    distance_limit=299792458.*guard/2
    propagation=[]
    for seed in seeds:
        for hs in (.5,1.5):
            cfg=Config(hs=hs);scene=SeaTruth(cfg,seed)
            times=np.arange(round(cfg.duration/cfg.dt)+1)*cfg.dt
            sea=scene.at(times)
            slow_bound=np.linalg.norm(cfg.mean_current)+np.linalg.norm(scene.slow_a,axis=-1).sum(-1)+np.linalg.norm(scene.slow_b,axis=-1).sum(-1)
            response=np.linalg.norm(scene.Hp,axis=-1);knorm=np.linalg.norm(scene.k,axis=-1)
            buoy_bound=slow_bound+(scene.w[None,:]+slow_bound[:,None]*knorm[None,:])@response
            slow_acc_bound=(scene.slow_w[None,:]*(np.linalg.norm(scene.slow_a,axis=-1)+np.linalg.norm(scene.slow_b,axis=-1))).sum(-1)
            encounter_bound=scene.w[None,:]+slow_bound[:,None]*knorm[None,:]
            buoy_acc_bound=slow_acc_bound+(encounter_bound**2+slow_acc_bound[:,None]*knorm[None,:])@response
            # All prescribed UAV trajectories have speed <=sqrt(15^2+2^2)<16.
            relative_bound=16.+float(buoy_bound.max())
            for kind in KINDS:
                p,v=original.reference.trajectory(times,kind,cfg)
                q=p[:,None]-sea['p'];range_center=np.linalg.norm(q,axis=-1)
                center_max=float(range_center.max())
                adc_bound=center_max+relative_bound*wave['half_span_s']
                min_range=float(range_center.min())-relative_bound*wave['half_span_s']
                radial_acc_bound=4.+float(buoy_acc_bound.max())+relative_bound**2/min_range
                radial_center=float(np.abs(np.sum(q*(v[:,None]-sea['v']),axis=-1)/range_center).max())
                radial_bound=radial_center+radial_acc_bound*wave['half_span_s']
                propagation.append(dict(seed=seed,hs=hs,trajectory=kind,max_center_range_m=center_max,
                    analytic_relative_speed_bound_mps=relative_bound,max_adc_range_bound_m=adc_bound,
                    fixed_tx_guard_s=guard,same_sweep_range_limit_m=distance_limit,
                    margin_m=distance_limit-adc_bound,model_domain_valid=adc_bound<distance_limit,
                    maximum_center_radial_speed_mps=radial_center,
                    minimum_cpi_range_bound_m=min_range,positive_range_bound=min_range>0,
                    analytic_radial_acceleration_bound_mps2=radial_acc_bound,
                    analytic_cpi_radial_speed_bound_mps=radial_bound,
                    unambiguous_radial_speed_mps=wave['velocity_unambiguous'],
                    doppler_domain_valid=radial_bound<wave['velocity_unambiguous']))
    write(OUT/'propagation_audit.csv',propagation)
    source=manifests[0]['source_sha256'];matching={p:Path(p).exists() and hashlib.sha256(Path(p).read_bytes()).hexdigest()==v for p,v in source.items()}
    if not all(m['source_sha256']==source for m in manifests):raise ValueError('Workers used different code')
    # These transitive trajectory definitions were not in the original worker
    # snapshot. Record the exact coverage gap rather than invent start hashes.
    additional_dependencies={}
    earliest_log_creation=min((p/'run.log').stat().st_ctime for p in chunks)
    for dep in (HERE.parent/'imm_dpc'/'run_experiment.py',HERE.parent/'gm_dpc_unified'/'experiment.py'):
        additional_dependencies[str(dep)]=dict(current_sha256=hashlib.sha256(dep.read_bytes()).hexdigest(),
            mtime_unix=dep.stat().st_mtime,mtime_before_formal_logs=dep.stat().st_mtime<earliest_log_creation,
            coverage='Current hash and filesystem mtime only; not included in worker start/end cryptographic snapshot. No edits made by this task.')
    manifest=dict(status='completed',rows=len(rows),sensor_rows=len(sensors),seeds=seeds,methods=list(METHODS),source_sha256=source,source_hash_match=matching,worker_elapsed_s=[m['elapsed_s'] for m in manifests],waveform=manifests[0]['waveform'],
        report_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        additional_dependency_record=additional_dependencies,
        physical_sweep=dict(tx_duration_s=wave['repetition_s'],adc_duration_s=wave['chirp_s'],unrecorded_guard_s=guard,
            total_tx_bandwidth_hz=wave['bandwidth_hz']*wave['repetition_s']/wave['chirp_s'],
            all_trials_same_sweep_valid=all(r['model_domain_valid'] for r in propagation),
            maximum_adc_range_bound_m=max(r['max_adc_range_bound_m'] for r in propagation),
            all_trials_unambiguous_velocity=all(r['doppler_domain_valid'] for r in propagation),
            all_minimum_range_bounds_positive=all(r['positive_range_bound'] for r in propagation),
            maximum_cpi_radial_speed_bound_mps=max(r['analytic_cpi_radial_speed_bound_mps'] for r in propagation)))
    (OUT/'manifest_all.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(dict(rows=len(rows),sensor_rows=len(sensors),seeds=len(seeds),all_sources_match=all(matching.values())),indent=2))

if __name__=='__main__':main()
