# 动态浮标FMCW距离率完整配对实验

研究问题：从相干复数回波得到的相对径向速度，能否改善UAV三轴速度、轨迹波动、预测置信度及发射相干？结果必须由完整实验判断，当前设计本身不证明有效。正式文献、数学推导和边界见[THEORY.md](THEORY.md)。

完整20种子实验已完成：32000行跟踪/相干指标、1600行测速指标。混合机动q=1时位置RMSE0.4389→0.2392m、速度0.8388→0.1495m/s，误差变化率1.7850→0.6151m/s；直线及爬升远段则退化，不能作为全场景最终方案。完整统计见[RESULTS.md](RESULTS.md)，判断和失效证据见[ANALYSIS.md](ANALYSIS.md)。图表均保留完整时间，无后平滑。

## 数据链与对照

- 原连续慢漂移＋JONSWAP海场、伪距/IMU/地磁及角距观测复用原种子；零杆臂，无真值初始化。
- 新增点体FMCW复数ADC：TX线性段50μs、每50μs重复，前10μs固定guard，后40μs采64点，共256段。77GHz为ADC起点瞬时频率；ADC有效扫频5MHz，TX总扫频6.25MHz。历史字段chirp_s=40μs实际代表ADC窗时长，原“40μs发射chirp”表述在此更正，信号数值不变。250m处ADC SNR=5dB，幅度r^-2；全二维频谱选峰，不用真实距离或速度挑峰。
- 新通道只供Doppler，旧角距传感器保留。因此是隔离测速贡献的补充通道实验，不是整套FMCW角距测速整机仿真。
- 五方法为q=1的无D/正确联合D、q=0.01的无D/正确联合D、q=1故意忽略浮标速度。全部收集相同回波并计相同半CPI等待，无D只是不用速度似然。
- 20种子、2海况、5轨迹、完整90s，32000行跟踪/相干指标及1600行测速指标。全部帧保留，无冒烟、小测试、测量门控、删帧、事后平滑或真值选参。
- 统计评价覆盖三轴位置/速度、误差变化率、NEES/覆盖率、未来差分距离、300MHz及1GHz相干、计算与CPI等待。通信载频与77GHz雷达载频分开。
- 同TX段解调要求往返延迟δ≤τ+10μs；统一充分距离条件约1498.962m。完整传播几何审计见outputs/propagation_audit.csv，检查全部种子及采样窗距离，不删帧、不用于选峰。快时间复I/Q的1918.672m距离频率周期不能取代这个传播时序条件；超过guard支持距离需要跨TX段接收模型。审计完成前不将该适用条件声称为已验证。
- 计算使用逐次线性化GNSS薄QR和观测稀疏列乘法，保留完整联合协方差及全部相关块。最初密集启动的源码与日志保存在initial_dense/及outputs/initial_dense/，完整正式批次在等价压缩后从头执行，不拼接部分结果。

## 完整运行

在四个独立PowerShell进程中运行不重叠的正式种子组。每组5种子、两海况，共8000行，不能用单组结果代替完整结论。bootstrap将BLAS线程默认设为1；如已有环境变量，请运行前统一设为1。

```powershell
cd D:\git\sim2
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
D:\deepseek\python.exe -u kf_dpc_coherent_5nodes/fmcw_doppler/run_full.py --trials 5 --seed-start 20261028 --output kf_dpc_coherent_5nodes/fmcw_doppler/outputs/chunks/20261028
```

```powershell
cd D:\git\sim2
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
D:\deepseek\python.exe -u kf_dpc_coherent_5nodes/fmcw_doppler/run_full.py --trials 5 --seed-start 20261033 --output kf_dpc_coherent_5nodes/fmcw_doppler/outputs/chunks/20261033
```

```powershell
cd D:\git\sim2
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
D:\deepseek\python.exe -u kf_dpc_coherent_5nodes/fmcw_doppler/run_full.py --trials 5 --seed-start 20261038 --output kf_dpc_coherent_5nodes/fmcw_doppler/outputs/chunks/20261038
```

```powershell
cd D:\git\sim2
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
D:\deepseek\python.exe -u kf_dpc_coherent_5nodes/fmcw_doppler/run_full.py --trials 5 --seed-start 20261043 --output kf_dpc_coherent_5nodes/fmcw_doppler/outputs/chunks/20261043
```

全部四组完成后生成完整合并报告和科学图：

```powershell
cd D:\git\sim2
D:\deepseek\python.exe kf_dpc_coherent_5nodes/fmcw_doppler/report.py
```

report.py检查全部种子、记录数及各worker源码hash，再汇总。outputs/chunks保留独立原始CSV、进度、manifest和首种子数据，最终输出RESULTS.md、summary.csv、paired_comparisons.csv、manifest_all.json及轨迹/误差/测速图。已有joint_gp和raw_sensor_maritime完整成果保留，不覆盖。

## 实现位置

- waveform.py：只在传感器前向仿真访问运动真值，输出全图FFT测速、固定0.15m/s似然尺度及独立评价记录。
- joint_filter.py：联合导航误差、WNOJ先验和相对距离率因子，保留交叉协方差及Gaussian历史到达信息。
- run_full.py：同原始角距和导航数据的五方法完整配对运行；所有版本计半CPI可用时间，分别报固定处理时延和测得CPU时延。
- report.py：完整汇总和绘图，不做曲线后平滑。

未建模旋翼、海杂波、多径、雷达互扰、相位噪声及完整网络负载。0.15m/s不是理论或实测已校准标准差；512点零填充不是物理速度分辨率翻倍。Nyquist、时间中心近似、点体与SNR设定边界及文献不能直接移植的内容均在THEORY.md列出。不能从一次较好的RMSE就宣布满足低波动、高置信度、实时性及1GHz全网相干。
