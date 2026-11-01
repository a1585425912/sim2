# 文献驱动的联合GP轨迹估计与相干发射

当前入口为run_joint.py，调用sequential.py：**持久高斯到达因子＋含加速度的GP动力学＋随机克隆固定滞后联合估计**。每次只加入当前原始观测一次，完整保留导航—目标及浮标间交叉相关。线性高斯特例与增广Kalman平滑等价，不能称为完全脱离KF的算法。

原始运动、GNSS伪距、IMU、地磁和雷达全部复用raw_sensor_maritime同种子仿真；无真值状态初始化、无真实慢快运动分量、无未来观测输入、零杆臂。当前UAV输出是在线末端，不是整段数据事后拟合结果。

## 正式完整实验

本轮另有完整20种子的T-Trend过程先验与Gaussian匹配协方差对照，见T_TREND.md。Gaussian主研究32000行，扩展12800行，最终报告合计44800行；两项均全部保留，并非只选较好的方案。

```powershell
cd D:\git\sim2
python -u kf_dpc_coherent_5nodes/joint_gp/run_joint.py --trials 20
python -u kf_dpc_coherent_5nodes/joint_gp/run_adaptive.py --trials 20
python kf_dpc_coherent_5nodes/joint_gp/merge_adaptive.py
python kf_dpc_coherent_5nodes/joint_gp/report.py
```

本次为缩短完整运行等待时间，使用4个相互独立、种子不重叠的正式进程；每个进程将OPENBLAS_NUM_THREADS/OMP_NUM_THREADS/MKL_NUM_THREADS设置为1：

```powershell
python -u kf_dpc_coherent_5nodes/joint_gp/run_joint.py --trials 5 --seed-start 20261028 --output kf_dpc_coherent_5nodes/joint_gp/outputs/chunks/20261028
python -u kf_dpc_coherent_5nodes/joint_gp/run_joint.py --trials 5 --seed-start 20261033 --output kf_dpc_coherent_5nodes/joint_gp/outputs/chunks/20261033
python -u kf_dpc_coherent_5nodes/joint_gp/run_joint.py --trials 5 --seed-start 20261038 --output kf_dpc_coherent_5nodes/joint_gp/outputs/chunks/20261038
python -u kf_dpc_coherent_5nodes/joint_gp/run_joint.py --trials 5 --seed-start 20261043 --output kf_dpc_coherent_5nodes/joint_gp/outputs/chunks/20261043
python kf_dpc_coherent_5nodes/joint_gp/merge.py
python -u kf_dpc_coherent_5nodes/joint_gp/run_adaptive.py --trials 10 --seed-start 20261028 --output kf_dpc_coherent_5nodes/joint_gp/adaptive_outputs/chunks/20261028
python -u kf_dpc_coherent_5nodes/joint_gp/run_adaptive.py --trials 10 --seed-start 20261038 --output kf_dpc_coherent_5nodes/joint_gp/adaptive_outputs/chunks/20261038
python kf_dpc_coherent_5nodes/joint_gp/merge_adaptive.py
python kf_dpc_coherent_5nodes/joint_gp/report.py
```

每组5种子、2海况、5轨迹、5方法、2载频、2种时延评价、8个统计时段，共8000行；完整合并32000行。全部0–90s为主结果，分段只解释，不删除启动/机动数据。不运行冒烟、小测试或门限测试。

## 对照的含义

|方法|估计|发射预测与功率|
|---|---|---|
|IMM_CI_linear|原始雷达IMM＋有限25轮CI，原样重跑|浮标直线预测，均匀功率|
|GP_independent_nav|持久GP轨迹，但近似各帧导航误差独立；同样集中观测|浮标参考导航直线预测，均匀功率|
|JointGP_linear|持久GP＋全部导航误差联合，保留完整交叉块|联合修正浮标直线预测，均匀功率|
|JointGP_accel|同一份联合后验|当前IMU比力、联合姿态/偏置恢复加速度，二阶预测，均匀功率|
|JointGP_power|同一份联合后验|同样惯性预测，加条件期望接收功率主特征向量分配|

后三种目标轨迹相同，只改变发射预测/功率策略。功率优化允许连续重分配每节点功率、总功率始终1W；若硬件要求各节点功率严格固定，应采用均匀功率主对照。没有节点门控或删帧，另报告有效节点数以防把功率集中当作全网高相干。

## 计算与网络边界

新估计是集中联合信息参考，不是已完成的分布式GP。全网观测共享与有限轮CI不等价，因此用GP_independent_nav控制集中观测带来的贡献，不能把全部变化归因于联合GP理论。

分别评价旧固定137ms起始时延，以及以实际求解＋导航平均耗时替代12ms处理假设的measured_cpu。通信仍假设125ms；并行电脑耗时不是专用硬件最坏时延，预测、编解码和带宽开销尚未完整计入。不得根据固定时延的更好结果声称已满足200ms帧内调度。

## 记录、理论与失败方案

- RESEARCH.md：正式文献、具体推导、实现与一致性边界。
- RESULTS.md：仅完整批次成功后生成。
- outputs/summary.csv、paired_comparisons.csv：种子级均值和配对区间。
- 最终有效合并表为outputs/per_seed_all.csv（44800行）。原始Gaussian outputs/per_seed.csv与各子进程文件保留执行当时版本；其中旧基线measured_cpu评估曾漏做预测CI，已按全部20种子重新计算并修正最终合并表，原始快照、修正程序和日志在outputs/corrections/。固定时延数据与全部估计轨迹不受该错误影响。
- outputs/manifest.json、report_audit.json：各组来源一致性、首种子旧基线复现差异、完整计数及图表。
- sequential.py：实际主估计器和联合几何预测。
- navigation.py：原始IMU转移和误差坐标参考，不把参考ESKF后验再次当作联合观测。
- solver.py：保留初版函数及当前复用的GP数学/参数助手；当前入口不调用其有限历史fit。

第一版仅2s有限历史、每次弱先验重置，正式批次首批结果显示误差变化率变大，故停止该结构并保留全部中间数据：源码在finite_history_v1/，4组CSV和日志在outputs/finite_history_v1/。更早固定时延启动的日志在outputs/initial_fixed_delay_run.log。这些均不冒充正式完整数据，也没有隐藏负面结果。当前修订没有通过真值调小q或美化曲线，而是持续保留历史信息。

会话跨日恢复后才汇总，manifest的overall_elapsed_s包含等待/会话暂停，不作为算法耗时；算法时间取每帧compute字段，研究进程时间取各worker_elapsed_s。checkpoint_prefix.py为曾拟定但未执行的调度助手，本次没有使用前缀截取，实际两项研究均按原完整分组自然结束。
