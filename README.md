# 动态浮标联合观测UAV与相干发射仿真

300MHz以下通信频率扫描与旧模型对照已完成，位于 `kf_dpc_coherent_5nodes/frequency_sweep/`：20种子、两海况、五轨迹、10策略、360频率，公共条件144万频率指标及1.6万轨迹分段指标。先看ASSESSMENT.md和MODEL_VALIDITY.md的模型有效性判断，再看LEGACY_RESULTS.md及频率曲线。旧模型低频下也能高相干，但目前没有一版已证明全场景轨迹稳定、机动准确、置信度可信并满足原实时约束。原始数据、所有失败场景与几何缓存保留；低频理想外推不等于真实可用/最优频段。

FMCW距离率补充实验位于 `kf_dpc_coherent_5nodes/fmcw_doppler/`，完整20种子批次已结束：从复数回波ADC全图选峰得到相对径向速度，比较联合浮标速度、忽略浮标速度、两档运动先验及无多普勒对照。THEORY.md给出正式来源和推导，RESULTS.md为全部统计，ANALYSIS.md解释收益与失败。悬停/圆周/混合改善明确，但直线/爬升远段异常造成退化，实时性和1GHz高相干仍未解决，不作为全场景最终默认方案。

此前连续轨迹研究位于 `kf_dpc_coherent_5nodes/joint_gp/`：连续时间GP动力学、持久联合到达因子、导航—目标交叉相关、T-Trend重尾过程及完整20种子对照。正式结果见其RESULTS.md，结论解读见ANALYSIS.md，具体推导与正式论文来源见RESEARCH.md、T_TREND.md。研究保留退化结果，不代表已同时解决抖动、机动与高相干问题，也未将新方法设为最终有效估计器。

传感器驱动基线位于 `kf_dpc_coherent_5nodes/raw_sensor_maritime/`：双时间尺度浮标运动、共享方向谱海浪、原始GNSS伪距/IMU/磁场导航解算、原始雷达CV/IMM、DPC/CI、未来双方运动预测与实际复场评价。暂不考虑杆臂，全部传感器与单个等效发射阵元共位。

此版本补齐传感器驱动仿真链路，不将真实导航状态加噪声作为导航解算结果，也不提供真实慢快运动分量给预测器。它仍是分层近似估计，不宣称全局最优联合后验或最大相干增益。正式结果以该目录完整批次的RESULTS.md为准。

## 运行

在本目录执行：

```powershell
python -m pip install -r requirements.txt
python -u kf_dpc_coherent_5nodes/raw_sensor_maritime/run_full.py --trials 20
python kf_dpc_coherent_5nodes/raw_sensor_maritime/report.py
```

重跑会更新对应实验outputs中的结果。最新模型推导见 `raw_sensor_maritime/THEORY.md`，接口与假设见其README.md。

## 保留内容

- `kf_dpc_coherent_5nodes/joint_gp/`：最新正式文献驱动研究、全部完整与失败实验记录；复现命令见其README.md。
- `kf_dpc_coherent_5nodes/raw_sensor_maritime/`：原始传感器驱动的新实验。
- `kf_dpc_coherent_5nodes/imm_ci/`：简化固定姿态条件下IMM+DPC/CI完整对照及结果。
- `kf_dpc_coherent_5nodes/imm_dpc/`：简化固定姿态条件下CV/CT/IMM完整对照及结果。
- `kf_dpc_coherent_5nodes/gm_dpc_unified/`：原一阶GM漂移＋六维CV基线、绘图、结果及原始数据。
- `kf_dpc_coherent_5nodes/pipeline.py`：有限轮DPC公共模块。
- `position_only_velocity_demo/model.py`：坐标转换与KF公共模块。
- `dfpc_experiment/sim_core.py`：共识矩阵与共识迭代公共模块。
- `docs/`：当前模型Markdown与Word文稿。

## 版本边界

旧GM统一实验采用六维CV，不能将旧文稿中的七维CT段落当作其代码说明。后来新增的CT/IMM和CI结果分别位于imm_dpc、imm_ci。新原始传感器实验采用20个单阵元浮标，旧统一实验每浮标2×2阵元，且导航假设不同；不可将新旧功率变化解释为单一算法的改善。

所有当前发射实验仍假设理想载波频相同步，尚未验证实际通信负载与硬件实时延迟。

历史代码与旧实验产物统一归档于 `已废弃/`，仅供追溯，不参与当前实验。归档清单与恢复说明见该目录。外部备份位于 `D:/git/sim2_cleanup_backup_20261030_164626`。
