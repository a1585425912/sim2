# IMM首版：近悬停S / CV / 未知转弯率CT

这是探索性新实验，保留原GM统一实验为基线。运行入口：

```powershell
python -B -m unittest discover -s kf_dpc_coherent_5nodes/imm_dpc -p test_imm.py -v
python -B kf_dpc_coherent_5nodes/imm_dpc/run_experiment.py --trials 20
python -B kf_dpc_coherent_5nodes/imm_dpc/report.py
```

工作目录为仓库根目录。输出见本目录outputs；结论、假设和限制见RESULTS.md。
完整理论设计见../../docs/UAV多运动模型_IMM与DPC接入设计.md。

三分支并行预测与更新，观测似然修正概率；交互混合保留模型间不确定度。
预测先按分支外推再按概率加权，随后有限轮DPC平均同一未来时刻的位置/速度。
对悬停、直线、圆周及连续混合运动做相同随机样本的配对比较。

## 实施契约

- 授权：用户明确要求实现首版并运行讲解。
- 决策：新增隔离实验，不重写旧模块；依赖只读复用。
- 范围：本目录5个实现/测试/说明文件及结果产物；理论设计与旧结果不修改。
- 门槛：零转弯率极限、协方差、概率归一、前缀因果性、无真值输入、旧CV基线复现。
- 回退：使用原gm_dpc_unified；新目录无全局配置或依赖变更。
- 首版固定参数，不以测试结果反复调参；不承诺IMM在所有工况优于单模型。
