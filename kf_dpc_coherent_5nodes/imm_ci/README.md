# 同一IMM后验：DPC与CI融合对照

本实验只改变融合规则，不调整IMM参数。CI指协方差交集的信息域凸融合。

同一时刻每节点以六维边际状态x与协方差P形成Y=P^-1、y=Yx。
25轮Metropolis邻居交换后，Y_i'=sum_j (W^25)_ij Y_j，y_i'=sum_j (W^25)_ij y_j。
P_i'=Y_i'^-1，x_i'=P_i' y_i'。权重和为1，不乘节点数，不做独立信息累加。

实现使用有限轮等效矩阵而非全网平均，单测核对与逐轮邻居交换一致。
固定拓扑权重，并非逐轮优化logdet/trace的CI。与归档ci.py的定义一致，但不导入历史代码。

三组：local（同样延迟的本地IMM）、dpc（均值共识）、ci（信息共识）。
同样20节点环网、25轮、137ms延迟及40ms时隙；无共识结果回灌。
对各IMM分支先传播均值/协方差，再做模型混合；保留模型间均值散布。
先边际化转弯率再取六维信息矩阵，不把转弯率已知作为条件。

```powershell
python -B -m unittest discover -s kf_dpc_coherent_5nodes/imm_ci -p test_ci.py -v
python -B kf_dpc_coherent_5nodes/imm_ci/run_experiment.py --trials 20
python -B kf_dpc_coherent_5nodes/imm_ci/report.py
```

## 实施契约与边界

- 用户明确授权CI实验；范围为新imm_ci目录和旧fit的可选分支协方差输出（无数值递推修改）。
- 旧结果不覆盖，首种子DPC均值和两载频功率曲线与旧IMM逐项对照。
- 新实验包含原四轨迹及平滑爬升。仅用观测、导航和本地后验融合，真值仅用于仿真及评价。
- CI一致性以输入协方差可信为前提。IMM/EKF及转换观测协方差是近似，不宣称严格误差覆盖。
- 未来177ms内冻结模型，沿用既有近似。固定时延比较隔离算法差异，尚未验证增加的协方差消息带宽与硬件耗时。
- 回退：继续使用旧IMM+DPC；没有改变全局配置或新增依赖。
- 新增算法、测试、运行、报告及说明约700行以内；不重写旧估计器。

理论参考：Julier & Uhlmann (1997), A non-divergent estimation algorithm in the presence of unknown correlations, https://doi.org/10.1109/ACC.1997.609105 。
MathWorks CI说明：https://www.mathworks.com/help/fusion/ref/fusecovint.html （该产品使用优化权重，与本实验固定权重不同）。
