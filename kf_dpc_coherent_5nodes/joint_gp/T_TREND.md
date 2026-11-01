# T-Trend过程先验：解决稳态波动与机动之间的矛盾

## Material Passport
- 在完整Gaussian联合GP批次继续运行的同时，新增一项完整20种子研究；不截取Gaussian成功场景，不隐藏退化。
- 正式依据：Aravkin, Burke & Pillonetto (2014), Robust and Trend-Following Student’s t Kalman Smoothers, SIAM Journal on Control and Optimization52(5):2891–2916, DOI https://doi.org/10.1137/130918861 。作者保存的正式SIAM排版全文：https://sites.math.washington.edu/~burke/papers/reprints/82-robust-trend-T.pdf ，首页有卷页、版权、DOI；不是预印本。
- 补充：Roth, Özkan & Gustafsson (2013), A Student’s t Filter for Heavy Tailed Process and Measurement Noise, ICASSP5770–5774, DOI https://doi.org/10.1109/ICASSP.2013.6638770 。正式IEEE会议版：https://dihana.cps.unizar.es/proceedings/ICASSP/2013/pdfs/0005770.pdf 。

## 1. 文献机制和本项目特化

SIAM论文区分measurement-t的T-Robust与process-t的T-Trend。前者抗观测异常，后者允许真实状态快速变化。这里雷达原本按Gaussian生成，因此不使用measurement-t去压低大雷达创新，不把真实机动当作观测异常。

本项目采用一个明确的因果特化：小背景WNOJ过程加Student-t加速度改变。不是声称复现论文完整批处理算法或其收敛定理。

\[
x_k=F x_{k-1}+G z_k+w_k,\quad
w_k\sim N(0,Q_{bg}),\quad
G=\begin{bmatrix}h^2I/2\\hI\\I\end{bmatrix}.
\]

x=[p,v,a]；G z表示区间起点加速度改变z对同一区间位移、速度和加速度的共同影响。位置/速度不独立抖动。背景q_bg=0.001m²/s⁵，沿用WNOJ积分Q公式；z的各轴Student-t scale标准差0.02m/s²，ν=3。所有参数在新增正式批次开始前固定，全部轨迹/种子/海况相同，无按真值阶段切换，也无参数扫描或留出集最优声明。

## 2. 为什么能强约束稳态又允许机动

采用Gaussian scale mixture：

\[
\lambda\sim Gamma(\nu/2,\nu/2),\qquad z|\lambda\sim N(0,S_0/\lambda).
\]

Gamma第二参数为rate。积分λ得到

\[
p(z)\propto[1+z^TS_0^{-1}z/\nu]^{-(\nu+3)/2},
\quad \rho(z)=\tfrac{\nu+3}{2}\log[1+z^TS_0^{-1}z/\nu].
\]

对平方残差δ=zᵀS0⁻¹z求导，Gaussian加权形式的精度为

\[
\bar\lambda=(\nu+3)/(\nu+\delta),\qquad S_{eff}=S_0/\bar\lambda.
\]

小过程偏离受到较强约束；大的真实加速度改变连续放宽。λ始终正，不使用硬阈值、状态标签、测量门控或删帧。

注意S0是scale矩阵，不是实际协方差。ν>2时Cov(z)=νS0/(ν−2)，本实验为3S0。因此设置Gaussian_matched_covariance对照，以相同的3S0协方差、相同背景Q、相同原始观测和计算结构，区分重尾形状与简单减小噪声的效果。

## 3. 保留到达因子，避免重复计数

仍维护全部浮标当前16维误差和UAV状态的联合Gaussian到达因子及全部交叉块。每帧先传播原始IMU，在一次线性化中加入当前伪距/磁场；然后只对当前Gaussian雷达似然和过程z因子进行迭代MAP。

给定λ，显式追加z作为3维潜变量。新目标块增加G S_eff Gᵀ，与z交叉块为G S_eff，z自身协方差S_eff。雷达虽然没有直接测z，仍能通过这些交叉块估计z。

每次IRLS都从同一个冻结导航条件到达因子重建候选过程先验，不能把当前雷达似然当成多份独立信息。计算后将z按协方差主子块边际化，历史信息持续保留。

Student-t损失非凸，采用两种固定尺度初值λ=1和0.001，每个三次IRLS；以到达因子Gaussian二次项＋Student-t过程损失＋Gaussian雷达残差的共同后验目标选较小者。选择依据只有观测和模型，没有真值、无误差达标门槛。初值是优化器多起点，不是独立试验或真值选参。

## 4. 等价消元和计算压缩

本扩展只保留前一与当前两个目标克隆，其他旧目标克隆先做Gaussian边际。未来原始观测只直接依赖当前目标和当前导航，且所有交叉块被保留，因此在同一个Gaussian近似下，删除不再被直接观测的旧克隆不改变当前后验。它不是把历史信息清零。

八颗GNSS卫星的线性化伪距只依赖位置3维和钟差1维。对各浮标8×4白化Jacobian作薄QR，将八条伪距的有效残差压缩成四个正交分量；被消掉的残差只给出与状态无关的常数。在线性化Gaussian导航因子内无信息损失；磁场三分量保留。此步骤没有减少卫星、选择好测量或门控样本。

初次未压缩执行的源码和日志在adaptive_outputs/initial_uncompressed，未作为正式完整结果。所有20种子在等价消元修订后从头执行，未拼接部分数据。

## 5. 置信度限制

输出协方差是给定估计λ的条件Gaussian近似，忽略λ的不确定度和完整非Gaussian多峰后验，仍可能过度自信。未来未知加速度改变也不能由当前观测预知。即使轨迹变平滑，也必须检查真实误差、机动滞后、覆盖率及实际相干增益，不能借Student-t名称保证可靠。

NAV一次线性化与原Gaussian联合GP中的三次合并线性化不同；因此T-Trend最严格的结构对照是其Gaussian_matched_covariance，而非把相对旧IMM的全部变化归因于重尾过程。

## 6. 完整实验及复现

20种子×2海况×5轨迹×2过程方法×2载频×2时延×8时段=12800行。与原32000行Gaussian研究组合44800行，原Gaussian与扩展CSV分别保留。两研究使用完全相同原始传感器数据与真值随机数，零杆臂、总功率1W、所有节点均匀发射。

```powershell
python -u kf_dpc_coherent_5nodes/joint_gp/run_adaptive.py --trials 20
```

本次以两组不重叠的10种子并行，seed-start20261028和20261038，输出adaptive_outputs/chunks/对应种子。之后依次merge.py、merge_adaptive.py、report.py，完整结果才生成RESULTS.md。
