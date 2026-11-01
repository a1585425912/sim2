# 仅位置观测通过KF增益修正速度：文献与推导

日期：2026-10-31。此次为文献核对与数学分析，没有新增数值实验、冒烟测试或门限筛选。以下推导针对明确的线性模型，不把推导出的公式冒称为未取得全文论文的原公式。

## 已正式发表的直接研究

1. B. Friedland, “Optimum Steady-State Position and Velocity Estimation Using Noisy Sampled Position Data,” IEEE Transactions on Aerospace and Electronic Systems, AES-9(6), 906–911, 1973. DOI: https://doi.org/10.1109/TAES.1973.309666 。作者所在NJIT的正式出版记录与摘要：https://researchwith.njit.edu/en/publications/optimum-steady-state-position-and-velocity-estimation-using-noisy/ 。直接研究带噪采样位置数据下的KF位置、速度估计精度，采用一维分段恒加速度模型。此次核对摘要与正式元数据，未取得该文完整正文；不采用未经核对的具体闭式系数。
2. P. R. Kalata, “The Tracking Index: A Generalized Parameter for α-β and α-β-γ Target Trackers,” IEEE Transactions on Aerospace and Electronic Systems, AES-20(2), 174–182, 1984. DOI: https://doi.org/10.1109/TAES.1984.310438 。经典稳态跟踪增益研究。出版社全文访问受限，不虚构定理号。
3. R. A. Singer, “Estimating Optimal Tracking Filter Performance for Manned Maneuvering Targets,” IEEE Transactions on Aerospace and Electronic Systems, AES-6(4), 473–483, 1970. DOI: https://doi.org/10.1109/TAES.1970.310128 。机动动态建模与跟踪性能的正式期刊基础。
4. G. Revach, N. Shlezinger, X. Ni, A. López Escoriza, R. J. G. van Sloun, and Y. C. Eldar, “KalmanNet: Neural Network Aided Kalman Filtering for Partially Known Dynamics,” IEEE Transactions on Signal Processing, 70, 1532–1547, 2022. DOI: https://doi.org/10.1109/TSP.2022.3158588 。已读取正式出版PDF，其封面明确标注Publisher’s PDF/Version of Record：https://pure.tue.nl/ws/files/201346390/KalmanNet_Neural_Network_Aided_Kalman_Filtering_for_Partially_Known_Dynamics.pdf 。该文把增益计算交给RNN学习，保留预测—创新修正结构，研究模型与噪声知识不完整时的状态估计。它不是专门验证本项目位置雷达速度估计或相干发射的论文，不能移植其性能数字。

本文不使用预印本作为论据。KF速度修正机制是成熟理论，本身不是可宣称的新算法。

## 一维推导：速度修正从何而来

设x=[p,v]^T，采样间隔T>0，F=[[1,T],[0,1]]，H=[1,0]。
位置观测z=p+n，位置观测噪声方差R>0，与预测误差独立。
预测误差协方差分块为P^-=[[Ppp,Ppv],[Pvp,Pvv]]。

S=H P^- H^T+R=Ppp+R。
K=P^- H^T S^-1=[Ppp/(Ppp+R), Pvp/(Ppp+R)]^T。

记位置创新r=z−p^-，则
p^+=p^-+Kp r，v^+=v^-+Kv r。

所以Δv=Kv r，而不是“速度等于卡尔曼增益”。Kv单位1/s，r单位m，乘积单位m/s。
关键是Pvp是预测速度误差与预测位置误差的交叉协方差，绝非额外观测到的真实速度。

即使上一帧Pvp^+=0，经预测仍有
Pvp^-=Pvp^+ + T Pvv^+ + Qvp。

当前连续白加速度CV离散模型的Qvp=q T²/2。因此只要速度存在不确定度，预测便把它传播到位置不确定度，使位置创新能够修正速度。这不需要预先知道真值速度。

## 三维对应

x=[p^T,v^T]^T，H=[I,0]，P^-按位置/速度分块，则

K_v=P_vp^- (P_pp^-+R)^-1，
v^+=v^-+K_v(z−p^-)。

三维Kv是3×3矩阵；协方差非对角时，各位置创新分量均可能影响某个速度分量。当前position_only_velocity_demo/model.py的_filter中，gain的后三行执行的就是这一更新；IMM中每个分支也保留相同的创新修正结构，只是预测模型及协方差不同。

## 可观测性与估计精度不是同一回事

一维可观测矩阵O=[H;HF]=[[1,0],[1,T]]，det O=T≠0。三维对应[[I,0],[I,T I]]满秩。
因此理想CV模型下，位置时间序列包含恢复速度的信息；不是单帧位置就能精确测速，也不是两帧带噪观测足以获得高精度。

在线性高斯且模型/噪声正确时，速度后验协方差满足

P_vv^+=P_vv^-−P_vp^- S^-1 P_pv^-。

被减项半正定，说明使用位置观测会减小该模型下的速度不确定度。它不意味着任意一次实际速度误差必然下降，也不保证模型失配时报告协方差覆盖真实误差。

## 与α–β跟踪器的联系

在一维定采样、稳态增益等条件下，令Kp=α、Kv=β/T，可写成

p^+=p^-+αr，v^+=v^-+(β/T)r。

这与经典α–β跟踪器相同形式。不是任何人为选择的α、β都等价于某个指定Q/R下的最优KF；需要满足对应模型的Riccati关系。

## 增益加大为什么不能消除模型失配

考虑稳定的固定增益α–β滤波器，真值有恒定加速度a，而预测仍用CV。
定义后验误差ep=p_true−p_hat，ev=v_true−v_hat。忽略测量噪声时：

r_k=ep_(k−1)+T ev_(k−1)+aT²/2，
ep_k=(1−α)r_k，
ev_k=ev_(k−1)+aT−(β/T)r_k。

若误差进入稳态，由第三式得到r*=aT²/β；再代入前两式：

ep*=(1−α)aT²/β，
ev*=aT(α/β−1/2)。

这表明持续加速时，即便滤波器稳定，位置与速度仍可能存在系统滞后。增益选择会改变滞后，却不自动把CV变为正确的加速度模型。以上是固定增益、恒加速度、稳定递推条件下的本项目推导，不是一般非线性IMM的精确误差公式。

## 仅位置测速的精度参照

对严格恒速、独立等方差位置噪声R、n个等间隔时刻t_i=iT，普通最小二乘速度斜率方差为

Var(v_hat)=R / Σ(t_i−t_bar)²=12R/[T² n(n²−1)]。

这是有限窗口恒速回归的解析参照，不能直接当作机动KF或当前IMM的误差下界。真实机动时延长窗口会增加模型失配；直接位置差分则有噪声方差2R/T²。速度估计精度必须结合时间跨度、采样、观测噪声和动态假设一起讨论。

## 对当前项目的结论

- 位置创新乘速度增益修正速度，具有标准统计理论基础，且存在直接研究同类问题的IEEE正式论文。
- 当前精度不足不能据此归咎于这个更新公式无效；应分析状态模型、Q/R、观测转换及预测时延。
- CA/Singer/IMM等路线主要改进预测及其不确定度，而不必抛弃创新修正机制。
- 学习增益是已有现代研究方向，但需要训练与独立评估；不因为有KalmanNet就假设它一定改善当前UAV跟踪，更不能直接套用它的误差或一致性结论。
- 论文表述宜为“基于位置观测的状态联合估计，利用位置—速度误差交叉协方差，通过创新修正未直接测量的速度”，而非“从卡尔曼增益中提取速度”。
