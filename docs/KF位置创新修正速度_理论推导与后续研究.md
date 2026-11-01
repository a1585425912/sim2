# 从位置创新修正速度：理论推导、局限与后续研究

日期：2026-10-30。研究对象是当前仿真的位置观测KF及其IMM分支。本文只做文献梳理和数学推导，不修改算法、不运行新的数值实验、不设置门限或筛选样本。

## Material Passport
- 阶段：理论分析。
- 核心问题：没有直接速度测量时，位置创新经卡尔曼增益修正速度是否成立，以及如何改进机动时的速度估计。
- 区分：下面可观测性、增益、误差下界和固定增益机动偏差是按明确假设展开的推导；项目现有实验只用于说明应用背景，不冒充这些结论的证明。
- 来源政策：仅引用正式发表期刊论文，不以预印本作为依据。部分IEEE/ASME正文受到访问限制，不编造原文定理号或声称逐式核对。

## 1. 速度修正量究竟是什么

状态与位置观测为

\[
x_k=\begin{bmatrix}p_k\\v_k\end{bmatrix},\quad
z_k=Hx_k+\epsilon_k,\quad H=[I_3\quad0],\quad
F=\begin{bmatrix}I_3&hI_3\\0&I_3\end{bmatrix}.
\]

预测误差协方差分块为

\[
P_k^-=\begin{bmatrix}P_{pp,k}^-&P_{pv,k}^-\\P_{vp,k}^-&P_{vv,k}^-\end{bmatrix}.
\]

创新与创新协方差：

\[
\nu_k=z_k-\hat p_k^-,\quad S_k=P_{pp,k}^-+R_k.
\]

于是

\[
K_k=P_k^-H^TS_k^{-1}=\begin{bmatrix}P_{pp,k}^-S_k^{-1}\\P_{vp,k}^-S_k^{-1}\end{bmatrix}.
\]

速度部分是

\[
\boxed{K_{v,k}=P_{vp,k}^-S_k^{-1}},\quad
\boxed{\hat v_k^+=\hat v_k^-+K_{v,k}\nu_k}.
\]

这里\(\Delta\hat v_k=K_{v,k}\nu_k\)是对预测速度估计的修正，不是真实速度、不是直接测速结果，也不是物理上的\(v_k-v_{k-1}\)。

## 2. 为什么速度没有测到，却能被修正

CV预测给出

\[
P_{vp,k}^-=P_{vp,k-1}^++hP_{vv,k-1}^++Q_{vp,k}.
\]

因为过去速度的不确定性会在一个采样间隔后造成位置不确定性，速度误差与下一帧的位置误差产生相关。位置创新因而包含关于速度误差的间接信息。

进一步令真实速度预测误差\(e_v=v-\hat v^-\)，寻找线性修正\(L\nu\)使

\[
J(L)=E[\|e_v-L\nu\|^2]
=\operatorname{tr}P_{vv}^- -2\operatorname{tr}(LP_{pv}^-)+\operatorname{tr}(LSL^T)
\]

最小。对L求导得

\[
\frac{\partial J}{\partial L}=-2P_{vp}^-+2LS=0,
\]

因此\(L=P_{vp}^-S^{-1}=K_v\)。这说明速度增益是位置创新对速度误差的最优线性回归系数。在正确线性高斯模型下，这也等于联合高斯条件均值修正。

相应的理论速度协方差更新为

\[
\boxed{P_{vv}^+=P_{vv}^- -P_{vp}^-S^{-1}P_{pv}^-}.
\]

减去的是半正定矩阵，说明有效位置观测可以降低理论速度不确定度。但当运动模型或噪声统计不匹配时，该理论协方差不等于真实误差二阶矩，不能据此保证实际速度RMSE下降。

## 3. 增益的量纲与启动抖动

\(P_{vp}\)量纲为m²/s，S为m²，因此\(K_v\)量纲为1/s。它不是0到1之间的无量纲权重，数值大于1并不说明实现错误。

单轴示例（不是新增实验或实际节点测量）：上帧\(P_{pp}^+=1\,m^2\)、\(P_{vv}^+=400\,(m/s)^2\)、\(P_{vp}^+=0\)，取h=0.2s、q=0.5m²/s³、R=1m²，则

\[
P_{vp}^-=0.2\times400+\tfrac12q h^2=80.01,
\]

\[
S=1+0.2^2\times400+qh^3/3+1\approx18.00133,
\quad K_v\approx4.445\,s^{-1}.
\]

所以1m位置创新可以造成约4.445m/s的速度修正。初速均值为零、标准差为20m/s是“允许速度很快被新信息修正”的先验，不是“确信UAV静止”。当前代码采用这类初速方差，启动时速度跳动有明确的统计来源；不能简单通过设初速方差为零来压平曲线。

## 4. 可观测不等于任意精度

对位置观测CV系统，前两帧可观测矩阵为

\[
\mathcal O_2=\begin{bmatrix}H\\HF\end{bmatrix}
=\begin{bmatrix}I_3&0\\I_3&hI_3\end{bmatrix}.
\]

只要h≠0，它的秩为6。这是在给定CV模型下，位置历史可以确定初始位置和速度的结构性结论，不代表两帧带噪观测就能精确测速，也不代表任意未建模加速度可由两帧恢复。

单轴、严格匀速、n帧等间隔、独立同方差位置噪声的简化模型为

\[
z_i=p_0+v t_i+\epsilon_i,\quad t_i=ih,\quad\epsilon_i\sim\mathcal N(0,\sigma_p^2).
\]

令\(\theta=[p_0,v]^T\)。Fisher信息矩阵及速度方差下界为

\[
\mathcal I(\theta)=\frac1{\sigma_p^2}
\begin{bmatrix}n&\sum_i t_i\\\sum_i t_i&\sum_i t_i^2\end{bmatrix},
\]

\[
\operatorname{Var}(\hat v)\ge[\mathcal I^{-1}]_{22}
=\frac{\sigma_p^2}{\sum_i(t_i-\bar t)^2}
=\boxed{\frac{12\sigma_p^2}{h^2n(n^2-1)}}.
\]

这个特殊模型下无先验加权的最小二乘斜率可以达到该界。机动、异方差、时间相关导航误差和转换观测偏置会改变结论，不能将该界直接当作本项目的性能保证。

它揭示的取舍是：更多时间历史可提高匀速测速精度，但历史过长时，真实加速或转弯会破坏匀速假设。

## 5. 机动时为什么存在滞后：从固定增益推导

考虑单轴、固定增益的CV校正，写成α–β形式：

\[
\hat p_k^+=\hat p_k^-+\alpha\nu_k,\quad
\hat v_k^+=\hat v_k^-+\frac\beta h\nu_k.
\]

定义误差为真实值减估计值\(e_p=p-\hat p^+,e_v=v-\hat v^+\)。假设真实加速度a恒定、测量噪声均值为零、滤波器稳定且已进入固定误差均值状态，则

\[
e_p^- =e_p+h e_v+\tfrac12 a h^2,
\quad e_p=(1-\alpha)e_p^-,
\]

\[
e_v=e_v+ah-\frac\beta h e_p^-.
\]

消元得到

\[
\boxed{E[e_p^+]=\frac{1-\alpha}{\beta}a h^2},\qquad
\boxed{E[e_v^+]=ah\left(\frac\alpha\beta-\frac12\right)}.
\]

这是明确假设下的机动偏差，不是所有时变KF的通式。但它直接说明：即使速度一直通过创新修正，错误的匀速预测结构仍可能保留持续速度偏差；不能认为“只要有K_v就不会滞后”。

增大增益可能减小响应滞后，也会将更多测量噪声带进速度；噪声部分通过\(K_v R K_v^T\)进入Joseph协方差更新。仅乘一个人工放大系数，无法兼顾所有工况，原协方差递推也需要相应重算。

## 6. 已发表研究推进了什么

### A. 加速度有记忆，而不是只用CV白噪声承受机动

Singer类模型将加速度作为状态，并给它有限相关时间。以单轴连续模型为例：

\[
\dot p=v,\quad \dot v=a,\quad
da=-\alpha a\,dt+\sqrt{2\alpha\sigma_a^2}\,dW.
\]

其条件均值转移直接积分得到

\[
E[a(t+h)]=e^{-\alpha h}a(t),
\]

\[
E[v(t+h)]=v(t)+\frac{1-e^{-\alpha h}}\alpha a(t),
\]

\[
E[p(t+h)]=p(t)+hv(t)+
\frac{\alpha h-1+e^{-\alpha h}}{\alpha^2}a(t).
\]

仍然使用位置创新，通过增益分别修正p、v、a，但速度预测已经包含估计加速度。恒加速度CA则是另一个明确假设的候选，不应把两者混为同一模型。

### B. 从数据辨识Q、R，进而改变K_v

K_v依赖P与R，P又依赖F与Q。自适应滤波及ALS研究通过创新的自协方差约束辨识噪声协方差；不是通过给速度增益任意加倍来实现自适应。

在正确线性高斯滤波模型下，创新满足

\[
E[\nu_k\nu_k^T]=HP_k^-H^T+R,\qquad
E[\nu_k\nu_{k-\ell}^T]=0\quad(\ell>0).
\]

辨识需要模型、噪声结构及可辨识性条件。只用零时延创新方差的一条关系，通常不足以同时唯一确定所有未知Q、R；持续机动造成的偏差也不能一律当作测量噪声去增大R。

Odelson等（Automatica 2006）研究了受半正定约束的ALS与可辨识条件；Kong（Automatica 2023）进一步讨论已知噪声作用矩阵结构下的辨识问题。它们的线性/稳定模型假设不能未经修改直接套在当前非线性、切换IMM上。

### C. 多模型与分布式多模型

各模式保留自己的F、Q、P和增益，用观测似然更新模式概率。进一步的分布式多模型方法保留模式条件后验及模式概率交换，不只是先合并为一个高斯再做CI。

因此，本地估计中的错误机动模型和下游融合不足是两个不同问题。已有CI实验表明，仅改变输出融合并没有消除加速/刹停段的问题，不能把这项结果理解成“位置创新无法估计速度”。

## 7. 对无人机相干合成意味着什么

仅看发射延迟h中的速度误差项，单程距离残差一阶为

\[
\delta d_i^{(v)}\approx h u_i^T\delta v_i,
\quad \delta\phi_i^{(v)}\approx\frac{2\pi h}{\lambda}u_i^T\delta v_i.
\]

真正影响阵元间合成的是\(\delta\phi_i-\delta\phi_j\)。沿某阵元LOS的速度误差若为0.5m/s，h=137ms，则该项距离误差为0.0685m；300MHz下约24.7°，1GHz下约82.3°。这只是单阵元速度项的条件数值例，不是全网相位RMSE实测。

所以研究目标应是：在有效预测时长内，降低LOS方向及节点间的速度预测误差，并校准其不确定度；不只是把\(\Delta\hat v\)变大或让速度曲线更平滑。

## 8. 针对当前项目的结论与优先级

1. 保留\(\Delta\hat v=K_v\nu\)这一有理论基础的校正步骤；目前没有证据表明需要放弃它。
2. 现阶段更值得改的是状态预测：增加CA或Singer机动分支，并区分水平/垂直模式，避免水平机动控制全部XYZ的过程噪声与模式权重。
3. 固定比例放大K_v不是优先方案；它会改变估计器，必须重算协方差，也不能消除错误模型带来的系统性偏差。
4. 从独立数据校准观测协方差与结构化Q，可利用创新相关性研究，但需要先满足相应可辨识假设，不能只看残差大就盲目改噪声。
5. 要定位现有速度误差，可在以后完整实验中保存\(\hat v^-\)、\(K_v\)、\(\nu\)、\(K_v\nu\)、\(\hat v^+\)和模型概率，将预测误差、修正不足及噪声放大区分开。本次没有新增实验或对性能作出新承诺。

## 正式期刊参考文献

- R. E. Kalman, “A New Approach to Linear Filtering and Prediction Problems,” Journal of Basic Engineering, 82(1), 35–45, 1960. [DOI](https://doi.org/10.1115/1.3662552)。基础滤波出处；本次出版社正文访问受限。
- R. A. Singer, “Estimating Optimal Tracking Filter Performance for Manned Maneuvering Targets,” IEEE Transactions on Aerospace and Electronic Systems, AES-6(4), 473–483, 1970. [DOI](https://doi.org/10.1109/TAES.1970.310128)。相关加速度模型的经典出处；未声称逐式核对原文。
- R. K. Mehra, “On the identification of variances and adaptive Kalman filtering,” IEEE Transactions on Automatic Control, 15(2), 175–184, 1970. [DOI](https://doi.org/10.1109/TAC.1970.1099422)。噪声辨识与自适应KF的经典出处。
- X. R. Li and V. P. Jilkov, “Survey of Maneuvering Target Tracking. Part I: Dynamic Models,” IEEE Transactions on Aerospace and Electronic Systems, 39(4), 1333–1364, 2003. [DOI](https://doi.org/10.1109/TAES.2003.1261132)。运动模型综述，不把其综述视为本仿真的性能证明。
- B. J. Odelson, M. R. Rajamani, and J. B. Rawlings, “A new autocovariance least-squares method for estimating noise covariances,” Automatica, 42(2), 303–308, 2006. [DOI](https://doi.org/10.1016/j.automatica.2005.09.006)。出版社索引提供摘要及方法/结论节选。
- H. Kong, “Noise covariance estimation via autocovariance least-squares with deadbeat filters,” Automatica, 153, 111049, 2023. [DOI](https://doi.org/10.1016/j.automatica.2023.111049)。依据正式出版摘要说明研究方向，不将未取得的完整证明移植到IMM。
- G. Battistelli, L. Chisci, C. Fantacci, A. Farina, and A. Graziano, “Consensus-based multiple-model Bayesian filtering for distributed tracking,” IET Radar, Sonar & Navigation, 9(4), 401–410, 2015. [正式全文](https://doi.org/10.1049/iet-rsn.2014.0071)。作为分布式多模型方向的补充文献，主理论依据仍为以上IEEE/Automatica正式论文。
