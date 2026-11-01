# 正式研究、数学推导与当前联合轨迹方案

## Material Passport
- 当前版本：持久高斯到达因子、WNOJ轨迹先验、随机克隆固定滞后联合MAP。
- 本页为针对性研究，不是系统性综述；只采用正式期刊/会议成果。
- 初版有限历史方案的完整文献笔记与源码在finite_history_v1/，其负面中间结果保留在outputs/finite_history_v1/。初版每次舍弃旧雷达信息导致端点波动，本版从结构上修正，不靠减小噪声参数或输出后平滑掩盖。

## 1. 已有研究如何处理这些问题

|问题|正式来源|采用的机制|不能借用的结论|
|---|---|---|---|
|匀速先验在加速下偏差|Tang, Yoon & Barfoot, IEEE RA-L 2019, 4(2):594–601, DOI10.1109/LRA.2019.2891492|WNOJ过程，显式包含加速度，限制其变化而非强迫匀速|本项目是R³特化，不是逐式复现SE(3)版本；不能保证任何q都更好|
|逐帧状态不形成稳定轨迹约束|Anderson等, Autonomous Robots 2015, 39:221–238, DOI10.1007/s10514-015-9455-y|连续时间SDE/GP运动因子、稀疏信息结构|线性高斯时与平滑等价，GP名称不是性能证明|
|目标与传感平台误差混淆|Taylor等, IPSN2006:27–33, DOI10.1145/1127777.1127785|联合概率状态、原始观测约束双方、保留交叉协方差|其传感节点静止，不能直接搬用动力学；其厘米级实测结果不属于我们|
|历史状态和当前状态被错误当独立|Roumeliotis & Burdick, ICRA2002:1788–1795, DOI10.1109/ROBOT.2002.1014801|随机克隆，保存不同时间状态及交叉协方差|不意味着任意非线性实现严格一致|
|非线性MAP与Gaussian到达因子的关系|Bell, SIAM J. Optimization1994, 4(3):626–636, DOI10.1137/0804035|迭代平滑与Gauss–Newton的联系；同一次迭代使用冻结先验|不能把它包装成与KF完全无关的新理论|
|边缘化造成虚假信息|Huai等, AAAI2021, 35(7):6084–6092, DOI10.1609/aaai.v35i7.16758|警惕边缘化/重线性化混用、重复计数和过度自信|本项目不是RI-FLS，不借用其VINS右不变一致性定理|
|位置更准仍不一定更相干|Ochiai等, IEEE TSP2005, 53(11):4110–4124, DOI10.1109/TSP.2005.857028|研究位置/相位误差对分布式波束的影响|其随机远场数组结论不等于当前动态全过程|

核实记录：RA-L作者正式发表目录https://asrl.utias.utoronto.ca/~tdb/pubs.html明确卷页和在ICRA2019展示；Autonomous Robots出版社页面可读；IPSN正式ACM会议版全文https://people.csail.mit.edu/jrb/papers/slat-ipsn-2006.pdf可读；AAAI正式页面https://ojs.aaai.org/index.php/AAAI/article/view/16758可读，正式PDF曾读到而后有机器人拦截；ICRA随机克隆由Caltech正式出版记录https://authors.library.caltech.edu/records/enrkd-zam24核实。不用同名MIT2005技术报告或arXiv替代正式出版来源。

## 2. WNOJ过程的逐式来源

定义x=[p,v,a]，p'=v、v'=a、a'=w，E[w(t)w(t')ᵀ]=qIδ(t−t')。积分给出

a(t+h)=a(t)+∫w，
v(t+h)=v(t)+ha(t)+∫(h−τ)w(τ)dτ，
p(t+h)=p(t)+hv(t)+h²a(t)/2+∫(h−τ)²w(τ)/2 dτ。

所以

\[
F(h)=\begin{bmatrix}I&hI&h^2I/2\\0&I&hI\\0&0&I\end{bmatrix},\quad
Q(h)=q\begin{bmatrix}h^5/20&h^4/8&h^3/6\\h^4/8&h^3/3&h^2/2\\h^3/6&h^2/2&h\end{bmatrix}\otimes I.
\]

例如Q_pp=q∫τ⁴/4 dτ=q h⁵/20；Q_pv=q∫τ³/2 dτ=q h⁴/8，其余同理。q单位m²/s⁵；本实验q=1、所有轴相同，没有锁Z、没有输入悬停标签。该先验允许非零加速度，并对其改变施加概率约束。固定探索性参数，没有独立留出标定，不宣称普适最优。

## 3. 持久联合到达因子与当前观测

维护状态χ=[窗口内的UAV p/v/a克隆，全部20浮标当前16维导航误差]。导航误差为[δp,δv,δθ,δb_a,δb_g,δc_clock]。保留整个均值与协方差，包括所有历史目标—当前导航、不同浮标之间的块。

每帧首先按GP与原始IMU推进得到到达因子N(χ^-,P^-)。仅加入当前的伪距、磁场和雷达似然，求

\[
\hat χ=\arg\min_χ\tfrac12\|χ-χ^-\|^2_{(P^-)^{-1}}
+\tfrac12\|z_k-h_k(χ)\|^2_{R_k^{-1}}.
\]

不是将含有窗口观测的上轮后验当先验后再重复添加窗口观测。每条原始观测只在其当前历元进入一次。

当前观测白化后R=I，在第j次线性化r_j=z−h(χ_j)、H_j=∂h/∂χ，令S_j=I+H_jP^-H_jᵀ、K_j=P^-H_jᵀS_j^-1，则Gauss–Newton等价更新为

\[
χ_{j+1}=χ^-+K_j[r_j+H_j(χ_j-χ^-)].
\]

三次迭代始终冻结同一个χ^-和P^-，不会将当前测量当三份独立信息。最后在线性化终点只更新一次协方差P^+=P^-−KHP^-，数值上对称化。该P是局部Gaussian/Laplace近似，不是严格非线性可信保证。

## 4. 随机克隆如何保留历史信息

追加目标状态x_new=F x_last+w，必须同时追加所有交叉块：P_new,new=F P_last,last Fᵀ+Q，P_new,rest=F P_last,rest。

目标窗口最多11个历元（2s），超窗只删除最旧目标变量。若联合后验是N([m_old,m_keep],P)，保留变量的边际恰为N(m_keep,P_keep,keep)。**取协方差主子块是边际化，不是将历史信息清零**。信息形式才需Schur补H_keep−H_keep,old H_old,old^-1 H_old,keep。

旧信息持续保留在到达因子中；旧非线性原始因子不再重新线性化。因此这是增广状态迭代MAP/随机克隆固定滞后估计，不是每次对全历史原始非线性因子重新优化，也不是RI-FLS。线性高斯特例与增广Kalman平滑等价，这一点必须明确承认。

## 5. 参考ESKF不是第二份导航观测

原始伪距/IMU/地磁ESKF只生成误差坐标参考。导航过程噪声从原始IMU白噪声与偏置随机游走累积得到，**不把参考ESKF当前后验P再当联合先验或过程噪声**。

令参考预测R_p、参考更新R_f=R_p Exp(d_θ)，新误差R_true=R_f Exp(δ_cur)，则

Log[Exp(d_θ)Exp(δ_cur)]≈d_θ+J_r(d_θ)^−1δ_cur。

令T姿态块为J_r^-1(d_θ)≈I+[d_θ]×/2+[d_θ]×²/12，其余为I，坐标转换为

\[
δ_{cur}=T^{-1}(F_{imu}δ_{prev}-d)+\xi,\qquad
Q_{cur}=T^{-1}Q_{imu}T^{-\mathsf T}.
\]

均值、协方差及所有交叉块用同一变换。当前参考状态依赖GNSS不等于重复似然：它仅是已知坐标重心；联合先验仍从上帧联合后验和IMU产生，再加当前原始GNSS/磁场似然一次。首帧伪距/磁场已用于初始化，联合首帧不再次加入这两类测量，只更新雷达。

原始IMU累计Q包含加速度密度平方乘Δt³/3、Δt²/2、Δt，角速率密度平方乘Δt，以及两类偏置、钟差随机游走。误差F用二阶矩阵近似，姿态坐标变换为小角近似；不能由此宣称原始非线性后验精确。

## 6. 原始观测如何约束双方

候选浮标p_b=p_ref+δp，R=R_ref Exp(δθ)。伪距预测||s−p_b||+c_ref+δc；其位置Jacobian为(p_b−s)ᵀ/d，钟差Jacobian为1。磁场预测m=Rᵀm_world，姿态Jacobian为[m]×J_r(δθ)。

雷达q=Rᵀ(p_u−p_b)，h=spherical(q)。令B=∂spherical(q)/∂q，原始残差对UAV和导航的Jacobian分别为

\[
H_u=[BR^T,0,0],\quad H_n=[-BR^T,0,B[q]_\times J_r(δθ),0,0,0].
\]

导航作为显式未知量后，雷达R只含原始测量噪声，不能再加J P_nav Jᵀ重复计入。当前联合更新会通过协方差交叉块修正浮标速度和偏置。

GP_independent_nav是受控近似对照：同样持久GP和全网雷达，但用R_eff=R_radar+J P_nav Jᵀ近似每帧独立导航误差。用于区分全网集中观测/运动先验的贡献与联合相关性的贡献。

## 7. 未来相对几何与相位不确定度

目标按p_u(h)=p_u+h v_u+h²a_u/2预测。浮标有受控对照：p_b+h v_b，或p_b+h v_b+h²a_b/2，其中a_b=R(f_m−b_a)+g，全部来自当前测量和联合估计，无真实波浪分解输入。

在联合后验上计算各浮标的相对[p,v,a]，保留完整N×9×N×9协方差Σ_rel。a_b对姿态/偏置的Jacobian为−R[f_m−b_a]×J_r(δθ)、−R。当前加速度测量噪声及未来未知运动仍近似处理，与同一IMU样本的状态相关性未完全建模，未来概率必须经验评价。

令u_i为预测LOS、G_i=[u_iᵀ,h u_iᵀ,h²u_iᵀ/2]，则距离协方差C_ij≈G_iΣ_rel,ij G_jᵀ，并增加未来共享UAV jerk噪声和独立浮标过程噪声。保留交叉项是为了计算Var(r_i−r_j)=C_ii+C_jj−2C_ij，而非只看单节点位置方差。

## 8. 相干目标与期望功率策略

均匀功率主对照完全保留固定总功率1W和所有20节点，距离命令d_hat=||p_u_hat−p_b_hat||。真实传播距离仅供评价，不进入命令。

附加功率策略在给定Gaussian距离误差和预测幅度时构造

\[
M_{ij}=\hat a_i\hat a_j\exp\{-k^2(C_{ii}+C_{jj}-2C_{ij})/2\}.
\]

令时隙均值为Mbar、||w||²=1，则预测期望接收功率≈P_tot wᵀMbar w，Rayleigh商最优w为主特征向量。M非负可取非负向量，整个时隙用一个功率分配。没有门控、删帧或真值选节点。

这只是在当前Gaussian近似及固定总功率下最优的期望接收功率，不保证真实全局最优，也不保证全网高相干。硬件若要求逐节点固定功率，必须使用均匀主对照。另报G=P_rx/P_noncoherent、统一均匀功率理想参考η和有效节点数1/Σw_i⁴，防止单节点功率集中被误称为全网相干。

## 9. 完整实验、波动与时延

20种子、两海况、五轨迹、五方法、两载频、两种时延评价，共32000行指标。主统计完整0–90s，无样本门控；启动、加速、刹停分段只解释现象。

增加位置误差变化率RMS(||(e_p,k−e_p,k−1)/dt||)、速度误差变化率，以及位置差分与相邻速度平均的一致性缺陷。前两者针对真实误差变化，不会只因真实轨迹有机动就判抖动。导数统计使用同阶段相邻点，是差分定义，不删原始位置样本。

固定时延采用12ms处理+125ms通信+40ms时隙；measured_cpu用实际求解加导航平均耗时替代12ms下限。它是当前电脑并行负载下的诊断，未计全预测、编解码与带宽，不是硬件实时保证。若末端超过200ms帧周期，则需额外调度设计，不删除这些样本。

新方法是集中联合信息参考，不宣称完成分布式GP。全网原始观测共享的代价未实测，不能把对有限轮CI的全部收益归于GP。参数固定、无留出标定、有限阶线性化和Gaussian到达因子均限制普适性。NEES与覆盖率据结果报告，不借用论文定理声称自动可信。

## 正式参考文献

1. Tang, T. Y., Yoon, D. J., & Barfoot, T. D. (2019). A White-Noise-on-Jerk Motion Prior for Continuous-Time Trajectory Estimation on SE(3). IEEE Robotics and Automation Letters, 4(2), 594–601. https://doi.org/10.1109/LRA.2019.2891492
2. Anderson, S., Barfoot, T. D., Tong, C. H., & Särkkä, S. (2015). Batch nonlinear continuous-time trajectory estimation as exactly sparse Gaussian process regression. Autonomous Robots, 39, 221–238. https://doi.org/10.1007/s10514-015-9455-y
3. Taylor, C., Rahimi, A., Bachrach, J., Shrobe, H., & Grue, A. (2006). Simultaneous localization, calibration, and tracking in an ad hoc sensor network. IPSN, 27–33. https://doi.org/10.1145/1127777.1127785
4. Roumeliotis, S. I., & Burdick, J. W. (2002). Stochastic cloning: a generalized framework for processing relative state measurements. ICRA, 1788–1795. https://doi.org/10.1109/ROBOT.2002.1014801
5. Bell, B. M. (1994). The Iterated Kalman Smoother as a Gauss–Newton Method. SIAM Journal on Optimization, 4(3), 626–636. https://doi.org/10.1137/0804035
6. Huai, J., Lin, Y., Zhuang, Y., & Shi, M. (2021). Consistent Right-Invariant Fixed-Lag Smoother with Application to Visual Inertial SLAM. AAAI, 35(7), 6084–6092. https://doi.org/10.1609/aaai.v35i7.16758
7. Ochiai, H., Mitran, P., Poor, H. V., & Tarokh, V. (2005). Collaborative beamforming for distributed wireless ad hoc sensor networks. IEEE Transactions on Signal Processing, 53(11), 4110–4124. https://doi.org/10.1109/TSP.2005.857028
