# 从正式研究到本项目方案：证据、推导与边界

## Material Passport
- 阶段：针对已观察到的轨迹抖动、机动滞后、导航相关误差和低相干收益，检索正式成果，设计并实现新的完整配对实验。
- 本页不是系统性综述，不声称覆盖全部研究。只采用正式期刊/会议文献，不引用预印本作为依据。
- 科学问题：能否通过连续时间运动先验和相关导航—目标联合后验，提高在线轨迹预测和相干合成，而非仅让曲线视觉平滑？

## 一、检索到的研究实际解决了什么

|正式研究|已核实内容|本项目采用|不能直接借用的结论|
|---|---|---|---|
|Tang, Yoon & Barfoot, RA-L 2019, 4(2):594–601, DOI 10.1109/LRA.2019.2891492|WNOA匀速先验在非零加速轨迹下会有偏差，提出WNOJ含加速度的连续时间先验；作者正式目录注明在ICRA2019展示|UAV状态从p/v扩展为p/v/a，采用欧氏空间WNOJ GP过程|不是直接复现论文SE(3)算法，也不能凭论文保证本项目更好；q太小仍会滞后|
|Anderson, Barfoot, Tong & Särkkä, Autonomous Robots 2015, 39:221–238, DOI 10.1007/s10514-015-9455-y|白噪声驱动SDE先验产生稀疏GP精度结构；非线性轨迹可联合优化；线性高斯时与经典平滑等价|因果有限窗口MAP，原始雷达残差加轨迹过程残差，Schur结构求解|GP名称不是优于KF的证明，不能使用未来数据后声称在线预测|
|Taylor等, IPSN 2006:27–33, DOI 10.1145/1127777.1127785|正式会议LaSLAT把传感器位置/校准与目标轨迹放入联合概率状态，MAP/Laplace近似，Schur消元|将浮标导航误差与UAV轨迹同时作为未知变量，保留交叉块|其传感器是静态节点，不能直接移植为动态浮标；其厘米级实测性能不是本项目证据|
|Huai等, AAAI 2021, 35(7):6084–6092, DOI 10.1609/aaai.v35i7.16758|边缘化和重新线性化可能产生虚假可观信息；在VINS中给出右不变FLS一致性分析|重视不重复计数、线性化参考和经验覆盖率；避免把同一窗口后验重复当先验|本项目有GNSS与地磁，不声称实现RI-FLS，也不直接套用其一致性定理|
|Ochiai等, IEEE TSP 2005, 53(11):4110–4124, DOI 10.1109/TSP.2005.857028|分布式阵列波束图、位置/相位误差的影响|用成对距离误差而非单独坐标RMSE评价合成|其随机远场阵列假设不等于当前动态近场全过程|

核实方式：RA-L由作者正式发表目录与正式DOI核对，未把作者稿冒称出版社排版全文；Autonomous Robots由出版社摘要核实；IPSN阅读全文并确认ACM版权及会议日期；AAAI由正式页面及可访问的正式论文核实（后续PDF请求有机器人拦截）。

主要原始记录：

- https://asrl.utias.utoronto.ca/~tdb/pubs.html （RA-L正式条目，非其arXiv条目）
- https://link.springer.com/article/10.1007/s10514-015-9455-y
- https://people.csail.mit.edu/jrb/papers/slat-ipsn-2006.pdf （正式IPSN2006版，不是MIT2005技术报告）
- https://ojs.aaai.org/index.php/AAAI/article/view/16758

## 二、为什么不是在旧输出后面再加平滑器

本项目直接求解

\[
\min_{X,E}\ \|x_s\|^2_{P_0^{-1}}
+\sum_{j=s}^{k-1}\|x_{j+1}-F x_j\|^2_{Q^{-1}}
+\sum_i\|E_i-m_i\|^2_{C_i^{-1}}
+\sum_{j=s}^{k}\sum_i\|z_{ij}-h(x_j,E_{ij})\|^2_{R_r^{-1}}.
\]

X是窗口内UAV的[p,v,a]；E_i是同一窗口浮标i的位置/姿态误差序列；导航高斯因子(m_i,C_i)由原始伪距/IMU/地磁形成，并包含窗口时间相关性。测量残差同时约束双方，轨迹过程项直接参与估计，不是结果端曲线美化。

窗口长度11个雷达历元=2s；始端位置/速度/加速度给很弱的通用先验，分别标准差10000m、100m/s、10m/s²。优化固定3次Gauss–Newton，随后在最终点重建一次法方程求Laplace协方差。固定参数，无真值标签、无真值选参、无剔除帧。

这是**有限历史MAP**：窗口外旧雷达信息被有意舍弃；上轮解只作数值初值，不作为含有旧测量的先验再加入。没有宣称执行全历史精确边缘化。导航历史则通过一次压缩因子保留。这样牺牲一部分旧目标信息，换取清晰的计数边界和独立节点Schur结构。

## 三、WNOJ公式从何而来

\[
\dot p=v,\quad\dot v=a,\quad\dot a=w,\qquad
E[w(t)w(t')^T]=qI\delta(t-t').
\]

积分得到a(t+h)=a(t)+∫w；v(t+h)=v(t)+ha(t)+∫(h−τ)w(τ)dτ；p(t+h)=p(t)+hv(t)+h²a(t)/2+∫(h−τ)²w(τ)/2 dτ。

因此

\[
F(h)=\begin{bmatrix}I&hI&h^2I/2\\0&I&hI\\0&0&I\end{bmatrix},
\quad
Q(h)=q\begin{bmatrix}h^5/20&h^4/8&h^3/6\\h^4/8&h^3/3&h^2/2\\h^3/6&h^2/2&h\end{bmatrix}\otimes I.
\]

例如Q_pp=q∫_0^h τ⁴/4 dτ=q h⁵/20；Q_pv=q∫τ³/2 dτ=q h⁴/8。其他块同理由积分得到。q单位m²/s⁵，本实验q=1，所有轴一致，没有锁死Z或注入静止标签。

物理含义：允许持续非零加速度，但要求其改变有统计依据。先验既不要求匀速，也不要求加速度永不改变。q为固定探索性参数，没有独立留出集标定，不能声称普适最优。

## 四、导航时间相关性如何保留

原始导航仍用传感器驱动ESKF产生一个线性化参考，但不将每帧后验当成独立观测。保存每历元更新前预测P^-、区间惯性误差转移F以及滤波P^+。线性高斯条件分布满足

\[
A_j=P_j^+F_j^T(P_{j+1}^-)^{-1},\qquad
D_j=P_j^+-A_jP_{j+1}^-A_j^T.
\]

若更新后的右姿态误差需转换到预测姿态切空间，令d为名义状态更新量，T的姿态块为J_r(d_θ)^{-1}≈I+[d_θ]×/2+[d_θ]×²/12。则固定滤波参考下的后向条件为

\[
e_j|e_{j+1},Z^{nav}_{0:k}\sim N(A_jd_{j+1}+B_je_{j+1},D_j),\quad B_j=A_jT_{j+1}.
\]

从当前e_k均值0、协方差P_k^+向窗口起点递推，可得到平滑均值和全跨时间协方差：C_{j,l}=B_j C_{j+1,l}。取位置/姿态6维边际，构成一个6W维联合导航因子，而不是W个互相独立的6维因子。

所有使用的导航观测时间都≤当前雷达时刻k。虽会修订窗口内历史，但保存的在线输出始终是当前末端，未来数据不会回写过去的已报告输出。

导航因素经过ESKF和右误差一阶近似压缩，因此仍非原始非线性全后验的精确表示。尤其首帧比力初始化、有限阶误差传播和重置会影响校准；实际NEES/覆盖率不能由公式外观代替。

## 五、雷达如何同时修正双方

浮标位姿在候选误差E=[δb,δθ]下为b=b_ref+δb，R=R_ref Exp([δθ]×)；雷达q=R^T(p_u−b)，h=spherical(q)。

令B=∂spherical(q)/∂q，则

\[
H_u=[BR^T,0,0],\qquad
H_E=[-BR^T,\ B[q]_\times J_r(δθ)].
\]

这是同一条原始测量对UAV和浮标误差的两个Jacobian。R_r只包含原始雷达测量噪声，因为导航不确定度已经作为显式变量和联合先验出现；不能再次加J P_nav J^T，否则重复计入。

## 六、Schur求解与相关后验

法方程按UAV序列与每个浮标误差序列分块，导航源在浮标间独立，因而H_NN为节点块对角。消元：

\[
S=H_{UU}-\sum_iH_{Ui}H_{ii}^{-1}H_{iU},\quad
b=g_U-\sum_iH_{Ui}H_{ii}^{-1}g_i,
\]
\[
δU=S^{-1}b,\quad δE_i=H_{ii}^{-1}(g_i-H_{iU}δU).
\]

不显式求全部全局大矩阵逆，采用Cholesky求解。最终P_UU≈S^{-1}，L_i=−H_ii^{-1}H_iU，则

\[
P_{iU}=L_iP_{UU},\quad
P_{ij}=\mathbf1_{i=j}H_{ii}^{-1}+L_iP_{UU}L_j^T.
\]

可见即使原始导航源独立，共同观测目标后浮标误差也会相关。代码保留这一结构，不用独立信息相加冒充联合后验。

速度、IMU偏置和钟差已在导航压缩中被边际化，但当前完整16维导航误差可由高斯条件恢复：R_c=C_{full,pose}C_pose^{-1}，E[e_full|E]=R_c(E−m_pose)，剩余条件协方差D_c=P_full−R_cC_{pose,full}。联合雷达更新会经这些交叉块修正浮标速度和偏置，无需读取它们的真值。

## 七、预测与相干发射

UAV末端[p,v,a]按F(h)预测。浮标有两个受控对照：线性p_b+h v_b；惯性二阶p_b+h v_b+h²a_b/2，其中a_b=R(f_m−b_a)+g，全部来自估计和当前IMU，不输入真实波浪分量。

距离d_i=||p_u−b_i||，一阶距离误差δd_i=u_i^T(δp_u−δb_i)。使用上述共享/条件分解，代码计算完整N×N预测距离协方差，包括目标—浮标和浮标间交叉项、未来UAV jerk噪声。未来姿态/IMU误差和模型失配仍有近似，不宣称严格概率校准。

固定总功率下先比较完全相同的均匀功率补偿，避免功率策略掩盖轨迹改善。附加一个明确区分的期望功率策略：在条件高斯距离误差和预测幅度下，

\[
M_{ij}=\hat a_i\hat a_j\exp\{-\tfrac12 k^2(C_{ii}+C_{jj}-2C_{ij})\},\quad
E[P|Z]\approx P_{tot} w^T\bar M w.
\]

\bar M为时隙内平均，||w||²=1，所以Rayleigh商最大值为λ_max(\bar M)，最优振幅权重是其主特征向量。M非负时可取非负主特征向量；无门控、无删帧、无真值选节点。该策略允许每节点功率连续变化，硬件若要求逐节点固定功率则不可采用。

注意：最大期望接收功率不等于所有浮标都有高相干增益。协方差很大时可能集中功率，故同时报告实际接收功率、相对非相干功率增益G、统一均匀功率理想参考下的η和有效节点数N_eff=1/Σw_i⁴。不能用单节点η_own接近1冒充全网相干成功。

## 八、公平比较与计算时延

基线IMM+CI原样重跑，首种子状态与旧完整数据逐项比较并保存最大差值。新增加GP_independent_nav，它同样使用全网原始雷达和WNOJ先验，但近似每帧导航误差独立，用于区分集中观测/运动先验与时间相关联合估计的贡献。JointGP_linear/accel/power共享完全相同联合轨迹后验，只改变发射预测或功率策略。

新方案是集中联合信息参考，不宣称实现有限轮分布式GP。仍假设125ms通信，未实测原始观测传输负载和拓扑聚合。固定延迟评价采用12ms处理+125ms通信+40ms时隙；另报告measured_cpu，以每帧实际求解时间加导航平均耗时替代12ms下限。后者是当前电脑并行负载下的求解/导航耗时诊断，不是硬件最坏时间保证，未计完整预测/通信编解码开销。时隙可能跨越200ms帧周期，不删这些样本，也不声称已经满足实时调度。

全时段和同一20种子、2海况、5轨迹、2载频、5方法、2种延迟评价，总32000行记录。新增位置误差变化率RMS(||(e_p,k−e_p,k−1)/dt||)、速度误差变化率和位置/速度一致性缺陷；它们对真实误差而非仅估计曲线粗糙度计分，避免将真实机动误判成抖动。边界差分使用同一阶段相邻点，是导数定义，不删位置样本。

## 正式参考文献

1. Tang, T. Y., Yoon, D. J., & Barfoot, T. D. (2019). A White-Noise-on-Jerk Motion Prior for Continuous-Time Trajectory Estimation on SE(3). IEEE Robotics and Automation Letters, 4(2), 594–601. https://doi.org/10.1109/LRA.2019.2891492
2. Anderson, S., Barfoot, T. D., Tong, C. H., & Särkkä, S. (2015). Batch nonlinear continuous-time trajectory estimation as exactly sparse Gaussian process regression. Autonomous Robots, 39, 221–238. https://doi.org/10.1007/s10514-015-9455-y
3. Taylor, C., Rahimi, A., Bachrach, J., Shrobe, H., & Grue, A. (2006). Simultaneous localization, calibration, and tracking in an ad hoc sensor network. IPSN, 27–33. https://doi.org/10.1145/1127777.1127785
4. Huai, J., Lin, Y., Zhuang, Y., & Shi, M. (2021). Consistent Right-Invariant Fixed-Lag Smoother with Application to Visual Inertial SLAM. AAAI, 35(7), 6084–6092. https://doi.org/10.1609/aaai.v35i7.16758
5. Ochiai, H., Mitran, P., Poor, H. V., & Tarokh, V. (2005). Collaborative beamforming for distributed wireless ad hoc sensor networks. IEEE Transactions on Signal Processing, 53(11), 4110–4124. https://doi.org/10.1109/TSP.2005.857028
