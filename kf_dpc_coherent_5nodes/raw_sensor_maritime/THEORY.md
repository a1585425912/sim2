# 零杆臂动态浮标—UAV仿真：数学模型与实现边界

## Material Passport
- 阶段：本目录实现的数学说明。数字结果以完整实验RESULTS.md为准。
- 文献仅引用正式出版的期刊/会议；本文给出的具体离散参数、工程接口和实验组合属于本项目设计，不冒称文献已经验证海上相干任务。
- 世界坐标为局部右手XYZ，Z向上；重力g=[0,0,-9.81]m/s²；R为机体系到世界系；右姿态误差R_true=R_hat Exp([δθ]×)。

## 1. 双时间尺度与共享海况

浮标参考点b_i(t)=c_i(t)+s_i(t)。慢速度选取有限谱随机过程：

v_c,i(t)=vbar+Σ_l[a_il cos(ν_l t)+b_il sin(ν_l t)]，

c_i(t)=c_i(0)+vbar t+Σ_l[a_il sin(ν_l t)/ν_l+b_il(1−cos(ν_l t))/ν_l]。

a,b各分量独立N(0,σ_v²/L)，所以每个水平速度分量Var(v_c−vbar)=σ_v²。ν_l=2π/[180,120,80]s；整个实验一次抽取系数，任意时间查询都沿同一条连续轨道。它有有限谱的周期性，不是OU过程，也不覆盖任意海流。

共享波场η(r,t)=Σ_m A_m cos(ω_m t−k_mᵀr+φ_m)。频率幅值来自JONSWAP形状并归一化Σ A_m²/2=Hs²/16。该式为波面谱方差，不保证带RAO的浮标升沉标准差也是Hs/4。

有限水深色散ω²=g k tanh(kh)，由Newton迭代求k。各浮标共享频率、波向、相位，采样位置不同，具有空间相关性。

平移、姿态采用复响应H_p,m、H_θ,m：s_i=Re Σ H_p,m exp(jψ_im)，角度亦然。升沉响应取H(ω)=1/[1−(ω/ω_n)²+2jζω/ω_n]；水平、倾角响应为显式简化参数。系数没有水动力标定，仅为轻量物理基线。

沿慢中心采样，ψ'=ω−kᵀv_c，ψ''=−kᵀa_c，因此

s'=Re Σ jψ'H_p exp(jψ)，

s''=Re Σ[−(ψ')²+jψ'']H_p exp(jψ)。

这保证位置/速度/加速度使用同一解析样本，没有独立位置抖动，也没有网格边界随机跳变。姿态由ZYX欧拉角构造，并将欧拉角导数转为机体系角速度，使Rdot=R[ω_B]×。

低频/波频分解的海洋运动背景见Fossen & Perez [1]。本项目的浮标响应不是文献船型参数的直接移植。

## 2. 原始观测，而非导航真值加噪声

GNSS第j颗卫星伪距：ρ_ij=||s_j−b_i||+c_clock,i+n_ρ,ij。

这里c_clock已乘光速，单位m；八颗合成可见卫星位置已知。接收机钟差真值为连续时间一阶有色过程在GNSS历元的精确离散样本；滤波用随机游走近似。未模拟大气、多径、载波模糊度或电波基带。

IMU测量f_m=Rᵀ(b''−g)+b_a+n_a，ω_m=ω_B+b_g+n_g。

偏置b_{k+1}=b_k+σ_rw sqrt(Δt)ξ；连续白噪声密度σ的区间平均读数标准差为σ/sqrt(Δt)。因此惯性过程协方差随Δt增加，而读数噪声方差随1/Δt增加，二者单位不能混用。以中点真值近似区间平均，包在区间末可用。

磁场m_m=Rᵀm_world+n_m，用已知非竖直地磁场提供航向参考。真实角度不是测量输入。磁干扰未模拟，存在磁干扰时不能借用当前结果。

雷达z=[r,az,el]+n，其中q=Rᵀ(p_u−b)，r=||q||，az=atan2(q_y,q_x)，el=atan2(q_z,sqrt(q_x²+q_y²))。无Doppler观测。

杆臂为零：GNSS/IMU/雷达/发射点共位。但R仍影响IMU、地磁和雷达测量；零杆臂不等于固定姿态。

## 3. 从观测初始化

位置与钟差联合求解min_{p,c} Σ_j(ρ_j−||s_j−p||−c)²/σ_ρ²。

残差对预测的Jacobian H_j=[(p−s_j)ᵀ/||p−s_j||,1]。每次Gauss–Newton步δ=(HᵀH)^−1 Hᵀr。六次迭代从局部坐标原点和零钟差开始；初始位置/钟差协方差近似σ_ρ²(HᵀH)^−1。

首个加速度观测归一化得到机体系向上方向u_B；磁向量去掉在u_B上的投影，得到北向n_B；e_B=u_B×n_B，R_0=[n_B,e_B,u_B]ᵀ。波浪加速度使该初始倾角存在误差，故给10°初始姿态标准差，不宣称首帧精确。初速度和两类偏置为零，配非零先验方差，未读取真值。

初始化已使用的首个GNSS/磁场观测不在首帧再重复更新。

## 4. 原始惯性导航与16维右误差ESKF

名义动力学p'=v，v'=R(f_m−b_a)+g，R'=R[ω_m−b_g]×。惯性运动方程与偏置模型的正式机器人学背景见Forster等[2]；海上MEMS捷联惯性导航与非线性估计的正式实船研究见Rogne等[8]。本项目实现普通误差状态递推，不声称复现这些论文的全部算法或实船性能。

右误差δx=[δp,δv,δθ,δb_a,δb_g,δc_clock]，维数16。将R_true=R_hat(I+[δθ]×)代入动力学并保留一阶，得

δp'=δv，
δv'=−R[f]×δθ−Rδb_a−Rn_a，
δθ'=−[ω]×δθ−δb_g−n_g。

这是代码A的非零块。F≈I+AΔt+A²Δt²/2；使用中点R更新平移。噪声Q保留加速度白噪声导致的pp=σ_a²Δt³/3、pv=σ_a²Δt²/2、vv=σ_a²Δt；姿态噪声σ_g²Δt，偏置与钟差噪声为各自随机游走密度平方乘Δt。该Q忽略部分高阶耦合，不称精确Van Loan全16维离散化。

伪距更新Jacobian仅含位置H_p=(p−s)ᵀ/d和钟差H_c=1。磁观测h=Rᵀm_world，h_true≈h+[h]×δθ，所以H_θ=[h]×。

K=PHᵀ(HPHᵀ+R_meas)^−1，δx=Kr；注入位置、速度、偏置、钟差，姿态R←R Exp([δθ]×)。P用Joseph更新，姿态注入后使用一阶误差重置J_θ≈I−[δθ]×/2。基本高斯线性最小均方更新见Kalman[3]；非线性和重置为本项目一阶近似。

## 5. 估计浮标位姿条件下的原始雷达IMM

目标模型沿用已有S/CV/CT及IMM交互[4]，不向其传入真实航迹标签或切换时刻。S速度回归零；CV含连续白加速度；CT确定性水平转弯加世界坐标离散加速度噪声，后者仍是明确的近似CT过程模型。

每分支直接预测h(x,nav_hat)，对原始球坐标残差更新，两个角度残差wrap到[−π,π)。不是先把整条未来曲线拟合后称为在线估计。

令B=∂spherical(q)/∂q，则H_x=[BRᵀ,0,0]。对导航右误差，δq≈−Rᵀδb+[q]×δθ，所以J_nav的位置和姿态块分别为−BRᵀ和B[q]×。

采用有效测量协方差R_eff=R_radar+J_nav P_nav J_navᵀ，含位置—姿态边际交叉块。它近似边际化当前导航误差，但忽略目标先验与导航误差的跨时间交叉相关。因此这不是完整联合MAP，也不自动保证置信度校准。NEES/覆盖率必须据结果报告。

首帧用测量与导航解算位置，速度先验0±20m/s；首帧雷达不重复更新。后续使用各分支创新似然更新模式概率，混合协方差含模式间均值散布。

## 6. 从估计导航分离慢、快分量

每轴状态ξ=[c,v_c,s,v_s]ᵀ，动力学A的慢块为CV，快块为阻尼振子；观测z=[p_nav,v_nav]ᵀ，H=[[1,0,1,0],[0,1,0,1]]。F=e^{Ah}，Q用Van Loan积分。慢速度弱扩散表示未知慢流，快噪声按名义波况标定。

该系统的慢、快可分性来自不同极点和时间观测，不直接提供真值分解。初始s=v_s=0，c=p_nav，v_c=v_nav；初始慢快误差设负交叉项，使总位置/速度方差与导航输入一致，避免首帧将未知波浪方差重复相加。

未来b_hat(h)=[1,0,1,0]e^{Ah}ξ_hat；线性对照为p_nav+h v_nav。两者都不读取未来运动或真实波浪相位。

重要限制：导航后验跨时间相关，二级滤波没有完整建模其时间相关噪声；其P不当作高可信联合后验使用。本实验用其均值预测，并直接评价真实未来位置误差。波浪预测不是必然优于短时直线外推。

## 7. 有限共识与发射

UAV每个分支预测至h，P_m(h)=F_m P_m F_mᵀ+Q_m；先矩匹配、取六维边际，再融合。短时预测冻结模式，未精确积分未来跳变模型。

DPC：x_i^f=Σ_j C_ij x_j，保守边际P_i^f=Σ_j C_ij P_j。CI：Y_i^f=Σ_j C_ij P_j^−1，y_i^f=Σ_j C_ij P_j^−1 x_j，P_i^f=(Y_i^f)^−1，x_i^f=P_i^f y_i^f[5]。C=W^25为固定环网Metropolis有限轮核，不用全网平均替代，不乘节点数，不优化真值权重。

所有方法在h=12ms+25×5ms+[0,40]ms发射。每浮标一个共位各向同性阵元，d_hat,i(h)=||p_hat_u,i(h)−b_hat_i(h)||。真实通道d_i(h)只供仿真器评价。

给定总功率P_tot，a_i=√(P_tot/N) λ/(4πd_i)，实际复场

E=Σ_i a_i exp{j(2π/λ)[d_hat,i−d_i]}，P_rx=|E|²，η=P_rx/(Σ_i a_i)²。

全时隙η使用各采样点功率和与理想功率和之比，不先对相位做平均。同步残差设为零，是理想载波同步假设，不是由轨迹估计解决了同步。

展开得到P_rx=Σ a_i²+2Σ_{i<j}a_i a_j cos[k(r_i−r_j)]。所以跨节点差分距离误差决定相干损失，位置RMSE不等价于相干目标[6,7]。本实验按点预测做相位补偿，尚未求解一般后验下最大期望相干功率的相位/功率联合优化。

## 8. 完整实验与统计边界

20独立种子×Hs=0.5/1.5m×5轨迹×4方法×2载频。每场景共用原始观测，成对比较。完整0–90s为主结果；启动、加速、转弯、刹停分段只解释现象，不删样本。无冒烟、门控或门限筛选。

RMSE=sqrt(mean ||e||²)先按种子计算，再报告种子均值；三轴指标分别报告。NEES=eᵀP^−1e，名义95%区域覆盖率为落入椭球eᵀP^−1e≤χ²_d(.95)的经验比例，不是发射门控，也不是独立样本假设检验。偏差和模型失配可使P失准，不能仅凭CI形式或小P声称可靠。

种子配对差Δ_s的95%区间为mean(Δ)±t_19(.975) std(Δ)/sqrt(20)，条件是固定轨迹/海况和独立种子，未做多重比较校正；不能用20种子证明普遍最优。时序和节点数据相关，不作为独立显著性样本。

## 9. 正式出版来源

1. Fossen, T. I., & Perez, T. (2009). Kalman Filtering for Positioning and Heading Control of Ships and Offshore Rigs. IEEE Control Systems Magazine, 29(6), 32–46. https://doi.org/10.1109/MCS.2009.934408
2. Forster, C., Carlone, L., Dellaert, F., & Scaramuzza, D. (2017). On-Manifold Preintegration for Real-Time Visual–Inertial Odometry. IEEE Transactions on Robotics, 33(1), 1–21. https://doi.org/10.1109/TRO.2016.2597321
3. Kalman, R. E. (1960). A New Approach to Linear Filtering and Prediction Problems. Journal of Basic Engineering, 82(1), 35–45. https://doi.org/10.1115/1.3662552
4. Blom, H. A. P., & Bar-Shalom, Y. (1988). The interacting multiple model algorithm for systems with Markovian switching coefficients. IEEE Transactions on Automatic Control, 33(8), 780–783. https://doi.org/10.1109/9.1299
5. Julier, S. J., & Uhlmann, J. K. (1997). A non-divergent estimation algorithm in the presence of unknown correlations. Proceedings of the American Control Conference, 2369–2373. https://doi.org/10.1109/ACC.1997.609105
6. Ochiai, H., Mitran, P., Poor, H. V., & Tarokh, V. (2005). Collaborative beamforming for distributed wireless ad hoc sensor networks. IEEE Transactions on Signal Processing, 53(11), 4110–4124. https://doi.org/10.1109/TSP.2005.857028
7. Mudumbai, R., Brown III, D. R., Madhow, U., & Poor, H. V. (2009). Distributed transmit beamforming: challenges and recent progress. IEEE Communications Magazine, 47(2), 102–110. https://doi.org/10.1109/MCOM.2009.4785387
8. Rogne, R. H., Bryne, T. H., Fossen, T. I., & Johansen, T. A. (2021). On the Usage of Low-Cost MEMS Sensors, Strapdown Inertial Navigation, and Nonlinear Estimation Techniques in Dynamic Positioning. IEEE Journal of Oceanic Engineering, 46(1), 24–39. https://doi.org/10.1109/JOE.2020.2967094

以上来源支持各自通用理论，不构成本项目参数已标定、置信度已校准或联合算法已最优的证据。完整实验以实际数据为准。
