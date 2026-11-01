# UAV多运动模型：IMM理论与DPC接入设计

版本：2026-10-30。状态：理论设计，尚未实现或验证新增算法；现有CV基线和20种子数据保持不变。

## Material Passport
- 阶段：方法设计与实验规划。
- 问题：单一CV滤波在悬停时产生速度噪声，在持续转弯时产生预测滞后，影响实际发射时刻的路径差与相干效率。
- 输入：当前gm_dpc_unified代码与已保存实验、下列公开研究。
- 输出：可实现的IMM模型、概率更新、发射预测与DPC接口、配对验证协议。
- 证据边界：文献支持IMM与CT跟踪方法，不证明本设计在当前海上协同发射场景一定优于CV。

## 1. 理论依据与设计归属

[1] Blom, H. A. P., & Bar-Shalom, Y. (1988). The interacting multiple model algorithm for systems with Markovian switching coefficients. IEEE Transactions on Automatic Control, 33(8), 780–783. https://doi.org/10.1109/9.1299
原始IMM理论出处；本次IEEE正文访问受到验证页面限制，不能声称已逐式核对该原文。

[2] Yuan, X., Lian, F., & Han, C. (2014). Models and Algorithms for Tracking Target with Coordinated Turn Motion. Mathematical Problems in Engineering, 2014, 649276. https://doi.org/10.1155/2014/649276
已读取公开正文。支持把未知转弯率扩充为状态、使用非线性滤波，并讨论位置观测下的CT跟踪。这里使用的是该类模型，不是复现文中的EM算法。

[3] Huang, H. A. J., Bar-Shalom, Y., Yang, R., & Ng, G. W. (2017). Tracking a Maneuvering Target Using Two Heterogeneous Passive Sensors on a Single Stationary Platform with IMM Estimation. Journal of Advances in Information Fusion, 12(1).
公开全文：https://isif.org/files/isif/2024-01/JAIF_Vol12_1_complete_170828_7.pdf
已读取正文。该研究提供多模型、观测似然与机动跟踪的实际应用依据；其被动传感器和延迟模型不同于当前雷达，不移植其性能数字。

[4] Lee, Y.-L. (2020). Using a New Circular Prediction Algorithm to Design an IMM Filter for Low Update Rate Radar System. Sensors, 20(18), 5035. https://doi.org/10.3390/s20185035
公开全文：https://pmc.ncbi.nlm.nih.gov/articles/PMC7570880/
已读取公开正文。说明线性和圆周预测的组合具有研究基础；不直接采用文中的两点几何圆周预测器。

下文的“近悬停分支、三模型组合、预测后DPC及相干指标验收”是针对本项目的工程设计，不冒称上述任一论文已经证明了这套组合。

## 2. 推荐方案及适用范围

采用本地IMM（近悬停S、近匀速CV、未知转弯率CT）→各分支预测→模型概率加权→有限轮DPC→几何相位命令。

它保留基于位置观测估计速度的路线，不要求多普勒或UAV反馈。模型概率仅来自观测历史，不读取真值轨迹类型、转弯率或切换时刻。

CT描述水平近匀速转弯，垂向采用CV；它不能覆盖所有三维机动。持续切向加减速应另以CA或Singer分支评估。先用三模型隔离悬停和转弯两个已发现问题，再依据加减速测试决定是否扩展，避免一次叠加过多不可辨识模型。

## 3. 模型与共同状态

各分支统一七维状态，避免混合三维、六维和七维协方差时遗漏状态映射：

x=[pX,pY,pZ,vX,vY,vZ,omega]^T，z=H x+e，H=[I3,0(3×4)]。

观测仍用当前雷达r/az/el经导航变换得到的位置及协方差R。坐标转换的一阶近似、导航误差假设与基线一致，避免混入其他改动。

### S：近悬停模型（本项目选用的软约束）

不写死世界坐标、不制造“观测速度=0”的无噪声伪测量。采用阻尼速度：

dp=v dt；dv=−beta v dt+sqrt(qS)dW。

一步条件均值：v'=rho v，p'=p+(1−rho)v/beta，rho=exp(−beta dt)。

该分支把近零速度作为软先验，允许小幅漂移，不把它误称为精确静止模型。beta和qS需要校准；过强阻尼会压制真实低速运动。

每轴连续白噪声强度qS对应精确离散协方差：
Qvv=qS(1−rho²)/(2 beta)
Qpv=qS(1−rho)²/(2 beta²)
Qpp=qS/beta²*[dt−2(1−rho)/beta+(1−rho²)/(2 beta)]。

小beta*dt时采用稳定展开或Van Loan积分，避免浮点消减。omega作为辅助状态采用趋零OU过程，参数与物理速度过程分开。

### CV：近匀速模型

p'=p+dt v，v'=v。p/v过程协方差沿用当前白加速度模型：
Qpv_block=qCV*[[dt³/3 I,dt²/2 I],[dt²/2 I,dt I]]。

omega为趋零辅助过程，不让CV预测调用转弯运动。初始速度不读取真值；沿用零均值与非零初始方差。

### CT：未知转弯率模型

令a=omega*dt，s=sin(a)，c=cos(a)：
pX'=pX+(s/omega)vX−((1−c)/omega)vY
pY'=pY+((1−c)/omega)vX+(s/omega)vY
vX'=c vX−s vY
vY'=s vX+c vY
pZ'=pZ+dt vZ；vZ'=vZ；omega'=omega。

omega由多帧位置观测间接估计，不使用仿真的0.075rad/s作为已知输入。omega→0时使用sinc等稳定极限及一致Jacobian，不能只把状态转移切成CV却遗漏对omega的导数。

首版推荐CT用EKF、S/CV用线性KF；这样相对现有KF变动较小。CT的Q按局部线性化连续运动模型离散，包含角加速度驱动；不能把六维CV的Q不加说明地复制为完整CT过程模型。以后可单独比较UKF，不预设UKF必然更优。

近悬停时omega不可充分观测，必须保留不确定度；不把此时的转弯率估计当可靠物理结论。

## 4. 每个节点的IMM递推

以下省略节点下标。r=3，pi_ab=P(m_k=b|m_(k−1)=a)，每行和为1。

1. 模型先验：c_b=sum_a pi_ab mu_a。
2. 交互混合权重：alpha_(a|b)=pi_ab mu_a/c_b。
3. 混合状态：x0_b=sum_a alpha_(a|b) x_a。
4. 混合协方差：P0_b=sum_a alpha_(a|b)[P_a+(x_a−x0_b)(x_a−x0_b)^T]。
5. 各模型以自己的f、Q独立预测，用同一份当前位置观测及R更新一次。
6. 模型创新：nu_b=z−H xminus_b，S_b=H Pminus_b H^T+R。
7. 对数似然：ell_b=−0.5[nu_b^T S_b^(−1)nu_b+logdet(S_b)+3log(2pi)]。
8. 后验概率：mu_b=softmax(log(c_b)+ell_b)。
9. 输出状态和协方差仍需包含模型间均值散布项，不是仅对各P加权。

似然必须包含logdet(S)，不能只用残差大小，也不能把后验残差拿来替代预测创新。用Cholesky/线性求解和log-sum-exp保证数值稳定，不显式求逆。

模型转移率根据独立校准数据的停留时间设定，不把测试轨迹切换表编码进矩阵。可用连续时间生成矩阵G，pi(dt)=expm(G*dt)，G_ab≥0(a≠b)，G_aa=−sum_(b≠a)G_ab。全部模型保持可达，避免概率被永久锁死。参数范围与初始模型概率在测试前固定。

## 5. 必须同步修改发射预测接口

当前发射命令使用p+tau*v；只替换local_kf而保留这一步，仍会把CT输出沿切线外推。

对时隙采样点s，实际预测跨度h=processing+L*per_round+s。每个节点从自己的IMM后验分支状态出发预测到同一个t_k+h，再按分支概率形成未来位置均值。

短时基线近似：p_hat_i(h)=sum_b mu_ib projection_p[f_b,h(x_ib)]。这里冻结预测窗内的模式，必须显式记录此假设。非线性分支的协方差用Jacobian传播，并计入分支均值散布。更严格版本在预测时传播马尔可夫模式概率和混合分布，作为独立对照。

先做模型内预测，再做模型间加权。不能先把CT/CV/S状态和omega简单平均，再套一个CT或CV转移：非线性转移一般不与加权交换。

节点得到共同时间的预测后，再对这些未来位置和速度做有限轮DPC。可逐时隙复用现有(T,N,6)的dpc接口。每帧初始化于独立本地结果，融合输出仍不回灌本地IMM。

实现上可用同一轮邻居消息携带各模型后验与概率，在本地线性传播邻居贡献；或打包全部时隙预测后并行共识。两者都须核算消息字节数。不能把9个时隙误算成只传一份状态的零成本增加；新增处理时间也必须进入tau。

浮标仍用当前GM条件均值预测到对应时刻，阵元位置由导航和杆臂计算。几何命令直接用预测目标与预测阵元的距离：

d_hat_im(h)=||p_DPC_i(h)−p_hat_elem_im(h)||，theta_im(h)=2pi*d_hat_im(h)/lambda。

这一接口不再额外施加一次tau*v，避免双重预测；每个时隙直接算距离也避免再以一次距离变化率线性外推替代CT曲线。

每个节点继续使用自身有限轮结果，不能偷偷换成集中式全网平均。不将相关的DPC后验假装成独立信息融合协方差。

## 6. 目标指标与可解释性

位置RMSE只是中间指标。主指标保持现有接收复电场功率与eta=P/Pideal。记录每个时隙的命令路径误差去除全网共同分量后的RMS，因为共同相位不影响功率。

由相位定义直接得到成对关系：delta_phi_ij=(2pi/lambda)*(delta_d_i−delta_d_j)。例如将成对差分相位RMS控制在30度，对应路径差RMS约lambda/12：300MHz约8.3cm，1GHz约2.5cm。这是工程参考门槛，不是保证eta达到某值的严格定理。

同时记录模型概率、混合熵、预测位置RMSE、速度RMSE、切换响应延迟、eta时间曲线和低分位数。模型权重不是“真实运动类型的确定标签”；低速CV和近悬停可能不可区分，只要预测合理即可。

IMM解决运动模型失配，不自动解决DPC未收敛、共模导航偏差、载波同步或多径。因此不承诺只换IMM即可解决1GHz相干损失。

## 7. 公平验证方案

固定GM浮标、导航误差、雷达噪声和种子，同一条真值供各算法配对比较。保留已有三条轨迹用于回归，新增连续的悬停→平滑加速→直线→转弯→平滑减速→悬停序列及S弯。位置与速度切换保持物理连续，滤波器不接收切换标记。

校准集与测试集分离：先在独立种子及轨迹参数上确定beta、各q、转移率和初始概率，再固定参数测试。测试包含未见过的速度、转弯率、加减速幅度与导航误差；20次仅是与旧结果对齐的起点，不代表统计精度自动足够。

对照组：
A. 原始CV+DPC，原命令方式；
B. 原始CV+DPC，新时隙距离计算方式，用来分离接口变化收益；
C. 单一CT；
D. S/CV双模型IMM；
E. S/CV/CT三模型IMM完整方案。

所有方法比较同一发射时刻；报告额外计算及消息开销。轮数扫描继续受帧周期约束。模型标签oracle仅可作为明确标注的离线上界，绝不混入主算法。

主报告对每个种子计算eta，再统计配对差值和区间；同时保留功率累加后的eta以兼容现有结果。悬停、直线、转弯、切换后短窗分别报告，不能只靠全程平均掩盖起动或刹停失败。

验收：悬停速度噪声与相干损失改善；转弯的未来位置与差分路径误差下降；直线不出现明显退化；混合轨迹无需人工重置便可恢复跟踪。改善幅度和不退化容限在校准后、测试前记录，失败时如实报告。

## 8. 实施边界

另建新实验目录实现imm_filter与预测接口，导入已有GM仿真器，不覆盖gm_dpc_unified的基线数据。先验证模式概率归一、协方差半正定、CT零转弯率极限、无噪声可恢复、前缀因果性及无真值接口，再跑配对性能实验。

先完成S/CV/CT并验证，若加减速段失败，再把CA或Singer作为单独可归因的扩展。不能把IMM本身当论文创新；可研究的问题是多模型延迟预测、有限轮共识和路径差误差之间的关系及其相干功率影响。
