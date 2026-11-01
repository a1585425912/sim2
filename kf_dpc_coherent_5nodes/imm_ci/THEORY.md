# CI融合与相干合成：正式文献及数学推导

日期：2026-10-30。本文将已发表研究、对本项目的推导、仿真实测三者分开。只列正式期刊与会议论文，不使用预印本作为依据。数学分析不是额外数值测试，也不设置筛选门槛。

## 1. 文献与证据边界

1. S. J. Julier and J. K. Uhlmann, “A non-divergent estimation algorithm in the presence of unknown correlations,” Proceedings of the American Control Conference (ACC), vol. 4, pp. 2369–2373, 1997. [正式DOI](https://doi.org/10.1109/ACC.1997.609105)。CI原始会议论文。本次未能取得出版社全文，不编造其定理号或页内公式编号。
2. G. Battistelli and L. Chisci, “Kullback–Leibler average, consensus on probability densities, and distributed state estimation with guaranteed stability,” Automatica, vol. 50, no. 3, pp. 707–718, 2014. [正式DOI](https://doi.org/10.1016/j.automatica.2013.11.042)。信息共识的期刊理论出处；出版社全文访问受限。
3. G. Battistelli, L. Chisci, G. Mugnai, A. Farina, and A. Graziano, “Consensus-Based Linear and Nonlinear Filtering,” IEEE Transactions on Automatic Control, vol. 60, no. 5, pp. 1410–1415, 2015. [正式DOI](https://doi.org/10.1109/TAC.2014.2357135)。用于区分信息共识、测量共识及混合方案；本实验不能继承论文在特定假设下的全部稳定性结论。
4. H. Ochiai, P. Mitran, H. V. Poor, and V. Tarokh, “Collaborative beamforming for distributed wireless ad hoc sensor networks,” IEEE Transactions on Signal Processing, vol. 53, no. 11, pp. 4110–4124, 2005. [正式DOI](https://doi.org/10.1109/TSP.2005.857028)。分布式相干发射及相位/位置误差的正式期刊依据；本次出版信息及研究范围由作者所在高校的出版记录核对，没有把作者稿当作正式全文。
5. G. Battistelli, L. Chisci, C. Fantacci, A. Farina, and A. Graziano, “Consensus-based multiple-model Bayesian filtering for distributed tracking,” IET Radar, Sonar & Navigation, vol. 9, no. 4, pp. 401–410, 2015. [正式出版全文](https://ietresearch.onlinelibrary.wiley.com/doi/10.1049/iet-rsn.2014.0071)。第2.4节给出KLA与归一化加权几何平均的关系，本文用作可核对的补充依据。其DIMM设计并不等于本项目“先合并本地IMM、再做输出CI”的实现。

下列公式按以上理论路线从定义展开，不声称逐字重现无法取得全文的论文，也不把本仿真的数值提升归因于论文已经证明。

## 2. 从KL目标推导高斯CI

设同一时刻、同一坐标系下，第j个节点的局部状态近似分布为

\[
p_j(x)=\mathcal N(x;\hat x_j,P_j),\quad P_j\succ0.
\]

权重非负且归一化：\(\omega_j\ge0,\ \sum_j\omega_j=1\)。考虑KLA的变分目标

\[
q^*=\arg\min_{q:\int q=1}\sum_j\omega_jD_{KL}(q\Vert p_j).
\]

展开并加入归一化乘子：

\[
\mathcal L(q)=\int q\log q\,dx-\int q\sum_j\omega_j\log p_j\,dx
+\lambda(\int q\,dx-1).
\]

一阶变分条件为

\[
\log q(x)+1-\sum_j\omega_j\log p_j(x)+\lambda=0,
\]

从而

\[
q^*(x)=\frac{\prod_jp_j(x)^{\omega_j}}{\int\prod_jp_j(\xi)^{\omega_j}\,d\xi}.
\]

令\(Y_j=P_j^{-1},\ y_j=Y_j\hat x_j\)。代入高斯对数密度，保留含x的项：

\[
\log q^*(x)=-\tfrac12x^T(\sum_j\omega_jY_j)x+x^T\sum_j\omega_jy_j+C.
\]

通过配方得到

\[
\boxed{P_{CI}^{-1}=\sum_j\omega_jP_j^{-1}},\qquad
\boxed{\hat x_{CI}=P_{CI}\sum_j\omega_jP_j^{-1}\hat x_j}.
\]

这解释了信息共识（Consensus on Information）与高斯协方差交集（Covariance Intersection）的联系。两者常使用同一个CI缩写，但不应据此把所有分布式滤波实现视作同一算法。

## 3. 为什么不乘N；何时与DPC完全相同

若输入完全重复：\(P_j=P,\hat x_j=\hat x\)，则

\[
P_{CI}^{-1}=(\sum_j\omega_j)P^{-1}=P^{-1},\quad \hat x_{CI}=\hat x.
\]

即重复信息不会使协方差凭空缩小。把凸平均替换成不加权信息求和会给出P/N，只有充分处理独立性和公共先验后才可能成立；不能作为未知相关轨迹的CI公式。

若只有协方差相同、均值不同，则

\[
\hat x_{CI}=P\sum_j\omega_jP^{-1}\hat x_j=\sum_j\omega_j\hat x_j.
\]

因此CI不天然优于DPC；输入不确定度一样时二者均值相同。本次不同结果来自局部后验的各向异性、不确定度差异及有限轮信息传播。

## 4. 未知交叉相关时的保守性：条件与证明

假定各节点误差\(e_j=\hat x_j-x\)满足二阶矩界

\[
E[e_je_j^T]\preceq P_j.
\]

这里是误差二阶矩，不仅是去掉偏置后的协方差。设\(P=(\sum_j\omega_jY_j)^{-1}\)，则

\[
e_{CI}=P\sum_j\omega_jY_je_j.
\]

对任意确定向量v，定义随机标量\(a_j=v^TPY_je_j\)。凸平方不等式逐样本成立：

\[
(\sum_j\omega_ja_j)^2\le\sum_j\omega_ja_j^2.
\]

因此，不需要知道任何跨节点交叉协方差，就有

\[
\begin{aligned}
E[(v^Te_{CI})^2]
&\le\sum_j\omega_jv^TPY_jE[e_je_j^T]Y_jPv\\
&\le v^TP(\sum_j\omega_jY_j)Pv=v^TPv.
\end{aligned}
\]

于是\(E[e_{CI}e_{CI}^T]\preceq P\)。这证明的是在输入二阶矩界有效时，输出仍保守；不是“CI总比算术平均更准确”。当前IMM的协方差来自近似模型与线性化，尤其机动时未必覆盖真实误差，不能无条件套用此结论。

## 5. 与本实验有限轮分布式实现对应

第l轮只交换邻居信息：

\[
Y_i^{(l+1)}=\sum_jw_{ij}Y_j^{(l)},\quad
y_i^{(l+1)}=\sum_jw_{ij}y_j^{(l)}.
\]

展开L轮后

\[
Y_i^{(L)}=\sum_j[W^L]_{ij}Y_j^{(0)},\quad
y_i^{(L)}=\sum_j[W^L]_{ij}y_j^{(0)}.
\]

本实验取20节点环网、Metropolis权重、L=25，\(W^L\)每行仍是非负归一化权重。它不是预设每节点拿到全网均值。代码预计算有限轮矩阵只是代数等价加速。

各本地IMM先在预测时刻h做分支传播，再做混合矩匹配：

\[
\hat x(h)=\sum_m\mu_m\hat x_m(h),
\]

\[
P(h)=\sum_m\mu_m\{P_m(h)+[\hat x_m(h)-\hat x(h)][\hat x_m(h)-\hat x(h)]^T\}.
\]

随后取位置/速度的六维边际子块，再求逆形成信息。不能先求七维信息矩阵再直接截取其六维块，因为后者通常对应对转弯率条件化。

本项目的近似边界：先把混合高斯矩匹配成单个高斯，再做CI；不同于直接对混合分布做几何融合，也不同于文献中保留模式概率、模式条件后验并参与递推的完整DIMM。共识输出不回灌本地IMM，因此论文闭环滤波稳定性定理不能直接移植。

## 6. 从状态误差推导到相干功率

以下将每个发射阵元统一记作i。真实距离\(d_i=\|p_u-b_i\|\)，算法命令距离\(\hat d_i\)，距离残差\(r_i=\hat d_i-d_i\)。实际复场及功率为

\[
E=\sum_i a_i e^{\mathrm jkr_i},\quad P=|E|^2,\quad k=2\pi/\lambda.
\]

直接展开平方得

\[
\boxed{P=\sum_i a_i^2+2\sum_{i<j}a_ia_j\cos[k(r_i-r_j)]}.
\]

因此功率依赖成对差分距离残差，不是三维坐标RMSE本身。对LOS单位向量\(u_i=(p_u-b_i)/d_i\)，一阶展开为

\[
r_i\approx u_i^T(\delta p_i-\delta b_i),\quad
\delta p_i(h)\approx\delta p_i(0)+h\delta v_i(0)+\text{运动模型失配项}.
\]

即使每个节点的位置误差范数变小，也不保证\(r_i-r_j\)变小，因为其方向、节点相关性和速度预测误差都会影响结果。共同位置误差并不必然抵消——各节点LOS不同；完全共同的距离残差才只对应全网公共相位而不改变P。

在阵元幅度固定、差分相位\(\Delta\phi_{ij}=k(r_i-r_j)\)近似高斯时，若其均值为\(m_{ij}\)、方差为\(s_{ij}^2\)，由高斯特征函数得到

\[
E[\cos\Delta\phi_{ij}]=e^{-s_{ij}^2/2}\cos m_{ij}.
\]

从而

\[
E[\eta]=\frac{\sum_i a_i^2+2\sum_{i<j}a_ia_j e^{-s_{ij}^2/2}\cos m_{ij}}{(\sum_i a_i)^2}.
\]

这只是额外高斯/固定幅度假设下的解释公式；实际实验不使用该近似计算功率，而是逐阵元叠加复场。它清楚说明：只压低绝对位置RMSE，并不能保证提升相干效率。

## 7. 与完整20种子结果对应

统计全部0–90s帧，没有门限或样本筛选：

- 混合运动三维位置RMSE：DPC 0.7977m → CI 0.6859m；300MHz相干效率0.7149 → 0.7322。
- 但同一混合运动起飞加速段：位置RMSE 2.484m → 2.115m，η却0.1498 → 0.1338。
- 刹停段：位置RMSE 0.876m → 1.064m，η从0.4947降到约0.342。
- 爬升高度RMSE约0.5558m → 0.5564m，未改善。

这些是本实验测量，不是文献结论。第6节给出了“位置更准却未必更相干”的数学原因；是否由哪种具体偏置或概率误判主导，还需额外针对性分析，不能仅凭现象直接认定根因。

综合判断：当前固定权重CI有全程平均融合收益，但不是运动模型失配、垂直模式混用或加减速滞后的完整修复，也不能宣称是完整分布式IMM的复现。
