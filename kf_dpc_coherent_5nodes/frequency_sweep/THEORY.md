# 300 MHz以下通信载频扫描：连续相干公式与适用边界

日期：2026-11-03。这里只研究通信载频，FMCW感知载频固定77 GHz。只改通信载频不会改变原始传感器、UAV/浮标联合状态估计、轨迹RMSE、覆盖率或既定处理时延。新结论只能由完整配对频率重评估得出，不以一次低频相干改善替代轨迹可靠性。

## 正式发表的研究依据

1. Ochiai H., Mitran P., Poor H. V., Tarokh V. (2005). Collaborative Beamforming for Distributed Wireless Ad Hoc Sensor Networks. IEEE Transactions on Signal Processing, 53(11), 4110–4124. DOI https://doi.org/10.1109/TSP.2005.857028 。正式机构发表记录：https://collaborate.princeton.edu/en/publications/collaborative-beamforming-for-distributed-wireless-ad-hoc-sensor-/ 。其研究位置/相位误差与随机阵列波束；本文用的是当前逐节点实际传播距离，不把其远场随机圆盘阵列定理直接套作本实验性能保证。作者网站公开beamform.pdf标DRAFT，不作为正式全文证据。
2. Friis H. T. (1946). A Note on a Simple Transmission Formula. Proceedings of the IRE, 34(5), 254–256. DOI https://doi.org/10.1109/JRPROC.1946.234568 。1946正式期刊扫描全文：https://capmimo.ece.wisc.edu/capmimo_papers/friis_original_1946.pdf 。支持自由空间、天线有效面积与功率传输关系；其成立条件不能因低频相干率更高而忽略。
3. Harrington R. F. (1960). Effect of Antenna Size on Gain, Bandwidth, and Efficiency. Journal of Research of the National Bureau of Standards—D. Radio Propagation, 64D(1), 1–12. DOI https://doi.org/10.6028/jres.064d.003 。NIST正式全文：https://nvlpubs.nist.gov/nistpubs/jres/64D/jresv64Dn1p1_A1b.pdf 。支持固定尺寸相对波长变化时，增益、带宽与损耗效率不能任意同时保持的限制。
4. Chu L. J. (1948). Physical Limitations of Omni-Directional Antennas. Journal of Applied Physics, 19(12), 1163–1175. DOI https://doi.org/10.1063/1.1715038 。引用正式期刊发表记录；MIT公开RLE Technical Report No.64不是同一份期刊排版全文，本文不以该报告替代正式论文。这里只用其电尺寸与带宽约束的定性机制，不借用未经本天线验证的数值界。

没有采用预印本作为研究发表依据。上述来源提供物理机制，不证明当前海面动态浮标系统在任何具体频段实际可部署。

## 1. 从预测距离得到连续频率响应

令k表示观测帧，s表示发射时隙内采样，i表示浮标。真实传播距离为d_ksi，预测命令为dhat_ksi，距离误差e_ksi=dhat_ksi−d_ksi。误差只用于仿真评价，发射命令不能访问真实距离。

单程通信残余相位为

\[
\phi_{ksi}(f)=\frac{2\pi f}{c}e_{ksi}.
\]

在本项目等功率、理想各向同性自由空间参考模型中，令a_ksi=1/d_ksi，则接收复场为

\[
E_{ks}(f)=\sqrt{P_{tot}/N}\frac{c}{4\pi f}
\sum_i a_{ksi}\exp[j2\pi f e_{ksi}/c].
\]

相干效率用同频、同距离、同功率的理想全相干场归一化，公共因子消去：

\[
\boxed{\eta(f)=\frac{\sum_{k,s}|\sum_i a_{ksi}e^{j2\pi f e_{ksi}/c}|^2}
{\sum_{k,s}(\sum_i a_{ksi})^2}}.
\]

相对非相干功率参考的增益为

\[
G(f)=\frac{\sum_{k,s}|\sum_i a_{ksi}e^{j2\pi f e_{ksi}/c}|^2}
{\sum_{k,s,i}a_{ksi}^2},\qquad G_{dB}=10\log_{10}G.
\]

固定N=20与等节点功率下，由Cauchy–Schwarz不等式，G≤N，理想上界约13.01 dB；不同距离使可达的等功率理想G可能略小于20。高η不是通过减少节点或集中功率得到。

展开分子可见真正影响相干的是差分距离误差：

\[
|\sum_i a_i e^{jke_i}|^2=\sum_i a_i^2+2\sum_{i<j}a_ia_j\cos[k(e_i-e_j)],\quad k=2\pi f/c.
\]

因此全体节点相同的距离误差不降低本相位模型中的相干。误差不变时降低f，会降低相位敏感性，但有限误差样本形成的余弦和不保证在全频段严格单调，不能用两个已有频点线性插值代替扫描。

## 2. 为什么低频趋近完美，不等于找到真实最优频段

有限距离误差下f→0+时每个相位因子趋于1，故在固定a_i假设内η(f)→1。用cos展开：

\[
\eta(f)=1-\left(\frac{2\pi f}{c}\right)^2
\frac{\sum_{k,s}\sum_{i<j}a_i a_j(e_i-e_j)^2}
{\sum_{k,s}(\sum_i a_i)^2}+O(f^4).
\]

这给出低频损失的二次敏感性，是当前距离误差的条件结论，不是低频硬件效率结论。它也说明降低通信载频不修复原轨迹异常：只减小那些异常投影到传播相位后的代价。

真实频率相关链路应至少写成

\[
P_r(f)=P_tG_t(f)G_r(f)\left(\frac{c}{4\pi f d}\right)^2
L^{-1}(f),
\]

其中Gt、Gr应含实际辐射效率、失配和方向图等约定，L含传播/其他损耗。固定天线尺寸时ka=2πaf/c随频率降低，不能默认增益、辐射效率与有效带宽保持不变；海面传播、多径及噪声也未由本次相位扫描补齐。因此不能把理想Pr∝f^-2外推成“越低频越适合实际系统”，更不能把f→0当作可用通信方案。

当波长与浮标间距相当或更大时，独立发射单元、忽略互耦及ΣPi=1W的功率约定也需重新核对。当前极低频曲线只是固定幅度叠加模型的数学极限，不证明同一组真实天线在固定总辐射功率下仍具有相同结果。

## 3. 频率范围如何解释

1 kHz–300 MHz数值网格是跨频率的数学敏感性分析；公式本身定义了全部f>0的连续曲线，有限网格不是枚举了每个实数频点。300 MHz保留作与原实验的边界对照。

30–300 MHz作为当前传播参考模型的主要展示区，低于30 MHz的LF/MF/HF结果明确标注为相位模型外推。30 MHz只是展示与解释的分区，不是已经证明天线/传播有效的硬分界；30–300 MHz本身也仍需实际天线、传播和链路预算验证。不给出许可频段、合法发射名单或硬件可用性保证。

当前是载频处的窄带场叠加，没有为各频率重建通信带宽、调制、编码、接收噪声或数据率。保留125ms通信开销仅用于受控比较相位效应，不证明低频链路能承载同样的数据量或相同实时协议。

## 4. 必须同时报告哪些指标

- 全程能量比η，以及逐帧η的均值与低分位。全程η按理想接收功率加权，近距离帧权重更高，不能掩盖远段异常。
- 本次频率扫描单列75–90秒与完整0–90秒；加速/刹停的轨迹精度和覆盖率沿用上一轮完整实验，通信频率不改变它们。几何缓存可用于进一步计算这些时段的相干曲线。
- 两种海况、五轨迹、全部20种子按配对统计，不把频率网格、时隙采样或节点当作新的独立重复。
- 轨迹误差、覆盖率与时延仍保留原判据。通信频率扫描不会使原来不可信的远段状态后验自动可信，也不会把超过200 ms的调度自动变为满足约束。
- 绝对接收功率只在明确的理想天线参考假设下解释；低频外推不用于宣称真实最高接收功率或最优可部署载频。

完整缓存应保留每帧时隙的预测距离误差及真实距离权重，从同一缓存评价全部f，避免为每个通信频率重新改变传感器或优化器参数。
