# 当前实验的方法来源与文献索引

所有条目为正式发表论文；以下区分文献机制与本项目适配。完整推导见 ../fmcw_doppler/THEORY.md、THEORY.md 和 ../joint_gp/RESEARCH.md。

## michalczyk2022

Michalczyk, Jan and Jung, Roland and Weiss, Stephan. *Tightly-Coupled EKF-Based Radar-Inertial Odometry*. 2022 IEEE/RSJ International Conference on Intelligent Robots and Systems (IROS), 2022. [DOI](https://doi.org/10.1109/IROS47612.2022.9981396)

采用范围：紧耦合FMCW距离/速度与IMU；原文估计雷达自身运动，当前外部移动目标因子为项目推导。

## cohen2023

Cohen, Itzik and Levanon, Nadav. *Bistatic Radar Positioning of a Moving Target, Using Only Range and Range-Rate Measurements*. IEEE Transactions on Aerospace and Electronic Systems, 2023. [DOI](https://doi.org/10.1109/TAES.2022.3202758)

采用范围：距离和距离率联合定位；原文双基地、特定运动参数，不直接证明本项目任意三维可观性。

## siddiq2017

Siddiq, Kashif and Hobden, Mervyn and Watson, Robert and Pennock, Steve and Martins, Steve. *On Phase Measurement in FMCW Radar Systems*. 2017 Sensor Signal Processing for Defence Conference (SSPD), 2017. [DOI](https://doi.org/10.1109/SSPD.2017.8233250)

采用范围：FMCW相位测量背景；非本项目完整测速器精度依据。

## tang2019

Tang, Timothy Y. and Yoon, David J. and Barfoot, Timothy D.. *A White-Noise-On-Jerk Motion Prior for Continuous-Time Trajectory Estimation on SE(3)*. IEEE Robotics and Automation Letters, 2019. [DOI](https://doi.org/10.1109/LRA.2019.2891492)

采用范围：WNOJ先验；当前采用欧氏位置/速度/加速度特例，非完整SE(3)复现。

## anderson2015

Anderson, Sean and Barfoot, Timothy D. and Tong, Chi Hay and Sarkka, Simo. *Batch Nonlinear Continuous-Time Trajectory Estimation as Exactly Sparse Gaussian Process Regression*. Autonomous Robots, 2015. [DOI](https://doi.org/10.1007/s10514-015-9455-y)

采用范围：随机微分方程与GP轨迹先验；当前不是完整批量非线性GP重优化。

## bell1994

Bell, Bradley M.. *The Iterated Kalman Smoother as a Gauss-Newton Method*. SIAM Journal on Optimization, 1994. [DOI](https://doi.org/10.1137/0804035)

采用范围：高斯状态估计与Gauss-Newton联系；不将三次迭代视为全局收敛证明。

## fossen2009

Fossen, Thor I. and Perez, Tristan. *Kalman Filtering for Positioning and Heading Control of Ships and Offshore Rigs*. IEEE Control Systems Magazine, 2009. [DOI](https://doi.org/10.1109/MCS.2009.934408)

采用范围：海洋运动低频/波频分解机制；具体浮标响应参数为仿真假设。

## ochiai2005

Ochiai, Hideki and Mitran, Patrick and Poor, H. Vincent and Tarokh, Vahid. *Collaborative Beamforming for Distributed Wireless Ad Hoc Sensor Networks*. IEEE Transactions on Signal Processing, 2005. [DOI](https://doi.org/10.1109/TSP.2005.857028)

采用范围：分布式相干叠加及相位/位置误差；本项目有限距离动态预测为场景构造。

## aravkin2014

Aravkin, Aleksandr Y. and Burke, James V. and Pillonetto, Gianluigi. *Robust and Trend-Following Student's t Kalman Smoothers*. SIAM Journal on Control and Optimization, 2014. [DOI](https://doi.org/10.1137/130918861)

采用范围：旧模型T-Trend过程重尾机制；当前为因果特化，条件协方差不包含尺度不确定性。

## friis1946

Friis, Harald T.. *A Note on a Simple Transmission Formula*. Proceedings of the IRE, 1946. [DOI](https://doi.org/10.1109/JRPROC.1946.234568)

采用范围：自由空间功率参考；不代表海面实际链路预算。

## harrington1960

Harrington, Roger F.. *Effect of Antenna Size on Gain, Bandwidth, and Efficiency*. Journal of Research of the National Bureau of Standards, Section D, 1960. [DOI](https://doi.org/10.6028/jres.064d.003)

采用范围：低频天线尺寸/带宽/效率限制。

## chu1948

Chu, Lan Jen. *Physical Limitations of Omni-Directional Antennas*. Journal of Applied Physics, 1948. [DOI](https://doi.org/10.1063/1.1715038)

采用范围：电小天线限制；非具体设备的已验证数值界。

## 正式公开全文与核对边界

- Fossen/Perez正式全文：https://www.fossen.biz/publications/2009%20Fossen%20and%20Perez%20IEEE%20CST.pdf
- Aravkin等正式SIAM全文：https://sites.math.washington.edu/~burke/papers/reprints/82-robust-trend-T.pdf
- Friis正式期刊扫描：https://capmimo.ece.wisc.edu/capmimo_papers/friis_original_1946.pdf
- Harrington正式NIST全文：https://nvlpubs.nist.gov/nistpubs/jres/64D/jresv64Dn1p1_A1b.pdf

IROS 2022、TAES 2023的机制核对包含正式摘要/发表记录，不声称本次已读其付费正式全文。作者标注Pre-print的IROS PDF、Ochiai作者网站的DRAFT及Chu技术报告不作为正式期刊全文替代。

77 GHz ADC参数、参考SNR、测速似然尺度、过程噪声q、迭代次数、通信时延均为本项目实验设置，不归因于上述论文。最新多普勒主模型集中保留联合协方差，没有使用旧版本CI。完整系统由上述已有机制加本项目相对运动因子和预测相位链路构成，不宣称复现某一篇完整算法或已经获得全场景最优性证明。
