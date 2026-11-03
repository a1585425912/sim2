# 300MHz以下通信载频完整敏感性实验

FMCW测速载频固定77GHz；只改变相干通信频率。保持完整20种子20261028–47、两海况、五轨迹和0–90s所有帧，比较无多普勒与联合多普勒q=1两个主估计器，其他参数不变。不做冒烟、小测试、测量门控或删帧。

主两方法576000行已完整生成，新增8个旧策略/受控预测对照的20种子重放也已完整结束。公共183.395ms条件合并10策略1440000频率行、16000轨迹分段行；同一估计器不同浮标预测不当作独立运动模型。先看ASSESSMENT.md与MODEL_VALIDITY.md判断模型，再看LEGACY_RESULTS.md及频率曲线；当前没有一版已证明全场景同时满足可靠轨迹和完整任务目标。

360频率点：1kHz–1MHz为61点对数网格，2–300MHz为整数MHz网格，300MHz是原参考端点。连续公式覆盖中间频率，有限网格不等于逐一枚举实数频率。主比较30–300MHz；更低频仅相位模型外推，不推断实际天线、海面传播或许可频段。

每个种子/海况/轨迹只重放原始导航、角距与FMCW数据并估计一次，再缓存两个时延策略的距离误差与真实距离，逐频率使用同一缓存，不为每个频率重新调参。缓存保留于每组geometry/，后续可在任意正频率重评估。

时延分别为原固定处理12ms（含半CPI后时隙末端183.395ms），以及上一轮同种子同场景的实测处理均值。后者在一条轨迹内恒定且跨频率一致，并非上一轮逐帧CPU时序精确重放。新的离线扫描耗时不算作每帧在线处理成本。

四组分别运行，每进程BLAS/OMP/MKL各1线程：

```powershell
cd D:\git\sim2
$env:OPENBLAS_NUM_THREADS='1'
$env:OMP_NUM_THREADS='1'
$env:MKL_NUM_THREADS='1'
D:\deepseek\python.exe -u kf_dpc_coherent_5nodes/frequency_sweep/run_full.py --trials 5 --seed-start 20261028 --output kf_dpc_coherent_5nodes/frequency_sweep/outputs/chunks/20261028
```

其他三组seed-start和输出尾目录分别为20261033、20261038、20261043。四组齐全后：

```powershell
D:\deepseek\python.exe kf_dpc_coherent_5nodes/frequency_sweep/report.py
```

全部完成后可直接重评估任意中间正频率，例如37.5、88.5、125MHz，不重新生成回波或重估轨迹：

```powershell
D:\deepseek\python.exe kf_dpc_coherent_5nodes/frequency_sweep/replay.py --mhz 37.5 88.5 125
```

旧模型全部完成后，在同一命令后加 `--include-legacy` 可同时重评估全部10个策略；旧模型只具有公共固定时序，工具不会补造其历史CPU均值时序。

旧模型完整重放同样分为20261028/33/38/43四组，每组trials=5，输出改用outputs/legacy_chunks/<seed>：

```powershell
D:\deepseek\python.exe -u kf_dpc_coherent_5nodes/frequency_sweep/run_legacy.py --trials 5 --seed-start 20261028 --output kf_dpc_coherent_5nodes/frequency_sweep/outputs/legacy_chunks/20261028
D:\deepseek\python.exe kf_dpc_coherent_5nodes/frequency_sweep/report_legacy.py
```

旧重放应生成1152000频率行和12800轨迹分段行；与主两方法的公共固定时序结果合并后1440000频率行、16000轨迹分段行。新增CV_CI_linear只改变原CV的浮标外推，以隔离波浪预测贡献。原137–177ms时序仅作300MHz完整复现，低频主比较全部统一143.395–183.395ms，不冒充真实CPU都满足此理想时序。早期不同导航/阵元假设的归档版本不混入公平排名。

应得到576000行：20种子×2海况×5轨迹×2估计器×2时延×360频率×2时段。时段为全部0–90s、75–90s。指标包括能量加权η、逐帧效率均值和5%分位、非相干参考增益；30MHz以上另报理想单位增益Friis参考功率，不是实地测量或硬件链路预算。低于30MHz绝对功率留空。

字段rx_ideal_model_w指“保留实际距离预测误差、但使用理想单位增益自由空间信道”的接收功率，不是完美相位下的接收功率上界。

本扫描未重建各频率下的通信带宽、调制、编码、噪声和数据率；125ms通信成本是受控比较假设，不证明极低频也能以相同速率传输原数据。

来源和推导见THEORY.md。降频不改变轨迹误差/覆盖率，也不消除200ms时限问题；不因η趋近1就宣称达到全部任务目标。结果在完整批次完成后由report.py生成，不以单组或首种子替代20种子结论。旧实验与其失败结果完整保留。

## 线上数据与文献

完整数据下载及恢复见 [PUBLICATION.md](PUBLICATION.md)，方法来源见 [LITERATURE.md](LITERATURE.md)，BibTeX见 [references.bib](references.bib)。
