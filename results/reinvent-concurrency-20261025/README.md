# RTX 4090 并发基准

目标是在同一张24 GB RTX 4090、10核CPU上测量实际吞吐和GPU利用率，确定适合现有
REINVENT + PLANET流程的并发数。并发数指独立RL任务数，CPU线程数是另一参数。

固定工作量：CP3A4、EGFR、JAK2、AA2AR各LibINVENT/LinkINVENT，共8任务，
每任务batch=100 ×10轮、seed=42，均从原始prior开始。沿用cap9 + size0.01奖励，
PLANET批量32并关闭CUDA缓存，生成训练进程使用默认分配器。

每个任务运行在独立目录、使用独立PLANET服务端口。任务队列顺序一致；同时运行的任务
数分别测试1、2、3、4，先固定每进程OMP/MKL/OpenBLAS线程数为2，再依据测量结果
决定是否测试其他配置。遇到失败会停止该轮并中止更高并发。

`src/benchmark_manual_planet_concurrency.py`负责隔离、调度和逐秒采样；
`run_benchmark.sh`记录首轮实际云端命令。计时包含worker启动、模型加载、采样、评分、
优化、checkpoint写入及任务切换；准备输入、独立结果验收和串行输出对照不计入耗时。
GPU利用率、显存、功率来自nvidia-smi；CPU忙碌比例来自/proc/stat，占整机CPU百分比。
`gpu_memory_peak_mib`是1秒采样峰值，可能低于短时真实峰值；不能单靠它证明不会OOM。

每档使用完全相同的8000次尝试，记录吞吐attempts/s及耗时。满并发窗口另行汇总GPU
利用率，避免将不同档的排队尾部当作纯计算效率。nvidia-smi利用率表示GPU忙碌程度，
不等价于CUDA核心计算吞吐或理论算力百分比。推荐并发以完整工作量吞吐、资源余量和
稳定性为依据，不以瞬时100%为唯一目标。

各成功任务进行原奖励验收，并与本次1路基线逐行比较SMILES及原始分数/训练奖励。
完整临时模型和输出保留在云端`/home/ubuntu/benchmarks/reinvent_concurrency_20261025/`。
本地原始证据放在已有忽略目录`results/server-evidence/reinvent_concurrency_20261025/`，
精简summary和逐秒资源记录保存在本目录。

本基准每任务独立启动一个PLANET服务，而正式串行流程通常每个靶点的两个模式复用服务；
因此不能直接把8任务耗时乘30/8当作全30任务的精确预测。
