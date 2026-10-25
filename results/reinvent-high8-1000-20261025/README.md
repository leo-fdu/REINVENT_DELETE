# high=8、每任务 1,000 次 RL 采样

本次在 high=9 对照基础上仅改变 sigmoid high=8；low=4、k=0.25、seed=42，
batch=100 × 10 轮，30 个任务各自从原始 prior 开始。无效和重复尝试保留。
复用同一云端项目和 Python 环境，输出位于独立 high8 目录，没有续训 high9 checkpoint。

- `all_generated.csv.gz`：完整 30,000 次尝试的压缩表。
- `experiment.json`、`run.json`：实际参数、输入与运行配置。
- `acceptance-report.json`：预算、编号、输入哈希、评分公式、惩罚、oracle 与 checkpoint 验收。
- `training-numerics.json`：最终权重及训练数值有限性、权重实际更新检查。
- `checksums.json`：云端完整 run 目录的逐文件哈希；相对原始镜像路径。
- `transfer-verification.json`：本地完整镜像的逐文件校验结果。
- `first-batch-baseline-check.json`：与 high9 首轮采样和原始分数比较。
- `environment.json`、`environment-freeze.txt`、`runtime-requirements.txt`：运行环境。
- `source-transfer-manifest.json`：复用且核验过的源码、输入和模型，加本轮配置及验收脚本。
- `pocket-validation-report.json`：15 个靶点的实际 PLANET 口袋检查。
- `run_high8_cloud.sh`：实际执行脚本；引用上轮云端项目及环境的绝对路径。

完整原始数据（含 30 个 checkpoint、CSV、日志和 TensorBoard）保留在：
`results/server-evidence/high8_1000_seed42_20261025/run-v2/`，沿用现有忽略规则。
分析与结论见：
`evaluation/descriptor/generation_analysis_reinvent/high8_1000_seed42_20261025/`。

## 显存问题与重跑记录

首次运行在 CP3A4 LinkINVENT 第 3 轮更新时 CUDA OOM。实际样本重放显示 PLANET
活跃张量约 25.83 MiB，而保留缓存达到 19,746 MiB；调用 empty_cache 后缓存降为 32 MiB。
正式采用 `planet_python_nocache.sh` 仅为 PLANET 设置
`PYTORCH_NO_CUDA_MEMORY_CACHING=1`，生成训练进程仍使用默认分配器。
评分批量 32、RL 批量 100、模型及数学计算均不改动。`memory-probe.json` 保留实测。
在 `run-v2/` 从 prior 完整重跑全部 30 任务，不拼接失败轨迹或续训失败 checkpoint。
初次 `run/` 保留作为失败证据，不计入最终 30,000 次采样统计。
实际重跑脚本为 `run_high8_cloud_v2.sh`；原启动脚本仍保留。
