# PLANET high=8 独立重跑：每任务 1,000 次采样

在 high=9 实验基础上，仅将 sigmoid high 改为 8；low=4、k=0.25、seed=42。
15 个靶点 × LibINVENT / LinkINVENT，每任务从原始 prior 开始，batch=100 × 10 轮。
30,000 次采样尝试保留无效和重复，不保证每任务 1,000 个有效唯一分子。ACES 仍跳过。
目标结构、片段、模型、DAP 参数、重复惩罚和云端运行环境均沿用 high=9。
既有 DRD3 设计配体问题和 PLANET 不读取辅因子的局限仍适用。

配置：`configs/manual_planet_rl/experiment_high8_1000.json`。
精简数据：`results/reinvent-high8-1000-20261025/`。
完整原始镜像：`results/server-evidence/high8_1000_seed42_20261025/run-v2/`。

```bash
python evaluation/descriptor/generation_analysis_reinvent/high8_1000_seed42_20261025/src/run_analysis.py
```

复用 `generation_analysis_reinvent_1000/src/run_analysis.py` 全部六类分析，严格校验
30 任务 × 10 轮 × 100 行和实际实验参数。输出原始 PLANET 分数、sigmoid 奖励、
重原子数、MW、批内重复数、无效分子数的曲线及数据；汇总包括 30 任务及 3 个均值行。
`analysis/comparison_high9_baseline.csv` 给出与 high=9 的逐任务、逐指标差值；
`analysis/provenance.json` 记录输入、共享分析脚本及基线总结表的 SHA-256。

high 不是原始分数的截断上限。不同 sigmoid 的奖励均值不能直接用来比较亲和力。
这是同预算、同种子的对照；单次结果不替代多种子验证和独立结合评估。

## 显存问题与重跑记录

首次运行在 CP3A4 LinkINVENT 第 3 轮更新时 CUDA OOM。实际样本重放显示 PLANET
活跃张量约 25.83 MiB，而保留缓存达到 19,746 MiB；调用 empty_cache 后缓存降为 32 MiB。
正式采用 `planet_python_nocache.sh` 仅为 PLANET 设置
`PYTORCH_NO_CUDA_MEMORY_CACHING=1`，生成训练进程仍使用默认分配器。
评分批量 32、RL 批量 100、模型及数学计算均不改动。`memory-probe.json` 保留实测。
在 `run-v2/` 从 prior 完整重跑全部 30 任务，不拼接失败轨迹或续训失败 checkpoint。
初次 `run/` 保留作为失败证据，不计入最终 30,000 次采样统计。
实际重跑脚本为 `run_high8_cloud_v2.sh`；原启动脚本仍保留。
