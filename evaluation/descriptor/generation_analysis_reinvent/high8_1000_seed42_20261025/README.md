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

## 本次结果与 high=9 对照

30/30 个任务完成，共 30,000 次尝试；29,637 条有效且评分成功，363 条无效（1.21%），
没有有效分子评分失败。任务内跨轮及批内重复共 3,463 条（11.5433%）。
30 个最终 checkpoint 均为有限参数且相对 prior 实际更新，训练 likelihood 与奖励全部有限。
预算、编号、输入哈希、sigmoid、重复惩罚与 oracle 身份通过独立验收。

以下连续指标均为 30 个任务等权平均，均使用 10 轮 / 1,000 次尝试预算：

| 指标 | high=9 | high=8 |
|---|---:|---:|
| 全过程原始 PLANET 均值 | 8.999649 | 8.751825 |
| 第 10 轮原始 PLANET 均值 | 9.933194 | 9.662877 |
| 第 10 轮 sigmoid 均值 | 0.935063 | 0.965130 |
| 第 10 轮重原子数 | 46.064178 | 45.363283 |
| 第 10 轮 MW | 651.842214 | 641.975361 |
| 无效比例 | 1.3067% | 1.2100% |
| 任务内全过程重复比例 | 10.1767% | 11.5433% |

high 从 9 降为 8 后，末轮平均重原子数与 MW 仅下降约 1.5%，全过程重复比例
反而增加约 1.37 个百分点。平均重原子数仍从首轮 30.46 增至末轮 45.36，
MW 从 432.09 增至 641.98：这次单种子试验没有显示尺寸增长问题已解决。
奖励曲线更容易饱和不代表真实亲和力更好；19,470 条 sigmoid 分数超过 0.95。
原始分数超过 8 的记录仍有 19,811 条，超过 9 的有 11,958 条，说明 high 不是硬上限。

默认缓存与关闭缓存的 11 个可比完整任务共 11,000 行逐字段完全一致，见
`results/reinvent-high8-1000-20261025/allocator-invariance-check.json`。
六类分析共 216 张 PNG（含分组副本）和 216 张曲线数据表，33 行总结表与
逐指标 high9 对照表均已生成。所有 PNG 已检查完整性，原始分数和 MW 总览图已目视抽查。

本轮与 high9 的首批共 3,000 个样本，SMILES 与原始 PLANET 分数完全一致，
见精简结果目录的 `first-batch-baseline-check.json`。

完整原始镜像约 2.5 GB，350 个文件均已通过云端 SHA-256 清单核对，
见 `results/reinvent-high8-1000-20261025/transfer-verification.json`。
