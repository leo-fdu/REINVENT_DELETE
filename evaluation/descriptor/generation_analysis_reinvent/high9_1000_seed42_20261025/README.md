# PLANET high=9 独立重跑：每任务 1,000 次采样

本次为独立的 30 任务 RL 实验，不是旧轨迹的前 1,000 行截取。

- 15 个靶点 × LibINVENT / LinkINVENT；ACES 仍跳过。
- 每任务从原始 prior 开始，batch=100、10 轮，共 1,000 次采样尝试。
- 全部 30,000 行均保留无效和重复样本，不保证 1,000 个有效唯一分子。
- sigmoid：low=4、high=9、k=0.25；seed=42。
- DAP sigma=128、rate=0.0001；重复惩罚乘数 0.5，其余参数沿用原协议。
- 输入采用现有手动设计片段与 JAK2/DRD3 受体/口袋修正；不改变设计配体。
- DRD3 设计配体身份问题、PLANET 不读取辅因子等既有局限仍适用。

配置：`configs/manual_planet_rl/experiment_high9_1000.json`。
原始结果：`results/reinvent-high9-1000-20261025/`。

```bash
python evaluation/descriptor/generation_analysis_reinvent/high9_1000_seed42_20261025/src/run_analysis.py
```

复用相邻 `generation_analysis_reinvent_1000/src/run_analysis.py` 的全部分析，
保持先导分子加载、有效性、重复、召回及均值口径一致。新入口在分析前严格校验
30 任务 × 10 轮 × 100 行，拒绝截断额外数据；校验实际实验参数。

产出包含 33 行总结表、6 类指标的 198 张图及对应数据表，以及 18 份分组图/表副本：
原始 PLANET 均分、sigmoid 奖励均分、重原子数、分子量、批内重复数、无效分子数。
`analysis/provenance.json` 记录输入和共享分析源码的 SHA-256。

与旧实验的完整 100 轮终点比较时，同时存在奖励转换和训练轮数差异。
比较 high 的影响优先使用旧实验同样的前 10 轮分析；即使如此仍为单随机种子实验。
