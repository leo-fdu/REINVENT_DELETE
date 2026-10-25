# high=9、原始分数截断至9、生成重原子扣分0.01

本次为30任务独立RL运行：15靶点 × LibINVENT/LinkINVENT，seed=42，batch=100 ×10轮。
所有任务从原始prior开始，每任务1000次尝试，保留无效和重复。

```text
sigmoid_reward = sigmoid(min(PLANET_raw,9); low=4, high=9, k=0.25)
base_reward = sigmoid_reward - 0.01 * generated_heavy_atoms
training_reward = base_reward * (0.5 if repeated else 1)
```

生成重原子数为完整产物重原子数减固定输入片段重原子数，不包括H和连接标记*。
允许负奖励；重复系数乘整个基础奖励。无效或不可评分样本沿用零奖励。
原始分数保留未截断值。sigmoid(9)约为0.94676，并不是1。

配置：`configs/manual_planet_rl/experiment_high9_cap9_size001_1000.json`。
精简结果：`results/reinvent-high9-cap9-size001-1000-20261025/`。
完整本地原始镜像：`results/server-evidence/high9_cap9_size001_1000_seed42_20261025/run/`。

```bash
python evaluation/descriptor/generation_analysis_reinvent/high9_cap9_size001_1000_seed42_20261025/src/run_analysis.py
```

复用已有1000次分析的全部六类指标，并额外输出三类逐任务/分组曲线：
生成重原子数（评分成功样本均值）、base_reward（全部尝试均值）、training_reward（全部尝试均值）。
共9类指标，各包含30任务、3组总览以及3份分组副本，共324张PNG及324份曲线CSV。
原有33行总结表继续保留；`comparison_high9_baseline.csv`为与无截断无尺寸惩罚high9的对照。
`provenance.json`记录输入、共享分析代码与基线表的哈希。

本次相对high9基线同时增加原始分数截断和尺寸扣分，不能将变化单独归因于尺寸项。
单种子观察不代表跨种子稳定性，也不代替独立结合评估。已有DRD3设计配体与辅因子局限仍适用。
云端沿用已验证的环境，PLANET进程关闭CUDA缓存，生成训练进程使用默认分配器。
