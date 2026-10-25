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

## 实测结果

正式30/30任务完成，共30,000次尝试；29,722条有效且PLANET评分成功，278条无效
（0.9267%），没有有效分子评分失败。全过程任务内重复3,954条（13.18%）。
300条训练奖励为负，均进入DAP，没有截零；sigmoid最大值0.94675976，符合cap=9。
全部30个checkpoint的参数有限且相对prior实际更新，训练likelihood及奖励均有限。

下表为30任务等权平均；两组均从相同prior开始，seed=42，batch=100×10轮：

| 指标 | high9原协议 | cap9 + size0.01 |
|---|---:|---:|
| 全过程原始PLANET均值 | 8.999649 | 8.630659 |
| 第10轮原始PLANET均值 | 9.933194 | 9.248891 |
| 第10轮sigmoid均值（扣尺寸前） | 0.935063 | 0.890226 |
| 第10轮重原子数 | 46.064178 | 41.210621 |
| 第10轮MW | 651.842214 | 583.698797 |
| 无效比例 | 1.3067% | 0.9267% |
| 全过程重复比例 | 10.1767% | 13.1800% |

末轮重原子数下降约10.54%，MW下降约10.46%；30任务中27个重原子数下降，26个MW下降。
平均重原子数仍从首轮30.46增至末轮41.21，MW从432.09增至583.70，所以增长减弱但未消失。
新增部分重原子数均值从11.59增至22.34；末轮基础奖励均值0.659394，
施加重复系数后的训练奖励均值0.639422（这两个均值都包含无效尝试的零奖励）。
11,119条原始PLANET分数仍超过9：截断的是奖励计算输入，不是模型原始预测值。

与high9首批3,000个样本的SMILES和原始分数完全一致，见精简结果目录的
`first-batch-baseline-check.json`。全部采样预算、公式各项、原子数及oracle身份均通过独立验收。
324张PNG已检查完整性，全部曲线为10轮，33行总结表及high9对照表已生成，
MW与实际训练奖励总览图已目视抽查。

完整原始镜像约2.5 GB，349个文件已逐一通过云端SHA-256清单核对，
详见精简结果目录的`transfer-verification.json`。云端驱动正常退出（exit=0）。
