# REINVENT 生成任务分析（generation_analysis_reinvent_1000）

对 REINVENT4（LibINVENT / LinkINVENT）在 15 个靶点上的 30 个生成任务做批量统计，
产出一张总结表和 6 类 batch 折线图（共 198 张图及其对应 CSV 数据表）。

本目录是已有生成轨迹的前 1,000 个分子分析，不是独立重新生成的实验。

## 数据源

- `results/reinvent-20260927/all_generated.csv.gz`（300,000 行 = 30 任务 × 100 batch × 100 分子）
- 先导分子（lead molecule）：`real-world_dataset/<target>/crystal.mol2`（共晶配体）

## 运行方式

```bash
python3 evaluation/descriptor/generation_analysis_reinvent/generation_analysis_reinvent_1000/src/run_analysis.py
```

依赖：`rdkit`、`matplotlib`、`numpy`（无随机过程，结果可完全复现）。

## 目录结构

```
generation_analysis_reinvent_1000/
├── README.md
├── src/run_analysis.py
└── analysis/
    ├── summary_table.csv                  # 33 行 × 11 指标列的总结表
    ├── average_oracle_score_origin/       # 以下 6 个指标文件夹结构相同
    │   ├── figure/                        #   33 张 PNG：30 单任务 + all_tasks + linkinvent + libinvent
    │   └── table/                         #   33 个 CSV：与图一一对应
    ├── average_oracle_score_sigmoid/      # {figure, table}
    ├── average_atom_number/               # {figure, table}
    ├── average_MW/                        # {figure, table}
    ├── repetition_count/                  # {figure, table}
    ├── noncanonical_count/                # {figure, table}
    ├── summary/    {figure, table}        # all_tasks 图/表的额外副本（30 条折线）
    ├── linkinvent/ {figure, table}        # linkinvent 图/表的额外副本（15 条折线）
    └── libinvent/  {figure, table}        # libinvent 图/表的额外副本（15 条折线）
```

命名约定：单任务文件 `<target>_<mode>.png/.csv`（如 `aa2ar_libinvent.png`）；
多任务图 `all_tasks`（30 线）、`linkinvent` / `libinvent`（各 15 线）。

## 通用口径

- 每任务按原始 CSV 出现顺序保留前 1,000 行，包含不合法和重复分子，不筛选或去重。
  已核对这等于 step 1–10，每步 100 行；总计 30,000 行。
- 所有均值、重复率和召回次数仅使用这些行；final 指 step 10。
- 不合法定义为 `smiles_state == 0`；state 2 是合法的批内重复分子。
- 重复计数仍为每个 SMILES 的 n−1；总结表全局重复率覆盖所选十批，曲线为单批计数。
- oracle 均值只计 `planet_status == "ok"`；重原子数及 MW 只计合法且可解析的分子。
- sigmoid 使用原始 `sigmoid_reward` 列。缺失均值留空。
- 单任务图 8×5 英寸、多任务图 12×7 英寸，300 dpi；横轴为整数 batch 1–10。

## 总结表（analysis/summary_table.csv）

33 行 = 30 个任务（`<target>_<mode>`，靶点按字母序，每靶点先 libinvent 后 linkinvent）
+ `libinvent_average` + `linkinvent_average` + `average`（最后三行）。
CSV 第一列 `task` 为行名，其余 11 个指标列依次为：

| 列 | 定义 |
|---|---|
| `oracle_score_original` | 任务内全部 batch 的平均原始 planet 分数（对所有打分成功的分子求均值） |
| `final_oracle_score_original` | 最后一个 batch 的平均原始 planet 分数 |
| `best_oracle_score_original` | 各 batch 平均原始 planet 分数的最大值（"最好 batch 的均分"，非单分子最高分） |
| `oracle_score_sigmoid` | 任务内全部 batch 的平均 sigmoid 分数 |
| `final_oracle_score_sigmoid` | 最后一个 batch 的平均 sigmoid 分数 |
| `best_oracle_score_sigmoid` | 各 batch 平均 sigmoid 分数的最大值 |
| `noncanonical_rate` | 全局不合法 SMILES 数 / 任务总生成分子数 |
| `final_noncanonical_rate` | 最后一个 batch 内不合法 SMILES 占比 |
| `repetition_rate` | 全局重复数（n−1 规则）/ 任务总生成分子数 |
| `final_repetition_rate` | 最后一个 batch 内重复数（n−1 规则）/ 该 batch 分子数 |
| `recall_count` | 先导分子被生成的次数（见下） |

三个平均行取对应 15/30 个任务值的**算术平均**（含 `recall_count`，即"平均召回次数"）。
数值统一保留 6 位小数。

## recall_count 与先导分子匹配协议

- 先导分子 = `real-world_dataset/<target>/crystal.mol2`，由
  `evaluation/descriptor/similarity_test/src/crystal_actives_similarity.py` 的
  `load_crystal_ligand()` 加载（该函数处理了这批 MOL2 的电荷归零、芳香性修复、
  kekulize，并从 3D 坐标指认立体化学），与本仓库其他分析保持同一转换口径。
- 匹配方式：**严格匹配**——生成分子与先导分子的 RDKit 异构 canonical SMILES
  完全相等（含立体化学）。计数覆盖全部可解析生成分子（含 `smiles_state == 2`
  的重复采样行，重复生成即多次计数），即"先导分子被生成的次数"。
- 本子集召回次数由原始分子重新计算，见总结表。

## 折线图（6 指标 × 33 张）

横坐标为 batch（1–10，10 点连续曲线），纵坐标指标：

| 指标 | 定义 |
|---|---|
| `average_oracle_score_origin` | batch 内平均原始 planet 分数（仅打分成功分子） |
| `average_oracle_score_sigmoid` | batch 内平均 sigmoid 分数（仅打分成功分子） |
| `average_atom_number` | batch 内生成分子平均重原子数（仅合法可解析分子） |
| `average_MW` | batch 内生成分子平均分子量（仅合法可解析分子） |
| `repetition_count` | batch 内重复 SMILES 数量（n−1 规则，batch 内计数） |
| `noncanonical_count` | batch 内不合法 SMILES 数量 |

每个指标 33 张图：30 个单任务图（1 条折线）+ `all_tasks`（30 线）+
`linkinvent`（15 线）+ `libinvent`（15 线）。
多线图的宽表 CSV：第一列 `batch`，其后每个任务一列（列名 = 任务名）。
