# REINVENT 生成任务分析（generation_analysis_reinvent）

对 REINVENT4（LibINVENT / LinkINVENT）在 15 个靶点上的 30 个生成任务做批量统计，
产出一张总结表和 6 类 batch 折线图（共 198 张图及其对应 CSV 数据表）。

## 数据源

- `results/reinvent-20260927/all_generated.csv.gz`（300,000 行 = 30 任务 × 100 batch × 100 分子）
- 先导分子（lead molecule）：`real-world_dataset/<target>/crystal.mol2`（共晶配体）

## 运行方式

```bash
python3 src/run_analysis.py
```

依赖：`rdkit`、`matplotlib`、`numpy`（无随机过程，结果可完全复现）。

## 目录结构

```
generation_analysis_reinvent/
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

- **batch** = 数据中的 `step` 列（1–100），每个 batch 含 100 个生成分子；
  **final** 一律指最后一个 batch（step 100）。
- **不合法（noncanonical/invalid）SMILES** = 数据 `smiles_state` 列为 `0` 的行。
  该列是 REINVENT4 的 `SmilesState` 枚举（见
  `REINVENT4/reinvent/models/model_factory/sample_batch.py`）：
  `0 = INVALID`（8,610 条，2.87%，无法解析/kekulize 失败）、
  `1 = VALID`（289,182 条）、
  `2 = DUPLICATE`（2,208 条，**合法分子**在采样批内的重复出现）。
  因此不合法统计只取 `state == 0`，不把重复分子计入不合法。
- **重复（repetition）**：同一 SMILES 字符串出现 n 次（n≥2）记重复 n−1 次。
  字符串以数据 `smiles` 列为准（REINVENT 规范化后的 SMILES；不合法行则是原始采样串）。
  注：REINVENT 自带的 `smiles_state == 2`（DUPLICATE）是采样批内去重标记，
  与本口径（批内所有字符串的 n−1 计数）接近但不必相等，本文一律使用上述 n−1 规则。
  - `repetition_rate`（总结表）：任务内**跨全部 100 个 batch 全局**计数，除以任务总分子数。
  - `final_repetition_rate`（总结表）：仅在 step 100 的 batch 内计数，除以该 batch 分子数。
  - `repetition_count`（折线图）：**单个 batch 内**的重复数量（跨 batch 的重复不计入）。
- **oracle 打分缺失处理**：8,610 行（2.87%）`planet_status` 为空（PLANET 打分失败，
  恰好对应 `smiles_state == 0` 的不合法分子）。所有分数均值**只对
  `planet_status == "ok"` 的行计算**（缺失行不计入分子也不计入分母），
  不会按 0 分处理。不合法 SMILES 的计数与占比不受影响（仍以全部生成分子为分母）。
- **重原子数 / 分子量**：RDKit 计算（`GetNumHeavyAtoms()` / `Descriptors.MolWt`，
  平均分子量），只统计合法（`smiles_state ∈ {1, 2}`）且可被 RDKit 解析的分子。
- **sigmoid 分数**：直接使用数据中的 `sigmoid_reward` 列
  （即本次 run 配置的 sigmoid 变换：low=4, high=16, k=0.25）。
- 若某 batch 无任何可计入的分子，对应 CSV 单元格留空、折线断开（实际数据中未发生）。

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
- 已验证：改用忽略立体化学的连接性匹配，30 个任务的匹配结果与严格匹配**完全相同**
  （凡被召回的分子立体化学均与晶体一致），故两种协议等价，采用表述更严格的前者。
- 试点观察：召回集中在 batch 1–2；`aa2ar_linkinvent` 召回 67 次显著高于其他任务
  （其切割的 linker 过小、先导容易被整体重拼出）。因此本表只统计召回**次数**，
  不再记录召回轮次（recall_batch），该列信息量低且易被切割方式混淆。
- 各任务召回次数（供核对）：`aa2ar_libinvent`=1、`aa2ar_linkinvent`=67、
  `dpp4_libinvent`=3、`egfr_linkinvent`=3、`gria2_libinvent`=4，其余任务为 0。

## 折线图（6 指标 × 33 张）

横坐标为 batch（1–100，100 点连续曲线），纵坐标指标：

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
