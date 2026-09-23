# 16 个靶点的二维任务设计

打开 [index.html](index.html) 浏览全部图；也可直接查看各靶点 PNG。每张图上半为 LibINVENT，下半为 LinkINVENT，共 32 个设计。这里只交付图、说明和复现绘图脚本，没有生成 `.smi`、SDF/MOL2 片段或运行配置。

## 输入规则与图例

- **LibINVENT**：一个连续骨架，在固定的出口放一个 `*`；本设计统一使用单出口，便于与 DELETE growing 比较。模型生成另一侧装饰片段，不将原装饰片段作为输入。
- **LinkINVENT**：两个独立且各带一个 `*` 的片段。实际 REINVENT4 记录用 `片段A|片段B` 的形式；不是点号分隔的 SMILES。图中竖线表达该分隔符。生成 linker 有两个出口，不能把一个连续骨架上的两个出口当作这两个输入片段。
- 蓝色/绿色为保留原子，橙色为原配体中待替换区域（参考答案）；红色表示切断键。所有切断均为非环单键，不破坏芳香环或切断双键。`*` 是 dummy attachment atom，不是真实碳原子、不计重原子数。
- 图中连接点不加化学同位素。切断键使用原 MOL2 的 **1-based atom ID**，图上在对应端点显示这些 ID。LibINVENT 内部重新编号、LinkINVENT 的片段顺序与 linker 出口配对由实际 sampler 处理，不能把图上的原子 ID 直接当作模型的连接点编号。

## 任务价值与比较范围

这些是按结构提出的初始任务，不是经蛋白接触分析确认的最佳药效团切分。没有运行生成、打分、对接或亲和力验证。每个靶点的理由见图与图册。

- **ACES**：保留完整多环骨架和末端甲基，替换乙基中的一个 CH2。模型格式与拼接验证通过，但甲基只有 1 个重原子，不能视为有充分结合贡献的独立片段。保留为边界案例，建议主 linking 汇总排除它并另报覆盖率；若强制要求两个较大的片段，应另选已知配体，而不是在此多环骨架上硬切环。
- **ABL1、CDK2、DPP4、EGFR、GRIA2**：参考 linker 为单个 NH、CH2 或羰基（GRIA2 的另一端是较小膦酸基）。格式有效，单独报告短 linker 层；不能与长链重建任务直接混为一种难度。
- **JAK2、PARP1**：参考 linker 含环，属于含环连接区/核心替换任务。LinkINVENT 支持其连接拓扑，但生成难度和固定药效团范围与柔性链替换不同，应单列结果。
- **DYR linking**：输入满足连接点和词表要求，参考答案可恢复全部键连接；但当前 REINVENT4 拼接实现丢失该切口对应的烯键 E/Z 信息。它不是立体化学完全还原的任务，需另加几何异构体约束或单列结果。LibINVENT 方案已选用远离烯键的切口，保留参考 E/Z。
- 其余任务同样记录实际保留和删除重原子数。跨模型比较必须让 DELETE 与相应 REINVENT4 任务使用完全相同的原子保留集合；DELETE 的片段坐标将来应从原晶体坐标提取，不能使用二维绘图坐标。
- 原配体恢复仅作完整性检查；正式基准应另外报告新颖性，避免把恢复已知配体当作创新。模型输入的二维结构并不自行施加蛋白结合位置或三维出口方向约束。

## 化学结构核对

原始 MOL2 没有显式氢，RDKit 直接读取有 7 个无法 sanitize，其他若干会出现异常自由基、价态或电荷。不能使用 `sanitize=False` 得到的图就宣称可作为模型输入。

本图册对 15 个配体采用 RCSB CCD 字典的完整配体参考结构，在**元素和重原子邻接图完全一致**的前提下映射回原 MOL2 原子 ID；用字典的键级、氢、形式电荷和立体化学绘图。CCD 理想坐标没有用于任务片段空间约束。参考字符串及原文件哈希嵌入脚本，复现不需要联网。原始数据未修改。

这不是对原 MOL2 进行无条件认可：例如 CP3A4、DYR、KIF11 等的键级与字典有差异；ADRB1 字典使用其记录中的亚胺/非芳香五元环互变异构表示；羧酸、膦酸与胺使用参考结构的电荷态。这里采用可解析的参考状态，并不等同于某一 pH 的主导质子化态。以后真正构建输入前，应确认两模型统一采用同一个化学标准化方案。

**DRD3 特例**：本地文件原子 1 是连接芳环的 C，3PBL 的 ETQ 字典相应位置为 Cl，元素图不能完整匹配。没有用 ETQ 替换本地配体。本图保留本地全部 23 个重原子的元素、键级和连接，补足正常价态的隐式氢，胺采用中性表示，手性取自本地三维坐标。图的任务设计对本地结构成立，但不可把它直接称作标准 ETQ 共晶配体。

## 验证

绘图脚本对全部 32 个方案执行：

1. 原 MOL2 SHA-256 与本次审阅版本一致；所有参考结构可 sanitize，无异常自由基。
2. 非环单键切分；LibINVENT 得到 2 个连续片段，固定端 1 个出口；LinkINVENT 得到 3 个连续片段，两个固定端各 1 个出口，待生成区 2 个出口。
3. 原子无遗漏或重复，dummy 均为度 1、单键连接。
4. 使用本仓库 `bond_maker.join_scaffolds_and_decorations` 和连接点编号逻辑，将参考答案拼回，比较去除 atom-map 后的 canonical isomeric SMILES；31 个方案完全一致，DYR linking 仅连接图一致、E/Z 丢失，已单独标注。单原子 linker 也经过实际拼接检查。
5. 使用 Transformer sampler 所用的标准化方法，确认条件片段化学结构不变；检查输入 token 和参考输出 token 全部存在于本地 prior 的词表中，输入 token 长度不超过 180。仅检查词表/结构，不加载网络执行生成。

当前验证环境：RDKit `2026.03.5`，PyTorch `2.12.0`。
- `libinvent`：metadata model_id `7c8e104e81124bafb5f62d14c64ec9f9`；max_sequence_length `180`。
- `linkinvent`：metadata model_id `9ef8cd169df24d58a16ee58b882d20af`；max_sequence_length `180`。

| 靶点 | Lib：保留 + 生成重原子 | Link：端 A + 端 B + 生成重原子 | Lib / Link 输入 token |
|---|---:|---:|---:|
| [aa2ar](aa2ar.png) | 16 + 9 | 7 + 16 + 2 | 30 / 43 |
| [abl1](abl1.png) | 21 + 8 | 8 + 20 + 1 | 40 / 55 |
| [aces](aces.png) | 19 + 2 | 1 + 19 + 1 | 37 / 40 |
| [adrb1](adrb1.png) | 17 + 4 | 12 + 5 + 4 | 30 / 34 |
| [akt1](akt1.png) | 18 + 12 | 7 + 18 + 5 | 34 / 47 |
| [cdk2](cdk2.png) | 23 + 7 | 14 + 15 + 1 | 38 / 50 |
| [cp3a4](cp3a4.png) | 43 + 7 | 29 + 18 + 3 | 77 / 83 |
| [dpp4](dpp4.png) | 18 + 13 | 16 + 13 + 2 | 37 / 57 |
| [drd3](drd3.png) | 15 + 8 | 14 + 7 + 2 | 29 / 39 |
| [dyr](dyr.png) | 22 + 2 | 8 + 11 + 5 | 42 / 34 |
| [egfr](egfr.png) | 25 + 8 | 7 + 25 + 1 | 40 / 53 |
| [gria2](gria2.png) | 21 + 6 | 22 + 4 + 1 | 45 / 56 |
| [jak2](jak2.png) | 21 + 12 | 11 + 12 + 10 | 41 / 46 |
| [kif11](kif11.png) | 25 + 8 | 8 + 21 + 4 | 48 / 51 |
| [parp1](parp1.png) | 19 + 6 | 12 + 6 + 7 | 38 / 34 |
| [ppara](ppara.png) | 23 + 10 | 16 + 13 + 4 | 46 / 58 |

运行方式（从仓库根目录；需已安装 RDKit、PyTorch、Pillow 和 REINVENT4 依赖）：

```bash
/opt/homebrew/Caskroom/miniforge/base/envs/reinvent4/bin/python real-world_dataset/task_design_astra/render_designs.py
```

字体路径 `FONT` 为本机 Arial Unicode；其他机器可改为支持中文的字体。脚本只重写本文件夹内的展示文件，所有候选模型输入字符串仅在内存中用于校验。

## 规范来源

以本地代码为准：

- `REINVENT4/reinvent/runmodes/samplers/libinvent.py`
- `REINVENT4/reinvent/runmodes/samplers/linkinvent.py`
- `REINVENT4/reinvent/chemistry/library_design/attachment_points.py`
- `REINVENT4/reinvent/chemistry/library_design/bond_maker.py`
- `REINVENT4/reinvent/models/transformer/core/vocabulary.py`
- [REINVENT4 官方文档](https://github.com/MolecularAI/REINVENT4/blob/main/contrib/reinvent-doc/index.md)
- [REINVENT4 论文：LibINVENT 骨架与 LinkINVENT 双片段条件](https://pmc.ncbi.nlm.nih.gov/articles/PMC10882833/)

RCSB 数据于 2026-09-22 查阅，CCD 原始下载地址形式为 `https://files.rcsb.org/ligands/download/<CCD>_ideal.sdf`。这些下载文件仅用于核对，没有作为交付输入文件保存。

- `aa2ar`：[3EML](https://www.rcsb.org/structure/3EML) / [ZMA](https://www.rcsb.org/ligand/ZMA)；原始 SHA-256 `fb4f772045e80d15e919fb0cc748f4e13429e73d9b6d284ab203f51233042cc7`。
- `abl1`：[2HZI](https://www.rcsb.org/structure/2HZI) / [JIN](https://www.rcsb.org/ligand/JIN)；原始 SHA-256 `351685615a31721c552492a3af427e4e5ac55f50061a67a41917457403efc0e6`。
- `aces`：[1E66](https://www.rcsb.org/structure/1E66) / [HUX](https://www.rcsb.org/ligand/HUX)；原始 SHA-256 `1703c50c57a2c24592b2e6fb571be4b7aa3bf1604e8d1ebd591b231cbfb2631e`。
- `adrb1`：[2VT4](https://www.rcsb.org/structure/2VT4) / [P32](https://www.rcsb.org/ligand/P32)；原始 SHA-256 `cf21c52fd56a44a837fc25989a185f07d3880e3e158e6a9e08a80c622ad0d0c7`。
- `akt1`：[4GV1](https://www.rcsb.org/structure/4GV1) / [0XZ](https://www.rcsb.org/ligand/0XZ)；原始 SHA-256 `4ad33e8a4308897c232f5b08213a1917a17b29a307c5e7fa6827338a9e6098ae`。
- `cdk2`：[1H00](https://www.rcsb.org/structure/1H00) / [FAP](https://www.rcsb.org/ligand/FAP)；原始 SHA-256 `721c460ebf2db564aae3b85953a5548971b476fa8e0c0863adfe1c2c67c57c6a`。
- `cp3a4`：[3NXU](https://www.rcsb.org/structure/3NXU) / [RIT](https://www.rcsb.org/ligand/RIT)；原始 SHA-256 `af38b4da66f376d216054a89849690ad6d90d6150d6e05ee99b7c65f11099bdb`。
- `dpp4`：[2I78](https://www.rcsb.org/structure/2I78) / [KIQ](https://www.rcsb.org/ligand/KIQ)；原始 SHA-256 `c66e797c433aa245a2ac2af30527e0f12178f38c08dcbe1aecf247d564b0a843`。
- `drd3`：[3PBL](https://www.rcsb.org/structure/3PBL) / [ETQ](https://www.rcsb.org/ligand/ETQ)；原始 SHA-256 `d6168ae6539f9de16c471c1dc2696757d6dd5269c65bfd122650aa93f80e2594`。
- `dyr`：[3NXO](https://www.rcsb.org/structure/3NXO) / [D2B](https://www.rcsb.org/ligand/D2B)；原始 SHA-256 `16136044517330b22775c8c859d1a39ecf431d97b5bb2f5173ab99d9eb68a142`。
- `egfr`：[2RGP](https://www.rcsb.org/structure/2RGP) / [HYZ](https://www.rcsb.org/ligand/HYZ)；原始 SHA-256 `c19cce10b32333f6a4a03b7e072865adc1c40ab5c83f8ad8bdc27882c3863bf8`。
- `gria2`：[3KGC](https://www.rcsb.org/structure/3KGC) / [ZK1](https://www.rcsb.org/ligand/ZK1)；原始 SHA-256 `1ef46e9db70370659bad4a1a32c36ca6f14594cc719653791534f54837a2cea1`。
- `jak2`：[3LPB](https://www.rcsb.org/structure/3LPB) / [NVB](https://www.rcsb.org/ligand/NVB)；原始 SHA-256 `cf672194d8e746a28a01968b8b30124b89423340d3b9c61e7ddfdad9dc3a80a6`。
- `kif11`：[3CJO](https://www.rcsb.org/structure/3CJO) / [K30](https://www.rcsb.org/ligand/K30)；原始 SHA-256 `70a9d35ee57ce2752f7e707d63102b5eefc9fb916ce59a871b5d2ddcd316b8fa`。
- `parp1`：[3L3M](https://www.rcsb.org/structure/3L3M) / [A92](https://www.rcsb.org/ligand/A92)；原始 SHA-256 `495f20545fafcd06d76c0edf88aec52d0ec85337ab6806ea54ee5e9d57143e8a`。
- `ppara`：[2P54](https://www.rcsb.org/structure/2P54) / [735](https://www.rcsb.org/ligand/735)；原始 SHA-256 `35db7f1d001b3d5723e2817f576394ac8436044a6890708689d2ebaedc29fd7a`。
