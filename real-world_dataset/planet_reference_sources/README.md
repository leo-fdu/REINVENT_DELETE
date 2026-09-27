# PLANET 口袋定位参考结构

这些文件是固定、版本控制的参考输入，不是生成分子或运行结果。
保留 DURIAN 原始 `receptor_out.pdb`、`crystal.mol2` 和现有设计 SMILES。
结构修正只在 `configs/manual_planet_rl/target_inputs.json` 指定；派生受体在 prepare 的新
运行目录生成，原子坐标保持不变。脚本不依赖网络，也不需要 Vina、PDBQT 或配体 3D 生成。

| 文件 | 官方来源 | SHA-256 |
| --- | --- | --- |
| 3LPB.pdb | https://files.rcsb.org/download/3LPB.pdb | 89cfb782bbed5c6113430b2c18d9a880d918bf48f2fcc6a4dec7cf0e4d1b945f |
| 3PBL.pdb | https://files.rcsb.org/download/3PBL.pdb | ad01099584fe8ac8b1d8eb97cd8f10a810f5607c0902921e7749abd07a213586 |

下载于 2026-09-27；未来官方修订可能改变字节内容，运行以本仓库固定文件和哈希为准。

## 选择依据

- [3LPB](https://www.rcsb.org/structure/3LPB)：人 JAK2 激酶域与 NVB 共晶。
  作者链 A 的 NVB 1133 共 33 个重原子，其坐标集合与现有 JAK2 `crystal.mol2` 完全一致。
  提取正确受体 A 链；不是对错误 JNK3 输入做平移或结构拟合。
- [3PBL](https://www.rcsb.org/structure/3PBL)：人 DRD3/T4 lysozyme 嵌合体与 ETQ 共晶，
  有 A/B 两个作者链单体。固定选择 A 的 ETQ 1200（23 个重原子）定位正构口袋；
  运行蛋白由现有正确的 DURIAN 3PBL 受体提取 A 链，保留原有预处理坐标。
  B/ETQ 1200 的中心为 (4.714435, 13.399043, -9.408696) Å，属于另一个晶体单体。
  原设计 MOL2 中心 (8.405391, 21.356087, 23.562652) Å 不能用作这两个 ETQ 口袋的中心。

非空口袋只能验证解析能力，不能证明靶点身份或口袋定位正确。来源、作者链、参考残基和
哈希是定位依据；实际 12 Å 口袋残基及张量指纹由 PLANET 验证/服务 manifest 记录。
