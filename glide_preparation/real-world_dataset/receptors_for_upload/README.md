# 工作站上传用：16 个靶点

将本目录整体上传到 Linux 工作站。16 个 `<target>.pdb` 均在本目录下，
每个文件包含受体、下表指定的保留组分，以及一个用于定位口袋的晶体参考配体。
**这是整理后的 Protein Preparation 输入，不是已完成化学预处理的受体，也不是 Glide grid。**

原始 `real-world_dataset/` 和旁边的 `pdb_complexes/` 均未修改。
10 个靶点沿用数据包蛋白；6 个特殊靶点选用已下载的参考晶体链，见下表。
所有保留原子坐标不变；只重排原子序号并同步整理连接记录。
没有加氢、补侧链、建模缺失环、选择替代构象或能量最小化。

| 文件 | 参考 PDB | 蛋白来源/保留链 | 建 grid 时选择并排除的参考配体 | 其他处理或保留组分 |
|---|---|---|---|---|
| aa2ar.pdb | 3EML | 数据包 A | ZMA A:401 | 保留融合构建体；不加入 STE、SO4 |
| abl1.pdb | 2HZI | 数据包 A/B | JIN A:600 | 保留原裁剪片段 |
| aces.pdb | 1E66 | 数据包 A | HUX A:803 | 保留原受体 |
| adrb1.pdb | 2VT4 | 数据包 B | P32 B:400 | 使用 B 链位点；保留原构建体 |
| akt1.pdb | 4GV1 | 参考晶体 A | 0XZ A:501 | 恢复 TPO308；替代构象待预处理选择 |
| cdk2.pdb | 1H00 | 参考晶体 A | FAP A:1300 | 恢复被裁掉的已有坐标，保留 ACE 端基 |
| cp3a4.pdb | 3NXU | 参考晶体 A | RIT A:600 | 保留 HEM A:508，包括 Fe |
| dpp4.pdb | 2I78 | 数据包 A/B/C/D | KIQ B:901 | 保留原多链环境，使用 B 链位点 |
| drd3.pdb | 3PBL | 数据包 A | ETQ A:1200 | 使用正确口袋参考，未改设计配体 |
| dyr.pdb | 3NXO | 参考晶体 A | D2B A:187 | 保留 NDP A:193，即 NADPH |
| egfr.pdb | 2RGP | 数据包 A | HYZ A:1 | 保留原受体 |
| gria2.pdb | 3KGC | 数据包 B | ZK1 B:263 | 保留本地编号及已调整侧链；不使用 GLU 位点 |
| jak2.pdb | 3LPB | 参考晶体 A | NVB A:1133 | 替换误用的 JNK3；保留 PTR1007/1008 |
| kif11.pdb | 3CJO | 参考晶体 A | K30 A:1 | 保留 ADP A:601、MG A:603 |
| parp1.pdb | 3L3M | 数据包 A | A92 A:351 | 保留原受体 |
| ppara.pdb | 2P54 | 数据包 A/B | 735 A:469 | 保留共激活因子肽链 B |

## 上传后在 Maestro 中怎么做

1. 导入一个 PDB，打开 **Protein Preparation Workflow**。
   保留参考配体参与准备；先确认链、参考配体和上表保留组分齐全。
2. 开启 **Preprocess**，勾选 **Fill in missing side chains**。
   在 **More Options** 检查 **Assign bond orders / Using CCD database**、
   **Replace hydrogens**、**Create zero-order bonds to metals** 和 **Disulfide bonds**。
   不批量开启 **Fill in missing loops**，不做额外结构对齐或编号变换。
   端基封端按实际断点检查，不对所有裁剪端点一律自动封端。
3. 开启 **Optimize H-bond Assignments** 和受限最小化。
   pH 按实验条件统一设置；若暂时没有实验 pH，可将 7.4 作为明确记录的起始假设。
   两种生成方法共用同一个最终准备结果和同一个 grid。
4. 作业成功后检查 Diagnostics 和日志，保存新生成的 prepared `.maegz`/`.mae`。
   对多构象/多状态输出选定一个合理状态；不要直接拿导入的 PDB 生成 grid。
   HEM、NADPH、ADP/Mg 和磷酸化残基须再次检查键级、质子化、金属配位及电荷。
   CP3A4 保留了原晶体 Fe–RIT 配位记录；排除 RIT 后必须确认 HEM/Fe 仍在受体中，
   不要把整个金属连接组一起删除。
5. 使用准备成功的结构打开 **Receptor Grid Generation**，按表选中参考配体定位口袋，
   让该配体从受体 grid 中排除。HEM、NADPH、ADP/Mg、修饰残基、共激活肽继续保留。
   先用相应晶体配体做 Glide SP 回对接，检查姿势和相互作用，再运行正式评估。

上述选项名称参照 [Schrödinger 官方对接教程](https://learn.schrodinger.com/public/TwS/gb-docking-ls.pdf)
和 [Protein Preparation Workflow 说明](https://www.schrodinger.com/life-science/learn/white-papers/protein-preparation-workflow/)。
不同版本界面可能有差异；本地没有 Schrödinger，因此尚未验证工作站上的预处理或 grid 是否通过。

## 仍需注意的边界

- **普通缺失侧链交给 Schrödinger 自动补全。** 本目录没有人工补建侧链。
- **CDK2 的实验未解析残基仍然缺失**：13–14、36–43、152–161。
  恢复参考晶体只能补回已有实验坐标，不能恢复这些残基。
  13–14 靠近 ATP 结合位置，应优先检查口袋影响；152–161 属于活化段，也需结合构象检查。
  不要求无差别补齐全部缺失环；它们的处理和回对接结果需要记录。
- **DRD3 原设计配体问题未消失**：这里添加 ETQ 只修正 docking 口袋定位。
  原始 `crystal.mol2`、设计 SMILES 和生成结果均未替换。
- **这是一套无水起始输入**，延续原数据包的处理，不代表所有晶体水均无作用。
  若回对接或关键氢键提示需要保留桥联水，须从参考晶体单独选取、验证并统一更新两种方法的受体。
- 保留的数据包裁剪、融合蛋白、物种差异及工程突变仍存在，不能视作野生型完整受体。
- 导入和坐标核对通过不等于价态、力场或 docking 已通过。
  如果仍出现 Lewis structure 错误，应针对 prepared 输出中报告的残基/原子排查，
  不要继续通过随意删除口袋组分来尝试消除错误。

## 可追溯性与复现

`manifest.json` 记录每个靶点的来源、来源与输出 SHA256、保留组分、参考配体中心以及处理原因。
15 个靶点的参考配体重原子元素及坐标与原 `crystal.mol2` 完全对应；DRD3 是已知例外。
在远程本目录运行 `sha256sum -c SHA256SUMS` 可检查 16 个 PDB 及 manifest 的传输完整性。

从仓库根目录重新生成到一个不存在的新目录：

```bash
python glide_preparation/prepare_upload_receptors.py --output-dir /tmp/glide_receptors_new
```

脚本仅依赖 Python 标准库，生成 PDB、manifest 和校验和；本说明为人工维护文件。
脚本不会覆盖已有输出目录。
