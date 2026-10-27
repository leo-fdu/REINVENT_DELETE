# Glide 靶点的原始 PDB 复合物

本目录保存 16 个靶点对应的 RCSB 原始 PDB 条目文件。文件按
`<target>_<PDB_ID>.pdb` 命名，下载字节不作修改。它们包含条目中所有已解析的链、
配体、水及其他组分；不是只含蛋白的文件，也不是已经完成预处理的 Glide 受体。
这是原始晶体条目坐标，不是另外生成的 biological assembly；“完整条目”不代表
完整蛋白序列都有实验坐标。未覆盖 DURIAN 原输入，也未修改生成任务。

`manifest.json` 保存每个文件的官方 URL、下载时间、HTTP Last-Modified、大小、
SHA-256，以及链、非水 HETATM 残基、水、替代构象和 REMARK 465 缺失残基记录。
这些是文本记录检查，不是化学结构或 docking 验收。配体/辅因子的作者链和编号
使用 PDB 文件中的标识。`SHA256SUMS` 用于上传后的完整性检查。

## 下载清单

| 靶点 | PDB ID | 后续重点 |
|---|---|---|
| AA2AR | 3EML | 识别 T4 lysozyme 融合部分，审查删减导致的断点 |
| ABL1 | 2HZI | 选配体对应结合位点；检查链、替代构象和与裁剪输入的差异 |
| ACES | 1E66 | 电鳐来源；检查替代构象；现有六批 REINVENT 导出不含 ACES |
| ADRB1 | 2VT4 | 火鸡来源、稳定化突变；选择对应晶体单体 |
| AKT1 | 4GV1 | 核查 TPO308 及其化学处理，不能自动改成普通 THR |
| CDK2 | 1H00 | 已解析部分可恢复；A:13–14、36–43、152–161 未解析，仍须评估 |
| CP3A4 | 3NXU | 保留目标单体 HEM 和 Fe，核查配位；RIT 为共晶抑制剂 |
| DPP4 | 2I78 | 检查目标结合位点及需要保留的链，不能机械删除所有其他链 |
| DRD3 | 3PBL | 作者链 A，ETQ A/1200 定位与重对接；原设计配体问题仍存在 |
| DYR | 3NXO | 保留 NDP/NADPH；不要在排除共晶抑制剂时一起排除 |
| EGFR | 2RGP | 常规缺失原子、质子化、结合位点和水检查 |
| GRIA2 | 3KGC | 大鼠来源；与本地编号/侧链不完全一致，先做结构映射 |
| JAK2 | 3LPB | 作者链 A，NVB A/1133；检查 PTR1007/1008 修饰；不再使用错误的原 JNK3 受体 |
| KIF11 | 3CJO | 含 K30 抑制剂、ADP 和 Mg；区分抑制剂位点与核苷酸位点，明确 ADP/Mg 保留方案 |
| PARP1 | 3L3M | 常规准备；可作为流程试跑候选 |
| PPARA | 2P54 | 核对并记录共激活肽处理；不能直接删除所有非主体链 |

靶点问题和选链依据来自仓库根目录 `AGENTS.md` 的结构核查记录。
除明确列出的选择外，本目录并未替各靶点完成最终选链或水处理决策。

本次文件级检查已确认：3NXU 有 HEM A/508 与 B/508；3NXO 有 NDP A/193；
4GV1 有 TPO A/308；3LPB 有 NVB A/1133 和 PTR A/1007、A/1008；
3PBL 有 ETQ A/1200；3CJO 有 K30、ADP 和 Mg。1H00 的 REMARK 465
仍列出上述 20 个未解析残基，重新下载并未为它们提供坐标。

## 本地上传与远程处理

本地无需安装分子建模软件。上传整个目录，保留原始文件作为来源备份。
在远程目录运行：

```bash
sha256sum -c SHA256SUMS
```

后续派生文件存到其他目录，例如 `../prepared_targets/<target>/`。
先确认远程 Schrödinger 版本及 Protein Preparation、Glide 的可用性；使用
LigPrep/Epik 或 Prime 时另确认对应功能和许可。

1. **选结构组分**：确认作者链、共晶抑制剂、辅因子、金属、相关肽段和水。
   处理替代构象，检查断链、缺失侧链、缺失主链和非标准残基。不要统一只留 A 链，
   不要统一删除所有 HETATM。若从本目录整套重新准备，记录与 DURIAN 生成输入的差异。
2. **Protein Preparation**：分配键级、检查金属/二硫键、加氢、补必要缺失原子，
   指定 pH、检查口袋质子化和 His 状态、优化氢键网络、受限最小化。
   pH 优先按实验条件；条件未知时预先声明统一基线。缺失环建模需另行验证，
   CDK2 不能依赖自动补环后直接放行。对 HEM/Fe、NADPH、TPO 单独检查化学状态。
3. **保存预处理结果**：保存 prepared complex 的 `.mae/.maegz`、参数、日志和
   结构选择记录；检查关键组分是否仍存在及口袋几何是否异常变化。
4. **Receptor Grid Generation**：用正确共晶配体定义结合位点，建网格时排除
   待对接抑制剂，保留受体环境中的 HEM/NADPH 等。设置并保存 inner/outer box、
   vdW 参数及约束；验证盒子适合预定生成分子的大小范围。保存 grid `.zip` 和输入文件。
   Vina PDBQT/网格 JSON 不能作为 Glide grid。
5. **重对接验收**：另制备共晶抑制剂做柔性 SP 重对接，不以刚性原位评分代替。
   保存未经优化的晶体 pose 作参考；在同一受体坐标系计算对称性校正的重原子 RMSD，
   不独立拟合配体。建议预先设定 Top-1 RMSD ≤2 Å 为初步目标，同时检查 Top-5、
   关键相互作用和碰撞。该标准是拟议验收线，不是本次已测结果。
   CP3A4 还要核查铁配位；重对接通过不等于亲和力排名已经可靠。
6. **小批量试跑，再批量 SP**：检查分子制备失败、超框、无 pose 等问题；冻结每个
   靶点的受体、grid 和评分参数，REINVENT/DELETE/不同批次共用。记录失败率，
   不静默丢弃失败项。正式配体沿用项目共同 SMILES 输入及统一制备策略，保留分子 ID。

当前状态：仅下载与文件级检查。没有选链输出、蛋白化学预处理、grid、重对接或 SP 结果。

## 参考

- 官方原始文件地址规则：`https://files.rcsb.org/download/<PDB_ID>.pdb`
- [Protein Preparation Workflow](https://www.schrodinger.com/life-science/learn/white-papers/protein-preparation-workflow/)
- [Glide inner/outer box](https://my.schrodinger.com/support/article/701)
- [Glide 对金属体系的处理边界](https://my.schrodinger.com/support/article/723)
- 用户提供的《Glide分子对接.pptx》第 5–7 页：复合物、Protein Preparation、grid；
  第 10–11 页的服务器 VSW 示例需要预先准备的配体和 grid，不执行蛋白修复。
