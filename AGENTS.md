# AGENTS.md

## Project Goal

This project benchmarks two ligand-generation approaches across 16 protein targets:

- **REINVENT4**
  - LibINVENT
  - LinkINVENT
- **DELETE**
  - Corresponding fragment growing / fragment linking tasks

The goal is to compare the generation performance of the two models under consistent experimental settings across all 16 targets.

## Scope

For REINVENT4, only work on:

- **LibINVENT**
- **LinkINVENT**

Do **not** introduce Mol2Mol or other REINVENT4 generation modes unless explicitly requested.

Experiments should use the same target structures and, as far as possible, equivalent ligand fragments / constraints for both models so that comparisons are fair.

## Evaluation Plan

Evaluation lives in `evaluation/` with two independent modules:

```
evaluation/
├── docking/      # AutoDock Vina scoring
└── descriptor/   # molecular descriptors
```

### Data flow

The single, shared contract for both modules is **SMILES strings as input**.

```
model outputs (SMILES / SDF / MOL2)
        │
        ▼
  convert → canonical SMILES table  (target, method, smiles)
        │
        ├──────────────► evaluation/docking/
        └──────────────► evaluation/dcriptor/
```

The docking module is implemented in `evaluation/docking/`
(`prepare_receptors.py` + `run_docking.py`; see its README). Input contract:
one CSV with columns `target,method,smiles` (`method` ∈ `reinvent`/`delete`);
output: `molecule_id,method,target,smiles,vina_score` plus error/summary side
files. Grid boxes and Vina parameters are frozen per target in
`evaluation/docking/configs/*.json`.

DELETE outputs are not SMILES by default; they **must** be converted to SMILES
before entering the evaluation pipeline. The conversion step is a dedicated,
reusable adapter so both DELETE and REINVENT4 produce the same intermediate
format. Raw model outputs are preserved; converted tables are written as new
derived files.

### `evaluation/docking/`

- **Inputs:** generated-molecule SMILES + target 3D structure (receptor).
- **Pipeline:** SMILES → 3D conformer generation (RDKit ETKDG + MMFF/UFF
  minimization) → ligand PDBQT → Vina docking against the prepared receptor
  PDBQT → binding affinity (kcal/mol) per molecule.
- Receptor preparation and the docking grid box are configured per target, one
  config per target shared by both models so the comparison stays fair.
- Outputs: per-molecule scores plus per-target summary statistics (best, mean,
  median, success count).

### `evaluation/descriptor/`

- **Inputs:** the same SMILES table (no 3D structure needed).
- **Pipeline:** SMILES → RDKit descriptors (e.g. MW, LogP, TPSA, HBD, HBA,
  rotatable bonds, ring counts, QED, SA score) and optionally fingerprints.
- Outputs: one descriptor row per molecule, plus per-target summary statistics
  for DELETE vs REINVENT4 comparison.

### Fairness rules

- Both modules consume the identical SMILES table for a given target, so any
  difference in the reported numbers comes from the generation model, not from
  the evaluation path.
- Target structures, grid boxes, and descriptor settings are defined once in
  `configs/` (or `evaluation/` configs) and reused across models.
- Evaluation scripts must be deterministic (fixed random seed for conformer
  generation) and reproducible from raw outputs.

## 靶点输入核查与方法限制（2026-09-27）

以下记录用于解释生成和 docking 结果，不代表所有靶点都已完成 docking 验收。
本节只记录分析；不能据此自动替换设计配体、修改 PLANET 或启动实验。

### 原始数据与已修正输入

- 已确认 DURIAN 官方数据包与本地 32 个原始受体/配体文件逐字节一致。
  这里说的“本地文件缺失/删除”是指数据包相对参考晶体少了内容，并非用户自行删除。
- JAK2：原始 `receptor_out.pdb` 错用了 JNK3，且有拼接、截断。
  修正输入使用官方 3LPB 的作者链 A；其 NVB A/1133 与现有设计配体坐标匹配。
  保留坐标系、原始文件和设计 SMILES。
- DRD3：受体对应 3PBL，但原始设计 MOL2 的中心选错口袋。
  口袋定位使用 3PBL 的 ETQ A/1200，受体使用作者链 A。
  “口袋定位参考”与“设计配体”必须分开；修正中心并不解决原设计配体不匹配的问题，
  不得静默替换设计 SMILES。
- PLANET 修正输入及来源记录见 `configs/manual_planet_rl/target_inputs.json`。
  当前 JAK2/DRD3 docking 配置也记录了正确受体来源和口袋中心；
  配置正确不能代替受体 PDBQT 检查和真实 docking 验证。

### PLANET 看不到血红素和 NADPH：程序无变化，不等于真实结合无影响

- 当前实现使用蛋白三维口袋的标准残基特征及 Cα 坐标，加上生成分子的二维图；
  不使用生成分子的三维结合姿势。蛋白读取只取 `ATOM`，HEM/NDP 的 `HETATM`
  不进入模型；TPO 等修饰残基也不属于其标准残基特征字典。
- 在蛋白坐标、口袋中心、其余输入及运行条件不变时，仅补回血红素或 NADPH，
  当前 PLANET 的评分输入不变，生成任务不会因此改变。
  如果修复过程同时改变蛋白坐标或口袋，则不能沿用这个结论。
- 可以在保持方法固定的比较中把这点记录为 PLANET 的方法局限；
  不能据此声称生成分子的真实结合不受辅因子影响。
  PLANET 可能漏掉分子与辅因子的碰撞或有利接触，影响大小尚未定量测量。
- CP3A4 / 3NXU 缺少 HEM 血红素：参考链 A 的 HEM 与现有配体最近重原子距离
  约 2.23 Å，参考配体直接与血红素铁结合，因此不是远处无关成分。
- DYR / 3NXO 缺少 NDP（NADPH）：与现有配体最近重原子距离约 3.29 Å，
  位于结合区域附近，不能视为无关成分。

### CDK2 缺失分两类，处理方式不同

- **已有实验坐标、但被数据包裁掉的 58 个残基**：最近重原子距现有配体
  约 23.59 Å，距当前口袋中心约 28.49 Å。只补回这些现成坐标后，
  当前 PLANET 选中的全部口袋 ATOM 记录不变；这类裁剪对当前生成评分的直接影响小。
- **原始 1H00 晶体就未解析的 20 个残基**：作者链 A 的 13–14、36–43、152–161。
  其中 13–14 位于靠近 ATP 结合位置的环，152–161 属于活化段并包含 Thr160。
  这部分值得重点检查；不能称为完全无关，也不能凭缺失数量断言生成任务失效。
  对生成评分及分子排名的影响尚未量化。
- “未解析”表示没有实验坐标，不能直接从原 PDB 补回；需要建模，或选用结构状态、
  结合口袋和配体条件合适的其他实验结构。补建结构仍需验证，不能默认更准确。
- CDK2 文件标记为 `POCKET`，说明经过口袋裁剪。裁剪通常用于减少输入规模，
  但上游没有交代具体删除规则，不能把这一常见目的当成每个残基被删的确定原因。

### 其他裁剪、修饰和结构来源

- ABL1 / 2HZI：裁剪后保留 A/B 两链片段；缺少的已解析残基最近距配体
  约 20.53 Å。只补回现成坐标，当前 PLANET 口袋 ATOM 记录不变。
- AKT1 / 4GV1：未保留磷酸化残基 TPO308，最近距现有配体约 17.70 Å，
  位于当前口袋选择范围之外。它有生物学意义，但补回它不会自动改变固定的蛋白构象。
  ILE186、ARG406 实际仍保留，不能因参考替代构象不同误判为删除。
- GRIA2 / 3KGC：裁剪、重编号并调整过侧链，参考来源为大鼠。
  缺少的已解析残基最近距配体约 20.40 Å；补回现成坐标不改变当前 PLANET 口袋。
  但参考 B:72 ASN（本地 B:128）的侧链调整跨过 12 Å 选择边界，
  原参考与本地结构的口袋选择存在差异；不能把侧链调整视为毫无影响。
- 本次临时补回已解析残基/辅因子的检查覆盖 ABL1、CDK2、AKT1、GRIA2、CP3A4、DYR：
  六个靶点的全部选中 ATOM 记录均不变。这是输入检查，不是实测亲和力误差验证；
  也不覆盖未解析残基的建模。距离均来自已有参考结构，不能给未解析残基编造坐标距离。
- AA2AR / 3EML 含 T4 lysozyme 融合部分；ADRB1 / 2VT4 来自火鸡且含稳定化突变；
  ACES / 1E66 来自电鳐；PPARA / 2P54 保留共激活因子肽段。
  这些不等于输入错误，但解释实验结果时应说明物种和构建体限制。
- DPP4 / 2I78、EGFR / 2RGP、KIF11 / 3CJO、PARP1 / 3L3M
  未发现同类受体/配体身份错配。能解析 PLANET 口袋，不等于 docking 已通过验证。

### Docking 的处理优先级

1. 使用正确的 JAK2 受体和 DRD3 口袋；单独保留并说明 DRD3 设计配体的问题。
2. CP3A4 和 DYR 的 docking 受体建议分别保留/补回正确位置的血红素和 NADPH，
   核对链、坐标、化学状态、原子类型及最终 PDBQT，避免制备时再次丢失辅因子。
   普通 Vina 对含铁血红素相互作用的准确性仍需验证，不能认为加上 HEM 就完全解决。
3. CDK2 可先恢复已有实验坐标的部分；对靠近口袋的未解析部分评估建模或换结构，
   不要盲目补齐后直接认定可靠。若换结构，必须同步核对坐标系、口袋及设计依据。
4. 用已知晶体配体重新对接，检查是否能回到已知结合位置，再决定是否用于正式评估。
   其余远处裁剪的优先级较低，但不能统一认定剩余靶点完全没有问题。
5. 任何处理都保留原始数据，输出为可追溯的派生文件；REINVENT 和 DELETE
   必须共用同一套受体、辅因子处理、网格和评分参数。
   如改变生成时的蛋白输入，记录该变化并评估与旧结果的可比性。

### 证据入口

- 本地实现：`src/planet_oracle/chemistry.py`、`src/planet_oracle/backend.py`、
  `src/prepare_manual_planet_rl.py`、`configs/manual_planet_rl/README.md`。
- 上游 PLANET：[蛋白和分子输入实现](https://github.com/ComputArtCMCG/PLANET/blob/main/chemutils.py)。
- 参考晶体：[CDK2 1H00](https://www.rcsb.org/structure/1H00)、
  [CP3A4 3NXU](https://www.rcsb.org/structure/3NXU)、
  [DYR 3NXO](https://www.rcsb.org/structure/3NXO)、
  [AKT1 4GV1](https://www.rcsb.org/structure/4GV1)、
  [GRIA2 3KGC](https://www.rcsb.org/structure/3KGC)、
  [ABL1 2HZI](https://www.rcsb.org/structure/2HZI)。
- 验证原则：[AutoDock Vina FAQ](https://autodock-vina.readthedocs.io/en/latest/faq.html)。

## Development Rules

- Keep scripts, configurations, intermediate data, and results organized and reproducible.
- Do not make unnecessary changes outside the current task.
- Prefer simple and transparent implementations over unnecessary abstractions.
- Preserve raw input data; derived or processed files should be stored separately.

## Git

**Every meaningful change must be committed to Git.**

After completing a coherent unit of work:

1. Check the changes with `git status` / `git diff`.
2. Before committing, audit for files that may be generated, local-only, sensitive, or otherwise inappropriate for version control but are not covered by `.gitignore`.
3. For every possible candidate, report its path, why it may need to be ignored, and the proposed ignore pattern to the user. Ask the user whether the pattern should be added; do not silently add an ignore rule, stage the candidate, or delete it.
4. After the user decides how to handle any candidates, commit the changes.
5. Use a short, descriptive commit message.

The pre-commit audit should include:

```bash
git status --short --untracked-files=all
git status --short --ignored --untracked-files=all
git ls-files --others --exclude-standard
```

Use `git check-ignore -v --no-index -- <path>` to verify whether a specific path is covered by an ignore rule. Review both the root `.gitignore` and the `.gitignore` files belonging to the `Delete` and `REINVENT4` submodules. Do not use broad rules that could hide raw inputs, configurations, reproducible scripts, or intentionally versioned results without the user's decision.

Do not leave completed work uncommitted.
