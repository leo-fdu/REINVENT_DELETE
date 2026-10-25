# ShallowBench 配体任务设计

本地页面：<http://127.0.0.1:8767>。复用原来的二维切键、保留组分选择、草稿和不可覆盖快照功能。
本轮仅准备任务设计，不执行 RL、PLANET 评分、DELETE 或 docking。

## 当前范围和已确认选择

- 七个 PDB，八个配体候选：4QYG/3DW、4U1O/FWF、4U1O/KAI、2ROY/P28、4N4F/已解析肽、3VHW/VHW、4O5K/0KX、4PWJ/30Z。
- 用户选择：4U1O 的两个候选都先设计，之后再选；4N4F 先展示，暂缓任务。
- 因而可编辑的是六个结构中的七个小分子候选；4N4F 可以查看、填写备注和标记跳过。
- 多副本结构暂取 A 链配体用于二维设计，不表示已证明该副本对应 ShallowBench 条目，也不表示受体必须只取 A 链。
- 所有 `pocket_mapping_confirmed` 均为 false。条目中心与当前候选晶体配体的质心相距约 13–73 Å，不能直接混用；原因（坐标变换、条目来源或口袋差异）尚未查明。

## 使用

在仓库根目录运行：

```bash
/opt/homebrew/Caskroom/miniforge/base/envs/reinvent4/bin/python \
  real-world_dataset/task_design_manual/app.py \
  --dataset shallowbench-7targets/task_design/dataset.json \
  --records shallowbench-7targets/task_design/records --port 8767
```

其他电脑使用安装了原界面所需依赖的 Python 即可，不必使用上述本机绝对路径。
浏览器访问同机地址。LibINVENT 选一条非环单键并保留一个组分；LinkINVENT 选两条并保留两个端组分。
点击“记录当前设计”会保存全部条目的当前状态，不要求全部完成；不适合的任务可标记跳过。
记录保存到 `records/<时间戳>/`，不会覆盖旧记录，也不会写入原来的 16 靶点记录目录。
浏览器草稿随同一浏览器/端口保留；请使用“记录”得到可追溯文件。
按用户选择，新的 `records/` **不添加忽略规则，设计快照也纳入版本管理**。
页面的“记录”负责写文件，不自动执行 Git 提交；新增快照在后续提交时一并保存。

页面中的数字是提取文件中稳定的 1-based 原子 ID，**不是 PDB serial**；二者关系在各候选的 `atoms.json`。
原页面不带 `--dataset` 启动时仍使用原来的 16 靶点。

## 文件与化学依据

- `ccd/*.cif`：从 RCSB 官方 Chemical Component Dictionary 下载的原始化学字典，URL 和 SHA-256 在 `dataset.json`。
- `ligands/<候选>/observed.pdb`：所选残基的原始重原子记录，晶体坐标不移动。
- `ligands/<候选>/atoms.json`：页面原子 ID、PDB serial、残基、原子名称、替代构象、原始坐标映射。
- `ligands/<候选>/ligand.sdf`：保留晶体重原子坐标，附 `authoritative_mapped_smiles` 图属性。
- `ligands/<候选>/design_graph.sdf`：CCD 声明的化学图及二维坐标，适合检查分子图，**不能当成晶体三维姿势**。
- `dataset.json`：提取来源、候选标签、警告、源哈希、口袋映射状态和编号化的权威 SMILES。
- `compatibility_audit.json`：枚举界面允许的切法，实际调用两个 prior 的词表/长度检查及标准化、拼接校验；不是用户的任务设计。

按 CCD 原子名建立键，验证重原子完整性、元素、连通性、价态及声明立体化学。
立体中心用 RDKit 新 CIP 标记器对照 CCD 的 R/S；VHW 上旧 CIP 标签算法会产生不同标签，不能据此翻转结构。
VHW 的三个双键立体声明均保留。0KX 的 CCD 没有声明磷中心手性，不从略有偏差的晶体几何凭空指定它。
因此普通软件重读三维 SDF 时可能从坐标推断出额外手性；本页面使用文件内的权威编号 SMILES，
而非再次从三维坐标猜测。对外交换化学图使用 `design_graph.sdf`；晶体坐标单独溯源。

4N4F 仅有 C/7–12 的 41 个已解析重原子；两端 `*` 是未解析肽链延续的显示符号，
没有实验原子坐标，其展示位置在映射中标为 `virtual_boundary`。未补建缺失肽段、未化学封端；不可用作正式生成输入。

## 复现和验证

已有冻结输入不需要重新提取。需重建时先审阅保存的设计，避免来源哈希变化：

```bash
python src/prepare_shallowbench_design.py --rebuild
python src/audit_shallowbench_design.py
python -m pytest tests/test_shallowbench_design.py real-world_dataset/task_design_manual/test_design.py -q
```

脚本读取归档 CCD，不依赖在线查询；提取不改变原始 PDB 或任何用户记录。
本次测试覆盖坐标/元素/原子映射、CCD 立体声明、权威图往返、肽片段阻止提交、真实 prior 校验、
不可覆盖记录、篡改拒绝、通用配体导出，以及原 16 靶点回归。

## 后续导出及边界

设计完成后，把未设计候选明确设为跳过，然后使用：

```bash
python src/export_manual_design.py \
  --designs shallowbench-7targets/task_design/records/<时间戳>/designs.json \
  --output configs/shallowbench_design
```

新快照/导出使用 `source_ligand`（SDF）而不是冒称 `source_mol2`，并携带完整 provenance。
**现有 PLANET 批量 prepare 仍限定旧 MOL2 路径，尚未接入这个新 manifest。**
正式生成前仍需解决口袋映射、受体输入、最终候选与化学状态、实验配置及真实模型小预算验收。
不能把本页面的设计/兼容性通过说成“七靶点可以直接运行”。
