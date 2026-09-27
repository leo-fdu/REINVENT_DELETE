# 手动设计任务的 PLANET RL 批量运行

本工作流把 `configs/manual_design/` 导出的手动设计输入（16 个靶点，30 个 designed 任务；
aces 两个模式已 skipped）批量运行为只用 PLANET 奖励的 RL 生成实验，覆盖 LibINVENT 与 LinkINVENT。
不修改 REINVENT4、DELETE、PLANET 上游：专用入口 `src/manual_planet_rl.py` 只在自己的进程内
替换这两个模式的 Transformer 学习类和 CSV 报告器，其余生成模式不受影响。

评分链路复用 `configs/planet_oracle/README.md` 的常驻服务：REINVENT `ExternalProcess`
经 `planet_oracle.client` 请求同机 `planet_oracle.server`，模型和口袋每个靶点只初始化一次。

## 固定实验协议

与 `configs/manual_planet_rl/experiment.json` 一致，正式批量不得逐项改动：

- 每个任务一次运行：batch_size 100 × max_steps 100 = **10000 次采样**，`termination = "null"` 不早停。
- 唯一评分组件为 PLANET 奖励（weight 1，arithmetic_mean）；原始亲和力经 REINVENT 内置
  sigmoid（low=4，high=16，k=0.25）转换为 0–1 奖励，客户端不做二次归一化。
- 重复完整分子（同批次内或跨批次）奖励 ×0.5（penalty_multiplier），并照常参与 agent 更新；
  无效分子奖励 0，不参与更新，也不进入重复记忆。批内重复由专用入口保证同样按 ×0.5 计奖。
- DAP：sigma=128，rate=0.0001。seed=42（同时写入子进程 `PYTHONHASHSEED`）。
  multinomial 采样，temperature=1.0；randomize_smiles 关，isomeric_smiles 开。
- 受体与口袋定位：默认使用原 `receptor_out.pdb` 和设计 `crystal.mol2` 的重原子质心；
  JAK2/DRD3 按 `target_inputs.json` 使用明确的受体链与晶体参考残基，见下节。
  prepare 固化来源、哈希、链和 center；不移动坐标，不改动设计 SMILES。
- 保留全部生成记录：每任务和合并的 `all_generated.csv` 收录每一次采样，含无效与重复。

## 可追溯的结构输入修正

上游 DURIAN 与本地原始文件已由用户确认逐字节一致；此处保留全部原始文件。
口袋能解析且非空不代表定位正确，DRD3 的旧中心正是这一类问题。
本修正仅改变 PLANET 的蛋白输入/口袋位置，不重新导出任何设计片段。

| 靶点 | 运行受体 | 独立口袋定位参考 | center（Å） |
| --- | --- | --- | --- |
| JAK2 | 官方 3LPB 作者链 A | 同结构 NVB / A / 1133（33 重原子） | 114.671000, 66.113515, 10.554818 |
| DRD3 | 原 DURIAN 3PBL 受体作者链 A | 官方 3PBL ETQ / A / 1200（23 重原子） | 0.085304, -14.828348, 10.432043 |

- JAK2 原受体错用了 JNK3，且有拼接/截断；现有设计配体坐标与 3LPB/A NVB 完全一致。
  使用正确蛋白，保留其坐标系、设计 MOL2 和 15 靶点的全部 30 份 SMILES。
- DRD3 原设计 MOL2 质心与两个 ETQ 晶体口袋没有残基交集。固定选择作者链 A（两个
  等价晶体单体中按链名的第一个），以该链的 ETQ 1200 定位正构口袋，排除 B 链混入。
  ETQ 只用来确定位置；它不替代 `configs/manual_design/` 中的设计配体或片段。
- 官方来源及固定 SHA-256 见 `target_inputs.json` 和
  `real-world_dataset/planet_reference_sources/README.md`。两个完整官方 PDB 是固定参考输入。
- `src/planet_target_inputs.py` 从源 PDB 提取选定链：保留 ATOM、MODRES 标明的聚合物
  HETATM 和所有选中原子坐标；删除水/游离配体/其他链。每个残基的替代构象按平均占有率
  最大、并列时按标签字典序选一个。不加氢、不修复缺失残基、不转换修饰残基。
  PLANET 上游仍只解析 ATOM；3LPB 的两处 PTR 保留为 HETATM，均不在 NVB 的 12 Å 口袋内。
- 仅 JAK2/DRD3 的派生 `receptor.pdb` 写入新运行目录，不在原始数据目录生成或覆盖受体。
  其余 designed 靶点保持原受体和原 MOL2 质心；ACES 仍跳过。
- `run.json` 的 `source_mol2` 始终表示设计来源；`pocket_reference` 独立记录定位参考。
  `input_correction` 记录来源/链/原因。源 PDB、修正配置加入 `input_hashes`；
  派生受体加入 `prepared_hashes`，runner 对两类文件都拒绝篡改。
- 这些选择限定于本 PLANET RL 工作流；Vina/docking 的旧配置未随本任务修改。
  后续 DELETE 或 docking 比较应显式复用相同正确受体/口袋，不能直接混用旧 JAK2/DRD3 设置。

## 前置条件（服务器）

- 两个 Python ≥ 3.11 环境：REINVENT 环境（REINVENT4 依赖及两个
  `*_transformer_pubchem.prior`）；PLANET 环境（Torch、RDKit、NumPy、SciPy、Pandas，
  以及实际 `PLANET/` 源码和 `PLANET/PLANET.param`——该目录被 Git 忽略，需单独同步）。
  prepare/run 脚本本身只用标准库，示例直接用 REINVENT 环境的绝对路径调用。
- `configs/manual_design/manifest.json` 和各 `.smi` 是 `src/export_manual_design.py` 的导出结果，
  不得改动；prepare 会核对 crystal.mol2 的 SHA-256、connectivity_match 和 SMILES 文件内容。
- run 目录必须是仓库外的**新**目录；仓库内或已存在的目录会被拒绝。

## 两步命令

```bash
export PROJECT=/srv/REINVENT_DELETE
export PLANET_PY=/srv/envs/planet/bin/python
export REINVENT_PY=/srv/envs/reinvent4/bin/python
export RUN=/srv/runs/manual-planet-rl-1

"$REINVENT_PY" "$PROJECT/src/prepare_manual_planet_rl.py" --run-dir "$RUN" \
  --reinvent-python "$REINVENT_PY" --planet-python "$PLANET_PY" \
  --device cuda:0 --port 8765
"$REINVENT_PY" "$PROJECT/src/run_manual_planet_rl.py" --run-dir "$RUN"
```

- `--inputs` / `--experiment` 默认指向 `configs/manual_design/manifest.json` 和
  `configs/manual_planet_rl/experiment.json`；`--device` 默认 cuda:0，`--port` 默认 8765。
- `--target-inputs` 可指定修正 JSON；默认使用项目内
  `configs/manual_planet_rl/target_inputs.json`（不存在时兼容原工作流）。配置中的路径相对
  项目根目录，必须在项目内，来源哈希、参考残基和受体链必须匹配，否则创建 run 目录前失败。
- prepare 一次性生成全部 designed 任务的配置并冻结所有输入/产物哈希到 `run.json`。
- run 默认执行全部靶点；`--targets aa2ar cdk2` 可选子集（同一目录仍只允许运行一次）。
- run 目录单次使用：启动即写 `.started` 标记并永久占用；重复运行直接报
  `FileExistsError`。失败或中断后不要复用目录，重新 prepare 一个新目录再跑。

## 运行方式与端口

- 按靶点串行：启动该靶点的 `planet_oracle.server`（所有靶点共用同一端口）→
  轮询 `/health` 就绪并核对 version/ready/target_id/oracle_id/配置指纹，把 oracle_id
  写入 client.json → 依次运行两个模式任务 → 停止服务 → 下一靶点。
- 端口探测使用 SO_REUSEADDR：上一靶点服务停止后遗留的 TIME_WAIT（数十秒）不会阻塞
  下一靶点绑定；端口上存在活动监听时探测仍失败，绝不连接或顶替已运行的服务。
- 每个子进程独立进程组；正常结束或 Ctrl-C/SIGTERM 后整组先 SIGTERM、15 秒后 SIGKILL，
  不留下孤儿评分服务。退出码：成功 0，失败 1，中断 130。
- `status.json` 实时记录整体状态（running/complete/failed/interrupted）和每个
  `<target>/<mode>` 任务状态；中断时未完成的任务记为 interrupted，已完成的产物保留。

## 校验与产物

run 启动时按 `run.json` 复核全部源输入与 prepared 文件哈希，任何改动拒绝启动。
每个任务结束后在导出前逐项校验，任一不符即整个 run 失败（已完成任务的产物仍保留）：

- 每个 step 的行数恰为 batch_size，总预算 10000 行不缺；CSV 列来自专用入口。
- 逐行核对 `training_reward = sigmoid_reward × (重复 ? 0.5 : 1)`，重复按全任务记忆重算。
- 已评分行的 target_id/oracle_id 必须与本靶点服务一致；无效行必须奖励 0、无亲和力、
  status 为 not_scored。

目录结构：

```
$RUN/
  run.json  experiment.json  status.json  .started  all_generated.csv
  <target>/server.json  client.json  oracle-manifest.json  server.log
    receptor.pdb  # 仅修正靶点，由 prepare 从冻结源结构提取
    <mode>/rl.toml  input.smi  summary_1.csv  all_generated.csv
             agent.chkpt  run.log  process.log  resolved.json  tensorboard_*/
```

合并的 `all_generated.csv` 只含本次选中的靶点，17 列：
`target, mode, seed, step, sample_index, input_smiles, generated_smiles, smiles,
smiles_state, is_repeat, planet_affinity, planet_status, planet_error,
sigmoid_reward, training_reward, planet_target_id, oracle_id`。
`smiles_state`：0=无效，1=有效，2=批内重复；`is_repeat` 由 runner 按全任务（跨批次）重算；
`planet_affinity` 是 PLANET 真实预测（未评分时为空）；`planet_status` 取
ok/invalid_input/not_scored；`training_reward` 是进入 agent 更新的最终奖励。
`agent.chkpt` 为任务结束时的最终 agent，runner 校验其存在且非空。

## 验证

本机自动化测试（fixture 子进程，不加载真实模型）：

```bash
cd "$PROJECT"
PYTHONPATH="$PROJECT/src:$PROJECT/REINVENT4" "$REINVENT_PY" -m pytest tests/manual_planet_rl -q
```

覆盖：预算/不早停/center 计算、批内与跨批次重复 ×0.5 且参与更新、CSV 列稳定、
输入篡改拒绝、预算不足拒绝导出、多靶点串行复用同一端口（TIME_WAIT 下仍能绑定）、
子进程清理与目录占用。

真实模型验收（不启动正式批量）：复制 experiment.json 到项目外，把 steps 改为 2，
保留 batch_size=100、其余奖励/采样参数，prepare 时传入这个副本并选择新 run 目录。
先检查全部 15 个将运行靶点，再只跑 JAK2/DRD3 的四个任务：

```bash
PYTHONPATH="$PROJECT/src:$PROJECT/REINVENT4" PYTHONHASHSEED=42 \
  "$PLANET_PY" "$PROJECT/src/validate_manual_planet_inputs.py" \
  --run-dir "$RUN" --report-dir /srv/runs/input-validation-1
"$REINVENT_PY" "$PROJECT/src/run_manual_planet_rl.py" \
  --run-dir "$RUN" --targets jak2 drd3
```

验证入口会先核对冻结文件，再逐靶点加载真实 PLANET、检查口袋、对 CCO 做有限值推理，
保存模型/口袋/源码指纹和报告。不会运行 RL，不会写 `.started`；任一失败非零退出。
四个 RL 冒烟任务应总计 800 条记录、各有非空最终 checkpoint。
仅小预算目录被实际运行；正式 30 任务配置可另行 prepare，不能以这次验收代替正式启动授权。

### 2026-09-27 服务器验收结果

Ubuntu 22.04 / RTX 4090 / Python 3.11.16 / Torch 2.12.0+cu126 / RDKit 2026.3.6：

- `tests/manual_planet_rl` 与 `tests/planet_oracle` 合计 **86 passed**。
- 全部 15 个运行靶点通过真实 CUDA PLANET 口袋解析及 CCO 有限值评分；
  JAK2 49、DRD3 63、GRIA2 57 个口袋残基。
- JAK2/DRD3 × LibINVENT/LinkINVENT 四个任务各跑 2 步 × 100 条，共 **800 条**。
  797 条有有限 PLANET 预测，3 条无效分子正确保留为未评分、奖励 0；
  逐行预算、重复奖励、靶点/oracle 身份检查通过，四份最终 checkpoint 均非空。
- 原始 32 份受体/配体、30 份设计 SMILES 及设计 manifest 共 63 个文件的 SHA-256
  与修正前一致；两个运行目录的派生受体逐字节符合最终提取代码。
  其余 13 个运行靶点的原受体路径和 MOL2 中心未改。
- 重新 prepare 的正式目录含 30 个任务（ACES 两项跳过），没有 `.started`，未运行正式批量。

服务器证据保存在 `/home/ubuntu/experiments/reinvent-20260927/input-fix-20260927/`：
`acceptance-report.json`、`pocket-validation-final/`、`smoke-rl/`、`prepared-full/`。
完整运行产物留在仓库外；本节记录的是小预算运行链路验收，不能作为正式实验效果结论。
