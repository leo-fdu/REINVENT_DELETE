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
- 口袋：`real-world_dataset/<target>/receptor_out.pdb`；center 为原 `crystal.mol2`
  重原子的几何质心，prepare 阶段计算并固化进 server.json，不改动任何坐标。
- 保留全部生成记录：每任务和合并的 `all_generated.csv` 收录每一次采样，含无效与重复。

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

真实模型小预算冒烟可在本机 CPU 进行：复制 experiment.json 改小 steps/batch_size，
prepare 时 `--device cpu`，run 时 `--targets` 选两个靶点。正式 16 靶点批量前，
先按 `configs/planet_oracle/README.md` 在服务器通过 smoke 与真实验收。
