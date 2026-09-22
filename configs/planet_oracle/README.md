# PLANET 接入 REINVENT4 RL

本接口只接入 LibINVENT / LinkINVENT，不修改 REINVENT4、DELETE、PLANET。
生成器先拼接完整分子，内置 `ExternalProcess` 把 SMILES 经 stdin 交给轻量客户端，
客户端请求同机常驻服务并返回原始亲和力，REINVENT 内置 sigmoid 将其转换为奖励。
模型、口袋和蛋白特征只初始化一次。

本机交付不建立 PLANET 环境，不运行真实 PLANET 或 RL；下面的真实验收在 Linux 执行。
示例片段和奖励参数只用于验证接线，未针对 16 个靶点校准，不是正式靶点任务设计。

## 环境和路径

服务器应分别准备能运行本地版本 PLANET 与 REINVENT4 的两个 Python 环境（Python >= 3.11）。
客户端只有标准库依赖；PLANET 服务使用上游依赖，包括 Torch、RDKit、NumPy、SciPy、Pandas。
不要求在 REINVENT 环境中安装 PLANET。不要用环境自动激活代替以下绝对 Python 路径。

`PLANET/` 被根仓库 Git 忽略，服务器必须单独具备实际源码和 `PLANET.param`；
上传也必须包含实际 `REINVENT4/` 目录及两个 Transformer prior，不能只上传根仓库子模块指针。
可用 `git -C REINVENT4 rev-parse HEAD` 记录上游版本；本地改动需要一起同步。
服务 manifest 记录 PLANET 顶层 Python 源码、适配器源码、权重和结构文件 SHA-256、
解析后的配置、依赖版本、设备、口袋残基、实际口袋张量指纹。

先填写下列变量；`RUN` 使用仓库外的新目录，避免覆盖原始输入或已有结果：

```bash
export PROJECT=/srv/REINVENT_DELETE
export PLANET_PY=/srv/envs/planet/bin/python
export REINVENT_PY=/srv/envs/reinvent4/bin/python
export RUN=/srv/runs/planet-demo
export PYTHONPATH="$PROJECT/src:$PROJECT/REINVENT4${PYTHONPATH:+:$PYTHONPATH}"
export PYTHONHASHSEED=42
mkdir -p "$(dirname "$RUN")"
mkdir "$RUN"
mkdir "$RUN/libinvent" "$RUN/linkinvent"
cp "$PROJECT/configs/planet_oracle/server.example.json" "$RUN/server.json"
cp "$PROJECT/configs/planet_oracle/client.example.json" "$RUN/client.json"
cp "$PROJECT/configs/planet_oracle/libinvent_rl.toml" "$RUN/libinvent.toml"
cp "$PROJECT/configs/planet_oracle/linkinvent_rl.toml" "$RUN/linkinvent.toml"
```

编辑运行目录内的四份配置，逐项确认：

| 文件 | 必须确认的路径或字段 |
| --- | --- |
| server.json | target_id、planet_root、checkpoint、protein_pdb、ligand_sdf 或 center、device、port |
| client.json | URL/端口、相同 target_id；启动后填写 oracle_id |
| 两份 TOML | device、prior_file、agent_file、smiles_file、params.executable、params.args、transform.low/high/k |
| 两份 TOML | tb_logdir、json_out_config、summary_csv_prefix、chkpt_file 均指向各自的新运行目录 |

`params.executable` 必须是 **REINVENT 环境 Python 的绝对路径**；
`params.args` 为 `-m planet_oracle.client --config /绝对路径/client.json`，路径有空格时在参数字符串内加引号。
TOML 不展开 shell 环境变量，因此不能直接填 `$PROJECT` 或 `$RUN`。
两种模式分别使用 `libinvent_transformer_pubchem.prior` 和 `linkinvent_transformer_pubchem.prior`，
每个模式 prior 和初始 agent 指向同一文件。
JSON 配置的相对文件路径相对于该 JSON 文件；REINVENT TOML 路径按其原有规则解析，示例使用绝对路径。

蛋白必须预先准备好。本接口不修复蛋白、不选择质子化态；服务遇到不完整口袋直接拒绝启动。
使用显式中心时删除 `ligand_sdf`，加入 `"center": [x, y, z]`，必须恰好三个有限数值。
省略 checkpoint 时默认使用 planet_root 下的 `PLANET.param`。
device 默认 cpu；示例明确写 cuda:0，CUDA 不可用时失败，不回退 CPU。
batch_size 默认 32，port 默认 8765，seed 默认 42。
上游使用集合遍历残基/替代构象，固定 `PYTHONHASHSEED` 有助复现；每次启动仍应重新核对 oracle_id。

## 服务启动和独立评分

```bash
"$PLANET_PY" -m planet_oracle.server --config "$RUN/server.json" \
  --manifest "$RUN/oracle-manifest.json" >"$RUN/server.stdout.log" 2>"$RUN/server.log" &
PLANET_PID=$!
printf '%s\n' "$PLANET_PID" > "$RUN/server.pid"
curl --fail --silent --show-error http://127.0.0.1:8765/health
```

模型加载和蛋白特征预计算结束后才监听，首次健康检查可能需稍后重试；查看 `server.log` 判断进度。
服务只监听 127.0.0.1，每次处理一个请求。`--manifest` 拒绝覆盖已有文件；重启用新的文件名。
日志含解析后的完整配置。以下命令确认 target_id 后把健康接口返回的 oracle_id 写入客户端配置：

```bash
"$REINVENT_PY" - "$RUN/client.json" <<'PY'
import json, pathlib, sys, urllib.request
p = pathlib.Path(sys.argv[1])
c = json.loads(p.read_text())
with urllib.request.urlopen(c['url'].rstrip('/') + '/health', timeout=120) as r:
    h = json.load(r)
assert h['version'] == 1 and h['ready'] is True
assert h['target_id'] == c['target_id']
print(json.dumps(h, indent=2))  # 核对设备、口袋和指纹
c['oracle_id'] = h['oracle_id']
p.write_text(json.dumps(c, indent=2) + '\n')
PY
printf 'CCO\ninvalid\nCCO\n' | "$REINVENT_PY" -m planet_oracle.client \
  --config "$RUN/client.json" > "$RUN/independent-score.json"
```

同一次 RL 运行不能更换服务评分设置或更新客户端指纹，否则 REINVENT 缓存可能混入旧分数；
更换口袋、权重、源码或环境后应建立新运行目录并重启 REINVENT。

客户端 stdout 只有 `{"version":1,"payload":{...}}` JSON；错误写 stderr 并非零退出。
stdin 每行一个完整 SMILES，内部空行保留位置，空 stdin 返回空数组；最后的换行仅结束当前行。
REINVENT 的行协议本身不能区分空列表与单个空字符串，也不能区分末尾空条目和行结束符；
正常 REINVENT 会预先屏蔽无效分子，本服务 HTTP 数组协议没有该歧义。
超时默认 120 秒，无自动重试；超时不代表服务端已经取消推理。

## RL 示例运行与停止

两份模板都是单阶段、seed 42、batch size 16、最多 5 步、DAP sigma=128/rate=0.0001，
启用 SMILES 随机化和 isomeric SMILES。只有 PLANET 评分、权重 1，无额外性质奖励或多样性过滤。
单评分采用 arithmetic_mean；原始亲和力由 REINVENT 内置 sigmoid 转换，不在客户端重复归一化。
输入分别为 `c1ccccc1*` 和 `c1ccccc1*|*c1ccccc1`。

```bash
"$REINVENT_PY" -m reinvent.Reinvent -l "$RUN/libinvent/run.log" "$RUN/libinvent.toml"
"$REINVENT_PY" -m reinvent.Reinvent -l "$RUN/linkinvent/run.log" "$RUN/linkinvent.toml"
kill -TERM "$(cat "$RUN/server.pid")"
```

前台启动服务时也可 Ctrl-C 停止。停止后重复独立评分命令必须非零退出，不应产生伪造零奖励。
日志、CSV、checkpoint、TensorBoard 和解析后配置写入两种模式各自的目录。
CSV 的 `PLANET reward` 为 REINVENT sigmoid 转换后的奖励；`PLANET reward (raw)` 为送入 sigmoid 的亲和力，
因此有效输入是 PLANET 原始预测，无效输入是有限 sentinel `-1000000.0`。
真实 nullable 预测在 `planet_affinity (PLANET reward)`，另有 status、error、target_id、oracle_id metadata 列。
JSON 无效预测为 null；当前 REINVENT CSV writer 会将它写为字符串 `None`。
被 REINVENT 自己屏蔽、未送到 oracle 的输入，其 metadata 由上游填充，不能当作 PLANET 已评分。

## 推理语义和协议

[PLANET 官方说明](https://github.com/ComputArtCMCG/PLANET)定义输入为蛋白口袋三维图及配体二维图。
本接口复用本地 `PLANET(300,8,300,300,3,10,1,device)`、`ProteinPocket` 的 12 Å 规则和
`mol_batch_to_graph(..., auto_detect=False)`；不逐分子 docking 或生成构象。
SMILES 经 RDKit sanitize 和 `Chem.AddHs`，保留输入电荷与立体信息；不做中和、互变异构枚举或质子化选择。
空分子、无效 SMILES、dummy atom、不连通结构、加氢后邻接度超过上游 MAX_NB 为 invalid_input。
当前原版图没有额外总原子数上限；显存不足属于推理故障，不是输入零分。
启动审计报告上游会跳过的残基，检查每个纳入口袋的残基恰有一个有效 Cα，并检查坐标与特征有限。

`GET /health`：version、ready、target_id、oracle_id、device、pocket。
`POST /score` 请求：

```json
{"version": 1, "target_id": "example_target", "smiles": ["CCO", "invalid", "CCO"]}
```

响应包含相同 version、target_id、oracle_id 和 results；每条结果具有原始 index、affinity、status、error。
成功为 `status="ok"`、有限 affinity、error=null；无效输入为 `status="invalid_input"`、affinity=null、错误代码字符串。
空数组返回空结果，重复输入保留位置，不排序也不丢弃。第一版无服务端分子缓存。
请求格式错误 HTTP 400，错误靶点 409，模型/批图构建异常及非有限输出 500，失败影响整个请求。
客户端核对版本、靶点、指纹、长度、顺序、索引、状态与数值；连接/超时/HTTP/校验错误都使 REINVENT 评分失败。

客户端 payload 的 `planet_affinity_for_scoring` 为 REINVENT endpoint 数值输入：有效分子使用原始 affinity，
无效分子使用固定有限 sentinel `-1000000.0`。`planet_affinity` 始终记录真实预测，无效时保持 null；
status 和 error 明确区分模型低预测与不可评分输入。模型/通信/协议故障不会转换为 sentinel，而是使评分失败。

REINVENT 对有效亲和力 p 计算
`1 / (1 + 10 ** (-10*k*(p-(low+high)/2)/(high-low)))`。模板要求 high > low、k > 0，
默认 low=4、high=10、k=0.5，p=7 对应 0.5；sentinel 在该 float32 sigmoid 下得到严格 0。
更换变换类型或使用极端奖励参数时必须重新验证无效输入仍为 0。

## 验证

已有 REINVENT 环境中的本机自动化测试（不安装环境、不运行真实模型/训练）：

```bash
cd "$PROJECT"
PYTHONPATH="$PROJECT/src:$PROJECT/REINVENT4" "$REINVENT_PY" -m pytest tests/planet_oracle -q
```

测试使用可注入假预测器，只在测试代码中定义；生产命令没有假模型模式。
覆盖真实 HTTP/客户端 subprocess/ExternalProcess、错误语义、REINVENT sigmoid 与 sentinel 零分、RDKit、口袋审计、
两份配置解析、真实片段拼接到评分器、CSV metadata 与 prior 词表兼容性；这些不代表真实 PLANET/CUDA 验收。

Linux CPU 真实验收必须先通过；CUDA 可用时再重复（输出目录须不存在）：

```bash
"$PLANET_PY" -m planet_oracle.smoke --planet-root "$PROJECT/PLANET" \
  --device cpu --output-dir "$RUN/smoke-cpu"
"$PLANET_PY" -m planet_oracle.smoke --planet-root "$PROJECT/PLANET" \
  --device cuda:0 --output-dir "$RUN/smoke-cuda"
```

smoke 使用 PLANET 自带 demo 蛋白/参考配体，从 mols.sdf 前五个有效分子导出相同 SMILES，
独立初始化官方 `PlanetEstimator`，通过 `VS_SMI_Dataset` 的官方 SMILES 路径推理，
与适配器的 HTTP 多批次、单分子、重复请求比较（atol=rtol=1e-4）。
同时检查顺序、无效输入/空批次、服务停止后失败，以及常驻模型/口袋/特征对象未重建。
保存 server.json、manifest.json、demo.smi、report.json；任何不一致非零退出，不宣称自动回退设备成功。

最后运行上面的两种 RL，各确认 5 步、有限奖励、CSV metadata、resolved.json 和非空 agent.chkpt。
这两次真实 RL 验收由服务器运行，不包含在本机测试内。后续 16 靶点正式实验仍需单独设计等价片段、
准备蛋白、校准奖励并控制双方实验条件。
