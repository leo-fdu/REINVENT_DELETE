# 手动设计 16 个靶点的 LibINVENT / LinkINVENT 任务

从仓库根目录运行：

```bash
/opt/homebrew/Caskroom/miniforge/base/envs/reinvent4/bin/python real-world_dataset/task_design_manual/app.py
```

在本机浏览器打开 `http://127.0.0.1:8766`。可用 `--port` 改端口。服务只监听 `127.0.0.1`。需要该 REINVENT4 环境中的 RDKit、PyTorch、Pillow 和本仓库的 REINVENT4 代码及 LibINVENT/LinkINVENT prior；无需新增依赖或联网。

## 操作

1. 选择靶点与任务。点击原配体二维图上的非环单键。LibINVENT 选 1 条，LinkINVENT 选 2 条；再次点击可取消。
2. 切分图出现后，点击要保留的组分。LibINVENT 保留 1 个组分，LinkINVENT 保留 2 个各有一个出口的端组分。其余部分是原配体中的参考待生成区域。
3. 可写任务设计理由。浏览器会把未记录的选择保存在本机 `localStorage`，重新打开页面可继续编辑。原 MOL2 的哈希变化时旧草稿不会载入。
4. 32 项全部完成后点击“记录全部设计”。后端重新核对全部选择和模型输入，再在本目录新建时间命名的记录目录。每次点击都创建新目录，已有记录绝不覆盖。

后端采用 [task_design_astra](../task_design_astra/README.md) 中核对的 16 个参考配体结构及原 MOL2 原子 ID，并在启动时和记录时核对原始文件 SHA-256。原始 MOL2 不会改动。仅接受非环单键切分；要求 LibINVENT 切出两个连续组分、LinkINVENT 切出三个连续组分，并检查保留端与待生成区域的连接点数量、原子覆盖、REINVENT4 实际拼接、模型标准化、prior 词表及 180 token 长度。连接图无法恢复时拒绝记录；立体化学无法完全恢复时，在记录中明确标为不一致。

## 记录格式

每次记录生成 `YYYY-MM-DDTHH-MM-SS-ffffff+0800/` 目录，包含：

- `designs.json`：完整结构化记录。每个靶点有 LibINVENT 和 LinkINVENT 两项，记录源文件哈希、切割键、各组分的原 MOL2 原子 ID、重原子数、连接点数、模型输入 SMILES、参考输出 SMILES、token 数、拼接检查和备注。LinkINVENT 的两个保留端按参考 linker 的出口顺序排列，可直接辨认 `片段A|片段B` 的输入次序。
- `index.html`：无需服务也能打开的可视化汇总报告。
- `figures/<靶点>_<任务>.svg`：32 张切分图。蓝/绿为保留部分，橙色为原配体参考待生成部分，红色为切割键，数字为 MOL2 原子 ID。

记录目录只保存在本地，不由 Git 追踪；界面代码、说明和 `.gitignore` 则由 Git 追踪。这里的参考输出来自原配体，仅用于记录与完整性校验，不代表模型实际生成结果、三维结合约束或已验证的药效团。对比 DELETE 时还需从同一原始结构提取相同的原子保留集合和晶体坐标。
