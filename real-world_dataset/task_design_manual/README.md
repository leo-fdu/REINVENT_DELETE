# 手动设计 16 个靶点的 LibINVENT / LinkINVENT 任务

这是一个在自己电脑上运行的网页界面。浏览器用于选择切割键、保留组分和填写备注；本地 Python 服务用 RDKit 切分配体、预览结构，并在点击“记录全部设计”时核对和保存结果。它只记录任务设计，不执行分子生成、对接或 DELETE 实验，也不修改原始配体。

## 使用前检查

- Python **3.11 或更新版本**，并在同一环境中安装 RDKit、PyTorch、Pillow 及 REINVENT4 所需依赖。依赖安装方式见 `REINVENT4/README.md`；本机已有可用的 `reinvent4` Conda 环境。
- 仓库根目录下须有可用的 `REINVENT4/` 代码，以及 `REINVENT4/priors/libinvent_transformer_pubchem.prior` 和 `REINVENT4/priors/linkinvent_transformer_pubchem.prior`。启动界面要导入 REINVENT4 代码；记录已设计的任务时要加载 prior 来检查模型词表和输入长度。
- `real-world_dataset/` 下的 16 个 `crystal.mol2` 和 `task_design_astra/` 须保留在原位置。程序会按图册记录的 SHA-256 核对原始 MOL2；文件内容变化时需要重新核对设计。

**刚从 GitHub 克隆此仓库的人，需要先补齐 REINVENT4。** 主仓库目前把 `REINVENT4` 记为 Git 子模块引用，却没有 `.gitmodules`，因此普通 `git clone` 不会自动取回其内容。本机所用的两个 prior 是 REINVENT4 中 `8562f4a` 提交新增的文件；推送主仓库不会一并分发该提交或文件。仅克隆官方 REINVENT4 代码也不能保证具备这两个 prior。请从有权使用的来源取得兼容的 REINVENT4 代码与 prior，并放在上述路径；在这些文件齐备前，克隆出的界面不能完整使用。

## macOS 启动

打开“终端”，进入**本仓库根目录**并激活已安装依赖的环境：

```bash
cd /path/to/REINVENT_DELETE
conda activate reinvent4
python real-world_dataset/task_design_manual/app.py
```

本机现有环境也可直接用以下命令启动，无须先激活 Conda：

```bash
/opt/homebrew/Caskroom/miniforge/base/envs/reinvent4/bin/python real-world_dataset/task_design_manual/app.py
```

## Windows 启动

打开 **Anaconda PowerShell Prompt**（或已配置 Conda 的 PowerShell），进入**本仓库根目录**：

```powershell
Set-Location C:\path\to\REINVENT_DELETE
conda activate reinvent4
python .\real-world_dataset\task_design_manual\app.py
```

Windows 上若在点击“记录”时出现 `ZoneInfoNotFoundError: 'Asia/Shanghai'`，在同一环境中执行 `python -m pip install tzdata`，然后重新启动界面。REINVENT4 的 Windows 支持在其文档中标为测试较少；请先确保该平台上的 REINVENT4 依赖和两个 prior 可正常加载。

两种系统启动成功后，终端会显示 `http://127.0.0.1:8766`。在**运行服务的同一台电脑**上用浏览器打开这个地址。保持终端进程运行；按 `Ctrl+C` 停止。关闭进程或重启电脑后，需重新运行启动命令。若 8766 端口被占用，可在命令末尾加 `--port 8767`，并改为打开 `http://127.0.0.1:8767`。

这是本机地址，不是 GitHub 网页：其他电脑不能通过你的 `127.0.0.1` 访问此服务。直接双击本目录的 `index.html` 也不能使用编辑和记录功能，因为它需要 Python 后端。

## 设计与记录

1. 选择靶点和任务。点击二维配体图中的**非环单键**；再次点击可取消。LibINVENT 选 1 条切割键，LinkINVENT 选 2 条。
2. 点击切出的保留组分。LibINVENT 保留 1 个组分，LinkINVENT 保留 2 个各有一个连接点的端组分。剩余部分是原配体中的参考待生成区域。可以填写选择理由。不适合设计的任务可勾选“暂不设计此任务”，并在备注中说明原因；取消勾选可以继续编辑。
3. 界面把未记录的选择保存在**当前浏览器、当前地址**的 `localStorage`。重开同一地址可继续编辑；更换浏览器或端口后不会自动带过去。原 MOL2 哈希变化时旧草稿不会载入。
4. **随时**点击“记录当前设计”，无需完成全部 32 项。记录会逐项标为“已设计”“跳过”“未完成”或“尚未开始”，保留跳过理由和未完成的切割/组分草稿。后端对“已设计”任务核对原子覆盖、连接点、REINVENT4 实际拼接、模型标准化、prior 词表及 180 token 长度；无法通过时会说明原因并拒绝将该项作为已设计任务记录。立体化学不能完全恢复时会在记录中明确标出。

每次成功记录都在 `task_design_manual/` 下创建一个**新的** `YYYY-MM-DDTHH-MM-SS-ffffff+0800/` 文件夹，不覆盖旧记录。文件夹包含：

- `designs.json`：列出全部 16 个靶点的 32 个任务槽位及各自状态。已设计任务记录源文件哈希、原 MOL2 原子 ID、切割键、保留和参考待生成组分、SMILES、token 数、验证结果及备注；跳过和未完成任务记录已有选择与备注。LinkINVENT 已设计任务的两个保留端按参考 linker 的出口顺序排列。
- `index.html`：可单独在浏览器打开的只读汇总报告；它与编辑界面的 `index.html` 用途不同。
- `figures/<靶点>_<任务>.svg`：32 张切分图，蓝/绿为保留组分，橙色为原配体参考待生成区域，红色为切割键。

时间记录文件夹由本目录的 `.gitignore` 忽略，**不会随 Git 提交或 push 自动分享**；要交给别人，请自行复制对应的整个时间文件夹。对比 DELETE 时，仍需从同一原始结构提取相同的保留原子集合和晶体坐标。
