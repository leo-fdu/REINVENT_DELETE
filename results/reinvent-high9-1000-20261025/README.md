# high=9、每任务 1,000 次 RL 采样

本目录保存本次独立重跑的压缩生成表、配置与验收证据，不是旧运行数据的截取。

- `all_generated.csv.gz`：30 任务 × 10 轮 × 100 行，共 30,000 次尝试，保留无效和重复。
- `experiment.json`：low=4、high=9、k=0.25，batch=100、steps=10、seed=42。
- `acceptance-report.json`：逐任务预算、样本编号、输入哈希、评分和奖励公式、重复惩罚、checkpoint 验收。
- `training-numerics.json`：30 个 checkpoint 均为有限参数且相对 prior 有更新，likelihood 有限。
- `transfer-verification.json`：完整原始镜像逐文件 SHA-256 核对结果。
- `first-batch-baseline-check.json`：与旧实验首轮 3,000 个样本逐一比较，SMILES 和原始分数完全一致。
- `checksums.json`：云端原始运行目录的文件清单；路径相对下述完整本地镜像，并非本精简目录。
- `environment.json`、`environment-freeze.txt`、`runtime-requirements.txt`：实际运行环境。
- `run.json`、`source-transfer-manifest.json`、`pocket-validation-report.json`：输入、模型及源码来源。
- `run_high9_cloud.sh`、`resume_high9_cloud.sh`：实际云端执行过程记录，含原服务器绝对路径。

完整原始镜像：`results/server-evidence/high9_1000_seed42_20261025/run/`。
该目录按仓库已有规则忽略，保留 30 个最终 checkpoint、所有原始 CSV、日志和 TensorBoard 文件。
同级目录保存部署及短测证据。初次短测因服务器缺少 libXrender 失败，补装 libxrender1/libsm6
后以新目录重新短测，6 个短测任务通过；正式 30 任务没有重试。

云端为 RTX 4090，Python 3.11.17、PyTorch 2.12.0+cu126、RDKit 2026.03.6、
NumPy 2.4.6、SciPy 1.17.1、pandas 3.0.6。原锁文件的下载 URL 失效，安装使用独立必需依赖清单
和镜像，核心科学计算版本按旧服务器记录固定；没有修改模型代码或权重。

本地全部分析位于：
`evaluation/descriptor/generation_analysis_reinvent/high9_1000_seed42_20261025/`。
其 README 给出实测结果、同预算对照与解释边界。
