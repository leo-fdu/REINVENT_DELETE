# 截断PLANET分数及生成尺寸惩罚实验

参数：low=4、high=9、k=0.25、affinity_cap=9、generated_heavy_atom_penalty=0.01；
seed=42，30任务 ×1000次采样（batch=100 ×10轮）。重复乘数0.5作用于整个扣尺寸后的奖励，允许负数。

- `all_generated.csv.gz`：保留所有尝试的压缩结果及完整奖励分解。
- `acceptance-report.json`：逐任务预算、原始输入、分数截断、sigmoid、生成原子数、尺寸扣分及重复系数验收。
- `training-numerics.json`：最终权重有限性、实际参数更新及训练数值检查。
- `checksums.json`：完整原始run目录的SHA-256清单；相对下述镜像，不是本精简目录。
- `transfer-verification.json`：完整本地镜像哈希验证。
- `first-batch-baseline-check.json`：与high9首轮SMILES及原始PLANET分数对照。
- `environment.json`、`environment-freeze.txt`、`runtime-requirements.txt`：真实运行环境。
- `source-transfer-manifest.json`、`run.json`、`experiment.json`：代码/模型/输入及运行配置。
- `smoke-validation.json`：正式运行前，真实模型四个任务共80次尝试的奖励公式核验。
- `run_cloud.sh`、`planet_python_nocache.sh`：实际云端运行与PLANET显存设置。

完整原始镜像保留在`results/server-evidence/high9_cap9_size001_1000_seed42_20261025/run/`，
包含全部checkpoint、原始CSV、日志和TensorBoard，沿用已有Git忽略规则。
分析目录：`evaluation/descriptor/generation_analysis_reinvent/high9_cap9_size001_1000_seed42_20261025/`。

正式运行30/30任务通过验收，共30,000次尝试、29,722条有效成功评分、278条无效。
300条负训练奖励均保留；全部30个最终模型参数有限、确实更新，训练数值无NaN/Inf。
与high9首批3,000个样本的SMILES和原始分数完全一致。
实现提交为7a03f8b，本地97项相关测试通过；云端真实短测4任务共80次尝试通过后启动正式运行。

完整原始镜像约2.5 GB，349个文件已逐一通过云端SHA-256清单核对，
详见精简结果目录的`transfer-verification.json`。云端驱动正常退出（exit=0）。
