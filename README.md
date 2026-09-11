# 2026 数学建模 A 题

四个小问的模型、数值求解、论文分节和结果工作簿已完成，并经过最终复算审查。总入口为 [最终审查记录](reports/final-review.md)。第二问当前逐秒输出覆盖表 3、表 4 的前三小时；原始模板没有限定总行数，相关范围解释见 [第二问](reports/problem2.md)。

| 问题 | 正文 | 结果 |
|---|---|---|
| 1：固定半径预热 | [第一问](reports/problem1.md) | [result1.xlsx](outputs/result1.xlsx) |
| 2：变物性耦合 | [第二问，含表 3、4](reports/problem2.md) | [result2.xlsx](outputs/result2.xlsx) |
| 3：固定半径终点 | [第三问](reports/problem3.md) | [result3.xlsx](outputs/result3.xlsx) |
| 4：收缩与移动边界 | [第四问](reports/problem4.md) | [result4.xlsx](outputs/result4.xlsx) |

下面历史命令中的 Python 路径可替换为本机 Python 3.11+。安装计算和审计依赖：`python -m pip install -e . openpyxl pytest`。

## 复现第一问

```powershell
$env:PYTHONPATH='src'
D:\miniconda3\python.exe -m pytest -q
D:\miniconda3\python.exe scripts/run_problem1.py
```

随后使用 Node.js 与 `@oai/artifact-tool`，把 `tmp/problem1_result.json` 写入题目提供的 `result1.xlsx` 模板：

```powershell
node tools/build_result1.mjs TEMPLATE_PATH tmp/problem1_result.json outputs/result1.xlsx
```

主要产物：

- `reports/problem1.md`：第一问模型、算法、结果、检验与参考文献；
- `outputs/result1.xlsx`：题目要求的温度和水分浓度完整结果；
- `reports/figures/`：径向分布与代表位置时间历程图；
- `docs/problem1-formula-audit.md`：题面公式与符号复核记录。

## 复现第二问

```powershell
$env:PYTHONPATH='src'
D:\miniconda3\python.exe scripts/run_problem2.py
node tools/build_result2.mjs TEMPLATE_PATH tmp/problem2_result.json outputs/result2.xlsx
```

第二问主要产物为 `outputs/result2.xlsx`、`reports/figures/problem2_*.png`、`reports/verification.md` 和 `docs/problem2-formula-audit.md`。

## 复现第三问

```powershell
$env:PYTHONPATH='src'
D:\miniconda3\python.exe scripts/run_problem3.py
node tools/build_result3.mjs TEMPLATE_PATH tmp/problem3_result.json outputs/result3.xlsx
```

第三问主要产物为 `outputs/result3.xlsx`、`reports/figures/problem3_*.png`、`reports/verification.md` 和 `docs/problem3-formula-audit.md`。基准终点采用附件 1 最后 1 h 平台均值，并通过长期边界敏感性方案核查。

## 第四问与统一审计

```powershell
$env:PYTHONPATH='src'
python scripts/run_problem4.py
python tools/build_result4.py TEMPLATE_PATH tmp/problem4_result.json outputs/result4.xlsx
python tools/audit_all_results.py --attachment-dir '原始A题/附件'
python -m pytest -q
```

统一审计前需运行四个 `scripts/run_problem*.py` 生成载荷。审计只读四份工作簿，与重新计算的载荷逐单元格比对，同时核验 CSV 与原始附件以及论文表 1–6。第四问更多数值实验见 [复核说明](reports/problem4-review-response.md)。
