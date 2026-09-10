# 2026 数学建模 A 题

当前已完成第一问固定半径预热模型，以及第二问变物性温湿耦合模型、数值求解和结果工作簿。

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
