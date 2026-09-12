# A题支撑材料

本包为论文的数值支撑材料。解压后在本目录运行命令。论文及AI说明的可编辑Word另行提供，参赛队应在提交前人工审定。

## 结果文件

- outputs/result1.xlsx：第一问每秒、每0.1 cm温湿分布。
- outputs/result2.xlsx：原第二问前三小时每秒温湿分布，保留原文件。
- outputs/result2_full_process.npz：第二问1—206940 s完整逐秒温湿分布，包含时间列和21个径向位置；无损压缩的是四位小数交付值，不是未舍入求解器状态。
- outputs/result3.xlsx：固定半径长期含水率、临界行和严格整分钟行。
- outputs/result4.xlsx：收缩模型含水率、实际表面及实体外空白。
- reports/data/problem2_full_process_audit.json：完整Excel在压缩前逐单元回读的记录。
- reports/data：网格、独立复算、情景敏感性和守恒诊断。
- AI工具使用详情.pdf：规定的AI使用详情说明。

第二问全程Excel约27.7 MB，不能直接塞入20 MB以内的支撑包。本包以NPZ保存全部表格数值，通过下列命令恢复成可编辑Excel。恢复结果行列与完整Excel一致，但ZIP内部元数据不保证逐字节相同。全程结果与第三问独立分钟积分最大舍入差为0.0001 kg/kg；前3小时与原第二问逐值一致。

```powershell
python -m pip install -r requirements.txt
python tools/restore_full_result2.py
```

恢复文件为 outputs/result2_full_process.xlsx；如正式提交要求该文件名必须为result2.xlsx，应在单独的提交目录中重命名恢复文件，不覆盖保留的前三小时历史文件。

## 从原始数值输入重新计算

Python 3.11以上。在本目录执行：

```powershell
$env:PYTHONPATH='src'
python -m pytest -q
python scripts/run_problem1.py
python scripts/run_problem2.py
python scripts/run_problem3.py
python scripts/run_problem4.py
python tools/rebuild_delivery_workbooks.py
python tools/extend_problem2_delivery.py
```

Linux/macOS使用 export PYTHONPATH=src。重算工作簿保存在 outputs/recomputed/，不覆盖交付值。全程扩展脚本计算2560/5120网格，运行较久，产生临时网格缓存及 deliverables/supporting_materials/outputs/result2.xlsx。文件可用 tools/compress_problem2_full_process.py 转成紧凑NPZ。

原始输入为 data/raw/attachment1.csv、attachment2.csv，已与题目附件逐值核对。95项自动化测试是程序检查，并非参赛队人工审定或真实设备实验验证。

## 解释边界

第三、四问4 h之后沿用附件1末一小时均值平台。采用一维中截面与有效表面传质；第四问经验密度与收缩存在相容性问题。临界值0.1500不等于严格达标，严格整分钟分别为206940 s和183960 s。所有情景范围不是统计置信区间。

目录中的既有分问报告用于保留数值来源，其中第二问报告描述前三小时版本；本包新增的全程输出以上述说明为准。论文、源码和AI说明均须由参赛队最终人工核验。
