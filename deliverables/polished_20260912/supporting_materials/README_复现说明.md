# A题支撑材料

本压缩包采用单层结构，所有文件均直接位于根目录，不含子文件夹。`result1.xlsx`、`result2.xlsx`、`result3.xlsx`、`result4.xlsx` 已填写计算结果，沿用题目指定文件名和现有表格格式；打包过程仅逐字节复制，没有用表格软件重新保存。

`attachment1.csv` 和 `attachment2.csv` 为数值输入。支撑材料共含 25 份 Python 源码：四问模型与独立验证模块、运行入口、测试程序以及工作簿构建和审计工具。`problem1.py` 至 `problem4.py` 为附录中展示的核心模型，完整源码均保存在本目录。安装依赖后可先运行全部测试，再依次运行四问入口：

```powershell
python -m pip install -r requirements.txt
python -m pytest -q -p no:cacheprovider
python run_problem1.py
python run_problem2.py
python run_problem3.py
python run_problem4.py
```

四个结果工作簿是已审定交付值，运行脚本不会覆盖它们。第二问全程逐秒结果以 `result2_full_process.npz` 保存；运行 `python restore_full_result2.py` 可在根目录生成 `result2_full_process.xlsx`，脚本发现同名文件时会拒绝覆盖。

`flat_support_audit.json` 记录四个结果工作簿的 SHA-256、工作表尺寸、非空单元格数和模板样式核验结果。第四问中药材收缩后位于实体外的固定位置按题意留空，不属于漏填。

`AI 工具使用详情.pdf` 为按竞赛规定提供的独立披露文件。本题核心建模与分析由参赛队主导并审定；AI 仅用于代码调试、重复性检查、文件核验和语言润色等辅助环节，所有输出均经队员审查。论文参考文献仅收录赛题与学术文献，不将 AI 工具或技能项目列入参考文献。
