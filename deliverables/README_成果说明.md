# A题论文与支撑材料

## 直接使用

- **A题论文_可编辑版.docx**：完整论文，含原生Word公式、可编辑表格和完整代码附录。优先在Microsoft Word中修改。
- **A题论文_提交预览.pdf**：与可编辑版对应的PDF，摘要第一页，正文随后，末尾为附录。用于预览及队伍人工审定后的提交。
- **A题支撑材料.zip**：代码、输入CSV、结果表、图像、数值审计、AI使用详情和复现说明。
- **AI工具使用详情.pdf**、**AI工具使用详情_可编辑版.docx**：已在支撑包内另放一份。
- **result2_full_process.xlsx**：第二问1—206940 s逐秒温湿分布。该表约27.7 MB，供本地直接查看；支撑包内用约3.6 MB的NPZ无损保存同一四位小数数据，可恢复Excel。

论文由既有四问模型及数值记录整理而成，保留一维中截面、长期环境平台、有效传质边界及密度—收缩相容性的明确限制。没有新增未经计算的实验精度结论。原有仓库文件未覆盖。

## 编辑和重新导出

直接修改Word后，在Word选择“另存为/导出PDF”。这不会自动同步Markdown；请选定Word或Markdown作为自己的主编辑源。

需要脚本重建时，在仓库根目录依次运行：

```powershell
node tools/render_paper.mjs
python tools/finalize_word_math.py
pwsh -File tools/export_paper_pdf.ps1
python tools/check_paper_layout.py
python tools/package_paper_delivery.py
```

Node需要docx包；公式转换需要Python的lxml、latex2mathml以及本机Microsoft Office的MML2OMML.XSL；PDF导出需要Windows、Word和PowerShell 7。版式检查需要PyMuPDF、Pillow。数值复现依赖另见支撑包requirements.txt。

paper_source.md为正文模板；paper_complete.md为填充表格和完整源码后的可编辑纯文本备份。equations.json为生成公式的中间索引。脚本重建会覆盖本目录生成的Word/PDF，因此不要用它覆盖尚未回填到Markdown中的手工Word修改。

## 提交前由队伍确认

自动化验证已完成，不替代参赛队对核心模型、AI输出和参考文献的人工审定。电子论文和支撑ZIP分别控制在20 MB以内；电子文件不含承诺书或编号页。不要把包外的大Excel与全部本目录文件再次合成一个提交包。

截至本次检查，GitHub原仓库的API返回private=false，与“仅队内观看”要求不一致。因此成果先保存在本地并准备Git提交，未公开上传本次论文。待原仓库设为私有后可继续推送，不需要重做论文。
