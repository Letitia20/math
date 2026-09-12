# 参考文献与公开内容相似性审查（2026-09-12）

## 结论摘要

- 参考文献 [2]、[3] 的 DOI、作者、题名、期刊、年份、卷期和页码均可由 Crossref 精确解析，未发现伪造题录或 DOI 不匹配。
- 参考文献 [1] 是本题官方赛题及附件的资料型引用，不是学术论文。其依据是本项目保存的题面、附件和竞赛规则文件；该类内部竞赛资料通常不在 Crossref 或 OpenAlex 建立论文记录，因此不能用学术数据库 DOI 方式验证。
- 对公开 GitHub 同题仓库进行了定向检索和代码/文本比对。发现两个高相关公开仓库，但未发现本文的结果数字、连续原句、函数实现片段或独特流程的实质性复制证据。
- 不能据此宣称“全网绝无雷同”：本审查覆盖的是可访问的公开索引和已下载候选仓库，结论应表述为“在给定检索范围内未检出实质性雷同”。

## 一、参考文献核验

### [1] 竞赛题面资料

题录：全国大学生数学建模竞赛组委会. 2026 年高教社杯全国大学生数学建模竞赛 A 题：药材的烘干问题及附件[Z]. 2026.

核验方式：与本项目使用的官方题面、附件 1/2 及题目规则文件逐项对照。该条目是本题输入资料的来源标识，不声称存在 DOI、期刊卷页或公开论文版本。它不应被误读为一篇可在 Crossref 检索的学术论文。

### [2] Journal of Food Engineering 论文

题录：da Silva WP, e Silva CMDPS, Gama FJA. *Estimation of thermo-physical properties of products with cylindrical shape during drying: The coupling between mass and heat*. Journal of Food Engineering, 2014, 141: 65–73. DOI: `10.1016/j.jfoodeng.2014.05.010`.

Crossref 查询（2026-09-12）：

`https://api.crossref.org/works/10.1016%2Fj.jfoodeng.2014.05.010`

返回字段与论文一致：publisher=Elsevier BV；type=journal-article；volume=141；page=65-73；published=2014-11；三位作者姓名与题名均一致。未发现 DOI 指向其他题名或作者的情况。

### [3] Foods 论文

题录：Adrover A, Venditti C, Brasiello A. *A non-isothermal moving-boundary model for continuous and intermittent drying of pears*. Foods, 2020, 9(11): 1577. DOI: `10.3390/foods9111577`.

Crossref 查询（2026-09-12）：

`https://api.crossref.org/works/10.3390%2Ffoods9111577`

返回字段与论文一致：publisher=MDPI AG；type=journal-article；volume=9；issue=11；page=1577；published-online=2020-10-30；三位作者姓名与题名均一致。OpenAlex 同 DOI 记录显示 `is_retracted=false`，主题为 Food Drying and Modeling；未发现撤稿或 DOI 错配迹象。

## 二、公开仓库检索范围

### 候选仓库

1. [Daxcq/cumcm-herbal-drying](https://github.com/Daxcq/cumcm-herbal-drying)，检索到的提交：`15b9c61b98e6491f559810eceec5104ec1272415`（2026-09-11）。
2. [asdfghjkoiuy/cumcm2026-a-herb-drying](https://github.com/asdfghjkoiuy/cumcm2026-a-herb-drying)，检索到的提交：`1025cc902c525c5737344c72063d69b0c63a86f6`（2026-09-11）。

GitHub 仓库主题搜索“药材 烘干 数学建模”返回 2 个结果；搜索“drug drying finite volume moving boundary”未返回结果。对候选仓库的 README、Markdown/文本、Python、JavaScript 和 HTML 文件进行了本地读取。

### 检查项目

- 论文中文段落与候选仓库文本的最长连续中文片段；
- 论文和核心代码中的结果数字（包括 57.4740、51.0920、129.8481、10.97% 等）是否在候选仓库出现；
- Python 代码 token 序列的最长连续匹配、函数名集合和模块结构；
- 公式/流程关键词：有限体积、BDF、Richardson、移动边界、归一化坐标、PCHIP、阈值事件、密度—收缩相容性；
- GitHub repository search 对上述独特结果数字和组合短语的定向查询。

## 三、比对结果与解释

两个候选仓库都独立处理同一年度、同一 A 题，因此共享“一维径向圆柱”“传热—传质 PDE”“有限体积离散”“收缩半径”等题目驱动的通用概念是预期现象，不能单凭这些术语认定抄袭。

本项目与候选仓库的可区分特征包括：

- 本项目以 `solve_ivp(method="BDF")`、对数含水率状态、稀疏 Jacobian 结构和事件函数处理四问；候选仓库主要采用节点中心全隐式/Crank–Nicolson 及 Picard 或显式分段推进，函数组织和求解器接口不同。
- 本项目有独立的 `problem4_independent.py` 单元中心复算、密度—收缩相容性诊断、跨问结果审计和 95 项回归测试；候选仓库未出现同名函数、同一测试集合或同一审计链。
- 对候选 Python 文件进行 token 序列比对时，未发现长度达到 12 个 token 的连续匹配块；未发现本文核心函数名组合在候选仓库中重复出现。
- 论文文本与候选文档未检出超过短语级的连续中文匹配；结果数字组合的 GitHub 仓库定向查询均返回 0 个结果。

因此，在本次可重复、可记录的公开检索范围内，没有证据表明论文代码、解题流程或文字是从上述公开仓库复制而来。论文仍使用了传热传质与移动边界领域的公开基础思想，并已通过 [2]、[3] 进行正常学术引用；“思想来源可追溯”与“代码/文字复制”是两件不同的事。

## 四、限制与提交前建议

本审查不是查重系统，也没有访问付费数据库、封闭仓库、搜索引擎全部网页缓存或比赛内部稿件。若学校或竞赛组委会提供指定查重系统，应在提交前用该系统复核；若队员继续修改模型、代码或论文，应重新运行本报告中的数字、文本和代码指纹检查，并保存新的提交版本记录。
