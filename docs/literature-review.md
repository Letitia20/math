# A 题同类干燥模型文献对比

检索日期：2026-09-10

## 1. 检索范围与筛选原则

本次检索围绕“中药材或植物物料热风干燥、圆柱内部温度和含水率分布、变物性、收缩与移动边界”展开。正文比较仅使用能够核实作者、题名、期刊、年份和 DOI 的研究论文。综述、聚合页面与二次引用只用于发现线索，不作为模型结论的证据。

这些文献研究的物料包括西洋参、香蕉、榅桲、马铃薯、茄子、梨及一般生物多孔介质。它们不能直接证明题中药材具有相同参数；可借鉴的是守恒方程、边界形式、坐标处理和验证方法。题目给出的经验物性公式与附件数据仍是本题参数的首要依据。

## 2. 文献中的主要模型形式

### 2.1 薄层经验模型

吴小华等针对西洋参分段热风干燥比较了多种经验动力学模型，采用 Modified Page 模型描述各阶段平均水分比，并通过阶段转换计算处理工况改变；论文报告其试验对比最大相对误差为 7.44%、平均相对误差为 1.78% [1]。

典型形式可写为

\[
MR=\frac{\bar C-C_e}{C_0-C_e}=\exp[-(kt)^n],
\]

分段工况下需要把上一阶段终态转换为下一阶段的等效初始时间或重新定义阶段时间。

**与 A 题的关系：** 该模型与“中药材、分阶段热风干燥”最接近，适合解释平均干燥曲线和做结果合理性对照。但它只给出平均含水率，不能直接产生题目要求的径向温度场与含水率场，也不能自然处理表面传质边界。因此不作为四问主模型，可作为第三问平均含水率曲线的辅助检验。

### 2.2 固定域 Fick–Fourier 模型

Silva 等采用圆柱一维表观液态扩散模型，在表面设置第三类传热传质边界，以全隐式有限体积法求解，并同时考虑变物性与收缩 [2]。Tzempelikos 等对圆柱形榅桲片建立一维非稳态热质传递模型：内部热量按导热处理，水分按 Fick 扩散处理，扩散系数使用 Arrhenius 温度关系，换热与传质系数由外部 CFD 计算 [3]。Kaya、Aydin 和 Dincer 也对圆柱湿物体使用隐式有限差分，并以 CFD 研究局部换热和传质系数 [4]。

与这些工作对应的本题基准形式为

\[
\rho c_p\frac{\partial T}{\partial t}
=\frac{1}{r}\frac{\partial}{\partial r}
\left(rk\frac{\partial T}{\partial r}\right),
\]

\[
\frac{\partial C}{\partial t}
=\frac{1}{r}\frac{\partial}{\partial r}
\left(rD\frac{\partial C}{\partial r}\right),
\]

中心采用对称条件，表面采用第三类边界。题目第二至第四问使 \(D=D(C,T)\)，并使 \(\rho,c_p,k\) 随 \(C\) 变化，因此两场通过系数构成双向或单向的非线性耦合。

**与 A 题的关系：** 这是第一至第三问最合适的主模型。题目明确给出了导热、扩散和对流系数或经验公式，同时明确要求径向分布；模型复杂度与数据可识别性匹配。全隐式有限体积法也适合圆柱中心奇点与长时积分。

### 2.3 二维轴对称模型

Hussain 和 Dincer 对有限圆柱建立了二维温湿传递模型，同时计算轴向和径向变化，并用显式有限差分求解 [5]。其模型说明，当圆柱端面传递或轴向梯度不可忽略时，仅做径向一维近似会丢失端部效应。

**与 A 题的关系：** 题中长径比为 \(L/(2R)=6.25\)，可以把中段近似为径向一维，但“长圆柱”并非无限长。二维模型适合作为少量工况的验证模型：比较几何中心及中截面的温湿结果，量化忽略端面的影响。若题目未给端面与侧面的不同边界参数，直接把二维模型作为主模型会增加计算量，却不会增加可识别信息。

### 2.4 Luikov 交叉耦合模型

Ribeiro、Cotta 和 Mikhailov 对 Luikov 毛细多孔介质热质传递方程给出了积分变换解和基准结果 [6]。Luikov 类模型除常规导热、扩散外，还可包括温度梯度驱动的水分迁移、含湿梯度引起的热流以及相变热效应，概念形式为

\[
T_t=a_T\nabla^2T+a_{TC}\nabla^2C,
\qquad
C_t=a_C\nabla^2C+a_{CT}\nabla^2T.
\]

**与 A 题的关系：** 它比有效 Fick–Fourier 模型具有更强的机理耦合，但需要热湿交叉系数、相变比例或相关热力学参数。题目没有提供这些量，也没有内部测温测湿数据供反演；把这些参数自行设定会使模型不可验证。因此不作为主模型，只在论文“模型改进”中说明。

### 2.5 热–湿–力多物理模型

王会林、卢韬、姜培学基于 Fick 扩散、Fourier 导热和热弹性建立生物多孔介质热–湿–力双向耦合模型，同时预测温度、含水率及应力应变，并考虑部分物性随干基含水率和温度变化 [7]。

**与 A 题的关系：** 该模型可以解释不均匀干燥导致的应力和变形，物理层次高于题目所需。然而本题只给整体半径随时间变化，没有弹性模量、泊松比、孔压、渗透率或局部位移数据，因此无法可靠识别力学子模型。第四问宜把实测半径作为已知移动边界，而不额外求应力场。

### 2.6 经验收缩与孔隙率修正

Aprajeeta、Gopirajah 和 Anandharamakrishnan 在马铃薯干燥模型中同时考虑收缩、孔隙率以及热质传递；其实验显示径向尺寸明显下降，且孔隙率在后期发生变化 [8]。Brasiello 等提出两种茄子收缩扩散模型：一种把结构变化吸收到随含水率变化的扩散系数中，另一种显式加入等效对流项，并在特定条件下证明两种形式等价 [9]。

**与 A 题的关系：** 这些研究说明“改变半径但仍在固定网格上使用原扩散式”缺少坐标运动项或等价修正。题目已经提供 \(R(t)\)，因此应把区域映射到固定坐标 \(\xi=r/R(t)\)，并从干物质和水分守恒推导方程。题目没有孔隙率数据，不能照搬孔隙率修正。

### 2.7 预测型移动边界模型

Adrover、Brasiello 和 Ponso 建立了预测食品体积收缩与表面运动的移动边界模型，把局部收缩速度与水分扩散通量通过收缩函数关联；该函数可从总体体积与含水率的标定曲线估计 [10]。Adrover、Venditti 和 Brasiello 随后把该思想扩展到非等温、连续与间歇干燥，同时加入热传递、表面蒸发和温度相关扩散系数 [11]。

其关键思想可概括为

\[
\boldsymbol v_s=-\alpha(\phi)\boldsymbol J_w,
\]

再由表面速度决定未知边界的位置。

**与 A 题的关系：** 这类模型适合在未知 \(R(t)\) 时由水分通量预测收缩。本题附件2直接给出 \(R(t)\)，如果同时强制使用文献中的收缩速度本构，就可能对同一边界重复施加条件。因此本题基准方案采用“实测 \(R(t)\)+守恒坐标变换”；文献模型作为扩展方案，用来讨论半径数据与含水量是否满足自洽收缩关系。

## 3. 模型形式横向比较

| 模型 | 主要未知量 | 可输出空间分布 | 收缩处理 | 额外参数需求 | 对本题的定位 |
|---|---|---:|---|---|---|
| Modified Page 等薄层模型 [1] | 平均水分比 | 否 | 通常不显式处理 | 经验常数 k、n | 平均曲线辅助校验 |
| 一维 Fick–Fourier [2–4] | T(r,t)、C(r,t) | 是 | 固定域或经验修正 | 题目已基本给出 | 第一至第三问主模型 |
| 二维轴对称 [5] | T(r,z,t)、C(r,z,t) | 是 | 可扩展 | 端面边界及更高算力 | 一维近似验证 |
| Luikov [6] | 温度、含湿势等 | 是 | 通常需另加 | 交叉系数、相变参数 | 理论扩展，不作为主模型 |
| 热–湿–力耦合 [7] | T、C、位移、应力 | 是 | 力学方程预测 | 大量力学与孔隙参数 | 数据不足，不采用 |
| 经验收缩/等效对流 [8–9] | T、C、尺寸或孔隙率 | 是 | 数据拟合或等效项 | 收缩、孔隙关系 | 支撑第四问坐标修正 |
| 预测型移动边界 [10–11] | T、C、未知边界 | 是 | 水分通量驱动 | 收缩函数及相平衡 | 第四问扩展与自洽性检查 |

## 4. 对现有计划的具体修正

### 4.1 第一至第三问

采用一维径向 Fick–Fourier 有效扩散模型，空间上使用守恒型有限体积，时间上使用隐式方法。该选择与圆柱一维研究 [2–4] 一致，也能直接生成题目要求的径向分布。

耦合的表述应写成“温湿相关物性的非线性耦合”。题目第二、三问中 \(T\) 通过 \(D(C,T)\) 影响水分输运，\(C\) 通过 \(\rho(C),c_p(C),k(C)\) 影响传热。由于未给交叉热湿系数，不采用完整 Luikov 模型 [6]。

热边界的基准形式为对流换热。水分边界暂按题设有效浓度差写第三类边界，但要声明：气相水分浓度与固相干基含水率虽都写成 kg/kg，其参考质量不同。更严格的文献模型会通过平衡分配关系或吸附等温线连接气固两相 [11]。本题缺少此关系，不能自行引入未经标定的分配系数。

### 4.2 一维近似的验证

保留一维径向模型作为主模型。二维轴对称抽样验证属于后续改进建议，本仓库尚未完成该项计算，不能写成已经验证了一维近似。后续可在明确端面边界假设后，比较若干代表时刻的几何中心和中截面径向结果 [5]。当前网格收敛只检验一维方程的数值解，不能量化被忽略的轴向效应。

### 4.3 第四问

附件2的 \(R(t)\) 是观测输入，基准模型不再建立另一条独立的半径预测方程。假定均匀径向收缩、长度不变，固体骨架速度取

\[
v_r(r,t)=\frac{R'(t)}{R(t)}r.
\]

若 \(C\) 是随材料运动的干基含水率，并额外假定当前干物质量体积密度在空间上均匀，则在物理坐标下写材料导数：

\[
\frac{\partial C}{\partial t}+v_r\frac{\partial C}{\partial r}
=\frac{1}{r}\frac{\partial}{\partial r}
\left(rD\frac{\partial C}{\partial r}\right).
\]

映射到 \(\xi=r/R(t)\) 后，网格运动项与均匀骨架对流相消，扩散项带 \(R(t)^{-2}\) 尺度。该推导与收缩文献中对等效对流项和移动边界的处理相符 [9–11]。若改用单位当前体积水分浓度，方程会出现不同的体积变化项；论文必须始终使用同一个水分变量定义。

### 4.4 终点与验证

薄层经验模型只校验平均曲线，不能代替全域阈值。第三、四问仍以 \(\max C(r,t)<0.15\) 判断完成。固定域和移动域都要做守恒检查、网格加密，并比较固定半径、实测收缩半径两种结果。文献中更复杂模型的拟合良好不构成本题精度证明，本题结论只能由附件输入、数值收敛和敏感性结果支撑。

## 5. 推荐论文中的文献综述表述

可在“模型准备”或“模型假设”前加入以下内容，正式写作时根据全文风格调整：

> 热风干燥模型可分为薄层经验模型和空间分布模型。针对西洋参分段干燥，Modified Page 模型能够描述平均水分比随阶段变化的规律，但不提供物料内部的温度和含水率分布 [1]。为获得圆柱物料的空间场，已有研究通常采用 Fourier 导热与 Fick 有效扩散方程，并在表面设置对流换热、传质边界，以有限差分或有限体积方法求解 [2–5]。Luikov 模型及热–湿–力模型能够进一步描述交叉输运和变形 [6–7]，但需要题目未提供的交叉系数和力学参数。考虑到本题给出了温湿相关物性及半径观测，本文采用一维径向 Fick–Fourier 模型；前三问使用固定圆柱域，第四问将实测半径引入归一化移动坐标，并通过网格收敛检验数值稳定性、通过情景敏感性检查输入解释的影响。二维轴向效应尚未量化，保留为模型局限。该处理吸收了收缩扩散及移动边界研究的基本思想 [8–11]，同时避免引入无法由附件识别的参数。

这段文字只说明模型选择依据。参赛队还需要结合最终计算结果人工修改并逐项核对，不应在未完成验证前写“结果准确”或“模型优于其他模型”。

## 6. 参考文献

[1] 吴小华, 马渊博, 宁旭丹, 等. 西洋参分段式热风干燥动力学模型构建[J]. 农业工程学报, 2020, 36(5): 318-324. DOI: [10.11975/j.issn.1002-6819.2020.05.037](https://doi.org/10.11975/j.issn.1002-6819.2020.05.037).

[2] SILVA W P, SILVA C M D P S, GAMA F J A. Estimation of thermo-physical properties of products with cylindrical shape during drying: The coupling between mass and heat[J]. Journal of Food Engineering, 2014, 141: 65-73. DOI: [10.1016/j.jfoodeng.2014.05.010](https://doi.org/10.1016/j.jfoodeng.2014.05.010).

[3] TZEMPELIKOS D A, MITRAKOS D, VOUROS A P, et al. Numerical modeling of heat and mass transfer during convective drying of cylindrical quince slices[J]. Journal of Food Engineering, 2015, 156: 10-21. DOI: [10.1016/j.jfoodeng.2015.01.017](https://doi.org/10.1016/j.jfoodeng.2015.01.017).

[4] KAYA A, AYDIN O, DINCER I. Numerical modeling of forced-convection drying of cylindrical moist objects[J]. Numerical Heat Transfer, Part A: Applications, 2007, 51(9): 843-854. DOI: [10.1080/10407780601112753](https://doi.org/10.1080/10407780601112753).

[5] HUSSAIN M M, DINCER I. Two-dimensional heat and moisture transfer analysis of a cylindrical moist object subjected to drying: A finite-difference approach[J]. International Journal of Heat and Mass Transfer, 2003, 46(21): 4033-4039. DOI: [10.1016/S0017-9310(03)00229-1](https://doi.org/10.1016/S0017-9310(03)00229-1).

[6] RIBEIRO J W, COTTA R M, MIKHAILOV M D. Integral transform solution of Luikov's equations for heat and mass transfer in capillary porous media[J]. International Journal of Heat and Mass Transfer, 1993, 36(18): 4467-4475. DOI: [10.1016/0017-9310(93)90131-O](https://doi.org/10.1016/0017-9310(93)90131-O).

[7] 王会林, 卢韬, 姜培学. 生物多孔介质热风干燥数学模型及数值模拟[J]. 农业工程学报, 2014, 30(20): 325-333. DOI: [10.3969/j.issn.1002-6819.2014.20.039](https://doi.org/10.3969/j.issn.1002-6819.2014.20.039).

[8] APRAJEETA J, GOPIRAJAH R, ANANDHARAMAKRISHNAN C. Shrinkage and porosity effects on heat and mass transfer during potato drying[J]. Journal of Food Engineering, 2015, 144: 119-128. DOI: [10.1016/j.jfoodeng.2014.08.004](https://doi.org/10.1016/j.jfoodeng.2014.08.004).

[9] BRASIELLO A, ADILETTA G, RUSSO P, et al. Mathematical modeling of eggplant drying: Shrinkage effect[J]. Journal of Food Engineering, 2013, 114(1): 99-105. DOI: [10.1016/j.jfoodeng.2012.07.031](https://doi.org/10.1016/j.jfoodeng.2012.07.031).

[10] ADROVER A, BRASIELLO A, PONSO G. A moving boundary model for food isothermal drying and shrinkage: General setting[J]. Journal of Food Engineering, 2019, 244: 178-191. DOI: [10.1016/j.jfoodeng.2018.09.018](https://doi.org/10.1016/j.jfoodeng.2018.09.018).

[11] ADROVER A, VENDITTI C, BRASIELLO A. A non-isothermal moving-boundary model for continuous and intermittent drying of pears[J]. Foods, 2020, 9(11): 1577. DOI: [10.3390/foods9111577](https://doi.org/10.3390/foods9111577).

## 7. 查阅入口

- 西洋参分段干燥论文：[农业工程学报全文页](https://www.aeeisp.com/nygcxb/article/doi/10.11975/j.issn.1002-6819.2020.05.037)
- 圆柱一维变物性模型：[Journal of Food Engineering 文章页](https://www.sciencedirect.com/science/article/pii/S0260877414002118)
- 圆柱榅桲片热质传递模型：[Journal of Food Engineering 文章页](https://www.sciencedirect.com/science/article/pii/S0260877415000357)
- 有限圆柱二维模型：[论文 PDF](https://electronicsandbooks.com/edt/manual/Magazine/I/International%20Journal%20of%20Heat%20and%20Mass%20Transfer/2003%20Volume%2046/21/09.pdf)
- 生物多孔介质热–湿–力模型：[农业工程学报全文页](https://www.aeeisp.com/nygcxb/en/article/doi/10.3969/j.issn.1002-6819.2014.20.039)
- 非等温移动边界模型：[Foods 全文页](https://www.mdpi.com/2304-8158/9/11/1577)
