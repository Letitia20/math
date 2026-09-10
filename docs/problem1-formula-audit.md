# 第一问公式复核记录

复核日期：2026-09-10。

## 题面与输入

- 两份 `A题.pdf` 的 SHA-256 完全相同，排除了下载版与临时版不一致。
- 逐页渲染并人工核对题面后，第一问水分扩散系数确定为
  \[
  D(C)=7\times10^{-9}\exp\!\left(-\frac{0.89}{C}\right)\;\mathrm{m^2/s}.
  \]
  其中指数前有负号，分母为局部水分浓度 \(C\)。
- 其余常数复核为：\(\rho=820\ \mathrm{kg/m^3}\)、\(c_p=2600\ \mathrm{J/(kg\,K)}\)、
  \(k=0.36\ \mathrm{W/(m\,K)}\)、\(h_T=25\ \mathrm{W/(m^2\,K)}\)、
  \(h_m=8\times10^{-7}\ \mathrm{m/s}\)。圆柱半径 \(R=0.02\ \mathrm m\)，初始温度
  \(28\ ^\circ\mathrm C\)，初始水分浓度 \(2.55\ \mathrm{kg/kg}\)。

## 方程与边界符号

采用轴对称、仅沿半径变化的一维圆柱模型：

\[
\rho c_p\frac{\partial T}{\partial t}
=\frac1r\frac{\partial}{\partial r}\left(kr\frac{\partial T}{\partial r}\right),
\qquad
\frac{\partial C}{\partial t}
=\frac1r\frac{\partial}{\partial r}\left(rD(C)\frac{\partial C}{\partial r}\right).
\]

中心对称条件为 \(T_r(0,t)=C_r(0,t)=0\)。外表面采用向外法向，故 Robin 条件为

\[
-kT_r(R,t)=h_T[T(R,t)-T_\infty(t)],\qquad
-D(C)C_r(R,t)=h_m[C(R,t)-C_\infty(t)].
\]

边界条件符号又通过“当环境温度更高时表面升温、当环境水分浓度更低时表面失水”的方向性测试复核。

## 数值实现复核

- 空间采用守恒的节点型环形有限体积法，界面扩散系数取调和平均。
- 时间积分采用刚性 BDF 法；水分方程以 \(\log C\) 为状态量，避免非物理负浓度。
- 附件 1 每 60 s 给出环境历史，计算时作分段线性插值。
- 640 与 1280 个半径区间在指定抽查点上的最大差异为：温度
  \(2.00\times10^{-6}\ ^\circ\mathrm C\)，水分浓度
  \(2.18\times10^{-5}\ \mathrm{kg/kg}\)。最终结果使用 1280 个半径区间。
