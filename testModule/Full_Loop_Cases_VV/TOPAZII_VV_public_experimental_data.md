# TOPAZ-II 热工水力、TFE 电输出与保温特性 V&V 公开实验文献调研

> 用途：为基于 Python 的 TOPAZ-II 系统分析程序提供公开 V&V 数据源。覆盖系统热工水力、TFE 电输出、单 TFE 热阻/保温、NaK 回路、管道/辐射器以及可作为背景的燃料/TFE 设计与试验历史。
> 版本：调研初稿，2026-07-01。
> 说明：本文优先收集“公开文献或公开报告中能追溯到实验工况、结构几何和可数字化曲线”的资料。TOPAZ-II 公开数据中，系统级数据以 V-71/TSET 电加热试验最可用；单 TFE 数据以 NPS Benke/Venable 学位论文和 AIP 会议论文最可用。公开点堆中子动力学实验数据未形成可直接 V&V 的完整开放数据集，建议在本论文中实事求是地将点堆动力学部分作为模型实现/设计计算，或另找临界/反应性扰动实验验证。

---

## 0. 结论性建议

### 0.1 建议优先用于 V&V 的数据组合

| 优先级 | 数据源                                                          | 适合验证的模型                                   | 数据可用性 | 备注                                                              |
| ---:| ------------------------------------------------------------ | ----------------------------------------- | ----- | --------------------------------------------------------------- |
| A   | Paramonov & El-Genk, V-71 unit tests                         | 系统热工水力、NaK 回路、体积补偿器、EM 泵-流量/压降、管道/辐射器温度分布 | 高     | 有系统几何、测量量、稳态功率扫描、压力/温度曲线、局部温度/欠热度、180 h 电输出趋势；但为非核电加热试验。        |
| A   | Benke, TOPAZ-II single-cell TFE test stand                   | 单 TFE 结构、保温/热阻、调压 He gap、冷却水热平衡、集电极套筒温度   | 高     | 有详细结构、材料、几何、热阻网络、TISA 功率修正、He gap 有效导热系数。                       |
| A   | Venable, TOPAZ-II single-cell TFE electrical characteristics | 单 TFE I-V 曲线、铯压优化、热输入-电输出效率               | 高     | 有 1.0-3.6 kWt 功率扫描、0.1-1.5 torr 铯压扫描、I-V 图、最优铯压范围和经验关系。         |
| B   | Luchau et al., TFE No.24 long-duration tests                 | TFE 长时稳定性、异常/边界工况趋势                       | 中     | 摘要级信息公开，部分图在 PDF 中；适合作长期稳定性背景和趋势校核。                             |
| B   | Thermal Power Tests of Single Cell TFEs and Systems          | 单 TFE 与系统热功率试验关联                          | 中     | 有 optimum power、station-keeping、high-output 操作背景；详细数值需读 PDF 图表。 |
| C   | Russian TOPAZ II System Test Program (1970-1989)             | 试验历史和资格鉴定背景                               | 中-低   | 有系统试验类别和约 28 套系统测试概述；详细核试验/动力学数据公开不足。                           |
| C   | UO2 fuel design/performance papers                           | TFE/燃料几何和材料背景                             | 中     | 适合说明燃料/TFE 结构；不是热工水力实验基准。                                       |
| C   | Wang et al. modified RELAP5                                  | 程序间对比、现代建模参考                              | 中     | 不是原始实验，但复现 V-71，给出温度/压力误差水平，可作论文对标。                             |

### 0.2 对模型验证边界的建议写法

1. **热工水力模块**：用 V-71/TSET 整机电加热数据进行系统级验证，包括堆芯进出口温度、系统压力、回路局部温度、NaK 欠热度、EM 泵流量趋势。
2. **TFE 电输出模块**：用 Venable 单 TFE 试验的 I-V 曲线、最优铯压、最大输出功率经验关系验证；系统级总电输出只作为弱验证，因为 V-71 的 TFE 铯压并非直接测量。
3. **保温/隔热模块**：用 Benke 单 TFE 试验台的调压 He gap、非调压 He gap、Al2O3 绝缘层、冷却水夹套和热电偶温度进行验证。
4. **管道/辐射器模块**：用 V-71 的 78 根辐射器管、上下集管、回路管路几何以及局部温度/压力图校核；若需要独立辐射器实验，公开资料不足，建议说明为“基于系统试验的子模型验证”。
5. **点堆中子动力学模块**：V-71、Benke、Venable 均为非核电加热试验，不能直接验证中子动力学。可以在论文中只验证热工水力/热电转换耦合，点堆动力学用公开反应性反馈系数、设计值或独立临界实验做一致性检查。

---

## 1. TOPAZ-II 系统设计背景参数

### 1.1 系统总体设计参数

主要来源：[S1] Voss, *TOPAZ II System Description*, LA-UR-94-4, 1994。

| 参数           | 数值/说明                         | 用途        |
| ------------ | -----------------------------:| --------- |
| 设计寿期         | 3 yr                          | 稳态/寿期分析背景 |
| TFE 端电功率     | 约 6 kWe                       | 电输出模块标称值  |
| 电压           | 27 V                          | 负载/母线边界   |
| 初寿期热功率       | 115 kWth                      | 额定热工况     |
| 最大热功率        | 135 kWth                      | 高功率包络     |
| 冷却剂          | NaK-78，约 22 wt% Na + 78 wt% K | 物性库设置     |
| 额定冷却剂质量流量    | 约 1.3 kg/s                    | 系统流量边界/验证 |
| 堆芯冷却剂入口/出口温度 | BOL 约 743/843 K；最大约 773/873 K | 堆芯热工水力设计值 |
| TFE 数量       | 37 根                          | 堆芯通道数     |
| 负载/泵供电分配     | 34 根 TFE 给外部负载；3 根 TFE 给 EM 泵 | 电路边界      |
| 堆芯高度         | 约 375 mm                      | 释热高度/轴向建模 |
| 堆芯直径         | 约 260 mm                      | 堆芯几何      |
| 一次系统材料       | 不锈钢                           | 管路/结构材料边界 |
| 燃料           | UO2，高富集                       | 燃料温度/物性背景 |
| 有效辐射器面积      | 约 7.2 m²                      | 辐射器边界     |
| 辐射器管数        | 78 根                          | 管道/辐射器模型  |
| 辐射器翅片        | 铜翅片，黑色搪瓷涂层                    | 辐射换热边界    |

### 1.2 堆芯和 TFE 结构背景

Voss 系统说明和 Hoth 燃料设计论文给出如下结构信息，可用于构造 TFE 的径向热阻和轴向释热模型：

- TOPAZ-II 堆芯含 37 根单电池 TFE。
- 燃料为环形 UO2 芯块。Hoth 等给出燃料芯块直径约 17 mm、高度约 9 mm、密度约 96% 理论密度，每根 TFE 约 40 个芯块；Voss 系统报告中还提到芯块中心孔直径随径向位置变化，约 4.5、6.0 或 8.0 mm。
- 发射极管为单晶 Mo-3%Nb，外表面有钨涂层以增强热离子发射。
- 集电极为多晶 Mo，外侧有 Al2O3 绝缘层。
- 发射极-集电极间隙为铯蒸气工作间隙，由 Sc2O3 绝缘/定位件保持。
- 集电极绝缘层外侧与冷却管之间存在 He gap，可作为热阻调节/保温特性模拟对象。

> 注：不同公开资料对燃料芯块高度存在 8-9 mm 量级表述差异。建议在论文中按引用来源说明，并在敏感性分析中处理。

---

## 2. 系统级 V&V：V-71/TSET 电加热整机试验

### 2.1 公开文献来源

核心来源：[S2] Paramonov & El-Genk, *Comparison of a TOPAZ-II Model with Experimental Data from the V-71 Unit Tests* / *Development and Comparison of a TOPAZ-II System Model with Experimental Data*, 1994；辅助来源：[S1] Voss 1994；现代程序对比：[S11] Wang et al. 2019。

### 2.2 实验介绍

V-71 是 TOPAZ-II 整机电加热试验单元，在美国 Albuquerque 的 Thermionic System Evaluation Test Facility (TSET) 进行。公开资料说明其试验时间覆盖 1992 年 11 月至 1993 年 5 月。该试验不是核加热，而是在每根 TFE 中用钨电阻加热器模拟堆芯释热。

关键实验特征：

- 堆芯中 37 根 TFE 均电加热。
- TFE 发射极有效长度约 0.375 m，其中中间 0.300 m 被钨加热器近似均匀加热。
- 34 根 TFE 接外部负载；3 根 TFE 与外部电源共同用于 EM 泵供电。
- 外部电源用于维持 EM 泵端电压，避免泵电压随 TFE 输出波动过大。
- 真空室壁面水冷，实验边界与空间真实辐射边界不同；因此辐射器验证应写成“地面 TSET 边界下的回路热工水力验证”。

### 2.3 实验段/系统结构参数

V-71 一次冷却回路几何参数如下，可直接转为一维系统程序的管段/控制体输入。

| 部件           | 数量  | 长度/mm | 内径或等效通道/mm | 备注                  |
| ------------ | ---:| -----:| ----------:| ------------------- |
| 堆芯 TFE 冷却流道  | 37  | 500   | 环隙约 0.7    | 堆芯流道，NaK 沿 TFE 外侧流动 |
| 堆芯出口至辐射器上集管  | 2   | 2500  | 30         | 并联管段                |
| 辐射器上集管       | 1   | 824   | 20         | 分配至 78 根管           |
| 辐射器管         | 78  | 1850  | 7          | 管道/辐射器模型核心对象        |
| 辐射器下集管       | 1   | 1346  | 20         | 汇流                  |
| 辐射器下集管至 EM 泵 | 2   | 3500  | 30         | 并联管段                |
| EM 泵至堆芯入口    | 6   | 1000  | 18         | 并联回流至堆芯入口           |

### 2.4 体积补偿器/系统压力边界

Paramonov & El-Genk 给出的体积补偿器和 NaK 装量信息适合用于系统压力验证。

| 参数              | 数值                                  |
| --------------- | -----------------------------------:|
| 初始气体温度          | 293 K                               |
| 初始气体压力          | 29.4 kPa                            |
| 初始气体体积          | \(8.5\times 10^{-3}\ \mathrm{m^3}\) |
| 补偿器波纹管 NaK 体积   | 0.4 L                               |
| 堆芯 NaK 体积       | 3.6 L                               |
| 进/回流管 NaK 体积    | 7.8 L                               |
| 辐射器集管和管束 NaK 体积 | 7.0 L                               |
| 总 NaK 体积        | 18.8 L                              |

可采用的压力边界关系：

$$
P_g = P_{g,0}\frac{V_{g,0}}{V_g}\frac{T_g}{T_{g,0}}
$$

其中 \(P_{g,0}\)、\(V_{g,0}\)、\(T_{g,0}\) 为初始气体压力、体积和温度。

### 2.5 典型工况和测量量

#### 功率扫描

- TFE 电加热总功率：从低功率到约 100-120 kW。
- 约 40 kW 附近接入负载并提高 EM 泵电压，泵电流从约 300 A 增加至约 700 A，导致堆芯冷却剂温度约下降 25 K。
- 约 106 kW 工况给出了局部冷却剂温度、局部饱和温度、局部压力和欠热度。
- 1993 年 5 月试验中，在约 104-120 kW 功率水平运行约 180 h，并给出负载电压、电流和电功率随时间变化。

#### 测量参数

| 测量量          | 用途                   |
| ------------ | -------------------- |
| 堆芯入口冷却剂温度    | 堆芯能量方程验证             |
| 堆芯出口冷却剂温度    | 堆芯能量方程验证             |
| EM 泵端电压      | 泵-流量/压降边界            |
| 铯储罐温度        | TFE 电输出边界；注意铯压不是直接测量 |
| 体积补偿器处冷却剂压力  | 系统压力验证               |
| 负载电压、电流、功率   | TFE 系统级输出弱验证         |
| 局部冷却剂温度/饱和温度 | 欠热度和安全裕量验证           |

### 2.6 可用于程序 V&V 的基本公式

堆芯能量守恒：

$$
Q_{\mathrm{core}} \simeq \dot m_{\mathrm{NaK}} c_{p,\mathrm{NaK}}
\left(T_{\mathrm{out}}-T_{\mathrm{in}}\right)
$$

堆芯温升：

$$
\Delta T_{\mathrm{core}} = T_{\mathrm{out}}-T_{\mathrm{in}}
$$

回路压降：

$$
\Delta p_{\mathrm{loop}}
=
\sum_j
\left(
 f_j \frac{L_j}{D_j}
 + \sum_i K_{i,j}
\right)
\frac{\rho_j u_j^2}{2}
$$

电磁泵和回路平衡：

$$
\Delta p_{\mathrm{EM}}(I_{\mathrm{pump}},T_{\mathrm{NaK}})
=
\Delta p_{\mathrm{loop}}(\dot m,T_{\mathrm{NaK}})
$$

欠热度：

$$
\Delta T_{\mathrm{sub}}(x)=T_{\mathrm{sat}}\left[p(x)\right]-T_{\mathrm{NaK}}(x)
$$

### 2.7 文献给出的误差水平

Paramonov & El-Genk 报告的模型-实验对比中，冷却剂温度计算值与实测值约在 15 K 内，系统压力约在 12% 内。Wang 等在 modified RELAP5 复现 V-71 时，报告冷却剂温度和系统压力最大相对误差分别约为 8% 和 10%。这些数值可作为你程序初版 V&V 的合理参考，但不宜直接写成绝对验收标准。

### 2.8 重要不确定度

- V-71 为电加热非核试验，不能验证裂变功率分布和点堆动力学。
- V-71 中 TFE 铯压并非直接测量，而是由铯节流阀开度标定曲线推算。因此系统级电输出误差可能包含铯压边界误差。
- 地面 TSET 水冷真空室边界不同于空间深冷辐射边界；若用于空间辐射器模型验证，应明确采用 TSET 边界。

---

## 3. 单 TFE 热阻/保温 V&V：Benke 单电池 TFE 试验台

### 3.1 公开文献来源

核心来源：[S3] Benke, *Operational Testing and Thermal Modeling of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand*, Naval Postgraduate School thesis, 1994。辅助来源：[S5] Benke & Venable AIP 论文。

### 3.2 实验介绍

该试验台使用未装燃料的 TOPAZ-II 单电池 TFE，并用 TISA 钨电加热器替代核燃料释热。试验台最初在俄罗斯 LUCH 研究所构建和测试，后来在 New Mexico Engineering Research Institute (NMERI) 重新组装，用于 TOPAZ International Program 的非核 TFE 地面实验。

适合验证的内容：

- 单 TFE 径向热阻网络；
- 发射极、铯间隙、集电极、Al2O3 绝缘层、He gap、冷却水夹套之间的传热；
- 调压 He gap 的低压稀薄气体导热/保温特性；
- 水冷夹套热平衡；
- 集电极套筒热电偶温度分布。

### 3.3 试验段结构设计

| 部件         | 材料/结构           | 几何尺寸/参数                                | 备注                    |
| ---------- | --------------- | --------------------------------------:| --------------------- |
| TISA 电加热器  | 钨电阻加热器          | 最大 29 VAC，170 A；最大约 4500 W；长度 300.0 mm | 模拟工作段释热               |
| TISA 内导体   | 钨               | 直径约 6.5 mm                             |                       |
| TISA 外导体   | 钨               | 外径约 7.0 mm，厚度约 0.4 mm                  |                       |
| 热离子工作段     | 同轴圆筒            | 长度 375.0 mm                            | 与 TOPAZ-II TFE 有效长度对应 |
| 发射极        | Mo-3%Nb 单晶，外覆 W | 厚度约 1.15 mm，外径约 19.6 mm，W 涂层约 0.1 mm   | 热离子发射表面               |
| 铯蒸气间隙      | Cs vapor        | 约 0.5 mm                               | 发射极-集电极间隙             |
| 定位/绝缘件     | Sc2O3           | 6 组 spacer                             | 保持间隙                  |
| 集电极        | 多晶 Mo           | 厚度约 1.4 mm，内径约 20.6 mm                 | 电子收集表面                |
| 集电极外绝缘层    | Al2O3           | 厚度约 0.15 mm                            | 电绝缘/热阻                |
| 非调压 He gap | He              | 约 0.05 mm；通常 200-300 torr              | 集电极绝缘层与套筒之间           |
| 集电极套筒      | 俄制 1X18H10T 不锈钢 | 厚度约 3.0 mm，外径约 29.9 mm                 | 对应真实 NaK 冷却边界位置       |
| 热电偶槽       | 套筒外侧槽           | 12 个槽，深约 2.0 mm                        | 测量轴向温度                |
| 调压 He gap  | He              | 约 0.5 mm；正常运行 1-10 torr                | 保温/热阻控制核心对象           |
| 水冷夹套       | 1X18H10T 不锈钢    | 内壁 2.5 mm，外壁 1.0 mm                    | 代替 NaK 冷却             |
| 冷却水螺旋流道    | 不锈钢管/夹套         | 螺距约 35 mm，约 6.5 圈                      | 水自下而上流动               |

### 3.4 功率定义与热平衡

TISA 输入功率：

$$
P_{\mathrm{TISA}} = V_h I_h
$$

有效工作段热输入。Benke 给出经验修正，约 88% 的 TISA 输入功率进入 active zone：

$$
P_{\mathrm{az}} = 0.88 P_{\mathrm{TISA}}
$$

冷却水热平衡：

$$
Q_w = \dot m_w c_{p,w}
\left(T_{w,\mathrm{out}}-T_{w,\mathrm{in}}\right)
$$

圆筒径向导热：

$$
q' = \frac{2\pi k L\left(T_1-T_2\right)}{
\ln\left(r_2/r_1\right)}
$$

同轴圆筒辐射换热：

$$
q'_{\mathrm{rad}}
=
\frac{2\pi r_1 L\sigma\left(T_1^4-T_2^4\right)}
{\frac{1}{\varepsilon_1}+\frac{r_1}{r_2}\left(\frac{1}{\varepsilon_2}-1\right)}
$$

总热阻网络可写作：

$$
R_{\mathrm{tot}} =
R_{\mathrm{Cs}}+R_{\mathrm{collector}}+R_{\mathrm{Al_2O_3}}
+R_{\mathrm{He,unreg}}+R_{\mathrm{sleeve}}
+R_{\mathrm{He,reg}}+R_{\mathrm{water}}
$$

### 3.5 典型保温/热阻工况

Benke 报告中适合直接用于保温特性验证的典型工况：

| 条件                | 数值/说明                 |
| ----------------- | ---------------------:|
| 有效 active zone 功率 | 约 3003 W              |
| 调压 He gap 压力      | 10 torr               |
| 调压 He gap 宽度      | 0.5 mm                |
| 反推有效导热系数          | 约 0.073-0.087 W/(m·K) |
| 连续介质估算 He 导热系数    | 约 0.276 W/(m·K)       |
| 冷却水流动状态           | Re 约 1480，过渡/层流附近     |
| 水侧换热系数            | 约 528-1012 W/(m²·K)   |

> 对保温模型特别重要：低压 He gap 的实际等效导热系数明显低于常压/连续介质估算值，应考虑稀薄气体效应、温度跳跃、接触/装配不确定性和辐射换热并联。

### 3.6 可用于验证的输出量

| 输出量                   | 推荐比较方式             |
| --------------------- | ------------------ |
| 集电极套筒 12 点轴向温度        | 曲线比较；误差可用 MAE/RMSE |
| 冷却水出口温度               | 热平衡误差              |
| 调压 He gap 等效导热系数      | 由实验反推值与模型值比较       |
| active-zone 输入功率与水侧吸热 | 能量闭合误差             |
| 轴向温度分布                | 图像数字化后比较趋势和峰值位置    |

---

## 4. 单 TFE 电输出 V&V：Venable 电特性试验

### 4.1 公开文献来源

核心来源：[S4] Venable, *Electrical Characteristics and Thermal Analysis of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand*, Naval Postgraduate School thesis, 1995。

### 4.2 实验介绍

Venable 使用与 Benke 相同或高度相近的 TOPAZ-II 单电池 TFE 试验台，重点测量 TFE 的电输出特性。实验采用 TISA 加热器提供热输入，通过改变铯储罐温度/铯压以及外部负载，获得 I-V 曲线、最大电功率和转换效率。

### 4.3 典型工况设置

| 参数          | 数值/说明                                                     |
| ----------- | ---------------------------------------------------------:|
| TISA 输入功率范围 | 1.0-3.6 kWt                                               |
| 功率步长        | 约 200 W                                                   |
| 对应全堆平均功率    | 115 kWth / 37 ≈ 3.11 kWt/TFE；135 kWth / 37 ≈ 3.65 kWt/TFE |
| 铯压力扫描       | 约 0.1-1.5 torr；0.4 torr 以下曲线可用性差                          |
| 每个铯压点       | 做 I-V 上扫/负载扫描                                             |
| 典型图         | 3200 Wt 热输入下，\(P_{Cs}=0.4,1.0,1.3\) torr 的 I-V 曲线         |

### 4.4 最优铯压范围

Venable 按有效热输入 \(Q_{az}\) 将最优铯压分为四类：

| 有效热输入 \(Q_{az}\) | 最优铯压     |
| ----------------:| --------:|
| 892-1405 Wt      | 0.4 torr |
| 1580-2112 Wt     | 0.5 torr |
| 2281-2637 Wt     | 0.8 torr |
| 2813-3162 Wt     | 1.0 torr |

经验优化线为：

$$
P_{\mathrm{out,max}}\ [\mathrm{W}]
=260\ P_{Cs}\ [\mathrm{torr}] - 68
$$

该式只建议在 Venable 实验范围内用于趋势验证，不建议外推。

### 4.5 电输出和效率定义

TFE 输出功率：

$$
P_{\mathrm{out}} = V_{\mathrm{TFE}} I_{\mathrm{TFE}}
$$

热离子转换效率：

$$
\eta_{\mathrm{TFE}} = \frac{P_{\mathrm{out}}}{Q_{az}}
$$

I-V 曲线可比较：

$$
I_{\mathrm{calc}}(V; Q_{az}, P_{Cs})
\quad \text{vs.}\quad
I_{\mathrm{exp}}(V; Q_{az}, P_{Cs})
$$

建议误差指标：

$$
\mathrm{RMSE}_I
=
\sqrt{\frac{1}{N}\sum_{i=1}^{N}
\left[I_{\mathrm{calc}}(V_i)-I_{\mathrm{exp}}(V_i)\right]^2}
$$

$$
\varepsilon_{P,\max}
=
\frac{P_{\mathrm{out,max,calc}}-P_{\mathrm{out,max,exp}}}
{P_{\mathrm{out,max,exp}}}
$$

### 4.5.1 Venable Table 7-1：最大输出功率与效率数据

Venable 论文 Table 7-1 给出的最大输出功率和效率可直接作为单 TFE 电输出模型的数值校核点。这里的 active-zone power 已采用论文中的有效热输入口径：

| Active-zone power / Wt | Maximum output power / W | Efficiency / % |
| ----------------------:| ------------------------:| --------------:|
| 892                    | 10.23                    | 1.15           |
| 1062                   | 17.80                    | 1.68           |
| 1237                   | 30.13                    | 2.44           |
| 1405                   | 45.00                    | 3.20           |
| 1580                   | 63.25                    | 4.01           |
| 1755                   | 77.28                    | 4.40           |
| 1933                   | 86.26                    | 4.46           |
| 2112                   | 103.97                   | 4.92           |
| 2281                   | 115.44                   | 5.06           |
| 2474                   | 129.87                   | 5.25           |
| 2637                   | 146.75                   | 5.57           |
| 2813                   | 167.06                   | 5.94           |
| 2999                   | 178.16                   | 5.94           |
| 3162                   | 192.46                   | 6.09           |

这些点适合先做一维标量验证，即比较最大功率和效率随热输入的单调趋势；随后再用 I-V 曲线做曲线级验证。

### 4.6 物理趋势说明

Venable 结果显示：在给定热输入下，铯压过低会导致电极表面铯吸附不足、发射极功函数未达到最优、空间电荷中和不足；铯压过高则会增加电子碰撞并抑制电子输运。因此存在最优铯压。3200 Wt 工况下，铯压从 0.4 torr 增加到 1.0 torr 时电输出增加，继续升至 1.3 torr 时电输出降低。

---

## 5. TFE 长时运行与边界工况：TFE No.24 与单 TFE 热功率试验

### 5.1 TFE No.24 长时热功率试验

来源：[S7] Luchau et al., *Output Power Characteristics and Performance of TOPAZ II Thermionic Fuel Element No. 24*, AIP Conf. Proc. 361, 1996。

公开摘要说明：TFE No.24 进行了超过 3000 h 的热功率试验，覆盖正常和不利条件，包括：

- 低功率运行；
- 高功率运行；
- 向电极间隙引入空气；
- 集电极温度优化；
- 热模型分析；
- 输出功率特性测量。

用途：

- 可作为 TFE 电输出长期稳定性和异常工况趋势验证；
- 适合写入 V&V 的“扩展验证/补充证据”；
- 若无法获得完整图表，不建议作为定量主基准。

归一化长期稳定性指标：

$$
\Pi_P(t)=\frac{P_{\mathrm{out}}(t)}{P_{\mathrm{out}}(t_0)}
$$

退化率：

$$
D_P = -\frac{1}{P_{\mathrm{out}}}\frac{dP_{\mathrm{out}}}{dt}
$$

### 5.2 单 TFE 与系统热功率试验关联

来源：[S6] Luchau et al., *Thermal Power Tests of Single Cell Thermionic Fuel Elements and Systems*, AIP Conf. Proc. 324, 1995。

公开摘要说明：作为 TOPAZ International Program 的一部分，单电池 TFE 与 TOPAZ-II 系统均进行了非核地面试验，实验关注：

- optimum power data；
- station-keeping mode 低功率运行；
- high output power performance；
- 将单 TFE 试验台结果和 TOPAZ-II 系统试验关联。

用途：

- 可用于说明为什么单 TFE 试验可以作为系统 TFE 电特性模型的验证基础；
- 可作为 Venable/Benke 与 V-71 之间的桥接文献。

---

## 6. 俄方 TOPAZ-II 系统试验历史与公开性边界

来源：[S9] Voss & Rodriguez, *Russian TOPAZ II System Test Program (1970-1989)*, AIP Conf. Proc. 301, 1994。

公开资料说明：俄方 TOPAZ-II 计划在约 20 年内开展了系统测试，约 28 套系统用于试验，试验类别包括：

1. 非核热物理试验；
2. 机械试验，包括静载、动态和冲击/振动；
3. 核地面试验；
4. 冷温试验，用于模拟发射前和发射环境低温条件。

对 V&V 的价值：

- 可作为试验历史和飞行鉴定背景；
- 对论文绪论/模型可信度背景有用；
- 但公开文献主要为概述，未提供足够详细的中子动力学、反应性扰动或瞬态实验数据，因此不建议作为点堆动力学定量验证主数据源。

---

## 7. 管道/辐射器模型可用资料

### 7.1 V-71 几何与系统温度/压力数据

对于你已有的管道辐射器热工水力模型，最直接的实验约束仍然来自 V-71：

- 78 根辐射器管，每根长 1850 mm、内径 7 mm；
- 上集管长 824 mm、内径 20 mm；
- 下集管长 1346 mm、内径 20 mm；
- 两根堆芯出口至上集管管路，长度 2500 mm、内径 30 mm；
- 两根下集管至 EM 泵管路，长度 3500 mm、内径 30 mm；
- 可数字化局部温度、局部饱和温度、局部压力和欠热度图。

### 7.2 辐射器详细模型文献

Paramonov & El-Genk 的 *A detailed thermal-hydraulic model of the TOPAZ-II radiator* 是模型方法参考而非独立实验数据源。公开摘要说明该模型可计算辐射器 78 根管和上下集管内的冷却剂温度、流量和压力分布。

建议用法：

- 将其作为辐射器节点划分和分配/汇流计算方法参考；
- 使用 V-71 的局部温度/压力图作为定量验证；
- 若论文主张“独立辐射器实验验证”，需另找 TSET 辐射器专门试验数据，目前本轮未发现足够公开资料。

---

## 8. 建议的 V&V 矩阵

| 模型模块        | 推荐数据                  | 输入量                                | 输出验证量                    | 建议误差指标                                       |
| ----------- | --------------------- | ---------------------------------- | ------------------------ | -------------------------------------------- |
| 堆芯热工水力      | V-71/TSET             | TFE 总加热功率、NaK 物性、37 流道几何、EM 泵电压/电流 | 堆芯入口/出口温度、温升、局部冷却剂温度     | \(\mathrm{MAE}_T\)、\(\max                    |
| 一次回路压降      | V-71/TSET             | 管路几何、泵特性、NaK 温度                    | 系统压力、泵压升、局部压力            | 相对误差、泵-回路交点误差                                |
| 体积补偿器       | V-71/TSET             | 初始气体体积/压力/温度、NaK 总体积、热膨胀           | 体积补偿器处系统压力               | \(\varepsilon_p\)                            |
| 辐射器/管网      | V-71/TSET             | 78 管几何、集管几何、TSET 边界                | 辐射器进出口温度、局部温度/压力、欠热度     | 曲线 RMSE，局部峰值误差                               |
| 单 TFE 热阻/保温 | Benke                 | TISA 功率、He gap 压力、水冷边界、材料/几何       | 套筒温度、He gap 有效导热系数、水侧热平衡 | \(\mathrm{MAE}_T\)、热平衡闭合误差                   |
| 单 TFE 电输出   | Venable               | \(Q_{az}\)、\(P_{Cs}\)、负载电阻/电压扫描    | I-V 曲线、最大功率、效率、最优铯压      | \(\mathrm{RMSE}_I\)、\(\varepsilon_{P,\max}\) |
| 长时稳定性       | TFE No.24             | 分段热功率、铯压/集电极温度                     | \(P_{out}(t)\)、归一化功率     | 漂移率、趋势一致性                                    |
| 点堆动力学       | 暂无直接 TOPAZ-II 开放实验主数据 | 设计反应性系数、控制鼓模型、功率反馈                 | 功率瞬态、反馈趋势                | 仅能做设计一致性或另找临界/反应性实验                          |

---

## 9. 推荐数字化图像清单

> 详细图像定位表见单独文件 `topaz_ii_vv_figure_locator.md`。这里给出论文 V&V 最建议优先数字化的图。

| 优先级 | 来源                       | 图/内容                               | 数字化用途          |
| ---:| ------------------------ | ---------------------------------- | -------------- |
| A   | Paramonov & El-Genk V-71 | 冷却剂入口/出口温度 vs TFE heater power     | 系统热工水力主验证      |
| A   | Paramonov & El-Genk V-71 | 系统压力 vs TFE heater power           | 压力/体积补偿器验证     |
| A   | Paramonov & El-Genk V-71 | 106 kW 工况局部冷却剂温度/饱和温度              | 局部热工安全裕量验证     |
| A   | Paramonov & El-Genk V-71 | 106 kW 工况局部欠热度                     | NaK 欠热度验证      |
| A   | Venable                  | Fig. 6-1，3200 Wt 下不同铯压 I-V 曲线      | 单 TFE 电输出主验证   |
| A   | Venable                  | Fig. 6-2，最优输出功率 vs 铯压              | 最优铯压/经验线验证     |
| A   | Venable                  | Fig. 6-5 至 Fig. 6-8，不同最优铯压下 I-V 曲线 | 多功率点 I-V 验证    |
| B   | Benke                    | TFE/test stand 横截面图                | 结构建模说明图        |
| B   | Benke                    | 典型轴向温度分布                           | 单 TFE 热阻模型验证   |
| B   | Benke                    | 冷却水壁温 vs regulated He gap 有效导热系数   | 保温/He gap 模型验证 |
| B   | Luchau TFE No.24         | 输出功率随时间或工况变化图                      | 长时稳定性/扩展验证     |

---

## 10. 参考文献与公开入口

[S1] Voss, S. S. **TOPAZ II System Description**. LA-UR-94-4, 1994.
公开入口：OSTI record `https://www.osti.gov/biblio/10120556`；PDF `https://www.osti.gov/servlets/purl/10120556`

[S2] Paramonov, D. V.; El-Genk, M. S. **Comparison of a TOPAZ-II Model with Experimental Data from the V-71 Unit Tests**. *AIP Conference Proceedings* 301, 829-843, 1994. DOI: `10.1063/1.2950277`.
同主题/相近版本题名：**Development and Comparison of a TOPAZ-II System Model with Experimental Data**，常见于 ResearchGate/预印本和 Nuclear Technology 引用中。
公开入口：ResearchGate PDF；AIP DOI 页面。

[S3] Benke, S. M. **Operational Testing and Thermal Modeling of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand**. Naval Postgraduate School thesis, 1994.
公开入口：CORE/DTIC/NPS Calhoun PDF。

[S4] Venable, J. R. **Electrical Characteristics and Thermal Analysis of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand**. Naval Postgraduate School thesis, 1995.
公开入口：NPS Calhoun；Internet Archive/Wikimedia PDF。

[S5] Benke, S. M.; Venable, J. R.; 等. **Operational Testing and Thermal Modeling of a TOPAZ-II Single Cell Thermionic Fuel Element Test Stand**. *AIP Conference Proceedings* 324, 1995. DOI 页面可由 AIP 检索。

[S6] Luchau, D. W.; Bruns, D. R.; Sinkevich, V. G.; 等. **Thermal Power Tests of Single Cell Thermionic Fuel Elements and Systems**. *AIP Conference Proceedings* 324, 189-194, 1995. DOI: `10.1063/1.47229`.

[S7] Luchau, D. W.; Bruns, D. R.; Izhvanov, O.; Androsov, V. **Output Power Characteristics and Performance of TOPAZ II Thermionic Fuel Element No. 24**. *AIP Conference Proceedings* 361, 1163-1168, 1996. DOI: `10.1063/1.50056`.

[S8] Wold, S. K.; Izhvanov, O. L.; Vibivanets, V. I.; Schmidt, G. L. **TOPAZ-II Thermionic Fuel Element Testing**. *AIP Conference Proceedings* 301, 1025-1030, 1994. DOI: `10.1063/1.2950098`.

[S9] Voss, S. S.; Rodriguez, E. A. **Russian TOPAZ II System Test Program (1970-1989)**. *AIP Conference Proceedings* 301, 803-812, 1994. DOI: `10.1063/1.2950273`.

[S10] Hoth, C. W.; Degaltsev, Y.; Gontar, A.; Rakitskaya, E. **Design and Performance of the UO2 Fuel for the TOPAZ-II Reactor**. *AIP Conference Proceedings* 301, 55-61, 1994. DOI: `10.1063/1.2950235`.

[S11] Wang, C.-L.; Liu, T.-C.; Tang, S.-M.; Tian, W.-X.; Qiu, S.-Z.; Su, G.-H. **Thermal-hydraulic analysis of space nuclear reactor TOPAZ-II with modified RELAP5**. *Nuclear Science and Techniques* 30, Article 12, 2019. DOI: `10.1007/s41365-018-0537-3`.

[S12] Paramonov, D. V.; El-Genk, M. S. **A detailed thermal-hydraulic model of the TOPAZ-II radiator**. 公开入口：ResearchGate/AIP 检索。用于辐射器模型方法参考。

---

## 11. Google Scholar 检索建议

为方便后续追踪全文和引用版本，建议在 Google Scholar 中用以下题名精确检索：

```text
"TOPAZ II System Description" Voss LA-UR-94-4
"Comparison of a TOPAZ-II Model with Experimental Data from the V-71 Unit Tests"
"Development and Comparison of a TOPAZ-II System Model with Experimental Data"
"Operational Testing and Thermal Modeling of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand"
"Electrical Characteristics and Thermal Analysis of a TOPAZ-II Single-Cell Thermionic Fuel Element Test Stand"
"Thermal Power Tests of Single Cell Thermionic Fuel Elements and Systems"
"Output Power Characteristics and Performance of TOPAZ II Thermionic Fuel Element No. 24"
"TOPAZ-II Thermionic Fuel Element Testing"
"Russian TOPAZ II System Test Program (1970-1989)"
"Design and Performance of the UO2 Fuel for the TOPAZ-II Reactor"
"Thermal-hydraulic analysis of space nuclear reactor TOPAZ-II with modified RELAP5"
"A detailed thermal-hydraulic model of the TOPAZ-II radiator"
```

---

## 12. 论文写作中可直接采用的表述模板

> 本研究的 V&V 采用分层验证策略。系统级热工水力模型使用 TOPAZ-II V-71 整机电加热试验进行验证，该试验在 TSET 设施中用 37 根 TFE 内部钨电阻加热器模拟堆芯释热，并记录堆芯进出口温度、系统压力、EM 泵电参数以及系统电输出。单 TFE 热阻和保温模型采用 Benke 的 TOPAZ-II 单电池 TFE 试验台数据进行验证，重点比较调压 He gap、冷却水边界和集电极套筒温度。TFE 电输出模型采用 Venable 的单 TFE I-V 曲线、铯压优化和最大输出功率数据进行验证。由于上述主要公开数据均为非核电加热试验，本文不将其用于点堆中子动力学模型的直接验证，而仅用于热工水力和热电转换耦合模型的 V&V。

---

## 13. 本次交付的图像包说明

随本文另附 `TOPAZII_VV_public_figures.zip`，共 50 张 PNG 图片和一个 `00_manifest.csv`。图片包只包含用于科研数字化的公开文献图页裁剪/渲染图，不包含全文 PDF。

图像包目录：

```text
TOPAZII_VV_public_figures/
├── 00_README.md
├── 00_manifest.csv
├── Paramonov_ElGenk_1994_V71_system/
│   ├── V-71 冷却回路 Table 1 几何参数图页
│   ├── Fig. 5 EM pump pressure-loss correction
│   ├── Fig. 8 core inlet/exit coolant temperature vs heater power
│   ├── Fig. 9 system pressure vs heater power
│   ├── Fig. 10 Reynolds number vs heater power
│   ├── Fig. 11 flowrate/pump voltage vs heater power
│   ├── Fig. 12 local coolant/saturation temperature at 106 kW
│   ├── Fig. 13 local subcooling at 106 kW
│   ├── Fig. 14 local pressure at 106 kW
│   └── Fig. 15-17 load voltage/current/power over 180 h
├── Venable_1995_single_cell_TFE_main/
│   ├── TFE test stand 结构图
│   ├── Fig. 6-1 至 Fig. 6-12 主 I-V/功率图
│   ├── Fig. 6-13 相对 TFE 功率分布
│   └── Table 7-1 最大效率表
└── Venable_1995_single_cell_TFE_appendixB_IV/
    └── Appendix B 中 0.4-1.5 torr 铯压下的更多 I-V 上扫曲线
```

Origin 数字化时建议把图片包内文件名、文献图号和数字化后的 CSV 文件名保持一一对应，便于论文 V&V 可追溯性审查。
