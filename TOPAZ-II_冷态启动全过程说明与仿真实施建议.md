# TOPAZ-II 冷态启动全过程说明与仿真实施建议

> 用途：TOPAZ-II 系统冷态启动仿真、控制逻辑开发、程序验证和论文撰写。  
> 说明：公开资料中最完整的启动时序来自 TITAM 仿真。原作者明确指出，该时序主要用于展示模型能力，不一定完全代表俄罗斯实际飞行启动规程。

---

## 1. 冷态启动的范围

本文所称冷态启动，是指 TOPAZ-II 从反应堆深度次临界、热离子燃料元件尚未输出电功率、一次回路处于低温状态，逐步过渡到反应堆临界、功率提升、辐射器投入、TFE 点火并最终达到稳态的全过程。

冷态启动涉及以下耦合子系统：

- 点堆动力学与温度反应性反馈；
- 3 个安全鼓和 9 个控制鼓；
- UO₂ 燃料、发射极、收集极和 ZrH 慢化剂传热；
- NaK-78 一次回路和电磁泵；
- 管翅式辐射器与遮热罩；
- 铯供应系统；
- TFE 热离子发电；
- 电池、泵用 TFE 与负载供电切换。

---

## 2. 两种初始条件口径

### 2.1 TITAM 基准工况

TITAM 启动论文假定反应堆、一次回路和辐射器初始温度均为：

\[
T_0=300\ \mathrm{K}
\]

并将总冷却剂质量流量固定为：

\[
\dot m=1.5\ \mathrm{kg/s}
\]

该配置适合用于复现论文中的功率、温度、鼓角和反应性曲线。

### 2.2 工程启动工况

TOPAZ-II 系统说明指出，发射前通过一次回路交流电加热器将系统预热至约：

\[
T_\mathrm{prelaunch}\approx373\ \mathrm{K}
\]

NaK-78 冻结温度约为 262 K。发射和轨道转移期间，遮热罩降低系统热损失，启动电池间歇驱动电磁泵，以防止辐射器细管和一次回路中的 NaK 冻结。

### 2.3 使用原则

| 配置 | 初始温度 | 流量 | 用途 |
|---|---:|---:|---|
| TITAM 复现 | 300 K | 固定 1.5 kg/s | 与原论文曲线比较 |
| 工程物理启动 | 373 K 或轨道保温终值 | 泵—回路耦合求解 | 系统级真实启动分析 |

300 K 是计算基准假设，373 K 是发射前工程预热值，不能混写为同一工况。

---

## 3. 启动前状态

### 3.1 反应堆

| 参数 | 初始状态 |
|---|---|
| 安全鼓 | 3 个全部内旋，0° |
| 控制鼓 | 9 个全部内旋，0° |
| 冷态反应性 | 约 -6.0 美元 |
| 有效增殖因子 | 约 0.952 |
| 裂变功率 | 中子源水平或极低功率 |
| 温度反馈 | 相对冷态参考温度取零 |

### 3.2 一次回路

- 冷却剂为 NaK-78；
- 遮热罩覆盖辐射器；
- 启动电池驱动或间歇驱动电磁泵；
- 体积补偿器维持静压并吸收热膨胀；
- TITAM 复现中流量固定为 1.5 kg/s；
- 工程模型中流量应由泵压升和回路压降共同决定。

### 3.3 TFE 和铯系统

- TFE 电极间隙初始充入氦气；
- 铯尚未达到正常工作压力；
- 热离子电流和电功率近似为零；
- 3 根泵用 TFE 尚未承担电磁泵供电；
- 铯储罐随冷却剂和结构逐步升温。

需要区分电极间隙中的启动氦气与 Al₂O₃/护套间隙中的传热氦气。后者在正常运行中仍存在，公开参数为 50～150 torr。

---

## 4. 冷态启动流程总览

```text
发射前预热至约373 K
        ↓
遮热罩覆盖，启动电池间歇驱动电磁泵
        ↓
收到启动许可
        ↓
3个安全鼓由0°外旋至180°
        ↓
9个控制鼓由0°向外旋转
        ↓
约125°附近首次临界
        ↓
控制鼓继续外旋至154°
        ↓
控制鼓内旋至145°
        ↓
保持至热功率达到5 kW
        ↓
以600 W/s升至35 kW
        ↓
以80 W/s升至110 kW
        ↓
堆芯入口达到400 K时抛离遮热罩
        ↓
临界后约1500 s启动TFE
        ↓
电极间隙He逐渐被约2 torr Cs替换
        ↓
热离子电流和电功率建立
        ↓
泵用TFE接替启动电池
        ↓
约40 min后达到稳态
```

文献中的多数时间以首次临界时刻为零点：

\[
t_\mathrm{crit}=0
\]

---

## 5. 阶段一：安全鼓退出

全部鼓内旋时：

\[
\rho\approx-6.0\ \$
\]

\[
k_\mathrm{eff}\approx0.952
\]

3 个安全鼓从 0°同步旋转到 180°。其总反应性价值约为 2.0 美元，单一转速为：

\[
\dot\theta_s=22.5^\circ/\mathrm{s}
\]

理论动作时间为：

\[
\Delta t_s=\frac{180}{22.5}=8\ \mathrm{s}
\]

安全鼓完全外旋后：

\[
\rho\approx-4.0\ \$
\]

\[
k_\mathrm{eff}\approx0.968
\]

该动作的目的不是使反应堆临界，而是解除发射安全状态，同时保留约 4 美元的次临界裕量。正常启动后，安全鼓保持完全外旋。

---

## 6. 阶段二：控制鼓接近临界

9 个控制鼓从 0°向外旋转，最大速度为：

\[
|\dot\theta_c|\le1.4^\circ/\mathrm{s}
\]

在安全鼓完全外旋的冷态条件下，控制鼓约在：

\[
\theta_c\approx125^\circ
\]

使反应堆首次临界。

TITAM 采用以下角度—反应性拟合式：

\[
\begin{aligned}
\rho_c(\theta)[\$]={}&-4.0-2.5\times10^{-3}\theta
+3.72\times10^{-4}\theta^2\\
&+2.21\times10^{-6}\theta^3
-3.57\times10^{-8}\theta^4
+9.41\times10^{-11}\theta^5
\end{aligned}
\]

程序应根据总反应性或 \(k_\mathrm{eff}\) 判断临界，而不是把 125°作为所有工况下的固定临界角。

按最大鼓速计算，控制鼓从 0°到 125°约需 89.3 s；加上安全鼓约 8 s 动作时间，理论最短的启动指令至首次临界时间约为 97 s，不包含控制确认和工程延迟。

---

## 7. 阶段三：初始超临界和低功率建立

### 7.1 外旋至 154°

首次临界后，控制鼓继续向外旋转至：

\[
\theta_c=154^\circ
\]

约在临界后 20 s，过剩反应性达到：

\[
\rho\approx+0.58\ \$
\]

\[
k_\mathrm{eff}\approx1.00424
\]

目的在于使中子功率脱离源区，建立可测量和可控制的裂变功率，并克服初期燃料和发射极的负温度反馈。

### 7.2 回调至 145°

随后控制鼓内旋至：

\[
\theta_c=145^\circ
\]

过剩反应性降低至：

\[
\rho\approx+0.45\ \$
\]

控制鼓保持在 145°，直到热功率达到：

\[
Q_\mathrm{th}=5\ \mathrm{kW}
\]

文献计算中，约在回调后 10 s 达到 5 kW。该回调用于限制初始功率上升和防止功率超调。

---

## 8. 阶段四：两阶段功率提升

### 8.1 快速功率斜坡

\[
5\rightarrow35\ \mathrm{kW}
\]

升功率速率：

\[
\frac{dQ}{dt}=600\ \mathrm{W/s}
\]

理论持续约 50 s。

### 8.2 慢速功率斜坡

\[
35\rightarrow110\ \mathrm{kW}
\]

升功率速率：

\[
\frac{dQ}{dt}=80\ \mathrm{W/s}
\]

理论持续约 937.5 s。TITAM 结果约在首次临界后 1070 s 达到 110 kW。此时控制鼓约位于 90°，过剩反应性仅约 0.01 美元。

升功率期间控制鼓可能总体向内移动。原因是前期已经插入过剩反应性，缓发中子和慢化剂正反馈仍推动功率上升，控制鼓内旋用于限制功率斜率，而不是立即使功率下降。

---

## 9. 阶段五：遮热罩抛离和辐射器投入

TITAM 假定遮热罩保持覆盖，直到堆芯入口冷却剂温度达到：

\[
T_\mathrm{core,in}=400\ \mathrm{K}
\]

该事件约发生在首次临界后 250 s。遮热罩抛离后，辐射器开始有效向空间排热。

主要影响：

1. 辐射器净排热快速增加；
2. 堆芯入口温度上升速率降低；
3. 一次回路和结构储能增长率下降；
4. 慢化剂因热惯性较大，短时间内温度变化有限；
5. 系统由储热主导逐步转入排热平衡主导。

程序中建议采用温度触发：

```python
if core_inlet_temperature >= 400.0:
    thermal_cover_jettisoned = True
```

250 s 是计算结果，不应作为脱离温度条件的固定事件时间。

---

## 10. 阶段六：温度反应性反馈建立

首次临界后约 5 s，燃料和发射极温度快速上升，产生负温度反应性。ZrH 慢化剂热惯性较大，约在临界后 200 s 才明显升温，并逐渐产生较强的正反馈。

总温度反应性先降低，在慢化剂开始升温后约 100 s 达到约：

\[
\rho_T\approx-0.23\ \$
\]

随后慢化剂正反馈逐渐超过燃料、电极、冷却剂和结构的负反馈，总反馈转为正值。

TITAM 对该启动方案估算的温度反馈时间常数为：

| 部件 | 时间常数 |
|---|---:|
| ZrH 慢化剂 | 337 s |
| UO₂ 燃料 | 486 s |

稳态总温度反应性反馈约为：

\[
\rho_T\approx+1.46\ \$
\]

最终控制鼓约在 88°附近提供相应负反应性，使总反应性接近零。

---

## 11. 阶段七：110 kW 后维持临界

达到 110 kW 后，控制鼓首先以最大允许速度向内旋转，直至：

\[
k_\mathrm{eff}=1
\]

随后根据温度反馈变化率微调鼓速，使：

\[
\rho_\mathrm{total}\approx0
\]

总反应性应包括：

\[
\rho_\mathrm{total}
=
\rho_\mathrm{drum}
+
\rho_\mathrm{fuel}
+
\rho_\mathrm{moderator}
+
\rho_\mathrm{electrode}
+
\rho_\mathrm{coolant}
+
\rho_\mathrm{structure}
+
\rho_\mathrm{burnup}
\]

TITAM 在部分阶段假设控制鼓能瞬时精确抵消其他反馈。原作者指出，该假设不可能由实际控制系统完全实现。工程模型应显式加入鼓速限制、反馈延迟、传感器滤波和控制器带宽。

---

## 12. 阶段八：TFE 点火和 He—Cs 置换

TITAM 假定在首次临界后约：

\[
t=1500\ \mathrm{s}
\]

启动 TFE。此时典型温度约为：

| 参数 | 数值 |
|---|---:|
| 发射极温度 | 约 1050 K |
| ZrH 慢化剂温度 | 约 795 K |
| 堆芯入口冷却剂温度 | 约 725 K |

电极间隙中的氦气逐渐被铯蒸气替换，最终铯压力约为：

\[
P_\mathrm{Cs}=2.0\ \mathrm{torr}
\]

由于铯蒸气导热系数低于氦气，置换后会出现：

- 燃料和发射极温度快速上升；
- 燃料和发射极储能增加；
- 向收集极和冷却剂的传热暂时下降；
- 收集极、冷却剂和慢化剂温度短暂下降；
- 辐射器排热短暂下降；
- 负温度反馈增强；
- 热离子电流和电功率逐渐建立；
- 控制鼓重新调整以维持临界。

不建议将铯压力瞬时阶跃到 2 torr。可采用平滑置换函数：

\[
x_\mathrm{Cs}(t)
=
\frac{1}{2}
\left[
1+\tanh\left(
\frac{t-t_\mathrm{Cs,start}}{\tau_\mathrm{Cs}}
\right)
\right]
\]

并使用时间和温度双重触发：

\[
t-t_\mathrm{crit}\ge1500\ \mathrm{s}
\]

且：

\[
T_e\ge1000\sim1050\ \mathrm{K}
\]

Voss 的系统说明还提到铯系统使用一次性穿刺阀，启动序列中由电池供电打开；公开文字称其在供电电流达到 40～60 A 时开启，但具体信号口径不够清楚，宜作为可配置参数而不是固化常数。

---

## 13. 阶段九：电磁泵和流量建立

TITAM 启动基准将总流量固定为 1.5 kg/s，并明确说明真实系统中流量会随反应堆热功率逐步增加。

Voss 系统说明给出的寿期初设计流量为：

\[
\dot m_\mathrm{BOL}=1.3\ \mathrm{kg/s}
\]

因此：

- 复现 TITAM 曲线时采用 1.5 kg/s；
- 建立 TOPAZ-II 设计模型时，建议由电磁泵与回路压降耦合，使稳态流量接近 1.3 kg/s。

物理模型应求解：

\[
\Delta P_\mathrm{pump}
(I_\mathrm{pump},T,\dot m)
=
\Delta P_\mathrm{loop}(T,\dot m)
\]

启动初期由电池驱动电磁泵，TFE 点火后，3 根泵用 TFE 逐步提供电流并接替电池。建议设置：

```text
BATTERY_PUMPING
    ↓
HYBRID_HANDOVER
    ↓
TFE_SELF_POWERED_PUMP
```

---

## 14. 阶段十：稳态建立

包括 TFE 点火在内，TITAM 基准约在首次临界后：

\[
t\approx2400\ \mathrm{s}
\]

即约 40 min 达到稳态。

| 参数 | TITAM 结果 |
|---|---:|
| 裂变功率 | 107 kW |
| 负载电功率 | 5.55 kW |
| 系统效率 | 约 5.2% |
| 单根 TFE 电流 | 190 A |
| 单根 TFE 电压 | 0.815 V |
| 控制鼓角度 | 约 88° |
| 总温度反应性反馈 | 约 +1.46 美元 |

107 kW 是该启动算例的最终热功率；115 kW 是系统 BOL 名义设计值，二者不属于同一个数据口径。

---

## 15. 冷态启动时序表

| 阶段 | 控制动作或事件 | 典型条件 | 目的 |
|---|---|---|---|
| 发射前预热 | 加热一次回路 | 约373 K | 防止NaK冻结 |
| 轨道保温 | 电池间歇驱动泵 | NaK保持液态 | 防止局部冻结 |
| 初始停堆 | 全部鼓内旋 | -6.0美元 | 深度次临界 |
| 安全解除 | 安全鼓0→180° | 约8 s | 反应性升至-4.0美元 |
| 接近临界 | 控制鼓0→约125° | ≤1.4°/s | 首次临界 |
| 初始超临界 | 控制鼓125→154° | 临界后约20 s | 建立功率 |
| 回调 | 控制鼓154→145° | +0.58→+0.45美元 | 限制超调 |
| 低功率保持 | 保持145° | 达到5 kW | 建立斜坡起点 |
| 快速斜坡 | 5→35 kW | 600 W/s | 快速升温 |
| 遮热罩抛离 | 辐射器投入 | 入口400 K，约250 s | 开始有效排热 |
| 慢速斜坡 | 35→110 kW | 80 W/s | 控制热应力 |
| 高功率保持 | 控制鼓调节 | 约1070 s | 维持临界 |
| TFE点火 | He→Cs | 约1500 s | 建立电功率 |
| 泵供电切换 | 电池→泵用TFE | 达到电流条件 | 建立自持循环 |
| 稳态 | 控制鼓约88° | 约2400 s | 温度功率稳定 |

---

## 16. 推荐状态机

```text
PRELAUNCH_HEATING
    ↓
ORBIT_THERMAL_HOLD
    ↓
SAFETY_DRUM_WITHDRAWAL
    ↓
CONTROL_DRUM_APPROACH
    ↓
INITIAL_SUPERCRITICAL_RAMP
    ↓
REACTIVITY_PULLBACK
    ↓
LOW_POWER_HOLD
    ↓
FAST_POWER_RAMP
    ↓
SLOW_POWER_RAMP
    ↓
CRITICAL_POWER_HOLD
    ↓
TFE_GAS_TRANSITION
    ↓
PUMP_POWER_HANDOVER
    ↓
STEADY_REGULATION
```

遮热罩抛离作为独立温度触发事件处理。

---

## 17. 推荐验证算例

### Case A：纯核动力学

- 安全鼓退出后反应性约 -4.0 美元；
- 控制鼓约 125°附近临界；
- 154°约 +0.58 美元；
- 145°约 +0.45 美元。

### Case B：核热反馈

- 临界后燃料和发射极先升温；
- 慢化剂约 200 s 后明显升温；
- 总温度反馈先负后正；
- 最低值接近 -0.23 美元。

### Case C：遮热罩事件

- 入口 400 K 时抛离；
- 辐射排热快速增加；
- 冷却剂升温速率降低；
- 能量守恒连续。

### Case D：功率斜坡

- 5→35 kW 平均 600 W/s；
- 35→110 kW 平均 80 W/s；
- 约 1070 s 达到 110 kW；
- 鼓速不超过 1.4°/s。

### Case E：TFE 点火

- 临界后约1500 s启动；
- 发射极约1050 K；
- 燃料和发射极快速升温；
- 收集极和冷却剂短暂降温；
- 电功率平滑建立。

### Case F：完整 TITAM 基准

```text
initial_temperature = 300 K
coolant_mass_flow_rate = 1.5 kg/s
cover_jettison_temperature = 400 K
TFE_start_time_after_critical = 1500 s
Cs_final_pressure = 2.0 torr
```

目标：

```text
Qth_steady ≈ 107 kW
Pe_steady ≈ 5.55 kW
eta ≈ 5.2%
TFE_current ≈ 190 A
TFE_voltage ≈ 0.815 V
control_drum_angle ≈ 88°
startup_duration ≈ 2400 s
```

### Case G：工程物理启动

```text
initial_temperature = 373 K
coolant_flow = pump-loop coupled
steady_flow_target ≈ 1.3 kg/s
thermal_cover = detailed radiation model
TFE_start = time-and-temperature triggered
pump_power = battery-to-TFE handover
```

---

## 18. 建模中需要避免的问题

1. 将 300 K 与 373 K 当成同一初始条件。
2. 将 TITAM 的 1.5 kg/s 与 BOL 设计值 1.3 kg/s 混用。
3. 把 125°固定为所有工况的临界角。
4. 忽略安全鼓和控制鼓的速度限制。
5. 固定在 250 s 抛离遮热罩而不检查入口温度。
6. 将铯压力瞬时阶跃至 2 torr。
7. 混淆电极间隙氦气和绝缘护套氦气。
8. TFE 点火时只打开电功率，不改变气体导热和能量输运。
9. 固定流量后仍声称完成泵—回路物理耦合。
10. 用控制鼓瞬时抵消全部温度反馈。
11. 用设计功率115 kW替换TITAM算例107 kW后仍逐点比较论文曲线。
12. 将TITAM假定时序表述为经过飞行验证的唯一实际启动规程。

---

## 19. 参考文献

[1] EL-GENK, M. S.; XUE, H.; PARAMONOV, D. V. **Start-up Simulation of a Thermionic Space Nuclear Reactor System**. *AIP Conference Proceedings*, 1993, 271: 935–950. DOI: 10.1063/1.43119.

[2] EL-GENK, M. S.; XUE, H.; PARAMONOV, D. V. **Transient Analysis and Startup Simulation of a Thermionic Space Nuclear Reactor System**. *Nuclear Technology*, 1994, 105(1): 70–86. DOI: 10.13182/NT94-A34912.

[3] VOSS, S. S. **TOPAZ II System Description**. Los Alamos National Laboratory, 1994. LA-UR-94-4; OSTI ID 10120556.

[4] EL-GENK, M. S.; PARAMONOV, D. V.; MARSHALL, A. C. **Startup Simulation of the TOPAZ-II Reactor System for Accident Conditions**. *AIP Conference Proceedings*, 1994, 301: 1059–1068. DOI: 10.1063/1.2950104.

[5] LUPPOV, A. N.; PRIKOT, K. N.; LISOCHKIN, G. A.; KWOK, K. S. **Control Drum Drive Mechanism and Regulation Characteristics of the TOPAZ II Reactor**. *AIP Conference Proceedings*, 1994, 301: 593–598. DOI: 10.1063/1.2950039.

---

## 20. 数据性质

| 数据 | 性质 |
|---|---|
| 初始温度300 K | TITAM基准假设 |
| 发射前373 K | 系统工程描述 |
| 固定1.5 kg/s | TITAM简化假设 |
| BOL流量1.3 kg/s | 系统设计值 |
| 400 K抛离遮热罩 | TITAM控制条件 |
| 约250 s抛离 | TITAM计算结果 |
| 125°临界 | 冷态TITAM近似结果 |
| 154°、145° | TITAM假定启动程序 |
| 600 W/s、80 W/s | TITAM假定功率斜坡 |
| 1500 s点火 | TITAM假定启动程序 |
| 2 torr Cs | TFE工作参数 |
| 约40 min稳态 | TITAM计算结果 |
| 实际飞行启动规程 | 公开资料不完整 |
