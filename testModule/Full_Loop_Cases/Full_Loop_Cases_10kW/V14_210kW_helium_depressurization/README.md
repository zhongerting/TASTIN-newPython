# V14 210 kW 全堆氦气瞬时完全失压事故算例说明

## 1. 目的与范围

本目录用于建立 V14 210 kW 稳态基础上的全堆氦气失压事故算例。

事故作用于 5 个代表性 TFE（`Center`、`Ring1`、`Ring2`、`Ring3`、`Ring4`），按
`1、6、9、18、24` 的倍率代表全部 58 根实际 TFE。所有 TFE 的接收极与内套管之间的
氦气隙同时失压。

本算例只研究氦气隙失压引起的热工和反应性反馈响应，不叠加外界反应性，不启用控制鼓，
不改变主回路流量，也不启用轨道外热流。

## 2. 初始状态

初始状态采用正常材料热容长期计算后的 210 kW 稳态：

```text
testModule/Full_Loop_Cases/Full_Loop_Cases_10kW/
  V14_210kW_fast_steady_temp/runs/
  physical_cp_plus2000s_from13164/checkpoint_t013864s.npz
```

该状态的绝对时间约为 `13864.2 s`。它是在早期 debug 状态基础上，恢复正常材料热容后
继续计算约 4700 s 得到的当前最终状态，不是热容缩放状态。

本目录的 `initial_state/` 中已保存该 restart 和对应 `run_config.json` 的独立副本，
避免原运行目录移动或清理后无法复现。

`initial_state/steady_restart_t013864s.npz` 的 SHA256 为
`B2576E69499234CE84DCA9ED22C2838C3FA8598CC79B426AB5ACDDAC0EA9C887`。

加载初始 restart 后：

1. 不再施加固定 210 kW 功率源；
2. 初始化点堆动力学，或在加载事故续算 restart 时恢复已保存的点堆状态；
3. 将当前温度反馈校准为零增量；
4. 记录一行相对时间 `t = 0` 的事故前状态；
5. 在第一个推进步之前触发氦气失压。

## 3. 事故模型

### 3.1 当前氦气隙

氦气位于每根 TFE 的接收极外表面和内套管内表面之间：

```text
Collector -> helium gap -> InnerClad -> NaK78 coolant channel
```

几何间隙为：

```text
r_collector_outer = 11.85 mm
r_inner_clad_inner = 11.90 mm
gap width = 0.05 mm
```

当前简化模型采用固定氦气等效换热系数：

```text
h_He,initial = 5678 W/(m2*K)
```

间隙总传热由氦气导热和表面对表面辐射并联组成。

### 3.2 失压定义

本算例采用保守的“瞬时完全失去气体导热”模型：

```text
h_He: 5678 -> 0 W/(m2*K)
```

全部 5 个代表性 TFE 在相对时间 `t = 0` 同时发生上述变化。间隙宽度、接收极和内套管
发射率以及表面对表面辐射传热保持不变。

当前代码没有氦气压力状态，也没有经过验证的压力—导热关系。因此本事故是“氦气导热
瞬时完全丧失”的保守代理模型，不声称能够解析真实泄压速率、稀薄气体传热或残余压力。

## 4. 实现边界

采用算例局部实现，不修改公共 `GapCouple2D`、`TFEUnit` 或现有反应性控制算例。

新 runner 在重建并加载系统后找到每个代表性 TFE 的 `collector_iclad_gap`，记录初始参数，
然后将耦合器的气体导热参数设为零。公共求解器继续负责间隙辐射、固体导热、流固换热、
水力、TEC 和点堆推进。

目录结构为：

```text
V14_210kW_helium_depressurization/
  README.md
  __init__.py
  run_v14_helium_depressurization.py
  test_v14_helium_depressurization.py
  initial_state/
    steady_restart_t013864s.npz
    run_config.json
  runs/
```

## 5. 反应性与功率模型

本算例保持：

```text
external reactivity = 0
control drum reactivity = 0
fixed power source = disabled
```

点堆使用当前 ReactorCore 温度反馈：

```text
rho_temperature
  = rho_fuel
  + rho_emitter
  + rho_collector
  + rho_moderator
  + rho_reflector
```

实际进入点堆的反应性为当前温度反馈相对于初始稳态参考值的增量。氦气失压改变接收极、
套管、冷却剂和堆芯温度后，反应性和裂变功率将由现有耦合链自动更新。衰变功率继续由
4 组衰变热状态根据裂变功率历史演化。

## 6. 首轮计算设置

```text
duration                  = 100 s
global time step          = 0.05 s
TEC update interval       = 0.05 s
record interval           = 0.1 s
checkpoint interval       = 10 s
solid conduction method   = implicit_euler
fluid-solid coupling      = local_implicit
total target flow         = 2.46 kg/s
external heat             = disabled
external reactivity       = 0
control drum              = disabled
```

首轮先运行 0.1 s smoke。只有 smoke 的状态连续性、反馈方向、有限值和求解收敛均通过后，
才执行 100 s 事故计算。

## 7. 温度与功率限值

每个已接受时间步后检查以下限值：

| 监视对象 | 限值 |
| --- | ---: |
| 冷却剂通道最高壁温 | 1058 K |
| 燃料芯块最高温度 | 2700 K |
| 接收极最高温度 | 1023 K |
| 慢化剂最高温度 | 930 K |
| 反射层最高温度 | 1000 K |

冷却剂通道最高壁温定义为全部代表性 TFE、全部轴向节点中内套管和外套管温度的最大值：

```text
T_wall,max = max(T_inner_clad, T_outer_clad)
```

慢化剂检查全部局部及全局慢化剂固体；反射层检查全局反射层。还应检查总功率不超过初始
功率的 2 倍，并拒绝任何非有限温度、功率或反应性。

若初始状态已经违反任一限值，事故计算不得启动。若推进后首次越限，立即停止后续推进并
保存紧急 restart。固定步长为 0.05 s，因此限值穿越时刻的分辨率为 0.05 s；第一版不增加
事件回溯或求根。

## 8. 参数记录

`history.csv` 在现有反应性控制诊断基础上至少记录：

- 绝对时间、事故相对时间、步长和事故激活状态；
- 氦气等效换热系数及剩余比例；
- 5 个代表性 TFE 的接收极平均温度和最高温度；
- 5 个代表性 TFE 的内套管平均温度和最高温度；
- 各代表性氦气隙的径向传热功率；
- 按 `1、6、9、18、24` 倍率汇总的全堆氦气隙传热功率；
- 气隙总热阻最小值和最大值；
- 通道最高壁温、芯块最高温度、慢化剂最高温度和反射层最高温度；
- 裂变功率、衰变功率和总功率；
- 各温度反馈分量、有效温度反馈和总反应性；
- 主回路流量、TEC 电参数和系统散热量。

输出目录包含：

```text
history.csv
accident_event.json
run_config.json
latest_state.json
run_summary.json
checkpoint_t*.npz
stage_01_restart.npz
```

越限时额外写入：

```text
emergency_restart.npz
limit_trip.json
```

`limit_trip.json` 记录触发对象、代表环或全局结构、轴向位置、限值、实际值和触发时间。

## 9. Restart 续算规则

气隙耦合器的气体导热参数当前不属于 `.npz` 状态。因此 runner 必须根据输入 restart 同目录
的 `run_config.json` 恢复事故状态：

- 初始稳态配置标记 `helium_accident_active = false`：先记录正常状态，再触发事故；
- 事故输出配置标记 `helium_accident_active = true`：加载后立即恢复零气体导热，不重复触发；
- 事故续算保留已保存的点堆状态和温度反馈参考值，不重新初始化或校准。

runner 必须拒绝缺少事故状态标记、点堆状态与配置不一致，或气隙数量和名称不符合预期的
restart，避免静默恢复到错误的氦气导热状态。

未越限事故 restart 的续算命令示例：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -u `
  testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\run_v14_helium_depressurization.py `
  --restart-in testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\runs\smoke_0p1s_final\stage_01_restart.npz `
  --output-dir testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\runs\restart_smoke_final `
  --duration 0.05 --dt 0.05 --tec-update-interval 0.05 `
  --record-interval 0.05 --checkpoint-interval 0.05
```

若输入事故 restart 已经越限，runner 在任何新时间步之前写入 `phase=restart_preflight`
的 `limit_trip.json` 并停止，不允许从已越限状态继续恶化。

## 10. 验证顺序与接受条件

1. 单元检查找到且只找到 5 个目标气隙，名称和倍率正确；
2. 检查事故前 `h_He = 5678 W/(m2*K)`，事故后 5 个气隙均为零；
3. 检查气体导热清零后间隙辐射通道仍保持有限；
4. 检查事故 restart 重新加载后自动恢复失压状态；
5. 使用真实稳态 restart 完成 0.1 s smoke；
6. smoke 通过后完成或按限值安全终止 100 s 计算；
7. 评估功率、反应性反馈、关键温度、氦气隙传热和数值收敛情况。

100 s 运行的接受条件不是必须算满 100 s，而是：计算过程有限且可复现；若没有越限则正常
结束并保存 restart；若越限则在首次检测到越限后停止，并完整保存触发证据和可恢复状态。

## 11. 运行命令与已完成验证

必须使用项目指定的 Python 3.12 Conda 环境。正式 100 s 计算命令为：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -u `
  testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\run_v14_helium_depressurization.py `
  --restart-in testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\initial_state\steady_restart_t013864s.npz `
  --output-dir testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\runs\accident_100s_final `
  --duration 100 --dt 0.05 --tec-update-interval 0.05 `
  --record-interval 0.1 --checkpoint-interval 10
```

2026-07-19 已完成以下验证：

1. 14 项事故开关、TEC 刷新、原始状态扫描、诊断、五类温限和 restart 规则单元测试全部通过；
2. 现有反应性控制算例 6 项回归测试全部通过；
3. 从独立初态完成 0.1 s 真实 smoke；
4. 从事故态 restart 完成 0.05 s 续算，事故绝对时刻保持不变且没有重复触发。

0.1 s smoke 的关键记录为：

| 事故相对时间 | 总功率 | 接收极最高温度 | 通道最高壁温 | 芯块最高温度 | 有效温度反馈 |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 s（事故前） | 210.000 kW | 894.192 K | 871.596 K | 2372.136 K | 0 |
| 0.1 s | 209.909 kW | 901.421 K | 866.310 K | 2372.135 K | -7.19863e-6 |

事故前按 58 根 TFE 倍率汇总的氦气隙传热约为 `198.498 kW`；完全失去气体导热后，
0.1 s 时仅由辐射承担的间隙传热约为 `4.126 kW`。该传热量按接收极流向内套管为正，
由两侧表面温差除以 `GapCouple2D.R_gap_total` 计算。

## 12. 首轮正式计算结果

正式结果目录为：

```text
runs/accident_100s_final/
```

计算没有推进满 100 s，而是在事故后 `1.75 s` 按预设限值安全终止。首次触发项是
`Ring1` 代表 TFE 的接收极，轴向位置 `z = 0.29874 m`，温度为 `1024.851 K`，
超过 `1023 K` 限值。上一个记录点 `1.70 s` 的接收极最高温度为 `1021.276 K`。
终止时的主要状态为：

| 参数 | 数值 |
| --- | ---: |
| 总功率 | 206.736 kW |
| 裂变功率 | 194.068 kW |
| 衰变功率 | 12.668 kW |
| 有效温度反馈/总反应性 | -1.17408e-4 |
| 通道最高壁温 | 832.473 K |
| 芯块最高温度 | 2371.841 K |
| 慢化剂最高温度 | 846.557 K |
| 反射层最高温度 | 794.368 K |
| 总流量 | 2.46 kg/s |

全部已记录步的流体求解和 TEC 电路计算均收敛；事故前初始记录也已按当前温度显式刷新
TEC。输出已包含 `emergency_restart.npz`、`stage_01_restart.npz`
和 `limit_trip.json`，两个终止 restart 的 SHA256 相同，说明它们对应同一已接受状态。

runner 对 TEC 调度阈值加入小于名义周期的浮点容差，并在 `run_config.json` 中同时记录
名义 `tec_update_interval_s` 和实际 `tec_scheduler_threshold_s`。正式 restart 中
末次 TEC 更新时间仅落后全局时间 0.05 s，符合每个时间步开始时刷新一次的执行顺序。

`runs/accident_100s/`、`runs/accident_100s_tec_dt/` 和 `runs/accident_100s_strict/`
是逐轮诊断时保留的对照计算，不作为推荐正式结果。最终使用 `dt = 0.025 s`、
`TEC update interval = 0.025 s` 的独立敏感性计算位于
`runs/sensitivity_dt0p025_final/`，在 1.725 s、同一代表元件和
同一轴向位置以 `1023.120 K` 触发限值。其末次 TEC 更新时间仅落后全局时间 0.025 s。
两组严格计算把越限时刻夹在
`1.725–1.75 s`，支持“瞬时完全失去氦气导热会很快触发接收极温限”的结论。

最终 runner 在预检和每个已接受时间步后直接扫描水力网络的 `T/P/h/rho/W` 原始数组及
所有注册固体的温度数组；这避免 `nanmin/nanmax` 汇总掩盖局部 NaN/Inf。事故 restart
预检前还会按当前温度显式刷新一次 TEC 电路，因此 `tec_main_converged` 表示一次真实的
当前状态求解，而不是刚重建 ThermoCalc 对象时的默认标志。
