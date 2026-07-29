# V14 210 kW 固定功率瞬时完全失冷剂事故（LOCA-1）

## 1. 算例定义

本算例从双轨外热流计算的下列 checkpoint 接管：

```text
V14_210kW_fixed_power_external_heat_2orbits/
  runs/two_orbits_from13864_20260720/checkpoint_t016864s.npz
```

checkpoint 保存的绝对时间即事故时刻，并保持原 N18 外热流周期、相位、TEC 电压、
查表库、导线电阻、辐射率和固体温度。事故后堆芯总功率固定为 `210000 W`，不启用
点堆功率演化；反应性反馈仅作为诊断记录，暂不考虑 NaK 排空/空泡反应性。

## 2. 事故模型

事故在 `t=0+` 瞬时完成：

1. 全回路 NaK78 冷却剂消失；
2. 不再求解压力、流量或流体能量方程；
3. 全部 `FluidSolidCouple` 从 `SystemManager` 移除，同时删除它们安装的
   `solid_bc` 和 `_local_implicit_flux_bc`；
4. 清零 `Q_wall/Q_vol/implicit_coeff` 和全部连接流量；
5. 五个代表性 TFE 的内外套管之间建立真空 `GapCouple2D`，气体导热为零，
   两侧发射率均为 `0.8`；
6. 普通管道和集流环的空管内壁绝热；
7. 集流环流体—壁面及流体—热管蒸发段换热停止。热管固体热容、内部状态、
   翅片、空间辐射和轨道外热仍保留，避免冻结或删除结构储能。

该模型是“瞬时完全排空后的热工后果包络”，不计算泄压、喷放、排空时间、液膜、
NaK 蒸气、两相流或局部液体滞留。

## 3. 接管顺序

```text
按 checkpoint 邻近 run_config 重建 V14
  -> 加载 checkpoint 并完成原耦合同步
  -> 按当前温度立即刷新 TEC
  -> 保存事故前节点快照
  -> 移除全部 FluidSolidCouple 及其边界条件
  -> 建立五组内外套管真空辐射耦合
  -> 清空流体源项与流量并冻结水力求解
  -> 刷新固体边界缓存
  -> 保持 210 kW 固定功率推进固体、TEC 和外热边界
```

## 4. 输出与记录频率

默认全局步长为 `0.05 s`，每 `0.2 s` 写一行 `history.csv`，同时写一份压缩节点快照：

```text
snapshot_pre_accident.npz
snapshot_tplus_00000.000s.npz
snapshot_tplus_00000.200s.npz
...
```

节点快照包括：

- 流体节点和连接名称、分类、事故前参考 `T/P/h/W/v`；
- 事故后冷却剂状态显式写为 `present=0`、`T/P/h=NaN`、`W/v=0`；
- 所有注册固体的名称、分类、原始形状和温度数组；
- 主 TEC 总电流、电压、功率和收敛标志；
- 五个代表性 TFE 的轴向电流密度、电极电势差、电子冷却/加热热流与功率、
  发射极和接收极焦耳热功率；
- 燃料、电极、慢化剂、反射层和总温度反应性反馈。

流体分类为 `core/ordinary_pipe/collector_ring`；固体分类为
`core_structure/pipe_wall/heat_pipe`。

## 5. 能量审计口径

事故前能量守恒沿用源双轨算例。事故瞬时删除 NaK 质量及其储能，跨越事故时刻不要求
能量闭合。事故后从新的固体参考状态重新累计相对显热，并记录：

```text
solid_relative_energy_J
solid_energy_rate_W
radiator_net_rejection_W
tec_main_electric_power_W
post_accident_energy_residual_W
```

当前固体相对能量采用 `rho * cp(T) * V * (T-T0)` 估算，只用于短期趋势和审计定位，
不替代材料焓积分。

## 6. 运行

短期默认计算 `0.4 s`：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -u `
  testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_fixed_power_LOCA_1\run_v14_210kw_fixed_power_loca_1.py `
  --output-dir testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_fixed_power_LOCA_1\runs\smoke_0p4s
```

后续发射率敏感性使用：

```text
--radiation-emissivity 0.5
--radiation-emissivity 0.2
```

本轮先验证 `0.8` 基线，不启动长时计算。

## 7. 2026-07-20 短期验证

`runs/smoke_0p2s_final/` 已从指定 checkpoint 完成 `0.2 s` 真实短测：

- 事故时刻为绝对时间 `16864.45 s`；
- 移除 `58` 个 `FluidSolidCouple`，新增 `5` 个套管真空辐射间隙；
- 水力求解关闭，事故后全部质量流量和流速为零；
- `history.csv` 仅记录 `0` 和 `0.2 s`，没有终点微小空步；
- TEC 收敛，主电流 `212.469 A`、端电压 `50.65 V`；
- `0.2 s` 时堆芯结构最高温度 `2372.133 K`，套管间总辐射传热约 `391.025 W`；
- 总温度反应性反馈相对事故时刻变化 `-2.22776e-6`；
- 节点快照包含 `300` 个流体节点、`315` 个连接、`4730` 个固体温度值和
  `5 x 37` 的 TEC 轴向电数据，固体及 TEC 数组均为有限值。

该短测只验证事故接管、边界切换、记录和数值推进，不构成长时事故结论。

## 8. 后处理 history 文件

从后续新运行开始，每个 `record_interval_s` 时刻同时写出以下四个独立的长表 CSV；原 `history.csv` 继续保留为运行概览和能量审计摘要。

### 8.1 `history_coolant.csv`

- `category`：`core/ordinary_pipe/collector_ring`。
- 控制体行（`entity_type=volume`）：冷却剂温度、压力、焓。
- 连接行（`entity_type=junction`）：质量流量、流速。
- 完全失冷剂事故后，实际 `temperature_K/pressure_Pa/enthalpy_J_kg` 为 `NaN`，`mass_flow_kg_s/velocity_m_s` 为零；`reference_*` 列保留事故前 checkpoint 数值。

### 8.2 `history_solids.csv`

- `category`：`core_structure/pipe_wall/heat_pipe`。
- 每行对应一个固体展平节点，记录 `solid_name`、原始 `solid_shape`、`flat_node_index` 和 `temperature_K`。

### 8.3 `history_electrical.csv`

- 每行对应一个代表性 TFE 的一个轴向节点。
- 重复记录主电路总电流、端电压、总电功率和收敛标志。
- 逐 TFE/轴向记录电流密度、发射极电势、接收极电势、电势差、电子冷却/加热热流密度及功率、发射极/接收极焦耳热轴向功率。

### 8.4 `history_reactivity.csv`

- 每个记录时刻一行。
- 记录燃料、电极、慢化剂、反射层、总绝对反应性反馈，以及相对事故时刻的总反馈变化。

`snapshot_tplus_*.npz` 继续保留同一时刻的完整数组，作为无 CSV 精度和表结构限制的复核数据。

## 9. 温度失效边界与发射率补充算例

后续 ε=0.2 和 ε=0.5 算例保持堆芯功率 `210000 W`，计算步长 `0.05 s`，每 `0.5 s` 写一次 history 和快照。每个计算步后检查：

- 接收极最高温度 ≥ `1500 K`；
- 发射极最高温度 ≥ `3000 K`；
- 存在冷却剂时，其最高温度 ≥ `1058 K`；
- 慢化剂最高温度 ≥ `930 K`；
- 反射层最高温度 ≥ `1000 K`。

任一条件满足即判定反应堆失效，保存触发时刻终态，并在 `history.csv/failure_reason` 和 `run_summary.json/stop_reason` 写入原因。由于本算例在 `t=0+` 完全排空 NaK，事故后冷却剂温度不存在，因此冷却剂判据为不适用；其 `coolant_max_T_K` 记录为 `NaN`。运行时长参数仅作为未触发失效时的安全上限。

## 10. ε=0.2/0.5 计算结果

两个算例并行运行并正常完成，输出目录分别为：

```text
runs/LOCA_1_eps020_until_failure_record0p5s
runs/LOCA_1_eps050_until_failure_record0p5s
```

- ε=0.2：事故后 `19.95 s` 接收极达到 `1500.747 K`，触发 `collector_temperature_limit`；发射极 `2183.396 K`，慢化剂 `847.623 K`，反射层 `796.365 K`，堆芯结构最高 `2408.613 K`。
- ε=0.5：事故后 `24.25 s` 接收极达到 `1500.512 K`，触发 `collector_temperature_limit`；发射极 `2198.147 K`，慢化剂 `850.579 K`，反射层 `796.364 K`，堆芯结构最高 `2425.339 K`。

两者功率均保持 `210000 W`，常规记录间隔为 `0.5 s`，并额外保存实际触发时刻；`stderr.log` 均为空。

## 11. 反应性反馈与 5 s 停堆六工况

启用 `--enable-reactivity-feedback` 后，事故初态以 210 kW 初始化点堆，后续功率由温度反馈和外加反应性共同决定，不再强制固定为 210 kW。停堆工况从事故后 5.0 s 起持续施加 `-2 $`；按当前 `β_eff=0.0079321` 换算为无量纲反应性 `-0.0158642`。

六个工况分别为 ε=0.2/0.5/0.8 的纯温度反馈，以及相同三种发射率下的 5 s、-2 $ 停堆。最长计算时间均为 `1000 s`，并保留第 9 节温度失效边界。

使用 `--staged-recording` 时，history 和 snapshot 采用：

- 0–20 s：每 0.5 s；
- 20–100 s：每 2 s；
- 100–400 s：每 5 s；
- 400–600 s：每 10 s；
- 600–1000 s：每 20 s；
- 温度边界触发或计算终点：额外保存准确终态。

`history_reactivity.csv` 额外记录点堆总功率、裂变功率、衰变热、外加反应性（无量纲和 dollars）、有效温度反馈及总反应性。

## 12. 停堆后 TEC 开路

5 s 停堆工况在主 TEC 电流降至 `0.01 A` 或以下时执行一次性开路切换：

- 所有代表性 TFE 的电流、电势、电子热流和发射极/接收极焦耳热清零；
- 保留 TEC 间隙的被动传热；
- 设置 `core.enable_tec_coupled=False`，后续时间步不再调用 ThermoCalc；
- 在 summary、snapshot 和 `history_reactivity.csv` 中记录开路标志及触发时间。

此前未使用该规则的三套停堆输出仅作诊断参考；正式停堆结果使用带 `tecopen001A` 后缀的新输出目录。
