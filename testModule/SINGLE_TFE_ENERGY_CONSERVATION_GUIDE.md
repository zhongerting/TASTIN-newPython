# 单 TFE 能量守恒算例指南

## 1. 定位

[`test_single_tfe_energy_conservation_v7.py`](./test_single_tfe_energy_conservation_v7.py) 是 v7 CaseA 中心通道单根 TFE 的专题诊断入口。它从全系统快照选择性映射中心通道状态，不调用全系统 `load_global_state()`，也不恢复外围慢化剂、点堆或全堆倍率。

算例分为两个阶段：

| 模式 | 用途 | 当前可运行性 |
| --- | --- | --- |
| `thermal-baseline` | 固定内热源、严格绝热外边界下的热工能量守恒基线 | 可直接运行 |
| `tec` | 在热工基线上增加 TEC 电子热流、焦耳热和端电功率 | 已完成首次 `1 s` 基线 |

## 2. 固定物理配置

| 项目 | 数值 |
| --- | ---: |
| 源快照 | `test_core_assemble_v7_caseA_faststeady_restart_t18800.npz` |
| 源快照时间 | `18800 s` |
| 代表元件 | `Center` |
| 冷却剂 | `Sodium` |
| 中心 TFE 内热源 | `3489.8792830760617 W` |
| 入口速度 | `0.7533237780835761 m/s` |
| 入口质量流量 | `0.03513560646568767 kg/s` |
| 出口定压 | `161961.33075474895 Pa` |
| TEC 目标电压 | `0.8 V` |
| 导线电阻 | `[0, 0, 0, 0] Ohm` |

入口速度使用源快照进口联箱状态下的 Na 密度和固定入口质量流量换算。出口压力继承源快照 `outlet_plenum`，因为单 TFE 出口边界直接连接流道出口。

CLI 不提供功率覆盖参数。改变功率时应增加独立工况，避免复用守恒基线后混淆审计口径。

## 3. 自动生成流程

脚本启动时依次执行：

1. 校验源快照 `Fluid/shape == [176, 179]`、时间为 `18800 s`、中心 TFE 轴向节点数为 `37`。
2. 按 CaseA 几何建立单根 TFE 和 37 节点非均匀 Na 流道。
3. 通过 `strict_adiabatic_single_tfe=True` 跳过 CO2 间隙和慢化剂链路，清空外套管右边界条件。
4. 导入中心 TFE 的 `Pellet`、`Emitter`、`Collector`、`InnerClad`、`OuterClad` 温度。
5. 导入中心流道 `P/T/h/W` 分布，并设置固定入口流量和出口压力。
6. 重置旧 `dTdt`、边界缓存、电子热流、焦耳热和电学缓存。
7. 重新施加固定中心 TFE 内热源。
8. 将单 TFE 相对时间重置为 `0 s`，并在输出 restart 中写入 `SingleTFE/source_snapshot_time_s=18800`。

## 4. 运行命令

从仓库根目录运行热工基线：

```powershell
python testModule\test_single_tfe_energy_conservation_v7.py `
  --mode thermal-baseline `
  --duration-s 10 `
  --max-dt-s 0.1
```

仅检查映射和严格绝热边界时可将 `--duration-s` 设为 `0`。

运行 TEC 阶段时必须使用 ABI 匹配的 Python 3.12 解释器：

```powershell
python testModule\test_single_tfe_energy_conservation_v7.py `
  --mode tec `
  --duration-s 1 `
  --max-dt-s 0.01
```

`tec` 模式在推进前执行能力检查。若 `bindings.cpp` 将 `sideAreaE/sideAreaC` 折叠为单根 TFE 标量，或 `SingleTEC` 未暴露 `phiE/phiC/Vd`，脚本会明确失败，不会退化为热工模式。2026-06-01 已用重建后的 `te_solver.cp312-win_amd64.pyd` 验证该检查通过。

## 5. 输出与审计

默认产物目录：

```text
testModule/single_tfe_energy_conservation_v7/
```

| 文件 | 内容 |
| --- | --- |
| `latest_restart.npz` | 单 TFE 相对时间 restart，附带源快照时间 |
| `latest_summary.json` | 固定配置、边界检查和最终审计摘要 |
| `energy_balance_history.csv` | 每个时间步的全局能量项 |
| `solid_energy_balance_latest.csv` | 每个固体的热源、边界净流入、有限差分储能、`cap*dTdt` 储能和残差 |
| `interface_balance_latest.csv` | 普通固固、TEC 极间隙、流固和流体源合并的逐节点闭合 |
| `fluid_volume_balance_latest.csv` | 每个流体控制体的焓流、实际应用源项、末态同步源项、储能和矩阵残差 |
| `fluid_profile_latest.csv` | 最新 37 节点流体剖面 |
| `solid_temperature_profile_latest.csv` | 最新五个固体域温度剖面 |
| `tec_node_balance_latest.csv` | TEC 模式下的节点电子热流、焦耳热和电势 |

全局瞬态残差定义：

```text
residual =
  nuclear_heat
  - coolant_enthalpy_pickup
  - electrical_output
  - solid_storage_rate
  - fluid_storage_rate
  - outer_clad_heat_loss
```

严格绝热模式下：

```text
outer_clad_heat_loss = 0 W
```

`thermal-baseline` 当前验收门槛：

- 外套管热损失逐节点严格等于 `0 W`；
- 普通固固接口累计残差绝对值 `< 1e-6 W`；
- 流固功率映射累计误差 `< 1e-6 W`；
- 流体能量矩阵有限控制体残差累计绝对值 `< 1e-6 W`；
- `10 s` 最后 `1 s` 全局平均相对残差 `< 1%`。

2026-06-01 正式 `10 s` 热工基线结果：最后 `1 s` 平均相对残差约 `0.279%`，普通固固累计误差约 `1.72e-10 W`，流固累计误差约 `2.82e-10 W`，流体矩阵有限控制体累计绝对残差约 `1.51e-9 W`。

TEC 首次 `1 s` 基线只记录误差，不固化硬阈值。当前方向定义下观测到：

```text
电子边界功率差 ~= 端电功率 + 电极焦耳热
```

2026-06-01 首次结果为：电子边界功率差约 `289.686 W`，端电功率约 `234.357 W`，焦耳热约 `55.050 W`，闭合误差约 `0.280 W`。若使用“端电功率减焦耳热”口径会产生约 `110.379 W` 差值，不应作为当前实现的验收式。

## 6.1 已修复的底层问题

`HeatConduction.step()` 曾把可变的 `self.T` 本体直接传给 `solve_ivp(y0=...)`，而 RHS 又会原地写回 `self.T`。试探点评估因此可能污染积分初值，造成 kW 到 MW 级伪储能。当前已改为传入 `self.T.copy()`。

## 7. 验证顺序

1. 用 `--duration-s 0` 检查快照映射、固定参数和严格零热流边界。
2. 用短时 `thermal-baseline` 运行检查储能项和全局残差。
3. 修复 ThermoCalc 后先运行 `--mode tec --duration-s 0` 执行能力检查。
4. 再运行短时 `tec` smoke，检查节点电子热流、焦耳热、端电功率和总残差。

## 8. 2026-06-02 FVM 一致焦耳热基线

TEC 模式当前使用 C++ `VcalcFVM()` 输出的逐轴向节点焦耳热功率：

```text
joulePowerE / joulePowerC [W]
```

Python 只按轴向列内控制体体积比例映射到二维电极网格。`tec_node_balance_latest.csv` 同时记录：

- C++ 权威节点焦耳热。
- 映射后二维焦耳热按轴向回聚合值。
- 映射值与 C++ 值的差。
- 旧节点梯度法焦耳热，仅用于修复前后对照。
- `电子边界功率差 - 端电功率 - C++焦耳热`。

2026-06-02 正式命令：

```powershell
python testModule\test_single_tfe_energy_conservation_v7.py `
  --mode tec `
  --duration-s 1 `
  --max-dt-s 0.01 `
  --output-dir testModule\single_tfe_energy_conservation_v7_tec_fvm_1s
```

结果：

| 指标 | 数值 |
| --- | ---: |
| 二维映射与 C++ 节点功率总差 | `0 W` |
| 旧梯度法相对 C++ 少计的焦耳热 | `0.254738 W` |
| TEC 转换闭合差 | `0.025015 W` |
| 最终全局残差 | `0.089127 W` |
| 最终全局相对残差 | `0.00255%` |
| 最后 `1 s` 平均相对残差 | `0.0814%` |
| 普通固固接口累计绝对误差 | `1.07e-10 W` |
| 流固接口累计绝对误差 | `2.97e-10 W` |

剩余 TEC 转换闭合差主要与外层电路电流停止条件有关。本轮没有收紧该阈值。
