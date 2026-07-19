# RadiatorPipeWithFin Detailed Intro

本文档说明 `Components/RadiatorPipeWithFin.py`。该组件用于 TOPAZ-II 原始 NaK 管翅式辐射器，不适用于钾热管辐射器。

## 1. 模型定位

`RadiatorPipeWithFin` 表示一根 NaK 辐射管及其焊接铜辐射带：

```text
NaK fluid channel
  -> FluidSolidCouple
  -> cylindrical tube wall HeatConduction2D
  -> bare tube dynamic radiation
  -> reduced-order copper fin branch
```

它参考 `HPwithFin` 的降维翅片思路，但不包含 `HeatPipe2D`、吸液芯、蒸汽区或热管相变等效模型。

## 2. 主要物理假设

- 管壁为圆柱坐标 `HeatConduction2D`，默认可用 `n_radial_wall=1` 表示薄壁。
- 管内流体换热由外部传入的 `fluid_channel` 和 `correlation_func` 决定。
- 裸管外表面通过 `DynamicRadiationResistanceBC` 向 `T_space` 辐射。
- 铜带按每个轴向节点独立求解一维准稳态导热。
- 铜带被视为管两侧对称半翅片：局部净宽度为 `fin_width(z) - tube_outer_diameter`，每侧翅高为该净宽度的一半。
- `fin_area_scale` 只保留为几何/标定面积缩放；V15 默认取 `1.0`，内外面遮挡由独立角系数计算，不再借此参数混合表达。

## 3. 关键参数

| 参数 | 含义 |
|---|---|
| `fluid_channel` | 与管壁内表面耦合的 NaK 流体通道 |
| `tube_inner_diameter` / `tube_outer_diameter` | 辐射管内外径 |
| `tube_length` / `n_axial` | 管长和轴向离散 |
| `n_radial_wall` | 管壁径向离散 |
| `fin_thickness` | 铜带厚度 |
| `fin_width_upper` / `fin_width_lower` | 截锥上/下端局部节距宽度 |
| `n_fin_width` | 每个轴向切片的翅片横向离散数 |
| `tube_emissivity` / `fin_emissivity` | 裸管和铜带发射率 |
| `fin_area_scale` | 翅片几何/标定面积缩放；V15 默认 `1.0` |
| `fin_conductivity` | 铜带导热系数，默认常数 |

## 4. 输出与诊断

`get_heat_exchange_breakdown()` 返回：

- `bare_radiation`：裸管辐射功率分布，单位 W；
- `fin_radiation`：翅片辐射功率分布，单位 W；
- `fin_net_from_root`：翅片从管壁根部抽取的净热量，单位 W；
- `fin_conductance`：翅片等效热导，单位 W/K；
- `fin_effective_temperature`：挂回管壁边界的等效外部温度；
- `fin_equivalent_resistance`：挂回管壁边界的等效热阻。

`get_fin_temperature_distribution()` 返回形状为 `(n_axial, n_fin_width)` 的翅片温度场。

2026-06-16 起，离散翅片求解使用上一时间步翅片温度作为 warm-start 初值；首步或状态无效时自动回退到根部温度常值初值。`get_heat_exchange_breakdown()` 同步返回 `fin_iteration_count`、`fin_max_delta` 和 `fin_used_warm_start`，用于性能诊断。

## 5. 使用原则

- 系统级 TOPAZ-II 管翅式辐射器应优先使用该组件，而不是把铜带折算成固定等效面积。
- 如果只做单点快速标定，可降低 `n_fin_width`；用于变工况预测时建议检查 `n_fin_width=5/10/20` 的网格敏感性。
- 该组件没有实现轴向铜带导热、焊缝详细热阻、支撑遮挡角系数或表面对表面辐射网络。

## 6. 推荐测试入口

显式翅片模型的稳态回归入口是 `CoolantLoop/run_topaz2_pipefin_steady_test.py`。默认工况为 `Tin=823 K`、总流量 `1.3 kg/s`、`eps_tube=eps_fin=0.80`、`n_axial=8`、`n_fin_width=12`、`duration=500 s`。该测试会生成 history、latest state 和 `*_steady_summary.json`，用于检查出口温度、末段温度斜率、能量残差比例、有效面积和翅片根端温差。

当前 `Tout=727 K` 是后续重新标定显式翅片模型的目标，不是该测试的硬通过条件。

## 7. 分布式翅片外热流

RadiatorPipeWithFin 可通过 set_fin_external_heat_source() 将外热流直接加入准稳态翅片方程。每根辐射管使用一个周向热流密度，该密度沿轴向各段相同；每段只乘本段面积。

翅片吸热采用单侧几何投影面积，并在 V15 中乘入 0.992 转换效率和 Kirchhoff 吸收率 alpha=epsilon_fin：

    A_abs_eff = fin_strip_width * fin_height * illumination_factor * 0.992 * epsilon_fin
    Q_abs = q_flux * A_abs_eff
    Q_net_from_root = Q_fin_radiation - Q_abs

fin_area_scale 和 fin_view_factor 仅修正双侧辐射几何面积，不参与单侧受照几何面积；吸收效率与吸收率通过 illuminated_area_scale 独立进入。get_external_heat_absorption_distribution() 返回管壁、翅片和总吸热功率分布，单位均为 W；get_heat_exchange_breakdown() 额外返回 fin_absorption。

能量审计中，gross_rejection 表示管壁与翅片实际向空间辐射的总功率；fin_net_from_root 和 net_rejection 用于描述翅片对管壁根部的净抽热。若审计式已单独加入外热输入，辐射器出口项必须使用 gross_rejection，不能再次用已经减去翅片吸热的 net_rejection。

## 8. 2026-06-25 runtime radiation background

`RadiatorPipeWithFin` now exposes `set_radiation_background_temperature(value)` for optional external boundary modifiers such as the V13 startup thermal shield. The value may be a scalar or one value per axial node. It updates:

- bare tube `DynamicRadiationResistanceBC` background temperature;
- the reduced-order fin radiation solve;
- `get_heat_exchange_breakdown()` tube/fin radiation diagnostics.

When the shield is disabled, the component keeps the historical default background `T_space`, so existing V12/V13 cases are unchanged. The current shield integration is intentionally a boundary modifier, not a new radiator solid or fluid component.

2026-06-16 warm-start 优化 100 s 对比：默认 `eps_tube=eps_fin=0.80`、`n_axial=8`、`n_fin_width=12` 下，墙钟时间由 `152.88 s` 降至 `143.04 s`，平均翅片迭代数由 `5.44` 降至 `2.33`；出口温度差约 `3.8e-8 K`，总排热相对差约 `6.0e-8 %`。


## 9. 2026-07-15 V15 内外面角系数与吸收率

V15 对管壁和翅片采用同一套四项角系数：外面上/下均为 1.0，内面上为 0.0、下为 0.17。一维等效过程保持各项可诊断：

    F_outer = (F_outer_up + F_outer_down) / 2 = 1.0
    F_inner = F_inner_up + F_inner_down = 0.17
    F_face = (F_outer + F_inner) / 2 = 0.585
    epsilon_eff = epsilon_surface * F_face

因此 V15 默认 epsilon_surface=0.8 时，管壁和翅片的 epsilon_eff 都为 0.468。管壁动态辐射边界、翅片准稳态方程、切线热导和排热诊断统一使用该有效发射率。未传入四项角系数的旧调用默认得到 F_face=1.0，保持原有组件行为。

外热仅作用于外侧受照面，不乘上述辐射角系数；管壁和翅片分别使用自身 alpha=epsilon_surface，并乘 0.992 吸收效率：Q_abs=0.992*alpha*q_incident*A_projected。0.992 不用于向外自发辐射。
