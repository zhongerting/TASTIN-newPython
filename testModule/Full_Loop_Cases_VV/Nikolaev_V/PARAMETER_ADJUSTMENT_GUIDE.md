# Nikolaev_V 参数调整边界说明

本文约束 Nikolaev 1995 `SPACE-R` 单电池 TFE mock-up 验证算例的参数调整。

## 验证目标

目标文献：

- Yuri V. Nikolaev 等，**A single-cell TFE mock-up of the thermionic nuclear power system "Space-R"**，AIP Conference Proceedings 324, 815 (1995)，DOI `10.1063/1.47120`。
- 本地 PDF：`e:\文献阅读\nikolaev1995.pdf`。

当前 `Nikolaev_V` 首版只做表格级验证：

- Table 1：TOPAZ-II 与 SPACE-R prototype 主要参数对比。
- Table 2：300 W 平均 TFE 的电压、电流、热输入、发射极温度和效率。
- Table 3：燃料最高温度随自由体积和径向非均匀因子的变化。
- Table 4：10 年寿命约束下最大 venting-system capillary diameter。

Figure 4 的 VAC 曲线当前没有数字化，暂不作为定量误差指标。

## 必须固定的量

| 参数/数据 | 固定原则 |
| --- | --- |
| Table 1-4 原始数值 | 必须按文献表格录入，不得用模型输出回填。 |
| 论文题名、DOI、出处 | 必须保持可追溯。 |
| mock-up heater length | 文献给出 tungsten heating element `350 mm`，作为几何锚点。 |
| interelectrode gap | 文献给出 `0.5 mm`，作为 TOPAZ-II 相似建模锚点。 |
| collector temperature for local VAC approximation | 文献说明局部 VAC 基于 `collector temperature 870 K`。 |

## 可以调整的量

论文没有给出足以唯一重建原始电-热耦合程序的边界条件，因此以下量可以作为全局参数调整：

| 参数 | 代码字段 | 说明 |
| --- | --- | --- |
| nominal output power | `nominal_output_power_w` | Table 2 三个电压点约为 300 W，可以小范围调整。 |
| voltage-dependent thermal input | `thermal_power_reference_kw`, `low_voltage_thermal_power_slope_kw_per_v` | 用于闭合 Table 2 的 `Q=4.1-4.2 kW`。 |
| emitter temperature closure | `emitter_temp_reference_k`, `emitter_temp_linear_k_per_v`, `emitter_temp_quadratic_k_per_v2` | 用二次函数重建 Table 2 的 1880/1890/1910 K。 |
| effective height / active core assumptions | `effective_height_m`, `active_core_length_m` | 文中说明 active core length 40 cm、effective height 50 cm，但 Table 1 OCR 行存在歧义。 |

所有调整必须全局作用于 Nikolaev 表格工况，不能逐点单独调参。

## 当前模型定位

`nikolaev_single_tfe_model.py` 是 TOPAZ-II 相似的紧凑积分模型，不是原始 Nikolaev/Davydov 节点电-热求解器。它的作用是：

1. 把论文中可追溯的数值锚点固化成可运行验证。
2. 暴露缺失闭合参数，便于后续调参和误差追踪。
3. 为后续 digitize Figure 4 后接入更高保真 TEC 计算提供基线。

## 后续升级路径

1. 数字化 Figure 4 的六条 VAC 曲线，加入 `experimental_data/figure4_vac_digitized.csv`。
2. 用当前 ThermoCalc 或 Venable_V 的单 TFE 电输出链路替代紧凑电流闭合。
3. 将 350 mm heater、0.5 mm IG、870 K collector boundary 明确映射到轴向节点模型。
4. 若高保真模型不能闭合，再按本文记录全局参数调整，不允许逐点拟合。

## 2026-07-03 ThermoCalc path status

The earlier table-level path in `nikolaev_single_tfe_model.py` is now explicitly treated as a reconstruction baseline only. It computes current as `nominal_output_power_w / voltage_v`, so the very small Table 2 current error is not an independent electrical prediction.

The preferred electrical-validation path is now `nikolaev_thermocalc_runner.py`: fixed Nikolaev geometry and heat input are mapped to emitter/collector temperature fields by `nikolaev_thermocalc_model.py`, then ThermoCalc solves the fixed-voltage current and reports `Iout`. This path does not use `I=P/V` to set current.

Current runs:

| Run | Key parameters | Current MAE | Max current error | Temperature status |
| --- | --- | ---:| ---:| --- |
| `20260703_nikolaev_thermocalc_baseline` | `R_EC=0.248 K/W`, `Tcs=610 K`, `Rwire=0 ohm` | 135.717 A | 306.731 A | Mean emitter temperature about 1886.8-1911.6 K, close to Table 2 scale, but low-voltage current is much too high. |
| `20260703_nikolaev_thermocalc_balanced_candidate` | `R_EC=0.260 K/W`, `Tcs=560 K`, `Rwire=0.0005 ohm` | 23.489 A | 45.109 A | Current improves, but mean emitter temperature rises to 1936-1962 K, above the Table 2 values. |
| Coarse unconstrained current fit | `R_EC=0.300 K/W`, `Tcs=560 K`, `Rwire=0.001 ohm` | 11.847 A | 15.834 A | Rejected for validation because mean emitter temperature is about 2100 K. |

Adjustment rule from this point: do not rank parameters by current error alone. A valid setting must also respect the published emitter-temperature scale. The present model is still prescribed-temperature/fixed thermal-network input to ThermoCalc; electron cooling, electron heat transport, and Joule heat feedback are not yet iterated back into the thermal network.

## 2026-07-04 thermoelectric closed-loop update

A local closed-loop path has been added for the Nikolaev single-TFE case:

- `nikolaev_thermoelectric_closed_loop.py` closes a two-node-per-axial-station thermal network around ThermoCalc.
- `nikolaev_closed_loop_runner.py` runs Table 2 with iterative thermal feedback and writes CSV/JSON/Markdown reports.
- The coupling follows the current `ReactorCore` convention: electron emitter cooling and collector heating use the `J/phiE/TE/UE/UC` fields, while electrode Joule heat uses ThermoCalc's authoritative `joulePowerE/C` arrays. The code does not reconstruct Joule heat from voltage gradients.

Closed-loop baseline:

- Run: `runs/20260704_nikolaev_closed_loop_baseline`
- Parameters: `R_EC=0.248 K/W`, `Tcs=610 K`, `Rwire=0 ohm`, `R_CB=0.010 K/W`, relaxation `0.25`.
- Result: all ThermoCalc and outer thermal iterations converged, but emitter temperature is too low; current MAE `85.292 A`, emitter-temperature MAE `163.277 K`.

Balanced closed-loop candidate:

- Run: `runs/20260704_nikolaev_closed_loop_balanced_candidate`
- Parameters: `R_EC=0.340 K/W`, `Tcs=560 K`, `Rwire=0.0005 ohm`, `R_CB=0.010 K/W`, relaxation `0.20`, tolerance `0.75 K`.
- Result: all three Table 2 points converge in 27 outer iterations.
- Metrics: current MAE `6.462 A`, max current error `12.175 A`, electric-power MAE `4.876 W`, emitter-temperature MAE `24.641 K`.
- Point comparison:
  - `0.7 V`: `I_calc=416.825 A` vs `429 A`, `Te_mean=1851.177 K` vs `1880 K`.
  - `0.8 V`: `I_calc=371.154 A` vs `375 A`, `Te_mean=1862.624 K` vs `1890 K`.
  - `0.9 V`: `I_calc=336.366 A` vs `333 A`, `Te_mean=1892.275 K` vs `1910 K`.

Interpretation: this is the first Nikolaev_V path where ThermoCalc electrical output, electron heat transport, and authoritative Joule heat are iterated back into the thermal network. It is still a simplified local network, not a full Core/TFEUnit/SystemManager model.

## 2026-07-04 physical thermal-hydraulic TFE loop

A more physical single-TFE validation path has been added:

- `nikolaev_physical_tfe_loop.py` models heater power, ThermoCalc electronic heat feedback, collector heat rejection, and coolant enthalpy rise with `m_dot * cp * dT`.
- `nikolaev_physical_loop_runner.py` runs Table 2 and records coupled electrical, thermal, coolant, and energy-balance metrics.
- This path does not prescribe collector temperature. Collector temperature follows local heat rejection to a flowing coolant stream.

Baseline physical-flow run:

- Run: `runs/20260704_nikolaev_physical_tfe_baseline`
- Parameters: `Tin=770 K`, `m_dot=0.040 kg/s`, `cp=1000 J/kg/K`, `h=8500 W/m2/K`, `R_EC=0.340 K/W`, `Tcs=560 K`, `Rwire=0.0005 ohm`.
- Result: energy balance residual below `1.2e-11 W`, but the `0.7 V` point did not converge within 40 outer iterations and the model remained too cold. Current MAE `25.023 A`; emitter-temperature MAE `64.603 K`.

Balanced physical-flow candidate:

- Run: `runs/20260704_nikolaev_physical_tfe_balanced_candidate`
- Parameters: `Tin=770 K`, `m_dot=0.035 kg/s`, `cp=1000 J/kg/K`, `h=8500 W/m2/K`, `R_EC=0.380 K/W`, `Tcs=560 K`, `Rwire=0.0005 ohm`, relaxation `0.20`, tolerance `0.75 K`.
- All three Table 2 points converged with max outer iterations `20`.
- Metrics: current MAE `7.272 A`, max current error `10.717 A`, electric-power MAE `5.840 W`, emitter-temperature MAE `4.208 K`, max coolant energy residual `5.912e-12 W`.
- Point comparison:
  - `0.7 V`: `I_calc=418.953 A` vs `429 A`, `Te_mean=1877.900 K` vs `1880 K`, `Tout=879.109 K`.
  - `0.8 V`: `I_calc=376.051 A` vs `375 A`, `Te_mean=1888.148 K` vs `1890 K`, `Tout=876.524 K`.
  - `0.9 V`: `I_calc=343.717 A` vs `333 A`, `Te_mean=1918.673 K` vs `1910 K`, `Tout=876.612 K`.

Parameter scan notes:

- Increasing `R_EC` raises emitter temperature and current.
- Reducing coolant flow raises coolant outlet and collector temperatures; it can improve temperature agreement but may also increase current.
- Increasing `Rwire` suppresses current, but in the local scan it pushed emitter temperature too high and worsened the combined score.
- The selected candidate was the best combined result in the scanned grid when current MAE, emitter-temperature MAE, and convergence were considered together.

Remaining limitation: this is a physical-flow reduced-order model, not yet a full `TFEUnit` + `SystemManager` 2D solid conduction solve. It explicitly models coolant enthalpy rise and heat rejection, but radial/axial solid conduction is still represented by lumped station-wise resistances.

## 2026-07-04 axial-conduction physical-flow refinement

The physical-flow single-TFE model now includes an optional reduced axial-conduction smoothing term:

- `axial_conduction_smoothing` applies explicit 1D diffusion with zero end heat flux to emitter and collector temperature profiles.
- `axial_conduction_passes` controls the number of diffusion passes.
- The smoothing preserves the mean temperature before clipping and does not change coolant heat gain; coolant energy remains controlled by `sum(Q_to_coolant) = m_dot * cp * (Tout - Tin)`.

Verification:

- Unit tests confirm smoothing reduces axial spread while preserving mean temperature and coolant heat gain.
- `py_compile` passes for `nikolaev_physical_tfe_loop.py` and `nikolaev_physical_loop_runner.py`.

Current axial-conduction candidate:

- Run: `runs/20260704_nikolaev_physical_tfe_axial_smoothing_candidate_v2`
- Parameters: `Tin=770 K`, `m_dot=0.035 kg/s`, `cp=1000 J/kg/K`, `h=8500 W/m2/K`, `R_EC=0.380 K/W`, `Tcs=560 K`, `Rwire=0.0005 ohm`, `axial_conduction_smoothing=0.1`, `axial_conduction_passes=3`, relaxation `0.20`, tolerance `0.75 K`.
- All three Table 2 points converged with max outer iterations `20`.
- Metrics: current MAE `7.056 A`, max current error `10.715 A`, electric-power MAE `5.625 W`, emitter-temperature MAE `3.884 K`, max coolant energy residual `1.273e-11 W`.
- Point comparison:
  - `0.7 V`: `I_calc=418.285 A` vs `429 A`, `Te_mean=1878.940 K` vs `1880 K`, `Tout=879.130 K`.
  - `0.8 V`: `I_calc=375.351 A` vs `375 A`, `Te_mean=1889.244 K` vs `1890 K`, `Tout=876.548 K`.
  - `0.9 V`: `I_calc=343.104 A` vs `333 A`, `Te_mean=1919.834 K` vs `1910 K`, `Tout=876.638 K`.

Compared with the no-axial-conduction physical-flow candidate, the axial-conduction candidate slightly improves current MAE (`7.272 A -> 7.056 A`) and emitter-temperature MAE (`4.208 K -> 3.884 K`). The remaining dominant mismatch is the high-voltage `0.9 V` current, which is still high by about `10.1 A`.

## 2026-07-06 physical axial-conduction replacement

The previous `axial_conduction_smoothing` / `axial_conduction_passes` pair is now treated as a historical compatibility option for reproducing `runs/20260704_nikolaev_physical_tfe_axial_smoothing_candidate_v2`. It should not be used as the formal axial-conduction model in new validation runs.

New physical axial-conduction path:

- Enable with `axial_conduction_enabled=True` or runner flag `--enable-axial-conduction`.
- Emitter axial conductance uses `MoNb().conductivity(T_face) * emitter_cross_area_m2 / dz`.
- Collector axial conductance uses `Molybdenum().conductivity(T_face) * collector_cross_area_m2 / dz`.
- Both axial ends use zero axial heat-flux boundary conditions.
- The collector-to-coolant heat rejected by convection is integrated through the coolant stream; the reported coolant energy balance remains `m_dot * cp * (Tout - Tin)`.

First real-axial run:

- Run: `runs/20260706_nikolaev_physical_tfe_real_axial_conduction`
- Parameters retained from the previous candidate except axial treatment: `Tin=770 K`, `m_dot=0.035 kg/s`, `cp=1000 J/kg/K`, `h=8500 W/m2/K`, `R_EC=0.380 K/W`, `Tcs=560 K`, `Rwire=0.0005 ohm`, `axial_shape_amplitude=0.18`, `axial_conduction_enabled=True`, `axial_conduction_smoothing=0.0`.
- Metrics: current MAE `14.664 A`, max current error `23.395 A`, electric-power MAE `11.140 W`, emitter-temperature MAE `30.010 K`, max coolant energy residual `3.541e-7 W`.

Adjustment implication: do not return to empirical smoothing just to recover the previous fit. Retune the physical uncertain parameters instead. The most direct next scan should vary `R_EC`, coolant mass flow, collector convective coefficient, cesium reservoir temperature, and possibly `axial_shape_amplitude` around the new real-axial baseline.

## 2026-07-06 recalibrated real-axial candidate

Wire resistance has now been included in the real-axial recalibration.

Scan outputs:

- `runs/20260706_real_axial_sensitivity_round1`
- `runs/20260706_real_axial_joint_scan_round2`
- `runs/20260706_real_axial_wire_scan_round3`
- `runs/20260706_real_axial_fine_scan_round4`

Best formal candidate:

- Run: `runs/20260706_nikolaev_real_axial_recalibrated_candidate`
- Parameters: `Tin=770 K`, `m_dot=0.035 kg/s`, `cp=1000 J/kg/K`, `h=8500 W/m2/K`, `R_EC=0.370 K/W`, `Tcs=567 K`, `Rwire=0.00050 ohm`, `axial_conduction_enabled=True`, `axial_conduction_smoothing=0.0`, relaxation `0.20`, tolerance `0.75 K`.
- Metrics: current MAE `4.049 A`, max current error `7.399 A`, electric-power MAE `3.348 W`, emitter-temperature MAE `5.296 K`, max coolant energy residual `3.612e-7 W`.

Point comparison:

| V | I exp A | I calc A | I err A | Te exp K | Te mean K |
| ---:| ---:| ---:| ---:| ---:| ---:|
| 0.7 | 429 | 424.857 | -4.143 | 1880 | 1874.979 |
| 0.8 | 375 | 374.395 | -0.605 | 1890 | 1890.300 |
| 0.9 | 333 | 340.399 | +7.399 | 1910 | 1920.568 |

Adjustment notes:

- Do not assume lower `Rwire` is better. In the real-axial model, `Rwire=0` and `0.0002 ohm` over-predict current and distort the thermal feedback.
- The current best value remains close to `0.00050 ohm`.
- `R_EC=0.370 K/W` and `Tcs=567 K` are the current best combined setting for the material/geometry axial-conduction model.
- `axial_shape_amplitude` is not an effective tuning parameter in the current physical-flow real-axial path because the explicit heat-source distribution comes from `centered_heater_power_profile`, not from the old prescribed-temperature shape model.
