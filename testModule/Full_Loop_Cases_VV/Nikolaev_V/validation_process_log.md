# Nikolaev_V validation process log

## 2026-07-03 initial implementation

- Read local PDF `e:\文献阅读\nikolaev1995.pdf`.
- Identified table-level validation targets from Nikolaev 1995.
- Created source data, compact single-TFE model, validation runner, parameter guide, and source extraction notes.
- Figure 4 VAC curves are not digitized in this implementation.

## 2026-07-03 baseline run

Run directory: `runs/20260703_nikolaev_table_validation_baseline/`.

Key metrics:

- Table 1 max abs error: `0`.
- Table 2 current MAE: `0.253968 A`.
- Table 2 electric power MAE: `0.2 W`.
- Table 2 emitter temperature MAE: `0 K`.
- Table 2 efficiency MAE: `0.0256678 percentage point`.
- Table 3 max fuel-temperature error: `0 K`.
- Table 4 max capillary-diameter error: `0 mm`.

## 2026-07-03 ThermoCalc fixed-voltage path check

- Root-cause check: the table-level path in `nikolaev_single_tfe_model.py` computes current as `nominal_output_power_w / voltage_v`; therefore the very small current error in the baseline report is table reconstruction, not independent prediction.
- Added ThermoCalc path tests to ensure the new runner records current from ThermoCalc `Iout`, not from the table power formula.
- Ran `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe -m unittest discover -s testModule\Full_Loop_Cases_VV\Nikolaev_V -p "test_*.py"`; result: `Ran 8 tests`, `OK`.
- Ran `nikolaev_thermocalc_runner.py --run-id 20260703_nikolaev_thermocalc_baseline`; result: current MAE `135.717 A`, max current error `306.731 A`, all three ThermoCalc solves finite and converged.
- Ran balanced candidate `R_EC=0.260 K/W`, `Tcs=560 K`, `Rwire=0.0005 ohm`; result: current MAE `23.489 A`, max current error `45.109 A`, but average emitter temperatures `1936-1962 K` are above the Table 2 values.
- A coarse unconstrained current fit found `R_EC=0.300 K/W`, `Tcs=560 K`, `Rwire=0.001 ohm` with current MAE `11.847 A`; this is rejected as a validation setting because average emitter temperature is about `2100 K`.
- Current conclusion: the ThermoCalc electrical path is now real, but the thermal model is still only prescribed-temperature/fixed-network input. It has not yet closed electron cooling, electron heat transport, or Joule heat feedback into the thermal network.

## 2026-07-04 closed thermoelectric feedback validation

- Added `nikolaev_thermoelectric_closed_loop.py` with a local thermal feedback loop around ThermoCalc.
- Added tests proving that closed-loop current is not table-derived, thermal iterations feed electronic heat terms back into temperatures, and Joule heat is taken from `joulePowerE/C` rather than reconstructed from `UE/UC`.
- Added `nikolaev_closed_loop_runner.py` for Table 2 closed-loop validation reports.
- Ran tests: `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe -m unittest discover -s testModule\Full_Loop_Cases_VV\Nikolaev_V -p "test_*.py"`; result: `Ran 10 tests`, `OK`.
- Baseline closed-loop run `20260704_nikolaev_closed_loop_baseline`: all three points converged, but final emitter temperatures were too low. Current MAE `85.292 A`; emitter-temperature MAE `163.277 K`.
- Balanced candidate run `20260704_nikolaev_closed_loop_balanced_candidate`: `R_EC=0.340 K/W`, `Tcs=560 K`, `Rwire=0.0005 ohm`, `R_CB=0.010 K/W`, relaxation `0.20`, tolerance `0.75 K`; all three points converged in 27 outer iterations. Current MAE `6.462 A`; max current error `12.175 A`; emitter-temperature MAE `24.641 K`.
- Remaining limitation: the network is local and reduced order. It captures ThermoCalc electron cooling/heating and Joule heat feedback, but it is not yet the full TFEUnit/Core/SystemManager solid conduction model.

## 2026-07-04 physical thermal-hydraulic TFE loop

- Added `nikolaev_physical_tfe_loop.py`, which advances from fixed collector-boundary closure to an explicit heater-to-coolant flow path.
- Added `nikolaev_physical_loop_runner.py` to produce Table 2 coupled thermal/electrical/coolant reports.
- Added tests proving the model tracks coolant inlet/outlet temperatures, coolant heat gain, energy balance, flow sensitivity, and ThermoCalc-derived current.
- Ran tests: `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe -m unittest discover -s testModule\Full_Loop_Cases_VV\Nikolaev_V -p "test_*.py"`; result: `Ran 12 tests`, `OK`.
- Baseline run `20260704_nikolaev_physical_tfe_baseline`: energy residual below `1.2e-11 W`, but current MAE `25.023 A`, emitter-temperature MAE `64.603 K`, and one point did not converge within 40 outer iterations.
- Coarse parameter scan over `R_EC`, `Tcs`, coolant mass flow, and collector heat-transfer coefficient identified `R_EC=0.380 K/W`, `Tcs=560 K`, `m_dot=0.035 kg/s`, `h=8500 W/m2/K`, `Rwire=0.0005 ohm` as the best combined candidate.
- Formal 50-node candidate run `20260704_nikolaev_physical_tfe_balanced_candidate`: all points converged; current MAE `7.272 A`; max current error `10.717 A`; electric-power MAE `5.840 W`; emitter-temperature MAE `4.208 K`; max coolant energy residual `5.912e-12 W`.
- Remaining limitation: the new path computes coolant heating explicitly, but it is still reduced-order radial/axial conduction rather than a full `TFEUnit`/`SystemManager` 2D solid conduction solve.

## 2026-07-04 axial-conduction physical-flow refinement

- Added optional axial-conduction smoothing to `nikolaev_physical_tfe_loop.py` and exposed it in `nikolaev_physical_loop_runner.py`.
- Added tests verifying axial smoothing preserves mean temperature, reduces axial spread, and does not alter coolant heat gain.
- Initial `0.1/3-pass` 50-node run `20260704_nikolaev_physical_tfe_axial_smoothing_candidate_v2` converged for all three Table 2 points.
- A local smoothing scan indicated weak smoothing is preferable; the formal 50-node `0.1/3-pass` candidate remains slightly better than the formal `0.05/2-pass` run by combined current/temperature metrics.
- Current best run: `runs/20260704_nikolaev_physical_tfe_axial_smoothing_candidate_v2`.
- Metrics: current MAE `7.056 A`; max current error `10.715 A`; electric-power MAE `5.625 W`; emitter-temperature MAE `3.884 K`; max coolant energy residual `1.273e-11 W`; all ThermoCalc and physical loops converged.
- Remaining mismatch: the `0.9 V` current is still high by about `10.1 A`, despite the average emitter temperature being within about `9.8 K` of Table 2.

## 2026-07-06 real material axial-conduction replacement

- Preserved the previous calibrated smoothing run `runs/20260704_nikolaev_physical_tfe_axial_smoothing_candidate_v2` unchanged.
- Added a real reduced 1D axial-conduction path to `nikolaev_physical_tfe_loop.py` behind `axial_conduction_enabled=True`.
- Emitter axial conductance now uses `Materials.Solids.MoNb.MoNb.conductivity(T) * A_emitter / dz` and `NikolaevTfeGeometry.emitter_cross_area_m2`.
- Collector axial conductance now uses `Materials.Solids.Molybdenum.Molybdenum.conductivity(T) * A_collector / dz` and `NikolaevTfeGeometry.collector_cross_area_m2`.
- Axial boundary condition is zero axial heat flux at both ends. Coolant heat gain is still computed from the convective heat removed by the collector and integrated as `m_dot * cp * (Tout - Tin)`.
- Verification: `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe -m unittest testModule.Full_Loop_Cases_VV.Nikolaev_V.test_nikolaev_model` passed with `Ran 12 tests`, `OK`; py_compile also passed for the modified runner, loop model, and tests.
- New run: `runs/20260706_nikolaev_physical_tfe_real_axial_conduction`.
- Parameters kept from the prior candidate except for axial treatment: `Tin=770 K`, `m_dot=0.035 kg/s`, `cp=1000 J/kg/K`, `h=8500 W/m2/K`, `R_EC=0.380 K/W`, `Tcs=560 K`, `Rwire=0.0005 ohm`, `axial_shape_amplitude=0.18`, `axial_conduction_enabled=True`, `axial_conduction_smoothing=0.0`.
- Metrics: current MAE `14.664 A`; max current error `23.395 A`; electric-power MAE `11.140 W`; emitter-temperature MAE `30.010 K`; max coolant energy residual `3.541e-7 W`; all ThermoCalc and physical loops converged.
- Interpretation: replacing the empirical smoothing with actual MoNb/Mo axial conduction is more physically defensible, but with the old calibrated parameters it over-raises the mean emitter temperatures and under-predicts current at all three Table 2 voltages. The next calibration should not re-enable smoothing; it should retune physical parameters such as `R_EC`, coolant mass flow, collector heat-transfer coefficient, cesium reservoir temperature, and possibly the axial power profile.

## 2026-07-06 real-axial recalibration with wire resistance

- User suggested including wire resistance in the real-axial recalibration. This was tested explicitly.
- Sensitivity round 1: `runs/20260706_real_axial_sensitivity_round1` varied `R_EC`, coolant flow, collector `h`, `Tcs`, and `axial_shape_amplitude` one at a time around the unretuned real-axial baseline.
  - `R_EC` and `Tcs` were dominant.
  - `axial_shape_amplitude` did not affect the physical-loop result in the current real-axial path because the explicit heat-source profile is supplied by `centered_heater_power_profile`; this parameter only affects the older ThermoCalc prescribed-temperature build path.
- Joint round 2: `runs/20260706_real_axial_joint_scan_round2` scanned `R_EC=0.38-0.41 K/W`, `Tcs=565/570/575 K`, and `m_dot=0.030/0.035 kg/s`.
  - Best combined point from that grid was `R_EC=0.380 K/W`, `Tcs=565 K`, `m_dot=0.035 kg/s`, with current MAE `5.356 A` and emitter-temperature MAE `17.045 K`.
- Wire round 3: `runs/20260706_real_axial_wire_scan_round3` scanned `R_EC=0.37/0.38/0.39 K/W`, `Tcs=560/565 K`, and `Rwire=0/0.0002/0.0005/0.0008 ohm`.
  - Very low wire resistance over-predicted current and worsened the coupled thermal feedback.
  - Best grid point was `R_EC=0.370 K/W`, `Tcs=565 K`, `Rwire=0.0005 ohm`, with current MAE `5.735 A` and emitter-temperature MAE `5.916 K`.
- Fine round 4: `runs/20260706_real_axial_fine_scan_round4` scanned around the wire-round optimum: `R_EC=0.365/0.370/0.375 K/W`, `Tcs=563/565/567 K`, and `Rwire=0.00045/0.00050/0.00055 ohm`.
  - Best combined score was `R_EC=0.370 K/W`, `Tcs=567 K`, `Rwire=0.00050 ohm`.
- Formal candidate run: `runs/20260706_nikolaev_real_axial_recalibrated_candidate`.
- Formal candidate parameters: `Tin=770 K`, `m_dot=0.035 kg/s`, `cp=1000 J/kg/K`, `h=8500 W/m2/K`, `R_EC=0.370 K/W`, `Tcs=567 K`, `Rwire=0.00050 ohm`, `axial_conduction_enabled=True`, `axial_conduction_smoothing=0.0`, relaxation `0.20`, tolerance `0.75 K`.
- Metrics: current MAE `4.049 A`; max current error `7.399 A`; electric-power MAE `3.348 W`; emitter-temperature MAE `5.296 K`; max coolant energy residual `3.612e-7 W`; all ThermoCalc and physical loops converged.
- Point comparison:
  - `0.7 V`: `I_calc=424.857 A` vs `429 A`, `Te_mean=1874.979 K` vs `1880 K`, `Tout=878.923 K`.
  - `0.8 V`: `I_calc=374.395 A` vs `375 A`, `Te_mean=1890.300 K` vs `1890 K`, `Tout=876.579 K`.
  - `0.9 V`: `I_calc=340.399 A` vs `333 A`, `Te_mean=1920.568 K` vs `1910 K`, `Tout=876.734 K`.
- Comparison with previous best smoothing candidate: current MAE improved from `7.056 A` to `4.049 A`, max current error from `10.715 A` to `7.399 A`, and electric-power MAE from `5.625 W` to `3.348 W`. Emitter-temperature MAE worsened slightly from `3.884 K` to `5.296 K`, but the model is now physically more defensible because axial conduction is material/geometry based.
