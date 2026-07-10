# Luchau V Single TFE Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a standalone full-length single thermionic fuel element model in `testModule/Full_Loop_Cases_VV/Luchau_V`, with a parameterized 30 cm central heater and single-TFE ThermoCalc fixed-voltage electrical coupling.

**Architecture:** The case creates one `TFEUnit` using the geometry, mesh, materials, gas gaps, coolant, and 37-node axial allocation from `testModule/Full_Loop_Cases/common_core_builder.py`. A separate helper configures `ThermoCalcModel(n_elements=1, n_nodes=37)` from the same axial faces and electrode geometry, using Cs pressure converted to `Tcs` and a caller-provided fixed voltage. Tests verify parameter validation, axial heating, geometry mapping, and ThermoCalc input setup through a fake model.

**Tech Stack:** Python 3.12 in `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`, `unittest`, TASTIN `TFEUnit`, `HydraulicNetwork`, `ThermoCalcWrapper`.

## Global Constraints

- Use `SodiumPotassium78` coolant.
- Coolant inlet temperature is `727 K`.
- Single-channel mass flow is `1.3 / 37 kg/s`.
- Full TFE length is `0.065 m + 0.377 m + 0.065 m = 0.507 m`.
- Axial nodes are `6 + 25 + 6 = 37`.
- Center heater length is `0.30 m`; total heater power is a required parameter.
- Cs pressure is `0.4 torr`; convert to `Tcs` using `Pcs_torr = 2.45e8 / sqrt(Tcs) * exp(-8910 / Tcs)`.
- Electrical mode is single-TFE `fixed_u`; target voltage is a required parameter.
- Do not wrap this in `ReactorCore` or add global moderator/barrel/reflector rings.

---

### Task 1: Single TFE Model API

**Files:**
- Create: `testModule/Full_Loop_Cases_VV/Luchau_V/luchau_single_tfe_model.py`
- Test: `testModule/Full_Loop_Cases_VV/Luchau_V/test_luchau_single_tfe_model.py`

**Interfaces:**
- Produces: `LuchauSingleTFEConfig`, `build_center_heater_profile(config)`, `build_luchau_single_tfe(config)`.
- Produces: build dictionary keys `system`, `network`, `tfe`, `channel`, `inlet`, `outlet`, `inlet_junction`, `outlet_junction`, `node_lengths_m`, `axial_faces_m`, `heater_profile`, `flow_area_m2`, `hydraulic_diam_m`.

- [ ] **Step 1: Write failing tests**

Tests import the not-yet-created module and assert:

```python
cfg = LuchauSingleTFEConfig(thermal_power_w=3000.0, target_voltage_v=1.0)
profile = build_center_heater_profile(cfg)
assert profile.shape == (37,)
assert abs(profile.sum() - 1.0) < 1e-12
assert np.count_nonzero(profile) > 0
assert abs(sum(length for length, weight in zip(build_node_lengths(), profile) if weight > 0.0) - 0.30) <= max(build_node_lengths())
```

- [ ] **Step 2: Run test to verify it fails**

Run: `& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_VV.Luchau_V.test_luchau_single_tfe_model`

Expected: fail with `ModuleNotFoundError` or missing symbols.

- [ ] **Step 3: Implement minimal model**

Create config, axial geometry helpers, material/gap setup, a fixed-flow inlet junction, outlet pressure boundary, and one `TFEUnit`. Use `tfe.update_neutronic_power(config.thermal_power_w)` after `SystemManager.initialize_system()` so the fuel `Q_source` is populated from the central heater profile.

- [ ] **Step 4: Run test to verify it passes**

Run the same unittest command. Expected: pass.

### Task 2: ThermoCalc Setup Helper

**Files:**
- Modify: `testModule/Full_Loop_Cases_VV/Luchau_V/luchau_single_tfe_model.py`
- Test: `testModule/Full_Loop_Cases_VV/Luchau_V/test_luchau_single_tfe_model.py`

**Interfaces:**
- Produces: `cesium_pressure_from_tcs(tcs_k)`, `tcs_from_cesium_pressure(pcs_torr)`, `configure_luchau_thermocalc(thermo_model, build, config)`.

- [ ] **Step 1: Write failing tests**

Tests use a fake ThermoCalc object with `_input_data`, `setup_circuit_mode`, `set_temperatures`, and `set_tcs`, then assert `fixed_u`, target voltage, 37-node arrays, `d_gap=0.5 mm`, electrode areas from TFE geometry, and finite `Tcs` for `0.4 torr`.

- [ ] **Step 2: Run test to verify it fails**

Run: `& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_VV.Luchau_V.test_luchau_single_tfe_model`

Expected: fail because `configure_luchau_thermocalc` is missing.

- [ ] **Step 3: Implement minimal helper**

Read emitter/collector mean axial temperatures from `tfe.solids["emitter"].T` and `tfe.solids["collector"].T`, set ThermoCalc geometry arrays from axial faces and radii, configure single fixed-voltage mode, and set uniform `Tcs`.

- [ ] **Step 4: Run test to verify it passes**

Run the same unittest command. Expected: pass.

### Task 3: Runner And Smoke

**Files:**
- Create: `testModule/Full_Loop_Cases_VV/Luchau_V/run_luchau_single_tfe.py`
- Modify: `testModule/Full_Loop_Cases_VV/Luchau_V/test_luchau_single_tfe_model.py`

**Interfaces:**
- Produces CLI args `--thermal-power-w`, `--target-voltage-v`, `--duration-s`, `--dt-s`, `--skip-thermocalc-calc`.

- [ ] **Step 1: Write failing test**

Test that the runner parser accepts required parameters and rejects missing `--thermal-power-w` or `--target-voltage-v`.

- [ ] **Step 2: Run test to verify it fails**

Run the unittest command. Expected: fail because runner is missing.

- [ ] **Step 3: Implement minimal runner**

Build the model, configure ThermoCalc, optionally call `calculate(verbose=False)`, optionally advance `SystemManager.step(dt)` for a short duration, and write a JSON summary under `Luchau_V/runs/<label>/summary.json`.

- [ ] **Step 4: Run tests and py_compile**

Run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_VV.Luchau_V.test_luchau_single_tfe_model
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile testModule\Full_Loop_Cases_VV\Luchau_V\luchau_single_tfe_model.py testModule\Full_Loop_Cases_VV\Luchau_V\run_luchau_single_tfe.py testModule\Full_Loop_Cases_VV\Luchau_V\test_luchau_single_tfe_model.py
```

Expected: both commands pass.

## Self-Review

- Spec coverage: covers full-length TFEUnit, 37-node axial mesh, center 30 cm heater, NaK inlet/flow, Cs pressure to ThermoCalc, single fixed-voltage circuit, and tests under `Luchau_V`.
- Placeholder scan: total heater power and target voltage are intentionally required runtime parameters per user instruction; no implementation step leaves behavior unspecified.
- Type consistency: all produced functions and dictionary keys are named before use.
