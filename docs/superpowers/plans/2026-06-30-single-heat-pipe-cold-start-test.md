# Single Heat Pipe Cold Start Test Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a focused single-heat-pipe cold-start test that validates the repository's pseudo-thermal-conductivity heat-pipe model against expected cold-start physics and the referenced PTC paper.

**Architecture:** Keep the production model unchanged. Add a `testModule` cold-start runner plus pytest-style checks around `HeatPipe2D`, `WickMaterial`, and `HPwithFin`; use short deterministic cases for CI-like validation and a longer benchmark case for literature comparison.

**Tech Stack:** Python 3.12 Conda environment at `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`, NumPy, existing `SystemManager`, `HeatPipe2D`, `HPwithFin`, `SodiumHP`, `SS316`, and CSV/NPZ outputs.

---

## Source Basis

User-provided paper:

```text
E:\文献阅读\会议论文\第十七届全国热管会议论文集\21-025 热管赝热导率模型及应用.pdf
```

Relevant points extracted with Ghostscript `txtwrite`:

- The paper defines a pseudo thermal conductivity (PTC) model that converts vapor flow plus phase-change heat transfer into an equivalent conduction process.
- The wick energy balance is written as internal energy change driven by wick axial conduction, pseudo-conduction heat flow, and external loss.
- The pseudo heat flow is mapped to a Fourier-like term; the resulting `k_pse` scales strongly with vapor pressure and latent heat and inversely with vapor viscosity and `T^3`.
- The paper stresses that high-temperature alkali-metal heat pipes during startup show a moving temperature front; an instant isothermal vapor-core assumption is not valid in the cold phase.
- The paper demonstrates startup simulation for a `0.6 m` nickel-alloy/sodium heat pipe under radiative heating, comparing PTC simulation with experiment and a detailed model.

Repository alignment:

- `Materials/Solids/WickMaterial.py` already implements structural conductivity plus pseudo-thermal axial conductivity.
- `Components/basicComponents/HeatPipe2D.py` can use anisotropic wick conductivity via `set_wick_conductivity_mode(True)`, making axial wick conductivity `conductivity_axial()` and radial wick conductivity `conductivity_radial()`.
- Existing historical cold-start prototype: `HeatPipe优化重构/test_single_hp_cold_start.py`.
- Existing reusable single-HP helper: `testModule/test_single_hp_fin_energy_conservation.py::build_hp_radiator()`.

## Test Strategy

Use three layers:

1. Material-level PTC validation: verify `SodiumHP` + `WickMaterial` produce finite, sharply temperature-dependent axial pseudo conductivity and sane structural/radial conductivity.
2. Minimal `HeatPipe2D` cold-start validation: isolate the heat pipe body without fins and verify temperature-front behavior, energy accounting, and numerical robustness.
3. Integrated `HPwithFin` cold-start validation: run the repository radiator-style single heat pipe with evaporator heating and condenser radiation/fin rejection.

Do not use the historical `HeatPipe优化重构/` script as the final entry point. Treat it as a reference prototype only.

## Proposed Files

- Create: `testModule/test_single_hp_cold_start_ptc.py`
  - Contains unit and integration tests for PTC material behavior, minimal heat-pipe startup, and short `HPwithFin` startup.
- Create: `testModule/run_single_hp_cold_start_ptc.py`
  - CLI runner for longer benchmark runs, CSV/NPZ history output, and optional paper-comparison output.
- Create: `testModule/SINGLE_HP_COLD_START_PTC_GUIDE.md`
  - Documents assumptions, commands, output fields, validation criteria, and paper mapping.
- Optional later: `testModule/reference_data/sodium_hp_cold_start_fig7_digitized.csv`
  - Digitized points from paper Fig. 7 if quantitative literature validation is required.

## Acceptance Criteria

- All short tests finish without non-finite temperatures, non-finite conductivity, or failed `SystemManager.step()`.
- Initial cold sodium starts below melt temperature, then passes through the apparent-heat-capacity melting interval around `371 K`.
- In the PTC-enabled case, a startup temperature front moves monotonically from evaporator toward condenser.
- PTC-enabled startup is materially faster and more uniform than a structural-conduction-only control case.
- Energy residual over each recorded interval remains bounded:
  - short isolated case: `abs(residual) <= max(1e-6 J, 0.005 * max(abs(Q_in) * dt, 1.0 J))`
  - integrated `HPwithFin` case: `abs(residual) <= max(5.0 J, 0.02 * max(abs(Q_in) * dt, 1.0 J))`
- Mesh/time sensitivity is acceptable:
  - halving `dt` changes selected probe temperatures by `< 3 K` after `20 s`
  - doubling axial nodes changes front location by `< 1 axial coarse cell` after `20 s`
- If Fig. 7 is digitized, simulated probe temperatures should match digitized experimental/model trends within:
  - temperature RMSE `< 30 K`, or
  - time-to-threshold error `< 15%` for thresholds such as `400 K`, `600 K`, `800 K`

---

### Task 1: Material-Level PTC Characterization

**Files:**
- Create/modify: `testModule/test_single_hp_cold_start_ptc.py`

- [ ] **Step 1: Add a test that samples sodium wick conductivity**

Test temperatures:

```python
temps = np.array([300.0, 360.0, 371.0, 450.0, 600.0, 800.0, 1000.0])
```

Assertions:

```python
assert np.all(np.isfinite(k_axial))
assert np.all(np.isfinite(k_radial))
assert np.all(k_axial >= 0.0)
assert np.all(k_radial >= 0.0)
assert k_axial[-1] > 10.0 * max(k_axial[0], 1.0)
assert np.all(k_axial >= k_radial)
```

- [ ] **Step 2: Verify against direct formula decomposition**

Use `WickMaterial.conductivity_structural(T)` and `WickMaterial.conductivity_pseudothermal(T)`.

Assertion:

```python
np.testing.assert_allclose(
    wick.conductivity_axial(temps),
    wick.conductivity_structural(temps) + wick.conductivity_pseudothermal(temps),
    rtol=1e-10,
    atol=1e-8,
)
```

- [ ] **Step 3: Run the material test**

Command:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m pytest testModule/test_single_hp_cold_start_ptc.py::test_sodium_wick_ptc_conductivity_temperature_response -q
```

Expected: pass in under a few seconds.

### Task 2: Minimal HeatPipe2D Cold-Start Case

**Files:**
- Create/modify: `testModule/test_single_hp_cold_start_ptc.py`

- [ ] **Step 1: Build a minimal `HeatPipe2D` fixture**

Use a no-fin, no-fluid-network heat pipe body:

- `L_eva = 0.06 m`
- `L_aba = 0.00 m`
- `L_con = 0.54 m` for a total `0.60 m` paper-style length
- `r_vapor = 8.5e-3 m`
- `r_in_wall = 9.0e-3 m`
- `r_out_wall = 11.0e-3 m`
- `porosity = 0.675`
- `n_wick = 2`
- `n_wall = 2`
- `n_eva = 6`
- `n_aba = 0`
- `n_con = 54`
- initial temperature `300 K`

Use `SS316` as the wall/wick structural proxy unless a nickel alloy material is added later. Record this as a model limitation in the guide.

- [ ] **Step 2: Configure startup boundaries**

Use radiative or power heating on `outer_eva`, not a near-Dirichlet boundary as the final validation condition.

Recommended deterministic short-test boundary:

```python
total_power_w = 150.0
weights = boundary.area / np.sum(boundary.area)
boundary.add_flux_condition(q_flux=total_power_w * weights)
```

Recommended paper-style benchmark boundary:

```python
boundary.add_dynamic_radiation_condition(
    emissivity=0.8,
    bare_area_array=boundary.area,
    T_env=864.0,
)
```

Use weak condenser radiation to `300 K` for short tests, then paper-style radiative boundary for benchmark runs.

- [ ] **Step 3: Enable the production PTC mode**

Use:

```python
hp.set_wick_conductivity_mode(True)
hp.set_face_conductance_mode("resistance_split_full")
hp.set_time_integrator("theta_implicit")
hp.set_theta_implicit_value(0.7)
hp.enable_frozen_property_correction = True
hp.max_outer_property_corrections = 3
hp.outer_property_tol = 1.0e-4
```

- [ ] **Step 4: Define front diagnostics**

For each recorded time:

- `T_outer_wall_z = T_2d[-1, :]`
- `T_wick_mean_z = mean(T_2d[:n_wick, :], axis=0)`
- `k_axial_wick_mean_z = mean(hp._k_axial_2d_view[:n_wick, :], axis=0)`
- `front_by_temperature = first z where T_wick_mean_z < 371.0 + 10.0`
- `front_by_ptc = first z where k_axial_wick_mean_z < 10.0 * k_radial_wick_mean_z`

Expected qualitative behavior:

- early time: only evaporator-side cells warm strongly
- mid time: front moves toward condenser
- late time: axial conductivity increases over a larger fraction of the pipe

- [ ] **Step 5: Add structural-only control**

For the control case, replace axial wick conductivity with structural conductivity only inside the fixture:

```python
hp.wick_mat.conductivity_axial = hp.wick_mat.conductivity_structural
hp.wick_mat.conductivity = hp.wick_mat.conductivity_structural
hp.wick_mat.invalidate_lookup_table()
hp._property_cache_initialized = False
```

Expected:

- PTC-enabled case reaches a downstream probe temperature threshold earlier than structural-only.
- PTC-enabled case has lower evaporator-to-condenser axial gradient after startup develops.

- [ ] **Step 6: Run short deterministic test**

Command:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m pytest testModule/test_single_hp_cold_start_ptc.py::test_minimal_heatpipe2d_cold_start_front_moves_with_ptc -q
```

Expected: pass in less than a minute.

### Task 3: Integrated HPwithFin Cold-Start Regression

**Files:**
- Create/modify: `testModule/test_single_hp_cold_start_ptc.py`
- Reuse: `testModule/test_single_hp_fin_energy_conservation.py`

- [ ] **Step 1: Build from existing helper**

Use `build_hp_radiator()` from `testModule/test_single_hp_fin_energy_conservation.py`, then reset it to cold initial state:

```python
hp_radiator.hp.T.fill(300.0)
hp_radiator.hp.current_time = 0.0
hp_radiator.last_fin_temperature.fill(300.0)
hp_radiator.hp.set_wick_conductivity_mode(True)
hp_radiator.hp.set_face_conductance_mode("resistance_split_full")
hp_radiator.hp.set_time_integrator("theta_implicit")
hp_radiator.hp.set_theta_implicit_value(0.7)
hp_radiator.hp.enable_frozen_property_correction = True
```

- [ ] **Step 2: Replace evaporator boundary**

Clear `outer_eva` and use a finite physical heating boundary:

```python
boundary = hp_radiator.hp.boundaries["outer_eva"]
boundary.conditions.clear()
boundary.clear_boundary_conditions()
area = boundary.area
boundary.add_flux_condition(q_flux=200.0 * area / np.sum(area))
```

- [ ] **Step 3: Run with dummy network through `SystemManager`**

Use `create_dummy_fluid_network()` from the existing energy-conservation test.

Short regression:

- `t_end = 20 s`
- `dt = 0.05 s`
- `inner_iter = 1`

Assertions:

- every `system.step()` returns successfully
- `T_hp_min >= 0 K`
- all heat rejection terms are finite
- `k_wick_max` increases from its initial value
- `T_eva_wall_mean > T_con_wall_mean` during early startup

- [ ] **Step 4: Run the integrated short test**

Command:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m pytest testModule/test_single_hp_cold_start_ptc.py::test_hpwithfin_cold_start_short_regression -q
```

Expected: pass in a few minutes.

### Task 4: Long Benchmark Runner

**Files:**
- Create: `testModule/run_single_hp_cold_start_ptc.py`

- [ ] **Step 1: Add CLI options**

Required options:

```text
--case minimal|hpwithfin
--heating-mode power|radiation
--total-power-w
--heater-temp-k
--initial-temp-k
--t-end
--dt
--record-interval
--out-dir
--structural-only-control
```

- [ ] **Step 2: Record CSV fields**

Minimum fields:

```text
t_s
Q_eva_in_w
Q_aba_out_w
Q_con_out_w
T_min_k
T_mean_k
T_max_k
T_eva_wall_mean_k
T_con_wall_mean_k
T_probe_25pct_k
T_probe_50pct_k
T_probe_75pct_k
k_wick_axial_min_w_m_k
k_wick_axial_mean_w_m_k
k_wick_axial_max_w_m_k
front_temperature_m
front_ptc_m
energy_residual_j
```

- [ ] **Step 3: Record NPZ arrays**

Minimum arrays:

```text
time
T_2d_history
T_outer_wall_history
k_wick_axial_history
k_wick_radial_history
q_eva_history
q_con_history
x_centers
y_centers
```

- [ ] **Step 4: Run smoke benchmark**

Command:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule/run_single_hp_cold_start_ptc.py --case minimal --heating-mode power --total-power-w 150 --initial-temp-k 300 --t-end 20 --dt 0.05 --record-interval 1 --out-dir testModule/single_hp_cold_start_ptc_smoke
```

Expected:

- writes CSV and NPZ
- no failed step
- front diagnostics finite or explicitly marked beyond pipe length

- [ ] **Step 5: Run paper-style benchmark**

Command:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule/run_single_hp_cold_start_ptc.py --case minimal --heating-mode radiation --heater-temp-k 864 --initial-temp-k 300 --t-end 600 --dt 0.1 --record-interval 10 --out-dir testModule/single_hp_cold_start_ptc_paper_style
```

Expected:

- reproduces moving startup front
- downstream probes show delayed heating
- final output suitable for comparing against Fig. 7 after digitization

### Task 5: Quantitative Literature Validation

**Files:**
- Optional create: `testModule/reference_data/sodium_hp_cold_start_fig7_digitized.csv`
- Modify: `testModule/run_single_hp_cold_start_ptc.py`
- Modify: `testModule/SINGLE_HP_COLD_START_PTC_GUIDE.md`

- [ ] **Step 1: Digitize Fig. 7**

Use the rendered PDF page or an external digitizer to extract time-temperature curves from Fig. 7.

CSV schema:

```text
curve_name,t_s,T_k
experiment_probe_1,0.0,300.0
experiment_probe_1,10.0,315.0
ptc_model_probe_1,0.0,300.0
```

Store only digitized numeric points if redistribution is allowed. Do not commit PDF images or copied paper figures.

- [ ] **Step 2: Add comparison postprocessor**

For each named probe, interpolate simulation to digitized times and report:

```text
rmse_k
max_abs_error_k
threshold_400k_time_error_pct
threshold_600k_time_error_pct
threshold_800k_time_error_pct
```

- [ ] **Step 3: Acceptance**

Before geometry/material calibration:

- qualitative front order and threshold timing are the main criteria

After geometry/material calibration:

- `rmse_k < 30 K` or threshold time error `< 15%`

### Task 6: Documentation and Verification Matrix

**Files:**
- Create: `testModule/SINGLE_HP_COLD_START_PTC_GUIDE.md`
- Modify if behavior changes: `Components/BASICCOMPONENTS_DETAILED_INTRO.md`

- [ ] **Step 1: Document model limitations**

Include:

- PTC is an equivalent conduction model, not an explicit vapor/liquid flow solver.
- `HeatPipe2D` cold-start behavior depends strongly on `WickMaterial` vapor pressure, viscosity, latent heat, and apparent heat capacity.
- `HPwithFin` fins are reduced-order quasi-steady branches, not explicit 2D fins.
- Paper-style benchmark geometry is only fully quantitative after Fig. 7 and all test fixture parameters are digitized.

- [ ] **Step 2: Document required commands**

Include the exact Conda Python command:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m pytest testModule/test_single_hp_cold_start_ptc.py -q
```

Include the long-run command from Task 4.

- [ ] **Step 3: Run full verification set**

Commands:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m pytest testModule/test_single_hp_cold_start_ptc.py -q
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule/run_single_hp_cold_start_ptc.py --case minimal --heating-mode power --total-power-w 150 --initial-temp-k 300 --t-end 20 --dt 0.05 --record-interval 1 --out-dir testModule/single_hp_cold_start_ptc_smoke
```

Expected:

- pytest passes
- runner writes CSV/NPZ
- guide documents where output artifacts are produced and says not to commit run-output directories by default

## Implementation Notes

- Prefer physical finite heat input for validation. Keep the old near-Dirichlet `T_ext=864 K, R_ext=1e-8` path only as a stress/smoke option.
- Use `BoundaryRegion.add_flux_condition()` with per-boundary-cell power in W, not W/m2.
- Use `HeatPipe2D.set_wick_conductivity_mode(True)` for production PTC cold-start tests.
- Use `resistance_split_full` for cold-start validation because boundary half-cell thermal resistance matters when the front is sharp.
- Keep outputs under `testModule/single_hp_cold_start_ptc_*`; do not commit CSV/NPZ/PNG run products unless explicitly requested.
