# V14 210 kW Fast-Steady Temporary Case Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build and run an isolated fixed-210-kW case that accelerates only solid thermal storage by a factor of 100, then restores physical heat capacities for confirmation.

**Architecture:** Add one runner and one small test inside the temporary case folder. The runner reuses the existing V14 debug assembly, wraps materials only on the returned in-memory solid objects, advances with the existing fixed-power step path, and records an objective temperature-window convergence result.

**Tech Stack:** Python 3.12, NumPy, unittest, existing TASTIN `SystemManager` and V14 10 kW builders.

## Global Constraints

- Use `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`.
- Do not modify files outside `testModule/Full_Loop_Cases_10kW/V14_210kW_fast_steady_temp/`.
- Keep external heat disabled, core power at 210000 W, NaK and potassium working-fluid heat capacities physical.
- Scale core solids, collector-ring solids, heat-pipe walls, and wick solid skeleton by 0.01.
- Write calculation outputs below this folder's `runs/` directory.

---

### Task 1: Isolated heat-capacity scaler and runner

**Files:**
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_fast_steady_temp/__init__.py`
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_fast_steady_temp/run_fast_steady.py`
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_fast_steady_temp/test_fast_steady.py`
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_fast_steady_temp/runs/.gitkeep`

**Interfaces:**
- Consumes: `build_debug_case(DebugRunConfig)`, `_apply_fixed_core_power(build, power_w)`, `collect_metrics(build, stage_index, dt_s)`, and `load_baseline_debug_config(runtime)`.
- Produces: `ScaledHeatCapacityMaterial`, `scale_system_solids(system, scale)`, and `run_fast_steady(config)`.

- [ ] **Step 1: Write the failing material-wrapper test**

```python
def test_scaled_material_changes_only_heat_capacity():
    base = FakeMaterial()
    scaled = ScaledHeatCapacityMaterial(base, 0.01)
    assert scaled.heat_capacity(800.0) == 0.01 * base.heat_capacity(800.0)
    assert scaled.conductivity(800.0) == base.conductivity(800.0)
    assert scaled.density(800.0) == base.density(800.0)
```

- [ ] **Step 2: Run the test and verify the missing implementation failure**

```powershell
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' -m unittest testModule.Full_Loop_Cases_10kW.V14_210kW_fast_steady_temp.test_fast_steady
```

Expected: import failure for `ScaledHeatCapacityMaterial`.

- [ ] **Step 3: Implement the scaler**

```python
class ScaledHeatCapacityMaterial:
    def __init__(self, base, scale):
        if not math.isfinite(scale) or scale <= 0.0:
            raise ValueError('heat-capacity scale must be finite and positive')
        self.base = base
        self.scale = float(scale)

    def __getattr__(self, name):
        return getattr(self.base, name)

    def heat_capacity(self, temperature):
        return self.scale * self.base.heat_capacity(temperature)


def scale_system_solids(system, scale):
    scaled_names = []
    for name, solid in system.solid_components.items():
        if isinstance(solid, HeatPipe2D):
            solid.wall_mat = ScaledHeatCapacityMaterial(solid.wall_mat, scale)
            solid.material = solid.wall_mat
            solid.wick_mat.solid = ScaledHeatCapacityMaterial(solid.wick_mat.solid, scale)
            solid._property_cache_initialized = False
        else:
            solid.material = ScaledHeatCapacityMaterial(solid.material, scale)
        solid._update_properties()
        scaled_names.append(name)
    return scaled_names
```

- [ ] **Step 4: Implement the minimal fixed-power run loop**

`run_fast_steady(config)` must load the sibling configuration, call `build_debug_case`, apply `scale_system_solids` when the scale differs from one, refresh couplers with a positive `dt`, call `_apply_fixed_core_power` before and after every `system.step`, capture all fluid and solid temperature arrays at the start and end of the final convergence window, and write `history.csv`, `run_config.json`, `run_summary.json`, and `stage_01_restart.npz`.

- [ ] **Step 5: Verify tests and syntax**

```powershell
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' -m unittest testModule.Full_Loop_Cases_10kW.V14_210kW_fast_steady_temp.test_fast_steady
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' -m py_compile testModule\Full_Loop_Cases_10kW\V14_210kW_fast_steady_temp\run_fast_steady.py
```

Expected: test passes and `py_compile` exits zero.

- [ ] **Step 6: Commit the isolated runner**

```powershell
git add -- testModule/Full_Loop_Cases_10kW/V14_210kW_fast_steady_temp
git commit -m 'add isolated 210kW fast-steady runner'
```

### Task 2: Accelerated and physical-capacity calculations

**Files:**
- Runtime outputs only: `testModule/Full_Loop_Cases_10kW/V14_210kW_fast_steady_temp/runs/`

**Interfaces:**
- Consumes: `run_fast_steady.py --restart-in --output-dir --duration --dt --heat-capacity-scale`.
- Produces: accelerated and physical-capacity restart chains plus convergence summaries.

- [ ] **Step 1: Run a 0.05 s accelerated smoke**

Use the current best 210 kW restart, `--heat-capacity-scale 0.01`, `--duration 0.05`, and `--dt 0.05`. Expected: finite temperatures, fixed 210000 W, external heat absent, and a restart written.

- [ ] **Step 2: Run accelerated stages**

Run 100 accelerated seconds with a 10 s convergence window and 0.05 K tolerance. If not converged, continue from the generated restart in further 100 s stages without changing parameters.

- [ ] **Step 3: Restore physical capacities**

Load the accepted accelerated restart with `--heat-capacity-scale 1.0` and run 100 s at fixed 210 kW. Require finite state and inspect the final 10 s maximum temperature change.

- [ ] **Step 4: Verify zero-reactivity handoff**

Use the physical-capacity confirmation restart with the existing `V14_210kW_reactivity_control` runner for 10 s. Require initial effective feedback within `1e-12`, zero external and drum reactivity, and less than 1% power drift.

- [ ] **Step 5: Report results**

Report accelerated-stage convergence, physical-capacity confirmation drift, zero-reactivity power drift, exact restart paths, and any TEC convergence warning. Do not commit runtime outputs.
