# V14 全堆氦气瞬时完全失压事故算例实施计划

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 在 V14 210 kW 正常热容稳态上建立全堆 58 根 TFE 氦气隙同时瞬时失去气体导热的事故算例，保留辐射、温度反应性反馈、点堆动力学和衰变热，并完成 0.1 s smoke 与首轮 100 s 计算。

**Architecture:** 新 runner 复用现有 `V14_210kW_reactivity_control` 的稳态加载、固定功率到点堆的交接和反应性诊断，不修改公共 `GapCouple2D`、`TFEUnit` 或已有算例。事故层只定位 5 个 `collector_iclad_gap`，在相对时间零点把其 `k_gas` 置零；事故状态写入输出 `run_config.json`，续算时据此重新施加零气体导热。

**Tech Stack:** Python 3.12、NumPy、标准库 `unittest/csv/json/dataclasses/pathlib`、现有 TASTIN `SystemManager/ReactorCore/GapCouple2D/PointReactor`。

## Global Constraints

- 所有 Python 命令必须使用 `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`。
- 新算例只能写入 `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/`；不得改变已有反应性控制算例和公共物理组件。
- 初始 restart 固定为正常热容长算状态 `checkpoint_t013864s.npz`，绝对时间约 `13864.2 s`。
- 全部 5 个代表性 TFE 同时失压，名称固定为 `Center/Ring1/Ring2/Ring3/Ring4`，倍率固定为 `1/6/9/18/24`。
- 事故模型固定为 `h_He: 5678 -> 0 W/(m2*K)`；只取消气体导热，保留间隙辐射。
- 外界反应性为 0，控制鼓关闭，固定功率源关闭，外热流关闭，主回路目标流量保持 `2.46 kg/s`。
- 首轮运行 `100 s`，固定步长 `0.05 s`，记录间隔 `0.1 s`，checkpoint 间隔 `10 s`。
- 温度限值：通道壁 `1058 K`、芯块 `2700 K`、接收极 `1023 K`、慢化剂 `930 K`、反射层 `1000 K`。
- 通道壁最高温度取全部 TFE 的 `inner_clad` 与 `outer_clad` 全部节点最大值。
- 不增加第三方依赖，不构造未经验证的氦气压力—导热关系，不实现限值事件回溯求根。

---

## 文件结构

- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/__init__.py`：事故算例包标记。
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/run_v14_helium_depressurization.py`：事故参数、状态切换、诊断、限值和运行入口。
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/test_v14_helium_depressurization.py`：纯函数和轻量 mock 测试。
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/initial_state/steady_restart_t013864s.npz`：初始稳态独立副本。
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/initial_state/run_config.json`：初态重建配置及 `helium_accident_active=false` 标记。
- Modify: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/README.md`：补充实际命令、文件校验值和完成后的运行结果。
- Modify: `testModule/AI_AGENT_TESTMODULE_ANALYSIS.md`：登记新事故入口、restart 语义和最小验证命令。
- Create at runtime: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/runs/smoke_0p1s/`。
- Create at runtime: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/runs/accident_100s/`。

---

### Task 1: 固化初始状态与包结构

**Files:**
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/__init__.py`
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/initial_state/steady_restart_t013864s.npz`
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/initial_state/run_config.json`
- Modify: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/README.md`

**Interfaces:**
- Consumes: 源 restart `V14_210kW_fast_steady_temp/runs/physical_cp_plus2000s_from13164/checkpoint_t013864s.npz` 及同目录 `run_config.json`。
- Produces: 可由现有 `load_baseline_debug_config()` 直接读取的独立初态目录；配置必须包含 `helium_accident_active: false`。

- [ ] **Step 1: 创建包标记**

```python
"""V14 210 kW 全堆氦气瞬时完全失压事故算例。"""
```

- [ ] **Step 2: 复制二进制 restart 并核对哈希**

Run:

```powershell
New-Item -ItemType Directory -Force testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\initial_state
Copy-Item testModule\Full_Loop_Cases_10kW\V14_210kW_fast_steady_temp\runs\physical_cp_plus2000s_from13164\checkpoint_t013864s.npz testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\initial_state\steady_restart_t013864s.npz
Get-FileHash testModule\Full_Loop_Cases_10kW\V14_210kW_fast_steady_temp\runs\physical_cp_plus2000s_from13164\checkpoint_t013864s.npz
Get-FileHash testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\initial_state\steady_restart_t013864s.npz
```

Expected: 两个 SHA256 完全相同，文件大小完全相同。

- [ ] **Step 3: 复制初始配置并添加明确事故标记**

先复制源配置，再用 `apply_patch` 在顶层 JSON 中加入以下键，其他物理配置逐字保留：

```json
{
  "helium_accident_active": false,
  "helium_accident_model": "instantaneous_total_loss_of_gas_conduction",
  "helium_h_initial_w_m2k": 5678.0,
  "helium_h_final_w_m2k": 0.0
}
```

Run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m json.tool testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\initial_state\run_config.json
```

Expected: exit code 0；`point_kinetics_enabled=false`、`external_heat_enabled=false`、`power_w=210000.0`、`target_flow_kg_s=2.46` 保持不变。

- [ ] **Step 4: 在 README 记录独立 restart 的 SHA256**

使用标准补丁工具加入 `initial_state/steady_restart_t013864s.npz SHA256:`，其后逐字填写 Step 2 命令返回的实际64位哈希。修改后再次运行 `Get-FileHash`，并用 `rg -n "SHA256:"` 确认 README 中只有这一条初态哈希记录且与命令输出一致。

- [ ] **Step 5: 提交初始状态包**

```powershell
git add -- testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/__init__.py testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/initial_state testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/README.md
git commit -m "add V14 helium accident initial state"
```

---

### Task 2: 氦气隙定位与事故状态切换

**Files:**
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/run_v14_helium_depressurization.py`
- Create: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/test_v14_helium_depressurization.py`

**Interfaces:**
- Consumes: `build['tfes']`、`build['ring_multipliers']`、每个 TFE 的 `couplers['collector_iclad_gap']`。
- Produces: `collect_helium_gaps(build) -> dict[str, tuple[Any, int]]`、`set_helium_h_eq(gaps, h_eq_w_m2k) -> None`、`read_source_accident_state(source_config) -> bool`。

- [ ] **Step 1: 写失败测试，固定名称、倍率和切换语义**

```python
import unittest

from testModule.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization import (
    run_v14_helium_depressurization as runner,
)


class Gap:
    def __init__(self):
        self.gap = 5.0e-5
        self.k_gas = 5678.0 * self.gap
        self.eps1 = 0.60
        self.eps2 = 0.80


class TFE:
    def __init__(self):
        self.couplers = {'collector_iclad_gap': Gap()}


class HeliumGapTests(unittest.TestCase):
    def make_build(self):
        names = runner.REPRESENTATIVE_NAMES
        return {
            'tfes': {name: TFE() for name in names},
            'ring_multipliers': dict(zip(names, runner.EXPECTED_MULTIPLIERS)),
        }

    def test_collects_exactly_five_expected_gaps(self):
        gaps = runner.collect_helium_gaps(self.make_build())
        self.assertEqual(tuple(gaps), runner.REPRESENTATIVE_NAMES)
        self.assertEqual(tuple(mult for _, mult in gaps.values()), runner.EXPECTED_MULTIPLIERS)

    def test_instantaneous_loss_only_clears_gas_conduction(self):
        gaps = runner.collect_helium_gaps(self.make_build())
        runner.set_helium_h_eq(gaps, 0.0)
        for gap, _ in gaps.values():
            self.assertEqual(gap.k_gas, 0.0)
            self.assertEqual(gap.eps1, 0.60)
            self.assertEqual(gap.eps2, 0.80)

    def test_source_accident_marker_is_required_and_boolean(self):
        self.assertFalse(runner.read_source_accident_state({'helium_accident_active': False}))
        self.assertTrue(runner.read_source_accident_state({'helium_accident_active': True}))
        with self.assertRaises(ValueError):
            runner.read_source_accident_state({})
        with self.assertRaises(ValueError):
            runner.read_source_accident_state({'helium_accident_active': 'false'})
```

- [ ] **Step 2: 运行测试并确认失败**

Run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.test_v14_helium_depressurization -v
```

Expected: FAIL，原因是 runner 或目标函数尚不存在。

- [ ] **Step 3: 写最小实现**

runner 顶部定义并实现：

```python
from __future__ import annotations

from pathlib import Path
import math
from typing import Any, Dict

REPRESENTATIVE_NAMES = ('Center', 'Ring1', 'Ring2', 'Ring3', 'Ring4')
EXPECTED_MULTIPLIERS = (1, 6, 9, 18, 24)
HELIUM_GAP_KEY = 'collector_iclad_gap'
HELIUM_H_INITIAL_W_M2K = 5678.0
HELIUM_H_FINAL_W_M2K = 0.0
HELIUM_GAP_WIDTH_M = 5.0e-5


def collect_helium_gaps(build: Dict[str, Any]) -> Dict[str, tuple[Any, int]]:
    tfes = build['tfes']
    multipliers = build['ring_multipliers']
    if tuple(tfes) != REPRESENTATIVE_NAMES:
        raise ValueError(f'unexpected TFE names/order: {tuple(tfes)}')
    actual_multipliers = tuple(int(multipliers[name]) for name in REPRESENTATIVE_NAMES)
    if actual_multipliers != EXPECTED_MULTIPLIERS:
        raise ValueError(f'unexpected TFE multipliers: {actual_multipliers}')
    result: Dict[str, tuple[Any, int]] = {}
    for name, multiplier in zip(REPRESENTATIVE_NAMES, EXPECTED_MULTIPLIERS):
        gap = tfes[name].couplers.get(HELIUM_GAP_KEY)
        if gap is None:
            raise ValueError(f'{name} missing {HELIUM_GAP_KEY}')
        if not math.isclose(float(gap.gap), HELIUM_GAP_WIDTH_M, rel_tol=0.0, abs_tol=1.0e-12):
            raise ValueError(f'{name} helium gap width is {float(gap.gap)} m')
        h_eq = float(gap.k_gas) / float(gap.gap)
        if not math.isclose(h_eq, HELIUM_H_INITIAL_W_M2K, rel_tol=0.0, abs_tol=1.0e-9):
            raise ValueError(f'{name} initial helium h_eq is {h_eq} W/(m2*K)')
        result[name] = (gap, multiplier)
    return result


def set_helium_h_eq(gaps: Dict[str, tuple[Any, int]], h_eq_w_m2k: float) -> None:
    h_eq = float(h_eq_w_m2k)
    if not math.isfinite(h_eq) or h_eq < 0.0:
        raise ValueError('helium h_eq must be finite and non-negative')
    for gap, _ in gaps.values():
        gap.k_gas = h_eq * float(gap.gap)


def read_source_accident_state(source_config: Dict[str, Any]) -> bool:
    value = source_config.get('helium_accident_active')
    if not isinstance(value, bool):
        raise ValueError('run_config.json must contain boolean helium_accident_active')
    return value
```

- [ ] **Step 4: 运行测试并确认通过**

Run: Task 2 Step 2 的同一命令。

Expected: 3 tests PASS。

- [ ] **Step 5: 提交事故切换核心**

```powershell
git add -- testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/run_v14_helium_depressurization.py testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/test_v14_helium_depressurization.py
git commit -m "add V14 helium gap accident switch"
```

---

### Task 3: 温度、气隙传热诊断与限值判定

**Files:**
- Modify: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/run_v14_helium_depressurization.py`
- Modify: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/test_v14_helium_depressurization.py`

**Interfaces:**
- Consumes: `core.tfes`、`core.mod_rings`、`core.reflector`、`solid.T`、`solid.mesh.geom_data.node_centers_y`、`gap.bc1.current_flux`、倍率字典。
- Produces: `collect_temperature_peaks(core) -> list[dict[str, Any]]`、`find_limit_trip(peaks) -> dict[str, Any] | None`、`collect_helium_metrics(build, gaps, accident_time_s, active) -> dict[str, Any]`。

- [ ] **Step 1: 写失败测试，覆盖五类限值和首次越限证据**

```python
import numpy as np


class Mesh:
    class Geom:
        node_centers_y = np.array([0.1, 0.2])
    geom_data = Geom()


class Solid:
    def __init__(self, values):
        self.T = np.asarray(values, dtype=float)
        self.mesh = Mesh()


class CoreForLimits:
    def __init__(self):
        self.tfes = {
            name: type('TFE', (), {'solids': {
                'inner_clad': Solid([900.0, 1000.0]),
                'outer_clad': Solid([890.0, 910.0]),
                'pellet': Solid([2600.0, 2690.0]),
                'collector': Solid([1000.0, 1024.0] if name == 'Ring3' else [1000.0, 1010.0]),
                'moderator': Solid([800.0, 850.0]),
            }})()
            for name in runner.REPRESENTATIVE_NAMES
        }
        self.mod_rings = [Solid([870.0, 900.0])]
        self.reflector = Solid([900.0, 950.0])


def test_collector_limit_reports_ring_and_axial_position(self):
    peaks = runner.collect_temperature_peaks(CoreForLimits())
    trip = runner.find_limit_trip(peaks)
    self.assertEqual(trip['component'], 'collector')
    self.assertEqual(trip['representative'], 'Ring3')
    self.assertEqual(trip['limit_k'], 1023.0)
    self.assertEqual(trip['actual_k'], 1024.0)
    self.assertEqual(trip['axial_position_m'], 0.2)
```

- [ ] **Step 2: 运行单个测试并确认失败**

Run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.test_v14_helium_depressurization.HeliumGapTests.test_collector_limit_reports_ring_and_axial_position -v
```

Expected: FAIL，原因是诊断函数不存在。

- [ ] **Step 3: 实现精确限值表和峰值收集**

```python
TEMPERATURE_LIMITS_K = {
    'channel_wall': 1058.0,
    'pellet': 2700.0,
    'collector': 1023.0,
    'moderator': 930.0,
    'reflector': 1000.0,
}


def _peak(solid: Any, *, component: str, representative: str) -> Dict[str, Any]:
    values = np.asarray(solid.T, dtype=float)
    index = int(np.nanargmax(values))
    axial = np.asarray(solid.mesh.geom_data.node_centers_y, dtype=float)
    return {
        'component': component,
        'representative': representative,
        'actual_k': float(values[index]),
        'limit_k': float(TEMPERATURE_LIMITS_K[component]),
        'axial_position_m': float(axial[index]),
    }


def collect_temperature_peaks(core: Any) -> list[Dict[str, Any]]:
    peaks = []
    for name in REPRESENTATIVE_NAMES:
        solids = core.tfes[name].solids
        wall_peaks = [
            _peak(solids['inner_clad'], component='channel_wall', representative=f'{name}:inner_clad'),
            _peak(solids['outer_clad'], component='channel_wall', representative=f'{name}:outer_clad'),
        ]
        peaks.extend(wall_peaks)
        peaks.append(_peak(solids['pellet'], component='pellet', representative=name))
        peaks.append(_peak(solids['collector'], component='collector', representative=name))
        if 'moderator' in solids:
            peaks.append(_peak(solids['moderator'], component='moderator', representative=name))
    for index, solid in enumerate(core.mod_rings):
        peaks.append(_peak(solid, component='moderator', representative=f'global_mod_ring_{index}'))
    if core.reflector is None:
        raise ValueError('V14 core missing global reflector')
    peaks.append(_peak(core.reflector, component='reflector', representative='global_reflector'))
    return peaks


def find_limit_trip(peaks: list[Dict[str, Any]]) -> Dict[str, Any] | None:
    violations = [peak for peak in peaks if peak['actual_k'] > peak['limit_k']]
    if not violations:
        return None
    return max(violations, key=lambda item: item['actual_k'] / item['limit_k'])
```

实现前必须 `import numpy as np`。如果 `solid.T` 含非有限值，单独返回 `component='nonfinite_temperature'` 的触发记录，不允许 `np.nanargmax` 静默忽略全部非有限数组。

- [ ] **Step 4: 实现氦气参数与传热诊断**

`collect_helium_metrics()` 对每个代表环记录：

```python
def collect_helium_metrics(build, gaps, accident_time_s, active):
    row = {
        'accident_elapsed_s': float(build['system'].global_time) - float(accident_time_s),
        'helium_accident_active': bool(active),
        'helium_h_eq_W_m2K': 0.0 if active else HELIUM_H_INITIAL_W_M2K,
        'helium_conduction_fraction': 0.0 if active else 1.0,
    }
    total_scaled = 0.0
    r_values = []
    for name in REPRESENTATIVE_NAMES:
        gap, multiplier = gaps[name]
        tfe = build['tfes'][name]
        collector = tfe.solids['collector']
        inner_clad = tfe.solids['inner_clad']
        q_out = -float(np.sum(np.asarray(gap.bc1.current_flux, dtype=float)))
        row[f'{name}_collector_mean_T_K'] = float(np.mean(collector.T))
        row[f'{name}_collector_max_T_K'] = float(np.max(collector.T))
        row[f'{name}_inner_clad_mean_T_K'] = float(np.mean(inner_clad.T))
        row[f'{name}_inner_clad_max_T_K'] = float(np.max(inner_clad.T))
        row[f'{name}_helium_gap_heat_out_W'] = q_out
        total_scaled += multiplier * q_out
        r_values.extend(np.asarray(gap.R_gap_total, dtype=float).ravel())
    row['helium_gap_heat_out_scaled_W'] = total_scaled
    row['helium_gap_R_total_min_K_W'] = float(np.min(r_values))
    row['helium_gap_R_total_max_K_W'] = float(np.max(r_values))
    return row
```

在每次采集前由运行器按现有 SystemManager 生命周期完成 coupler 同步；不得为统计再次推进物理状态。

- [ ] **Step 5: 运行全部轻量测试**

Run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.test_v14_helium_depressurization -v
```

Expected: 全部 PASS。

- [ ] **Step 6: 提交诊断和限值逻辑**

```powershell
git add -- testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/run_v14_helium_depressurization.py testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/test_v14_helium_depressurization.py
git commit -m "add helium accident diagnostics and limits"
```

---

### Task 4: 运行编排、CSV、restart 恢复与 CLI

**Files:**
- Modify: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/run_v14_helium_depressurization.py`
- Modify: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/test_v14_helium_depressurization.py`

**Interfaces:**
- Consumes: 现有反应性 runner 的 `ReactivityControlRunConfig`、`load_baseline_debug_config()`、`prepare_reactivity_control()`、`collect_reactivity_metrics()`、`_next_step_dt()`、`_write_json()`，以及 debug runner 的 `build_debug_case()`。
- Produces: `HeliumAccidentRunConfig`、`run_helium_accident(config) -> dict[str, Any]`、命令行入口。

- [ ] **Step 1: 写失败测试，固定默认参数和事故续算恢复**

```python
from pathlib import Path


def test_default_run_config_matches_approved_case(self):
    config = runner.HeliumAccidentRunConfig(restart_in=Path('steady.npz'))
    self.assertEqual(config.duration_s, 100.0)
    self.assertEqual(config.dt_s, 0.05)
    self.assertEqual(config.record_interval_s, 0.1)
    self.assertEqual(config.checkpoint_interval_s, 10.0)
    self.assertEqual(config.wall_limit_k, 1058.0)
    self.assertEqual(config.pellet_limit_k, 2700.0)
    self.assertEqual(config.collector_limit_k, 1023.0)
    self.assertEqual(config.moderator_limit_k, 930.0)
    self.assertEqual(config.reflector_limit_k, 1000.0)


def test_accident_restart_is_reapplied_without_retrigger(self):
    gaps = runner.collect_helium_gaps(self.make_build())
    event = runner.restore_or_trigger_accident(
        gaps,
        source_config={
            'helium_accident_active': True,
            'helium_accident_time_absolute_s': 13864.2,
        },
        current_time_s=13874.2,
    )
    self.assertFalse(event['triggered_now'])
    self.assertEqual(event['accident_time_absolute_s'], 13864.2)
    self.assertTrue(all(gap.k_gas == 0.0 for gap, _ in gaps.values()))
```

- [ ] **Step 2: 运行新增测试并确认失败**

Run: Task 3 Step 5 的命令。

Expected: FAIL，原因是运行配置和恢复函数不存在。

- [ ] **Step 3: 实现配置与事故恢复**

```python
from dataclasses import dataclass


@dataclass(frozen=True)
class HeliumAccidentRunConfig:
    restart_in: Path
    output_dir: Path = Path(__file__).resolve().parent / 'runs' / 'default'
    duration_s: float = 100.0
    dt_s: float = 0.05
    record_interval_s: float = 0.1
    checkpoint_interval_s: float = 10.0
    min_fluid_temperature_stop_k: float = 500.0
    max_power_factor: float = 2.0
    wall_limit_k: float = 1058.0
    pellet_limit_k: float = 2700.0
    collector_limit_k: float = 1023.0
    moderator_limit_k: float = 930.0
    reflector_limit_k: float = 1000.0


def restore_or_trigger_accident(gaps, source_config, current_time_s):
    active = read_source_accident_state(source_config)
    if active:
        if 'helium_accident_time_absolute_s' not in source_config:
            raise ValueError('active helium accident missing absolute event time')
        event_time = float(source_config['helium_accident_time_absolute_s'])
        triggered_now = False
    else:
        event_time = float(current_time_s)
        triggered_now = True
    set_helium_h_eq(gaps, HELIUM_H_FINAL_W_M2K)
    return {
        'helium_accident_active': True,
        'triggered_now': triggered_now,
        'accident_time_absolute_s': event_time,
        'h_before_W_m2K': HELIUM_H_INITIAL_W_M2K,
        'h_after_W_m2K': HELIUM_H_FINAL_W_M2K,
        'affected_representatives': list(REPRESENTATIVE_NAMES),
        'physical_tfe_count': sum(EXPECTED_MULTIPLIERS),
    }
```

运行时限值来自 `HeliumAccidentRunConfig`；构建候选峰值时将配置值传入，不在两个位置维护第二份可变限值。

- [ ] **Step 4: 实现 `run_helium_accident()` 主流程**

主流程必须按以下精确顺序实现：

```python
runtime = ReactivityControlRunConfig(
    restart_in=config.restart_in,
    output_dir=config.output_dir,
    duration_s=config.duration_s,
    dt_s=config.dt_s,
    record_interval_s=config.record_interval_s,
    checkpoint_interval_s=config.checkpoint_interval_s,
    min_fluid_temperature_stop_k=config.min_fluid_temperature_stop_k,
    max_power_factor=config.max_power_factor,
)
debug, source_config = load_baseline_debug_config(runtime)
build = build_debug_case(debug, apply_fixed_power=False)
system = build['system']
core = build['core']
handoff_type = prepare_reactivity_control(
    core,
    source_point_kinetics_enabled=bool(source_config['point_kinetics_enabled']),
    expected_power_w=float(debug.power_w),
)
gaps = collect_helium_gaps(build)
source_active = read_source_accident_state(source_config)

# 初始稳态只在事故触发前记录一次；事故续算的首行已经是失压状态。
if not source_active:
    pre_event = collect_all_metrics(build, gaps, handoff_type, debug.power_w, system.global_time, False, 0.0)
    preflight_trip = find_limit_trip(collect_temperature_peaks(core))
    if preflight_trip is not None:
        raise RuntimeError(f'initial state violates limit: {preflight_trip}')

event = restore_or_trigger_accident(gaps, source_config, system.global_time)
```

然后写入输出 `run_config.json`，必须包含：

```python
{
    'point_kinetics_enabled': True,
    'reactivity_control_mode': 'temperature_feedback_only',
    'external_reactivity': 0.0,
    'control_drum_enabled': False,
    'external_heat_enabled': False,
    'helium_accident_active': True,
    'helium_accident_model': 'instantaneous_total_loss_of_gas_conduction',
    'helium_accident_time_absolute_s': event['accident_time_absolute_s'],
    'helium_h_initial_w_m2k': 5678.0,
    'helium_h_final_w_m2k': 0.0,
    'temperature_limits_k': {
        'channel_wall': config.wall_limit_k,
        'pellet': config.pellet_limit_k,
        'collector': config.collector_limit_k,
        'moderator': config.moderator_limit_k,
        'reflector': config.reflector_limit_k,
    },
}
```

推进循环继续调用：

```python
system.step(
    dt,
    inner_iter=int(debug.inner_iter),
    fail_on_fluid_nonconvergence=False,
    fluid_max_iter=int(debug.fluid_max_iter),
    reactivity_control=0.0,
)
```

每步之后依次检查非有限值、5类温度限值、2倍功率上限和最低流体温度。任何触发都写
`limit_trip.json`、保存 `emergency_restart.npz` 并退出循环；正常结束保存 `stage_01_restart.npz`。

- [ ] **Step 5: 实现 CSV 和 JSON 输出**

CSV 字段由现有 `REACTIVITY_HISTORY_FIELDS` 加事故字段组成；使用标准库：

```python
def append_history(path: Path, fields: list[str], row: Dict[str, Any]) -> None:
    write_header = not path.exists()
    with path.open('a', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        if write_header:
            writer.writeheader()
        writer.writerow(row)
```

`accident_event.json` 在首次触发和事故续算两种情况下都写出，`triggered_now` 用于区分本次运行是否发生跳变。`latest_state.json` 和 `run_summary.json` 必须包含 `stop_reason`、最新指标和最新 restart 路径。

- [ ] **Step 6: 实现 CLI**

CLI 参数固定为：

```text
--restart-in
--output-dir
--duration
--dt
--record-interval
--checkpoint-interval
--min-fluid-temperature-stop
--max-power-factor
--wall-limit-k
--pellet-limit-k
--collector-limit-k
--moderator-limit-k
--reflector-limit-k
```

默认值必须与 `HeliumAccidentRunConfig` 一致。`--restart-in` 必填，输出目录已存在 `history.csv` 时拒绝覆盖。

- [ ] **Step 7: 运行轻量测试和语法检查**

Run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.test_v14_helium_depressurization -v
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m py_compile testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\run_v14_helium_depressurization.py
```

Expected: 全部测试 PASS；`py_compile` exit code 0。

- [ ] **Step 8: 提交完整 runner**

```powershell
git add -- testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/run_v14_helium_depressurization.py testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/test_v14_helium_depressurization.py
git commit -m "add V14 helium depressurization runner"
```

---

### Task 5: 真实 restart smoke 与 restart 续算验证

**Files:**
- Modify if required by observed defect: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/run_v14_helium_depressurization.py`
- Modify if required by observed defect: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/test_v14_helium_depressurization.py`
- Create at runtime: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/runs/smoke_0p1s/`
- Create at runtime: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/runs/restart_smoke_0p05s/`

**Interfaces:**
- Consumes: Task 1 初态和 Task 4 CLI。
- Produces: 可审计的 0.1 s 事故输出，以及证明事故 restart 不恢复氦气导热的 0.05 s 续算输出。

- [ ] **Step 1: 运行 0.1 s 真实 smoke**

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\run_v14_helium_depressurization.py --restart-in testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\initial_state\steady_restart_t013864s.npz --output-dir testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\runs\smoke_0p1s --duration 0.1 --dt 0.05 --record-interval 0.05 --checkpoint-interval 0.1
```

Expected:

- exit code 0；
- `history.csv` 第一行是 `helium_accident_active=false`、`h=5678`；
- 后续行为 `helium_accident_active=true`、`h=0`；
- `accident_event.json` 的 `physical_tfe_count=58`；
- 裂变、衰变、总功率和反应性均为有限值；
- 5类初始温度没有在事故前越限；
- 5个气隙 `R_gap_total` 有限，证明辐射仍存在；
- `stage_01_restart.npz` 存在。

- [ ] **Step 2: 从 smoke 输出续算 0.05 s**

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\run_v14_helium_depressurization.py --restart-in testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\runs\smoke_0p1s\stage_01_restart.npz --output-dir testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\runs\restart_smoke_0p05s --duration 0.05 --dt 0.05 --record-interval 0.05 --checkpoint-interval 0.05
```

Expected:

- 首行已经是 `helium_accident_active=true` 和 `h=0`；
- `accident_event.json.triggered_now=false`；
- 事故绝对发生时间与第一段完全相同；
- 点堆为 `reactivity_continuation`，没有重新初始化或校准。

- [ ] **Step 3: 若 smoke 暴露缺陷，先写对应失败测试再做最小修复**

只允许修复 smoke 直接证明的根因。修复后重建新的空输出目录重新执行 Step 1 和 Step 2，不覆盖失败目录；失败目录保留用于比较。

- [ ] **Step 4: 提交 smoke 验证后的必要修复**

如果没有代码改动则跳过提交；如果有改动：

```powershell
git add -- testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/run_v14_helium_depressurization.py testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/test_v14_helium_depressurization.py
git commit -m "fix V14 helium accident restart handoff"
```

---

### Task 6: 文档和模块导航同步

**Files:**
- Modify: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/README.md`
- Modify: `testModule/AI_AGENT_TESTMODULE_ANALYSIS.md`

**Interfaces:**
- Consumes: smoke 的实际命令、实际初态哈希和实际输出字段。
- Produces: 新会话可直接定位、运行和续算该事故的中文说明。

- [ ] **Step 1: 更新事故 README**

加入：

- 初次运行的完整 PowerShell 命令；
- 从事故 restart 续算的完整命令；
- 初始 restart 实际 SHA256；
- `history.csv` 字段解释和热流正方向；
- 五类限值的监视对象映射；
- smoke 的起止时间、功率、关键温度、反馈和气隙传热摘要；
- 明确运行产物默认不提交，初始 restart 因交付需求作为算例输入保留。

- [ ] **Step 2: 更新 testModule 手册**

在 V14 事故算例部分登记：

```text
入口：Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/run_v14_helium_depressurization.py
初态：initial_state/steady_restart_t013864s.npz
模型：全部58根TFE的collector-inner-clad氦气导热瞬时归零，辐射保留
反应性：外界0、控制鼓关闭、仅温度反馈增量
续算：必须读取run_config.json恢复helium_accident_active
```

- [ ] **Step 3: 检查文档无占位符并提交**

```powershell
rg -n "T[B]D|T[O]DO|待定|place[holder]" testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\README.md testModule\AI_AGENT_TESTMODULE_ANALYSIS.md
git add -- testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/README.md testModule/AI_AGENT_TESTMODULE_ANALYSIS.md
git commit -m "document V14 helium accident operation"
```

Expected: `rg` 无匹配；提交只包含两份文档。

---

### Task 7: 首轮 100 s 事故计算与结果评估

**Files:**
- Create at runtime: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/runs/accident_100s/`
- Modify after run: `testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/README.md`

**Interfaces:**
- Consumes: Task 1 初态、通过 smoke 的 runner。
- Produces: 100 s 历史、每10 s checkpoint、正常或紧急 restart、限值证据和评估结论。

- [ ] **Step 1: 从原始稳态独立启动 100 s**

不得从 0.1 s smoke 续算，避免把验证段重复混入正式结果：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\run_v14_helium_depressurization.py --restart-in testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\initial_state\steady_restart_t013864s.npz --output-dir testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\runs\accident_100s --duration 100 --dt 0.05 --record-interval 0.1 --checkpoint-interval 10
```

- [ ] **Step 2: 持续监测并确认终止类型**

运行期间至少每60秒向用户更新一次：当前事故相对时间、总/裂变/衰变功率、有效温度反馈、
通道壁/芯块/接收极/慢化剂/反射层最高温度、最接近限值的裕量、流体收敛状态。

Expected final state is exactly one of:

- `stop_reason=completed` 且存在 `stage_01_restart.npz`；
- 温度或功率限值触发，存在 `emergency_restart.npz` 和 `limit_trip.json`。

- [ ] **Step 3: 做结果完整性检查**

Run:

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -c "import csv,json,math,pathlib; p=pathlib.Path(r'testModule\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization\runs\accident_100s'); rows=list(csv.DictReader((p/'history.csv').open(encoding='utf-8'))); assert rows; assert all(all(math.isfinite(float(r[k])) for k in ('core_total_power_W','fission_power_W','decay_power_W','effective_temperature_feedback')) for r in rows); cfg=json.loads((p/'run_config.json').read_text(encoding='utf-8')); assert cfg['helium_accident_active'] is True; assert cfg['helium_h_final_w_m2k']==0.0; print(len(rows), rows[-1]['time_s'])"
```

Expected: 约1001个记录点；若提前越限则记录数更少，但所有检查仍通过。

- [ ] **Step 4: 评估并写入 README**

从 `history.csv` 和终止 JSON 写入实际数值：

- 起止绝对时间和事故持续时间；
- 起止及峰值总功率、裂变功率和衰变功率；
- 有效温度反馈和各反馈分量的起止及极值；
- 5类受限温度的起止、峰值、峰值位置和最小裕量；
- 全堆氦气隙传热功率变化；
- 主回路流量、流体收敛和 TEC 状态；
- 是否触发限值、触发时间及可续算 restart 路径。

- [ ] **Step 5: 最终回归和提交评估文档**

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_10kW.V14_210kW_helium_depressurization.test_v14_helium_depressurization -v
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases_10kW.V14_210kW_reactivity_control.test_v14_210kw_reactivity_control -v
git diff --check
git add -- testModule/Full_Loop_Cases_10kW/V14_210kW_helium_depressurization/README.md
git commit -m "record V14 helium accident result"
```

Expected: 两组测试全部 PASS；`git diff --check` 无输出；运行产物保留在本地 `runs/`，不批量加入提交。

