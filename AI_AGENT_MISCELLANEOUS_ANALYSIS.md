# 辅助模块 — AI Agent 速查手册

> 覆盖项目根目录及辅助目录中未纳入其他文档的源码文件。
> 包括: MathSolvers/, profiler.py, main.py, 重构目录

---

## 1. MathSolvers/ 数学求解器

### 1.1 solver_module.py — NuclearODESolver

```python
class NuclearODESolver:
    """scipy solve_ivp 封装 (替代 Fortran GEAR.FOR)"""
    def __init__(self, method='BDF', rtol=1e-6, atol=1e-8)
    def solve(ode_func, t_span, y0, jac=None, events=None) → ODEResult
```

- 封装 `scipy.integrate.solve_ivp`
- 默认 BDF (后向差分公式), 适合刚性方程
- 也支持 'Radau', 'LSODA'
- 用于 PointReactor 等 ODE 求解

### 1.2 optimization_utils.py — FluidJacobianBlockLayout

```python
class FluidJacobianBlockLayout:
    """流体系统块状雅可比矩阵构建器 (2N+M 维)"""
    def __init__(self, volumes, junctions)
    # 状态向量布局: [P₀...P_{N-1}, h₀...h_{N-1}, W₀...W_{M-1}]
    # 提供 p_row/p_col/p_h_block, h_row/h_col/h_block 等块索引
```

- 用于 HydraulicNetwork 的稀疏雅可比矩阵构建
- 适配 SystemManager 的块状变量布局

---

## 2. profiler.py — 性能剖析器

```python
class TEASAProfiler:
    """轻量级性能剖析器 (装饰器模式)"""
    stats: Dict[str, {'count': int, 'time': float}]  # 类级别统计

    @classmethod
    def profile(cls, func) → wrapper:
        # 装饰器: 统计 func 调用次数 + 累积时间
        # 使用 functools.wraps + time.perf_counter()
```

**使用方式:**
```python
@TEASAProfiler.profile
def some_heavy_method(self, ...):
    ...
```
- 被 Couplers.py 中的 `FluidSolidCouple.execute()` 和 `GapCouple2D.sync()` 使用
- 类级别 `stats` 字典全局累积

---

## 3. main.py

- 空文件 (0 行), 可能是未来入口预留

---

## 4. 重构目录 (可选参考)

| 目录 | 文件 | 说明 |
|------|------|------|
| `SystemManager优化重构/` | `validate_inner_iter_picard.py` | Picard 内迭代验证 |
| `HeatPipe优化重构/` | `benchmark_ringhp_v5.py`, `run_ringhp_v5_fin_tangent_1800s_dt025.py`, `test_single_hp_cold_start.py` | RingHP v5 性能测试 |
| `HeatConduction优化重构/` | `benchmark_systemmanager_heatconduction.py` | 导热求解器性能基准 |
| `HydraulicNetwork优化重构/` | `benchmark_open_loop_complex_300cv.py` | 水力网络性能基准 |

这些目录中的 .py 文件属于开发过程中的性能测试脚本，非核心功能代码。

---

## 5. inputs/ 目录

| 文件 | 说明 |
|------|------|
| `CoreInput.txt` | 堆芯输入参数 |
| `GloInput.txt` | 全局输入参数 |
| `HPInput.txt` | 热管输入 |
| `PipeInput.txt` | 管道输入 |
| `PipeNetInput.txt` | 管网输入 |
| `PumpInput.txt` | 泵输入 |
| `RadiatorInput.txt` | 散热器输入 |
| `steady.txt` | 稳态参数 |
| `transient.txt` | 瞬态参数 |
| `trip.txt` | 触发条件 |
| `sheildtemp.txt` | 屏蔽层温度 |
| `density1.txt`, `density2.txt`, `density3.txt` | 密度数据 |

这些是 Fortran 时代的输入文件，Python 版可能已不再直接使用。

---

## 6. 文件覆盖完整检查清单

| 目录 | .py文件 | 已覆盖? | AI_AGENT文档 |
|------|---------|---------|-------------|
| Components/ | 16 | ✅ | AI_AGENT_BASICCOMPONENTS_ANALYSIS.md |
| Solvers/ | 9 | ✅ | AI_AGENT_SOLVERS_ANALYSIS.md |
| ThermoCalc/ | 5(.py)+6(.cpp/.h) | ✅ | AI_AGENT_THERMOCALC_ANALYSIS.md |
| Materials/ | 20 | ✅ | AI_AGENT_MATERIALS_ANALYSIS.md |
| Correlations/ | 1 | ✅ | AI_AGENT_CORRELATIONS_ANALYSIS.md |
| CoolantLoop/ | 13 | ✅ | AI_AGENT_COOLANTLOOP_ANALYSIS.md |
| testModule/ | 47 | ✅ | AI_AGENT_TESTMODULE_ANALYSIS.md |
| MathSolvers/ | 2 | ✅ | 本文档 §1 |
| profiler.py | 1 | ✅ | 本文档 §2 |
| main.py | 1 | ✅ | 本文档 §3 |
| 重构目录 | 5 | ✅ | 本文档 §4 |
| inputs/ | 14(.txt) | ⬜ | 本文档 §5 (Fortran遗留) |