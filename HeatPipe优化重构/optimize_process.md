# HeatPipe 优化过程记录

## 1. 优化目标

- 目标模块：
  - `Components/HPwithFin.py`
  - `Components/basicComponents/HeatPipe2D.py`
  - `Components/RingHP.py`
- 目标：在不改变 `CoolantLoop/test_coolant_loop_v5.py` 当前建模口径的前提下，建立可重复 benchmark，先完成静态 review 和热点定位，再按优先级推进修复与性能优化。

## 2. 基线脚本

- 脚本：`HeatPipe优化重构/benchmark_ringhp_v5.py`
- 设计思路：
  - 复用 `CoolantLoop/test_coolant_loop_v5.py` 的真实装配与参数；
  - 用固定步长 `dt * n_steps` 替代原脚本的长时运行，得到更稳定、可复现的优化基线；
  - 在运行时对 `RingHP / HPwithFin / HeatPipe2D / SingleVolumeProxy` 的关键方法挂接 `TEASAProfiler.profile`；
  - 保留 `FluidSolidCouple.execute`、`BaseHeatConduction.step` 等已有 profiler 数据，便于判断热点是否真的集中在热管模块。

## 3. 当前基线配置

- Python：`E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`
- 运行配置：
  - `dt = 0.05 s`
  - `n_steps = 20`
  - `inner_iter = 1`
  - `print_every = 5`
- 输出文件：`HeatPipe优化重构/ringhp_v5_baseline_history.csv`

## 4. Round 0：建基线 + 静态分析

- 代码改动：
  - 新建 `HeatPipe优化重构/benchmark_ringhp_v5.py`
  - 新建 `HeatPipe优化重构/TODO-list.md`
  - 新建 `HeatPipe优化重构/optimize_process.md`
- 语法检查：
  - `python -m py_compile Components\HPwithFin.py Components\basicComponents\HeatPipe2D.py Components\RingHP.py`
  - `python -m py_compile HeatPipe优化重构\benchmark_ringhp_v5.py`
- benchmark 结果：
  - `20` 步总 wall time：`6.181586 s`
  - 末步 `T_out`：`844.678 K`
  - 末步 `Q_total`：`49665.981 W`

### 4.1 Profiler 热点

按热管模块内部函数排序：

| 函数 | 调用次数 | 总耗时 |
| --- | ---: | ---: |
| `HeatPipe2D._update_properties` | 8769 | `2.111064 s` |
| `HeatPipe2D._compute_fluxes` | 8286 | `1.603796 s` |
| `HeatPipe2D._update_boundaries_state` | 8769 | `0.734677 s` |
| `RingHP.pre_step` | 40 | `0.692757 s` |
| `HPwithFin.pre_step` | 460 | `0.692228 s` |
| `HPwithFin._solve_fin_quasi_steady` | 920 | `0.666283 s` |
| `SingleVolumeProxy.add_coupling_source_distribution` | 483 | `0.000349 s` |

按全系统 profiler 排序时，热管模块也已经是主要热点来源：

| 函数 | 调用次数 | 总耗时 |
| --- | ---: | ---: |
| `SystemManager.step` | 20 | `6.155332 s` |
| `BaseHeatConduction.step` | 500 | `5.231320 s` |
| `BaseHeatConduction.get_derivatives` | 8248 | `4.291870 s` |
| `HeatPipe2D._update_properties` | 8769 | `2.111064 s` |
| `HeatPipe2D._compute_fluxes` | 8286 | `1.603796 s` |
| `HeatPipe2D._update_boundaries_state` | 8769 | `0.734677 s` |
| `RingHP.pre_step` | 40 | `0.692757 s` |
| `HPwithFin.pre_step` | 460 | `0.692228 s` |
| `HPwithFin._solve_fin_quasi_steady` | 920 | `0.666283 s` |

### 4.2 结论

- 当前热点高度集中在 `HeatPipe2D` 的属性更新、热流拼装和边界状态更新。
- `HPwithFin.pre_step()` 与 `_solve_fin_quasi_steady()` 已经形成第二层明显热点，其中“双求解估算等效导热率”是最直接的优化切入点。
- `RingHP` 宏观装配层本身耗时不高，但静态 review 暴露出两个需要先修的正确性问题：
  - `RingHP.py:215` 的闭包绑定错误；
  - `HPwithFin.py:297` 的翅片导热系数硬编码。

## 5. 下一轮建议

1. 先做正确性修复：闭包绑定、材料硬编码、总量统计口径。
2. 再做 `HeatPipe2D` buffer 化：`_update_properties()`、`_compute_fluxes()`、`_update_boundaries_state()`。
3. 最后压缩 `HPwithFin` 翅片求解成本，重点处理双求解和 Thomas 临时数组。

## 6. 运行命令

### 6.1 语法检查

```powershell
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' -m py_compile 'Components\HPwithFin.py' 'Components\basicComponents\HeatPipe2D.py' 'Components\RingHP.py'
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' -m py_compile 'HeatPipe优化重构\benchmark_ringhp_v5.py'
```

### 6.2 基线 benchmark

```powershell
$env:BENCH_DT='0.05'
$env:BENCH_N_STEPS='20'
$env:BENCH_INNER_ITER='1'
$env:BENCH_PRINT_EVERY='5'
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HeatPipe优化重构\benchmark_ringhp_v5.py'
```

### 6.3 快速 smoke test

```powershell
$env:BENCH_DT='0.05'
$env:BENCH_N_STEPS='2'
$env:BENCH_INNER_ITER='1'
$env:BENCH_PRINT_EVERY='1'
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HeatPipe优化重构\benchmark_ringhp_v5.py'
```
## 7. Round 1-3：WickMaterial + HeatPipe2D 稳定优化

### 7.1 代码改动

- `Materials/Solids/WickMaterial.py`
  - 新增 `k(T)` 查表插值路径，默认启用；
  - 保留原始解析计算 `_conductivity_direct()` 作为回归基准；
  - 记录高非线性温区，用于 `HeatPipe2D` 的保守局部更新。
- `Components/basicComponents/HeatPipe2D.py`
  - `_update_properties()` 改为复用 2D view，并按温度变化阈值进行局部物性更新；
  - `_update_boundaries_state()` 改为复用边界热阻 buffer；
  - `_compute_fluxes()` 改为复用 `Q_net_2d_buffer`、`_flux_x_buffer`、`_flux_y_buffer`；
  - 实现了冻结物性外层修正框架，但系统级回归中出现 `Required step size is less than spacing between numbers`，因此默认关闭 `enable_frozen_property_correction`，保留为后续实验开关。

### 7.2 WickMaterial 查表误差

- 检查温区：`650 K ~ 900 K`
- 全区最大绝对误差：`3.659982e+02`
- 全区最大相对误差：`3.665578e-04`
- 陡变区（`716 K ~ 800 K`）最大绝对误差：`3.659982e+02`
- 陡变区最大相对误差：`3.665578e-04`
- 自动识别的高非线性温区：`[273.074 K, 799.933 K]`

### 7.3 回归结果

- `testModule/test_single_hp_fin_energy_conservation.py --t-end 0.01 --dt 0.01 --no-csv --no-restart`
  - `Final Q_con_out = 867.988121 W`
  - `Final Tmin = 799.476195 K`
  - 与优化前基线一致
- `testModule/test_ringhp_node_coupling_energy_conservation.py --t-end 0.01 --inner-iter 1 --print-every 1 --no-csv --no-restart`
  - `coupling_valid = True`
  - 该 1-step smoke 主要用于确认耦合守恒未破坏，不作为总能量收敛判据
- `HeatPipe优化重构/benchmark_ringhp_v5.py`
  - 末步 `T_out = 844.678 K`
  - 末步 `Q_total = 49665.942 W`
  - 与优化前主基线一致

### 7.4 耗时对比

| 指标 | 优化前 | 优化后 |
| --- | ---: | ---: |
| benchmark wall time | `6.841774 s` | `4.860568 s` |
| `HeatPipe2D._update_properties` | `2.207253 s` | `0.600521 s` |
| `HeatPipe2D._compute_fluxes` | `1.689430 s` | `1.623842 s` |
| `HeatPipe2D._update_boundaries_state` | `0.776220 s` | `0.555204 s` |

结论：本轮稳定优化主要收益来自 `WickMaterial` 查表化和 `HeatPipe2D` 的局部/缓冲更新；`_update_properties()` 耗时下降最明显，系统级 benchmark 总耗时下降约 `28.96%`。
