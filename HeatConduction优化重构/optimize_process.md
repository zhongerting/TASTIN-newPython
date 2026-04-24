# HeatConduction 优化过程记录

## 1. 优化目标

- 目标文件：`Solvers/HeatConduction/HeatConduction.py`
- 核心目标：
  - 在不牺牲物理准确性的前提下提升导热求解速度
  - 保持 `SystemManager`、`FluidSolidCouple`、`SolidSolidCouple2D` 现有耦合接口兼容
  - 保留边界半节点长度处理、动态辐射边界、多层固固耦合的真实性

## 2. Benchmark 算例

- 基准脚本：`HeatConduction优化重构/benchmark_systemmanager_heatconduction.py`
- 调度器：`SystemManager`
- 拓扑来源：`testModule/test_parallel_channels.py`
- 几何/工况量级来源：`testModule/test_core_assemble_v6.py`
- 轴向功率形状来源：`testModule/test_core_assemble_v5.py`

### 2.1 主结构

- `4` 根并联流道
- 每根流道对应 `4` 个 `HeatConduction2D` 固体层：
  - `inner_liner`
  - `transition_shell`
  - `heated_shell`
  - `outer_shield`
- 层间 `SolidSolidCouple2D`
- 内壁 `FluidSolidCouple`
- 外壁动态辐射边界
- 发热层采用稳态功率 + 周期扰动

### 2.2 物理代表性

- 非均匀轴向网格：下缓冲段 + 活性段 + 上缓冲段
- 非均匀径向网格：各层内外两端加密
- 边界半节点显式进入：
  - 内壁流固换热
  - 外壁动态辐射
  - 顶/底绝热热阻路径
  - 层间固固耦合

## 3. 默认配置

### `smoke`

- `BENCH_MODE=smoke`
- 默认 `3` 步

### `baseline`

- `BENCH_MODE=baseline`
- 默认 `8` 步

### `stress`

- `BENCH_MODE=stress`
- 默认 `12` 步

## 4. 运行命令

### 4.1 语法检查

```powershell
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' -m py_compile 'HeatConduction优化重构\benchmark_systemmanager_heatconduction.py'
```

### 4.2 `smoke`

```powershell
$env:BENCH_MODE='smoke'
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HeatConduction优化重构\benchmark_systemmanager_heatconduction.py'
```

### 4.3 `baseline`

```powershell
$env:BENCH_MODE='baseline'
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HeatConduction优化重构\benchmark_systemmanager_heatconduction.py'
```

### 4.4 `stability`

```powershell
$env:BENCH_MODE='stability'
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HeatConduction优化重构\benchmark_systemmanager_heatconduction.py'
```

### 4.5 `stress`

```powershell
$env:BENCH_MODE='stress'
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HeatConduction优化重构\benchmark_systemmanager_heatconduction.py'
```

## 5. 当前有效基线（回退后）

- 当前代码状态：
  - `Solvers/HeatConduction/Boundary.py` 已回退到未优化状态
  - `Solvers/Couplers.py` 已回退到未优化状态
  - 当前所有后续优化，都以本节结果为唯一对照组
- 验证结果：
  - `py_compile`：通过
  - `smoke`：通过
  - `baseline`：通过
  - `stability(200步)`：通过

### 5.1 `smoke` 结果

- 配置：
  - `BENCH_MODE=smoke`
  - `n_steps = 3`
  - `dt = 0.02 s`
- 结果：
  - 初始化时间：`0.004374 s`
  - 步进总时间：`0.442687 s`
  - 平均单步时间：`0.147562 s/step`
  - 最大发热层温度末值：`795.778 K`
  - 平均出口流体温度末值：`757.937 K`
  - 外辐射功率末值：`3482.442 W`
- profiler 主热点：
  - `BaseHeatConduction.get_derivatives`：`0.274005 s`
  - `HeatConduction2D._compute_fluxes`：`0.199580 s`
  - `HeatConduction2D._compute_internal_resistance`：`0.027640 s`
  - `FluidSolidCouple.execute`：`0.001366 s`

### 5.2 `baseline` 结果

- 配置：
  - `BENCH_MODE=baseline`
  - `n_steps = 8`
  - `dt = 0.03 s`
- 结果：
  - 初始化时间：`0.004955 s`
  - 步进总时间：`1.923487 s`
  - 平均单步时间：`0.240436 s/step`
  - 最大发热层温度末值：`793.184 K`
  - 平均出口流体温度末值：`763.062 K`
  - 外辐射功率末值：`3484.121 W`
- 通道末态摘要：
  - `ChanA`：`W = 0.035054 kg/s`，`T_out = 755.777 K`
  - `ChanB`：`W = 0.034879 kg/s`，`T_out = 760.630 K`
  - `ChanC`：`W = 0.034699 kg/s`，`T_out = 765.488 K`
  - `ChanD`：`W = 0.034506 kg/s`，`T_out = 770.352 K`
- profiler 主热点：
  - `BaseHeatConduction.get_derivatives`：`1.040417 s`
  - `HeatConduction2D._compute_fluxes`：`0.728811 s`
  - `HeatConduction2D._compute_internal_resistance`：`0.109997 s`
  - `BaseHeatConduction._update_properties`：`0.070394 s`
  - `FluidSolidCouple.execute`：`0.004500 s`
  - `HydraulicNetwork.step_Picard`：`0.015506 s`
- 结论：
  - 短 `baseline` 继续承担快速 profiler 对比职责
  - 当前 benchmark 的主耗时仍集中在 `HeatConduction` 路径，而不是 `HydraulicNetwork`
  - 在多层耦合、动态辐射和并联通道同时存在时，`_compute_fluxes` 与 `get_derivatives` 仍是首要观察对象

### 5.3 `stability` 结果

- 配置：
  - `BENCH_MODE=stability`
  - `n_steps = 200`
  - `dt = 0.03 s`
- 结果：
  - 初始化时间：`0.005355 s`
  - 步进总时间：`48.042741 s`
  - 平均单步时间：`0.240214 s/step`
  - 最大发热层温度末值：`793.430 K`
  - 平均出口流体温度末值：`769.184 K`
  - 外辐射功率末值：`3463.403 W`
- 历史尾段检查：
  - `time_s` 尾段：`[5.94, 5.97, 6.00]`
  - `heated_layer_max_k` 尾段：`[793.400, 793.192, 793.430]`
  - `fluid_outlet_mean_k` 尾段：`[769.181, 769.229, 769.184]`
  - `radiation_power_w` 尾段：`[3463.097, 3463.020, 3463.403]`
- profiler 主热点：
  - `BaseHeatConduction.get_derivatives`：`26.540871 s`
  - `HeatConduction2D._compute_fluxes`：`18.608792 s`
  - `HeatConduction2D._compute_internal_resistance`：`2.804830 s`
  - `BaseHeatConduction._update_properties`：`1.792278 s`
  - `FluidSolidCouple.execute`：`0.105426 s`
  - `HydraulicNetwork.step_Picard`：`0.406605 s`
- 结论：
  - `200` 步长回归可以稳定完成，没有出现中途失稳或后段发散
  - 尾段波动平稳，适合作为后续单步优化的长期稳定性门槛
  - 后续优化必须同时通过短 `baseline` 和 `stability(200步)`，才允许保留

### 5.4 边界与网格审计摘录

- 轴向半节点长度：
  - `axial_dz_min = 0.010833 m`
  - `axial_dz_max = 0.015080 m`
  - `axial_dz_ratio = 1.392`
- 径向最强非均匀层：
  - `outer_shield.dr_ratio = 6.314`
  - `heated_shell.dr_ratio = 5.027`
- 说明：
  - benchmark 已覆盖明显非均匀径向网格
  - `left/right/top/bottom` 半节点长度均已进入审计输出

## 6. 单步优化执行规则

- 每一轮只允许一个技术思路进入代码
- 每一轮只允许一个主代码文件发生真实实现变更
- 每一轮都必须先过 `py_compile`，再过 `smoke`，再过短 `baseline`，最后过 `stability(200步)`
- 任一步若出现 `nan/inf`、奇异矩阵或关键物理量异常漂移，立即只回退该步
- 禁止把 `jac_sparsity`、返回缓冲、边界热阻缓存、内部热阻缓存打包到同一轮

## 7. 后续 Round 模板

### Round N

- 修改目标：
- 涉及文件：
- 修改内容：
- 风险点：
- 验证命令：
- benchmark 结果：
- profiler 变化：
- 边界/能量审计变化：
- 结论：

## 8. Round 1：`get_jac_sparsity()` 缓存

- 修改目标：
  - 只对 `HeatConduction2D.get_jac_sparsity()` 加缓存
  - 不修改热流、边界、物性、时间步进路径
- 涉及文件：
  - `Solvers/HeatConduction/HeatConduction.py`
- 修改内容：
  - 在 `HeatConduction2D.__init__()` 中新增 `_jac_sparsity_cache`
  - `get_jac_sparsity()` 首次生成后缓存，后续直接返回
- 风险点：
  - 仅依赖网格形状，不依赖温度或物性，物理风险极低
  - 网格若未来变为动态重构，需要额外失效机制；当前模型中网格固定，因此可安全缓存
- 验证命令：
  - `& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' -m py_compile 'Solvers\HeatConduction\HeatConduction.py' 'HeatConduction优化重构\benchmark_systemmanager_heatconduction.py'`
  - `BENCH_MODE=baseline`
  - `BENCH_MODE=stability`

### 8.1 结果对比

- `baseline` 结果对比：
  - 优化前：`0.219800 s/step`
  - 优化后：`0.219600 s/step`
  - 变化：`-0.000200 s/step`，约 `-0.09%`
- `stability(200步)` 结果对比：
  - 优化前：`0.202924 s/step`
  - 优化后：`0.201988 s/step`
  - 变化：`-0.000937 s/step`，约 `-0.46%`
- 关键物理量对比：
  - `baseline` 末态：`heated_layer_max_k`、`fluid_outlet_mean_k`、`radiation_power_w` 无变化
  - `stability` 尾段：`heated_layer_max_k`、`fluid_outlet_mean_k`、`radiation_power_w` 无变化
  - 结论：本轮缓存没有改变计算结果

### 8.2 profiler 对比

- `baseline`：
  - `BaseHeatConduction.get_derivatives`：`0.954559 s -> 0.958544 s`
  - `HeatConduction2D._compute_fluxes`：`0.669620 s -> 0.671780 s`
  - `HeatConduction2D._compute_internal_resistance`：`0.101821 s -> 0.102620 s`
- `stability(200步)`：
  - `BaseHeatConduction.get_derivatives`：`22.093438 s -> 22.186831 s`
  - `HeatConduction2D._compute_fluxes`：`15.501587 s -> 15.579570 s`
  - `HeatConduction2D._compute_internal_resistance`：`2.338184 s -> 2.340469 s`
- 结论：
  - 本轮优化的理论方向正确，但实际收益极小
  - 当前改善量接近 benchmark 波动水平，尚不足以视为实质性优化

### 8.3 结论

- 本轮修改安全，结果一致，可保留
- 但收益不显著，`get_jac_sparsity()` 不是当前主瓶颈
- 下一轮应继续按计划，只单独验证 `_compute_fluxes()` 末尾的 `flatten()` 替换

## 9. Round 2：`_compute_fluxes()` 末尾 `flatten()` 替换

- 修改目标：
  - 只替换 `_compute_fluxes()` 末尾 `Q_net_2d.flatten()`
  - 使用独立一维返回缓冲，避免每次 RHS 调用都重新分配
- 涉及文件：
  - `Solvers/HeatConduction/HeatConduction.py`
- 修改内容：
  - 在 `HeatConduction2D.__init__()` 中新增 `Q_net_flat_buffer`
  - 在 `_compute_fluxes()` 末尾用 `np.copyto(..., Q_net_2d.reshape(-1))` 返回一维缓冲
- 风险点：
  - 理论上物理无关，只影响返回数组分配方式
  - 若收益不稳定或更慢，不应保留这一步修改
- 验证结果：
  - `py_compile`：通过
  - `baseline`：通过
  - `stability(200步)`：通过

### 9.1 结果对比

- `baseline`：
  - 优化前：`0.219800 s/step`
  - 优化后：`0.221914 s/step`
  - 变化：`+0.002114 s/step`，约 `+0.96%`
- `stability(200步)`：
  - 优化前：`0.202924 s/step`
  - 优化后：`0.203000 s/step`
  - 变化：`+0.000075 s/step`，约 `+0.04%`
- 关键物理量对比：
  - `baseline` 末态的 `heated_layer_max_k`、`fluid_outlet_mean_k`、`radiation_power_w` 无变化
  - `stability` 尾段历史无变化
  - 结论：计算结果保持一致，没有引入数值错误

### 9.2 结论

- 这一步没有带来可确认的性能收益，且短 `baseline` 单次测试略慢
- 本轮修改不保留，代码已回退到修改前状态
- 下一轮不再继续围绕 `flatten()` 做深挖，转入 `_update_boundaries_state()` 的边界热阻临时数组复用
## 10. Round 3：`_update_boundaries_state()` 边界热阻临时数组复用

- 修改目标：
  - 只优化 `HeatConduction2D._update_boundaries_state()`
  - 尝试复用边界 `R_int` 临时数组，减少每次边界更新中的瞬时分配
- 涉及文件：
  - `Solvers/HeatConduction/HeatConduction.py`
- 实施内容：
  - 在 `HeatConduction2D` 中为 `left/right/top/bottom` 预分配 `R_int` scratch buffer
  - 使用 `np.multiply` / `np.maximum` / `np.divide` 原地写入边界热阻
- 风险点：
  - 若 `BoundaryRegion.update_internal_state()` 对输入数组存在延迟引用，scratch buffer 可能引入别名污染
  - 若收益不稳定，不应保留该复杂度
- 验证结果：
  - `py_compile`：通过
  - `smoke`：通过
  - `baseline`：通过
  - `stability(200步)`：通过

### 10.1 结果对比

- `baseline`
  - 优化前：`0.219350 s/step`
  - 优化后：`0.219003 s/step`
  - 变化：`-0.000348 s/step`，约 `-0.16%`
- `stability(200步)`
  - 优化前：`0.196740 s/step`
  - 优化后：`0.196945 s/step`
  - 变化：`+0.000205 s/step`，约 `+0.10%`
- 关键物理量：
  - `baseline` 末态 `heated_layer_max_k`、`fluid_outlet_mean_k`、`radiation_power_w` 不变
  - `stability` 尾段历史 `heated_layer_max_k`、`fluid_outlet_mean_k`、`radiation_power_w` 不变
  - 结论：本轮没有引入数值错误，也没有改变物理结果

### 10.2 profiler 对比

- `baseline`
  - `BaseHeatConduction.get_derivatives`：`0.956098 s -> 0.958809 s`
  - `HeatConduction2D._compute_fluxes`：`0.671392 s -> 0.672151 s`
  - `HeatConduction2D._compute_internal_resistance`：`0.101348 s -> 0.101425 s`
- `stability(200步)`
  - `BaseHeatConduction.get_derivatives`：`21.562895 s -> 21.539707 s`
  - `HeatConduction2D._compute_fluxes`：`15.169220 s -> 15.089627 s`
  - `HeatConduction2D._compute_internal_resistance`：`2.261565 s -> 2.266338 s`

### 10.3 结论

- 这一步物理上安全，但总耗时收益不成立
- 长回归 `stability(200步)` 略慢，不满足“可保留优化”的门槛
- 本轮修改不保留，代码已回退到修改前状态
- 下一轮转入 `_compute_internal_resistance()` 的中间数组复用

## 11. Round 4: `_compute_internal_resistance()` scratch reuse

- Scope:
  - Only changed `HeatConduction2D._compute_internal_resistance()`
  - Reused preallocated conductance buffers and work buffers
  - Did not change boundary handling, flux assembly, or solver flow
- Files:
  - `Solvers/HeatConduction/HeatConduction.py`
- Validation:
  - `py_compile`: pass
  - `smoke`: pass
  - `baseline`: pass
  - `stability(200 steps)`: pass

### 11.1 Result comparison

- `baseline`
  - before: `0.230001 s/step`
  - after: `0.217434 s/step`
  - delta: `-0.012567 s/step` (`-5.46%`)
- `stability(200 steps)`
  - before: `0.223109 s/step`
  - after: `0.216020 s/step`
  - delta: `-0.007088 s/step` (`-3.18%`)

### 11.2 Physics check

- Final `heated_layer_max_k`: unchanged
- Final `fluid_outlet_mean_k`: unchanged
- Final `radiation_power_w`: unchanged
- `stability` tail history: unchanged
- Conclusion: no numerical error and no physical drift introduced

### 11.3 Profiler comparison

- `baseline`
  - `BaseHeatConduction.get_derivatives`: `1.004269 s -> 0.946155 s`
  - `HeatConduction2D._compute_internal_resistance`: `0.107277 s -> 0.092356 s`
  - `HeatConduction2D._compute_fluxes`: `0.702890 s -> 0.669926 s`
- `stability(200 steps)`
  - `BaseHeatConduction.get_derivatives`: `24.639062 s -> 23.811306 s`
  - `HeatConduction2D._compute_internal_resistance`: `2.600637 s -> 2.349183 s`
  - `HeatConduction2D._compute_fluxes`: `17.293228 s -> 16.824543 s`

### 11.4 Decision

- Keep this round
- This is the first post-benchmark `HeatConduction.py` optimization with clear speedup under both short and long regression
- Next round should move to `BoundaryRegion.compute_net_flux_for_solver()` persistent buffers

## 12. Round 5: `BoundaryRegion.compute_net_flux_for_solver()` persistent buffers

- Scope:
  - Only changed `BoundaryRegion.compute_net_flux_for_solver()`
  - Reused persistent aggregation buffers and work arrays
  - Replaced repeated temporary array allocation with in-place operations
  - Did not change coupling order, boundary formulas, or radiation direction
- Files:
  - `Solvers/HeatConduction/Boundary.py`
- Validation:
  - `py_compile`: pass
  - `smoke`: pass
  - `baseline`: pass
  - `stability(200 steps)`: pass

### 12.1 Result comparison

- `baseline`
  - before: `0.225267 s/step`
  - after: `0.212206 s/step`
  - delta: `-0.013061 s/step` (`-5.80%`)
- `stability(200 steps)`
  - before: `0.230699 s/step`
  - after: `0.194828 s/step`
  - delta: `-0.035871 s/step` (`-15.55%`)

### 12.2 Physics check

- Final `heated_layer_max_k`: unchanged
- Final `fluid_outlet_mean_k`: unchanged
- Final `radiation_power_w`: unchanged
- `stability` tail history: unchanged
- Boundary audit values remain unchanged
- Conclusion: no numerical drift and no physical regression detected

### 12.3 Profiler comparison

- `baseline`
  - `BaseHeatConduction.get_derivatives`: `0.977132 s -> 0.917378 s`
  - `HeatConduction2D._compute_fluxes`: `0.689375 s -> 0.642439 s`
  - `HeatConduction2D._compute_internal_resistance`: `0.095923 s -> 0.091424 s`
- `stability(200 steps)`
  - `BaseHeatConduction.get_derivatives`: `25.259944 s -> 21.062128 s`
  - `HeatConduction2D._compute_fluxes`: `17.843292 s -> 14.731873 s`
  - `HeatConduction2D._compute_internal_resistance`: `2.468504 s -> 2.094852 s`

### 12.4 Decision

- Keep this round
- This round delivers clear speedup under both short and long regression
- The next hotspot is still `HeatConduction2D._compute_fluxes()`

## 13. Round 6: `_compute_fluxes()` flux buffer reuse

- Scope:
  - Only changed `HeatConduction2D._compute_fluxes()`
  - Added persistent `flux_x` / `flux_y` scratch buffers
  - Replaced per-call temporary flux arrays with in-place `np.subtract` and `np.multiply`
  - Did not change boundary flux assembly or return layout
- Files:
  - `Solvers/HeatConduction/HeatConduction.py`
- Validation:
  - `py_compile`: pass
  - `smoke`: pass
  - `baseline`: pass
  - `stability(200 steps)`: pass

### 13.1 Result comparison

- `baseline`
  - before: `0.214535 s/step`
  - after: `0.212206 s/step`
  - delta: `-0.002329 s/step` (`-1.09%`)
- `stability(200 steps)`
  - before: `0.204165 s/step`
  - after: `0.194828 s/step`
  - delta: `-0.009337 s/step` (`-4.57%`)

### 13.2 Physics check

- Final `heated_layer_max_k`: unchanged
- Final `fluid_outlet_mean_k`: unchanged
- Final `radiation_power_w`: unchanged
- `stability` tail history: unchanged
- Conclusion: no physical drift and no numerical regression

### 13.3 Profiler comparison

- `baseline`
  - `BaseHeatConduction.get_derivatives`: `0.927025 s -> 0.917378 s`
  - `HeatConduction2D._compute_fluxes`: `0.647464 s -> 0.642439 s`
  - `HeatConduction2D._compute_internal_resistance`: `0.092234 s -> 0.091424 s`
- `stability(200 steps)`
  - `BaseHeatConduction.get_derivatives`: `22.123591 s -> 21.062128 s`
  - `HeatConduction2D._compute_fluxes`: `15.462611 s -> 14.731873 s`
  - `HeatConduction2D._compute_internal_resistance`: `2.193400 s -> 2.094852 s`

### 13.4 Decision

- Keep this round
- The gain is moderate but consistent under both short and long regression
- The next optimization target should move out of `HeatConduction2D` inner arithmetic and into coupling overhead
