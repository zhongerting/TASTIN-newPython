# HydraulicNetwork 优化过程记录

## 1. 优化目标

- 优化目标文件：`Solvers/Hydrodynamics/HydraulicNetwork.py`
- 优化重点：提升单相流动换热求解的运行速度
- 约束条件：
  - 保持现有 `SystemManager` 通用接口兼容
  - 不引入固体导热、热管、RingHP 等额外物理模块干扰基准
  - 优化前后优先保证数值稳定性和物理量量级合理

## 2. 基准算例

- 基准脚本：`HydraulicNetwork优化重构/benchmark_open_loop_complex_300cv.py`
- 运行环境：`E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`
- 模型类型：开式回路、纯流体、使用 `SystemManager`
- 拓扑特征：
  - 入口边界 + 入口绝热缓冲段
  - 热端主干 + 分配腔
  - 4 个并联复杂支路
  - 每支路由 `8 + 16 + 16 + 8` 个控制体组成
  - 汇流腔 + 两段等效冷却主通道 + 回流段 + 出口绝热缓冲段 + 定压出口
- 规模：
  - 流体节点数：`302`
  - 流动连接数：`304`
- 热工加载：
  - 支路上下集流段施加人工 `Q_vol`
  - `Q_vol` 采用多频周期型热源，不再使用简单定值热源
  - 冷却主通道施加隐式冷却项 `Q = hA * T_sink - hA * T_fluid`

## 3. 当前基线配置

- 时间步长：`dt = 0.02 s`
- 步数：`n_steps = 2000`
- 内迭代：`inner_iter = 1`
- 入口温度：`968 K`
- 初始温度：`863 K`
- 出口压力：`1.61e5 Pa`
- 总质量流量：`2.2 kg/s`
- 重力：`0.0 m/s^2`
- 基础体积热源均值：`124000 W`
- 周期热源形式：多频正弦叠加
- 主调制幅值：`±22%`
- 次调制幅值：`±7%`
- 总冷却导热系数：`900 W/K`

## 4. 当前基线结果

### 4.1 运行结果

- 初始化耗时：`0.292870 s`
- 2000 步总耗时：`27.020079 s`
- 平均单步耗时：`0.013510 s/step`

### 4.2 最终状态摘要

- 入口质量流量：`2.200000 kg/s`
- 出口质量流量：`2.196360 kg/s`
- 支路流量分配：`[26.88, 25.50, 24.25, 23.38] %`
- 热端主干出口温度：`968.001 K`
- 冷却段 1 出口温度：`986.660 K`
- 冷却段 2 出口温度：`960.336 K`
- 出口缓冲段出口温度：`959.048 K`
- 最终步瞬时热源：`113162.994 W`
- 估算有效冷却功率：`141757.201 W`

### 4.3 Profiler 热点

按累计耗时排序：

| 函数 | 调用次数 | 总耗时 |
| --- | ---: | ---: |
| `HydraulicNetwork._calc_momentum_coeffs` | 3176 | `20.578968 s` |
| `HydraulicNetwork._assemble_pressure_system` | 3176 | `3.044892 s` |
| `HydraulicNetwork._step_energy_implicit` | 2000 | `1.014521 s` |
| `HydraulicNetwork._update_fluid_properties` | 2000 | `0.961917 s` |
| `HydraulicNetwork._update_flow_rates` | 3176 | `0.733231 s` |
| `HydraulicNetwork._solve_linear_system` | 3176 | `0.392596 s` |
| `SyntheticFluidSourceCoupler.execute` | 2000 | `0.132657 s` |

结论：

- 当前主要瓶颈明确集中在 `HydraulicNetwork` 内部
- 其中 `_calc_momentum_coeffs` 是绝对主热点
- `_assemble_pressure_system` 和 `_update_flow_rates` 也已经进入优先优化区

## 5. 已完成事项

- [x] 明确 `HydraulicNetwork.py` 静态热点分析方向
- [x] 创建纯流体 benchmark 算例
- [x] 保留 `SystemManager` 作为统一调度接口
- [x] 接入 `TEASAProfiler` 统计调用次数和累计耗时
- [x] 跑通 3 步 smoke test
- [x] 将定值热源调整为周期型热源
- [x] 跑通 2000 步 baseline benchmark

## 6. 当前优化顺序

第一阶段建议优先项：

1. 延后对象同步，减少 Picard/初始化阶段的对象回写
2. 补完整拓扑缓存，减少热路径中的对象跳转和属性访问
3. 优化 `_update_flow_rates()`
4. 优化 `_calc_momentum_coeffs()`

第二阶段候选项：

1. 压力矩阵装配结构复用
2. 降低数组重复分配与 `copy()`
3. 评估热膨胀源项和能量方程中的重复计算

## 7. 变更记录模板

后续每轮优化按下面格式追加：

### Round N

- 修改目标：
- 涉及文件：
- 修改内容：
- 风险点：
- 验证命令：
- benchmark 结果：
- profiler 变化：
- 结论：

## 8. 运行命令

### 8.1 语法检查

```powershell
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' -m py_compile 'HydraulicNetwork优化重构\benchmark_open_loop_complex_300cv.py'
```

### 8.2 默认基准运行

```powershell
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HydraulicNetwork优化重构\benchmark_open_loop_complex_300cv.py'
```

说明：

- 当前默认配置为 `2000` 步
- 当前默认热源为多频周期型热源

### 8.3 快速 smoke test

```powershell
$env:BENCH_N_STEPS='3'
$env:BENCH_PRINT_EVERY='1'
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HydraulicNetwork优化重构\benchmark_open_loop_complex_300cv.py'
```

## 9. 备注

- 当前 benchmark 的目的不是复现完整系统物理，而是稳定放大 `HydraulicNetwork` 的计算热点
- 后续每次优化后，默认先用同一 benchmark 配置回归，避免“优化了别的模块而不是 HydraulicNetwork”

## 10. Round 2

- 修改目标：`Group C`
- 主要改动：
  - 为压力方程预构建固定拓扑下的 CSR 非零结构，并缓存 data/RHS 映射
  - `_assemble_pressure_system()` 改为只更新矩阵 `data` 和右端项 `B`
  - 增加定压边界目标压力缓存刷新
  - `step_Picard()` / `step_Picard_over()` 在单个时间步内冻结 `S_thermal`
  - `_calc_enthalpy_time_derivative_explicit()` 改成基于 `np.add.at` 的向量散加
- 验证：
  - `py_compile` 通过
  - `3` 步 smoke test 通过
  - `inner_iter = 2` smoke test 通过
  - 默认 `2000` 步 benchmark 通过

### 10.1 Benchmark 对比（相对 Round 1）

- 初始化耗时：`0.043040 s -> 0.013804 s`，`-67.93%`
- 2000 步总耗时：`6.490997 s -> 3.340302 s`，`-48.54%`
- 平均单步耗时：`0.003245 s/step -> 0.001670 s/step`，`-48.54%`

### 10.2 Profiler 对比（相对 Round 1）

- `HydraulicNetwork.step_Picard`：`6.283436 s -> 3.148679 s`，`-49.89%`
- `HydraulicNetwork._assemble_pressure_system`：`2.958539 s -> 0.101819 s`，`-96.56%`
- `HydraulicNetwork._calc_momentum_coeffs`：`0.670127 s -> 0.523187 s`，`-21.93%`
- `HydraulicNetwork._solve_linear_system`：`0.335671 s -> 0.270914 s`，`-19.29%`
- `HydraulicNetwork._update_flow_rates`：`0.025572 s -> 0.015206 s`，`-40.54%`
- `HydraulicNetwork._step_energy_implicit`：`0.926639 s -> 0.908983 s`，`-1.91%`
- `HydraulicNetwork._update_fluid_properties`：`0.994081 s -> 0.931701 s`，`-6.28%`

### 10.3 结论

- 第二轮已经把压力矩阵组装从主要热点中基本移除
- 当前新的主要耗时已经转移到：
  - `_update_fluid_properties`
  - `_step_energy_implicit`
  - `_calc_momentum_coeffs`
- 下一轮建议进入 `Group D`

## 11. Benchmark Case Update：周期入口流量

- 修改目标：在现有“周期热源”基准上，进一步引入入口流量周期扰动，避免回路过快进入近稳态响应，增强 `HydraulicNetwork` 在非定常流动条件下的 benchmark 压力。
- 修改文件：`HydraulicNetwork优化重构/benchmark_open_loop_complex_300cv.py`
- 修改内容：
  - 入口不再保持固定 `target_W`
  - 通过 `SyntheticFluidSourceCoupler` 在每个时间步按 `SystemManager.global_time` 周期性更新 `InletJunction.target_W`
  - benchmark 输出增加 `W_target`，用于直接对比入口实际流量和目标流量

### 11.1 默认周期入口流量参数

- `BENCH_FLOW_PRIMARY_AMP = 0.16`
- `BENCH_FLOW_SECONDARY_AMP = 0.05`
- `BENCH_FLOW_PRIMARY_PERIOD = 6.5 s`
- `BENCH_FLOW_SECONDARY_PERIOD = 2.4 s`
- `BENCH_FLOW_MIN_SCALE = 0.72`
- `BENCH_FLOW_MAX_SCALE = 1.28`

说明：
- 周期入口流量采用双频叠加后再截断到最小/最大倍率范围
- 基准总步数继续保持 `2000` 步
- 周期热源配置保持不变，当前 benchmark 为“周期热源 + 周期入口流量”联合扰动

### 11.2 测试命令

语法检查：

```powershell
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' -m py_compile 'HydraulicNetwork优化重构\benchmark_open_loop_complex_300cv.py'
```

`3` 步 smoke test（`inner_iter = 1`）：

```powershell
$env:BENCH_N_STEPS='3'
$env:BENCH_PRINT_EVERY='1'
Remove-Item Env:BENCH_INNER_ITER -ErrorAction SilentlyContinue
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HydraulicNetwork优化重构\benchmark_open_loop_complex_300cv.py'
```

`3` 步 smoke test（`inner_iter = 2`）：

```powershell
$env:BENCH_N_STEPS='3'
$env:BENCH_PRINT_EVERY='1'
$env:BENCH_INNER_ITER='2'
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HydraulicNetwork优化重构\benchmark_open_loop_complex_300cv.py'
```

默认 `2000` 步 benchmark：

```powershell
Remove-Item Env:BENCH_N_STEPS -ErrorAction SilentlyContinue
Remove-Item Env:BENCH_PRINT_EVERY -ErrorAction SilentlyContinue
Remove-Item Env:BENCH_INNER_ITER -ErrorAction SilentlyContinue
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' 'HydraulicNetwork优化重构\benchmark_open_loop_complex_300cv.py'
```

### 11.3 测试结果

- `py_compile`：通过
- `3` 步 smoke test（`inner_iter = 1`）：通过
- `3` 步 smoke test（`inner_iter = 2`）：通过
- 默认 `2000` 步 benchmark：通过

`2000` 步 benchmark 结果：
- 初始化时间：`0.013693 s`
- 总耗时：`3.737493 s`
- 平均单步时间：`0.001869 s/step`

入口流量跟踪检查：
- `step 200`：`W_in = 1.8795 kg/s`，`W_target = 1.8713 kg/s`
- `step 400`：`W_in = 2.6466 kg/s`，`W_target = 2.6447 kg/s`
- `step 2000`：`W_in = 2.3935 kg/s`，`W_target = 2.3944 kg/s`

最终摘要：
- `Inlet mass flow = 2.393517 kg/s`
- `Inlet target flow = 2.394428 kg/s`
- `Outlet mass flow = 2.393839 kg/s`

### 11.4 Profiler 摘要

- `SystemManager.step`：`3.732047 s`
- `HydraulicNetwork.step_Picard`：`3.540540 s`
- `HydraulicNetwork._update_fluid_properties`：`0.939479 s`
- `HydraulicNetwork._step_energy_implicit`：`0.890159 s`
- `HydraulicNetwork._calc_momentum_coeffs`：`0.763366 s`
- `HydraulicNetwork._solve_linear_system`：`0.379334 s`
- `HydraulicNetwork._assemble_pressure_system`：`0.147024 s`
- `SyntheticFluidSourceCoupler.execute`：`0.125631 s`
- `HydraulicNetwork._update_flow_rates`：`0.022319 s`

### 11.5 结论

- benchmark 已从“周期热源”升级为“周期热源 + 周期入口流量”联合扰动
- 入口流量对目标流量的跟踪正常，说明当前液压迭代路径在周期流量边界下工作稳定
- 在第二轮优化后的代码基础上，默认 `2000` 步 benchmark 仍能保持约 `1.87 ms/step`
- 当前热点排序没有本质变化，`Group D` 仍应优先处理：
  - `_update_fluid_properties`
  - `_step_energy_implicit`
  - 剩余 `_calc_momentum_coeffs`

## 12. Round 3

- 修改目标：`Group D`
- 主要改动：
  - 为能量方程新增固定拓扑下的 CSR 结构缓存，`_step_energy_implicit()` 不再每步 `list -> coo -> csr`
  - 为热工路径增加 scratch buffer，复用 `mass/cp/lam_h/Q_mod/RHS/T_old/h_old`
  - `_calc_enthalpy_time_derivative_explicit()` 改为复用焓导数、流量通量、质量和热源缓冲区
  - `_calc_thermal_expansion_source()` 改为基于复用缓冲区计算，避免重复分配
  - `_update_fluid_properties()` 去掉 `vol.P/T` 的步内对象回写，仅保留 `rho/mu` 兼容同步和流体源项提取
  - 将物性函数能力检测缓存为 `primary_material` 相关标志，减少 `hasattr(...)` 热路径开销
- 验证：
  - `py_compile` 通过
  - `3` 步 smoke test（`inner_iter = 1`）通过
  - `3` 步 smoke test（`inner_iter = 2`）通过
  - 默认 `2000` 步 benchmark 通过

### 12.1 Benchmark 对比（相对 Group D 前基线）

- 初始化耗时：`0.013427 s -> 0.013477 s`，`+0.37%`
- `2000` 步总耗时：`3.802042 s -> 2.956531 s`，`-22.24%`
- 平均单步耗时：`0.001901 s/step -> 0.001478 s/step`，`-22.25%`

### 12.2 Profiler 对比（相对 Group D 前基线）

- `HydraulicNetwork.step_Picard`：`3.602237 s -> 2.766202 s`，`-23.21%`
- `HydraulicNetwork._update_fluid_properties`：`0.953852 s -> 0.859642 s`，`-9.88%`
- `HydraulicNetwork._step_energy_implicit`：`0.907468 s -> 0.246120 s`，`-72.88%`
- `HydraulicNetwork._calc_momentum_coeffs`：`0.778078 s -> 0.753152 s`，`-3.20%`
- `HydraulicNetwork._solve_linear_system`：`0.378691 s -> 0.356913 s`，`-5.75%`
- `HydraulicNetwork._assemble_pressure_system`：`0.150408 s -> 0.143945 s`，`-4.30%`
- `HydraulicNetwork._update_flow_rates`：`0.022750 s -> 0.022138 s`，`-2.69%`
- `SyntheticFluidSourceCoupler.execute`：`0.128168 s -> 0.121553 s`，`-5.16%`

### 12.3 当前基准结果

- 初始化时间：`0.013477 s`
- `2000` 步总耗时：`2.956531 s`
- 平均单步时间：`0.001478 s/step`

关键状态检查：
- `Inlet mass flow = 2.393517 kg/s`
- `Inlet target flow = 2.394428 kg/s`
- `Outlet mass flow = 2.393839 kg/s`
- `step 200`：`W_in = 1.8795 kg/s`，`W_target = 1.8713 kg/s`
- `step 400`：`W_in = 2.6466 kg/s`，`W_target = 2.6447 kg/s`
- `step 2000`：`W_in = 2.3935 kg/s`，`W_target = 2.3944 kg/s`

### 12.4 结论

- `Group D` 的核心收益已经兑现，最大提速来自 `_step_energy_implicit()` 的矩阵结构复用
- 在周期热源和周期入口流量同时存在的 benchmark 下，物理结果与 Group D 前保持一致，未出现流量跟踪或温度场回归
- 当前新的主要热点已经转移为：
  - `_update_fluid_properties`
  - `_calc_momentum_coeffs`
  - `_solve_linear_system`
- 下一轮可以优先评估：
  - 是否继续压缩 `_update_fluid_properties` 的源项提取与物性更新开销
  - 是否进一步瘦身 `_calc_momentum_coeffs`
  - 是否进入 `Group E` 评估求解器层复用空间

## 13. Round 4

- 修改目标：`Group E`
- 工作内容：
  - 评估压力方程求解器替换或复用空间
  - 在不改变物理路径的前提下，对以下方案做 benchmark 级实测：
    - `spsolve` 基线
    - `splu(..., permc_spec='COLAMD')`
    - `splu(..., permc_spec='NATURAL')`
    - `bicgstab`（初始化阶段直接接管）
    - `bicgstab + 严格残差检查 + spsolve 回退`（仅瞬态阶段）
    - `gmres + 严格残差检查 + spsolve 回退`（仅瞬态阶段）
    - 原生 `CSC` 压力矩阵缓存 + `spsolve`
    - dense solve
- 代码改动：
  - 保持 `_solve_linear_system()` 的 `spsolve` 路径不变
  - 清理 `_solve_linear_system()` 的过时说明，明确 Group E 的结论

### 13.1 实测结果

候选求解器实测结论：
- `splu(COLAMD)`：比当前 `spsolve` 更慢
- `splu(NATURAL)`：比当前 `spsolve` 更慢
- `bicgstab` 直接接管：会破坏初始化收敛，不可接受
- `bicgstab + 回退`：即使只在瞬态阶段启用，也明显更慢
- `gmres + 回退`：显著更慢
- `CSC` 原生缓存 + `spsolve`：与当前 `CSR` 缓存 + `spsolve` 基本持平，无实质收益
- dense solve：明显更慢

### 13.2 结论

- 在当前固定拓扑、约 `300` 节点、周期热源 + 周期入口流量 benchmark 下，`spsolve` 仍然是最快且最稳定的压力求解路径
- `Group E` 没有产出值得并入主代码的 solver 替换方案
- 本轮的正确工程决策是：
  - 保持当前 `spsolve` 直解路径
  - 不引入实验性迭代求解器分支
  - 把后续优化重心重新放回 `_update_fluid_properties` 和 `_calc_momentum_coeffs`

## 14. Restart 边界集合一致性检查

- 修改目标：补齐 restart 场景下的防御性校验，避免“定压边界成员变化”在加载后静默错算
- 适用范围：只检查定压边界集合的一致性，不影响已有定压边界节点的 `target_P` 数值恢复
- 代码改动：
  - 在 `HydraulicNetwork.get_state_dict()` 中写入 `fixed_pressure_idx` 指纹
  - 在 `HydraulicNetwork.load_state_dict()` 中校验 restart 文件与当前模型的定压边界集合是否一致
  - 若集合不一致，抛出明确的 `ValueError`
  - 若旧 restart 文件没有该指纹，仅给出 warning，保留向后兼容
  - 在恢复 `target_P / target_W` 后，立即刷新缓存的边界目标量

### 14.1 设计结论

- 仍然支持在运行过程中修改“已有定压边界”的 `target_P`
- 不支持在 restart 前后修改“哪些节点属于定压边界”
- 这样做的原因是压力矩阵和能量矩阵的缓存结构默认假定定压边界集合固定

### 14.2 验证

- `py_compile` 通过
- 同一定压边界集合的 restart 恢复通过
- 人工清空 `fixed_pressure_indices` 后再加载 restart，能被明确拦截
- `3` 步 benchmark smoke test 通过，说明新增检查未破坏标准运行路径

### 14.3 风险收敛效果

- 把原先“边界集合变化后可能静默错算”的风险，收敛为“启动恢复时立即失败并报清晰错误”
- 对历史 restart 文件保持兼容，不会因为缺少新字段而直接无法读取

### 14.4 运行时保护

- `fixed_pressure_indices` 已改为受保护集合
- 初始化完成后，任何运行时成员修改都会在修改当下直接报错
- 以下路径均已验证会被拦截：
  - `fixed_pressure_indices.clear()`
  - `fixed_pressure_indices.add(...)`
  - `fixed_pressure_indices.remove(...)`
  - `fixed_pressure_indices.update(...)`
  - `network.fixed_pressure_indices = set(...)`
- 标准 benchmark 路径不受影响，`3` 步 smoke test 通过
