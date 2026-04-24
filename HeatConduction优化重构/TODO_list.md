# HeatConduction 安全实施 TODO List

## 1. 目标

- 目标文件：`Solvers/HeatConduction/HeatConduction.py`
- 基准脚本：`HeatConduction优化重构/benchmark_systemmanager_heatconduction.py`
- 记录文件：`HeatConduction优化重构/optimize_process.md`
- 核心约束：
  - 绝对优先保证物理计算准确
  - 优先压缩 `HeatConduction2D`、边界更新、层间耦合、流固耦合的总耗时
  - 不允许为了提速而弱化半节点热阻、边界热流、动态辐射或耦合路径

## 2. Benchmark 固定结构

- `SystemManager` 统一调度
- `4` 根并联流道，拓扑参考 `testModule/test_parallel_channels.py`
- 每根流道外侧 `4` 层二维柱坐标固体：
  - `inner_liner`
  - `transition_shell`
  - `heated_shell`
  - `outer_shield`
- 层间使用 `SolidSolidCouple2D`
- 最内层与流道使用 `FluidSolidCouple`
- 最外层边界使用动态辐射边界
- 发热层使用稳态功率 + 周期扰动
- 轴向分段、温度和流量量级参考 `testModule/test_core_assemble_v6.py`

## 3. 当前热点判断

预期高优先级热点：

1. `BaseHeatConduction.get_derivatives`
2. `HeatConduction2D._compute_internal_resistance`
3. `HeatConduction2D._update_boundaries_state`
4. `HeatConduction2D._compute_fluxes`
5. `FluidSolidCouple.execute`
6. `SolidSolidCouple2D.sync`
7. `SystemManager.step`

说明：

- benchmark 中流体网络故意保持最小真实复杂度，避免把主瓶颈重新拉回 `HydraulicNetwork`
- 如果 profiler 显示流体侧重新主导，需要先判断是 benchmark 参数失衡还是导热优化尚未进入主热点

## 4. 强制实施规则

### 规则 A：先冻结 benchmark，再改求解器

- 在没有固定 benchmark 配置、边界审计和能量审计之前，不进入正式性能优化
- 每轮优化必须使用同一 benchmark 配置回归

### 规则 A+：每轮只改一个思路

- 每轮只允许一个优化思路进入代码
- 每轮只允许一个主代码文件发生功能性修改
- 同一轮内禁止同时改 `HeatConduction2D` 内部热阻、边界更新和返回缓冲

### 规则 B：边界半节点处理不能动成“近似写法”

- 不能为了提速跳过 `left/right/top/bottom` 的 `R_internal`
- 不能把动态辐射简化为常数 `h_rad`
- 不能把层间耦合改成只传温度不传热阻

### 规则 C：优化必须按链路成组推进

- `_compute_internal_resistance`、`_update_boundaries_state`、`_compute_fluxes` 三者不能孤立修改后长期停留
- 若修改边界缓存策略，必须同步检查 `FluidSolidCouple` 和 `SolidSolidCouple2D`

## 5. 建议分组

### Group A：Benchmark 与审计冻结

- [x] 建立 `SystemManager` 主 benchmark
- [x] 固定并联流道、多层固体、动态辐射、发热层扰动
- [x] 输出 profiler、边界审计、能量审计
- [x] 完成短 `baseline` 结果登记
- [x] 建立 `stability(200步)` 长回归档并完成首版结果登记

### Group B：HeatConduction2D 内部热阻路径

- [ ] 评估 `k_face`、`G_x_inner`、`G_y_inner` 的重复分配
- [ ] 评估 `_compute_internal_resistance()` 中可复用的几何项
- [x] 检查 `reshape`、临时数组、`flatten()` 的残余开销
- [ ] 将 `HeatConduction.py` 的优化拆成单项小步回归，禁止成组修改

### Group C：边界状态与边界热流路径

- [ ] 评估 `_update_boundaries_state()` 中重复读取几何矩阵的开销
- [ ] 评估 `BoundaryRegion.compute_net_flux_for_solver()` 的重复数组构造
- [ ] 仅在独立一轮中尝试 `BoundaryRegion` 持久缓冲
- [ ] 检查动态辐射路径是否存在可缓存但不破坏正确性的中间量

### Group D：耦合层

- [ ] 评估 `FluidSolidCouple.execute()` 的物性访问和边界温度读取开销
- [ ] 评估 `SolidSolidCouple2D.sync()` 的重复参数更新开销
- [ ] 检查 `SystemManager._sync_solid_boundaries_for_coupling()` 是否存在重复工作

### Group E：求解器层

- [x] 评估 `HeatConduction2D` 的 `jac_sparsity` 是否需要缓存复用
- [ ] 评估多固体逐个 `solve_ivp(BDF)` 的总开销是否值得做更深层重构
- [ ] 在确保物理与接口不变前，不提前尝试实验性替代求解路径

## 6. 不允许单独推进的事项

- [ ] 只改 `_compute_fluxes()`，不检查边界与热阻路径
- [ ] 只改 `FluidSolidCouple.execute()`，不验证壁面温度与边界热阻读写
- [ ] 只改 `SolidSolidCouple2D.sync()`，不验证层间热流连续性
- [ ] 在 benchmark 尚未冻结前直接大改 `HeatConduction.py`

## 7. 每轮优化后必须检查

- [ ] `py_compile` 通过
- [ ] `smoke` benchmark 通过
- [ ] `baseline` benchmark 通过
- [ ] `stability(200步)` benchmark 通过
- [ ] profiler 主热点仍然集中在导热与耦合链路
- [ ] 边界 `R_internal` 审计无异常跳变
- [ ] 动态辐射热流方向正确
- [ ] 发热层、外壁、流体出口温度量级合理
- [ ] 200 步后段历史无明显失稳或异常漂移
- [ ] 能量残差未显著恶化

## 8. 当前进度

### 已完成

- [x] 建立 `HeatConduction优化重构/benchmark_systemmanager_heatconduction.py`
- [x] 建立 `HeatConduction优化重构/TODO_list.md`
- [x] 建立 `HeatConduction优化重构/optimize_process.md`
- [x] 完成回退后基线重建
- [x] 固定“每轮只改一个思路”的实施规则
- [x] 完成 `stability(200步)` 长回归配置与首轮验证
- [x] 完成 `get_jac_sparsity()` 缓存单步验证
- [x] 完成 `_compute_fluxes()` 末尾 `flatten()` 替换单步验证，并判定不保留

### 下一轮优先级

1. 单独验证 `_update_boundaries_state()` 的边界热阻临时数组复用
2. 单独验证 `_compute_internal_resistance()` 的中间数组复用
3. 单独验证 `BoundaryRegion.compute_net_flux_for_solver()` 的持久缓冲
### Round 3 记录

- [x] 完成 `_update_boundaries_state()` 的边界热阻临时数组复用单步验证，并判定不保留

### 当前下一轮优先级（更新）

1. 单独验证 `_compute_internal_resistance()` 的中间数组复用
2. 单独验证 `BoundaryRegion.compute_net_flux_for_solver()` 的持久缓冲
3. 单独验证 `FluidSolidCouple.execute()` 中的重复物性访问精简

### Round 4

- [x] Validate `_compute_internal_resistance()` scratch reuse
- [x] Keep Round 4 in main code path

### Next Priority

1. Validate `BoundaryRegion.compute_net_flux_for_solver()` persistent buffers
2. Validate `FluidSolidCouple.execute()` repeated property access cleanup
3. Re-check whether `SystemManager` has duplicated solid-state pre-sync work

### Round 5

- [x] Validate `BoundaryRegion.compute_net_flux_for_solver()` persistent buffers
- [x] Keep Round 5 in main code path

### Next Priority

1. Validate `HeatConduction2D._compute_fluxes()` temporary flux-array reuse
2. Validate `FluidSolidCouple.execute()` repeated property access cleanup
3. Re-check whether `SystemManager` still does duplicated solid pre-sync work

### Round 6

- [x] Validate `_compute_fluxes()` flux buffer reuse
- [x] Keep Round 6 in main code path

### Next Priority

1. Validate `FluidSolidCouple.execute()` repeated property access cleanup
2. Re-check whether `SystemManager` still does duplicated solid pre-sync work
3. Evaluate whether any residual boundary-state updates are still duplicated
