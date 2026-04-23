# HydraulicNetwork 安全实施 TODO List

## 1. 目标

- 目标文件：`Solvers/Hydrodynamics/HydraulicNetwork.py`
- 优化目标：在不破坏现有 `SystemManager` 接口和数值稳定性的前提下，优先降低 `HydraulicNetwork` 总运行耗时
- 当前基准：`HydraulicNetwork优化重构/benchmark_open_loop_complex_300cv.py`
- 当前默认基准配置：`dt = 0.02 s`，`n_steps = 2000`，`inner_iter = 1`

## 2. 总体判断

- 当前任务之间没有“方向性冲突”
- 主要风险不在于目标互相矛盾，而在于实现时出现“半改状态”
- 真正需要控制的有两条主线：
  - 状态同步策略
  - 稀疏矩阵结构复用策略

换句话说：

- 可以全部做
- 但不能随意拆开做

## 3. 两条强制实施规则

### 规则 A：迭代内部必须以向量状态为唯一真源

在 Picard 内部和初始化迭代内部，必须统一以下原则：

- 真状态：`P_vec / W_vec / T_vec / h_vec`
- 对象状态：`vol.P / vol.T / vol.h / junc.W / junc.vel`
- 若采用“延后对象同步”策略，则对象状态只能在时间步结束后统一刷新

否则会出现以下冲突：

- 一部分函数读向量状态
- 另一部分函数仍读对象状态
- 同一时间步内出现“半旧半新”混用

### 规则 B：矩阵结构复用只允许在固定拓扑下启用

只有满足下面条件，才允许预构建并复用稀疏矩阵结构：

- 节点数固定
- junction 数固定
- 边界节点集合固定
- 矩阵非零模式不变

当前 benchmark 满足这些条件，因此可以安全做结构复用。

## 4. 当前热点排序

按最新 baseline profiler 排序：

1. `HydraulicNetwork._calc_momentum_coeffs`
2. `HydraulicNetwork._assemble_pressure_system`
3. `HydraulicNetwork._step_energy_implicit`
4. `HydraulicNetwork._update_fluid_properties`
5. `HydraulicNetwork._update_flow_rates`
6. `HydraulicNetwork._solve_linear_system`

## 5. 强绑定关系

下面这些任务不能拆开独立推进。

### 强绑定组 1：状态同步与液压迭代

核心绑定任务：

- `3. _update_flow_rates`
- `4. step_Picard`
- `7. initialize_hydraulics`

同轮必须做兼容性审查的任务：

- `5. _update_fluid_properties`
- `8. _step_energy_implicit`

原因：

- `3/4/7` 决定“迭代内部谁是真状态”
- 一旦改成延后对象同步，`5/8` 里所有对象写回路径都必须同步检查

结论：

- `3` 不能单做
- `4` 不能单做
- `7` 不能单做

### 强绑定组 2：压力方程与热膨胀源项

核心绑定任务：

- `2. _assemble_pressure_system`
- `6. _calc_thermal_expansion_source`

桥接任务：

- `9. _calc_enthalpy_time_derivative_explicit`

原因：

- `2` 每轮都会调用热膨胀源项
- 若 `6` 改成时间步内冻结，而 `9` 仍保留旧的对象扫描或旧的数据路径，就会形成新的不一致

结论：

- `2` 和 `6` 必须捆绑推进
- 如果 `9` 仍依赖对象扫描，就必须在同轮至少做最小兼容修正

### 强绑定组 3：物性更新与能量方程

核心绑定任务：

- `5. _update_fluid_properties`
- `8. _step_energy_implicit`

原因：

- 这两个函数都直接碰 `T_vec / h_vec / vol.T / vol.h`
- 它们必须遵守同一套对象同步策略

结论：

- 若 Group B 已经改为延后对象同步，则 `5` 和 `8` 必须按同样规则改

### 强绑定组 4：求解器层

任务：

- `10. _solve_linear_system`

原因：

- 它依赖前面的矩阵结构、边界处理和状态路径都已经稳定

结论：

- 必须最后做

## 6. 可安全实施的分组版本

下面是建议采用的正式实施分组。

### Group A：基础缓存层

包含任务：

- `1. _calc_momentum_coeffs` 的前置缓存补全部分

本组目标：

- 补齐所有后续向量化要用到的拓扑和几何缓存
- 不改变求解流程
- 不改变对象同步策略

本组可单独实施：

- 是

本组建议内容：

- [ ] 增加 `idx_from_arr`
- [ ] 增加 `idx_to_arr`
- [ ] 增加 `is_inlet_junction_mask`
- [ ] 增加 `target_W_arr`
- [ ] 增加 `D_up_arr`
- [ ] 增加 `D_down_arr`
- [ ] 增加 `L_up_arr`
- [ ] 增加 `L_down_arr`
- [ ] 增加 `A_in_node_arr`
- [ ] 增加 `A_out_node_arr`
- [ ] 统一 `z_vec` 的使用口径

进入条件：

- 无

完成条件：

- 不改变 benchmark 结果趋势
- 所有下游函数都可以只依赖缓存而不是频繁跳对象

### Group B：状态同步与液压迭代层

包含任务：

- `3. _update_flow_rates`
- `4. step_Picard`
- `7. initialize_hydraulics`

本组必须同时审查：

- `5. _update_fluid_properties`
- `8. _step_energy_implicit`

本组目标：

- 统一“向量状态是真状态”的规则
- 将对象同步从迭代内部移到时间步结束
- 消除每轮 Picard 的大量对象回写

本组禁止拆分：

- [ ] 不允许只改 `_update_flow_rates`
- [ ] 不允许只改 `step_Picard`
- [ ] 不允许只改 `initialize_hydraulics`

本组建议内容：

- [ ] `_update_flow_rates` 改成向量计算
- [ ] 增加 `_sync_state_to_objects()`
- [ ] `step_Picard` 内部只操作向量
- [ ] `initialize_hydraulics` 与 `step_Picard` 统一策略
- [ ] 去掉每轮 Picard 内的对象回写
- [ ] 用 scratch buffer 替代重复 `copy()`
- [ ] 审查 `_update_fluid_properties` 是否仍隐式依赖对象最新状态
- [ ] 审查 `_step_energy_implicit` 是否仍在步内提前写对象

进入条件：

- Group A 完成

完成条件：

- Picard 和初始化阶段内部只依赖向量状态
- 时间步结束后对象状态与向量状态一致

### Group C：压力方程层

包含任务：

- `2. _assemble_pressure_system`
- `6. _calc_thermal_expansion_source`

可能一起带上的桥接任务：

- `9. _calc_enthalpy_time_derivative_explicit`

本组目标：

- 复用压力矩阵非零结构
- 冻结时间步内热膨胀源项
- 避免每轮 Picard 重复组装完整结构

本组禁止拆分：

- [ ] 不允许只改 `_assemble_pressure_system` 而不处理热膨胀源项调用策略
- [ ] 不允许只改 `_calc_thermal_expansion_source` 而不改其上游调用方式

本组建议内容：

- [ ] 预构建压力矩阵非零结构
- [ ] 只更新 `data` 和 `B`
- [ ] 将 `fixed_pressure_indices` 相关逻辑改成 mask + 预缓存位置
- [ ] 在 `step_Picard` 中冻结 `S_thermal`
- [ ] 若 `_calc_enthalpy_time_derivative_explicit` 仍对象扫描，则一起做最小向量化兼容

进入条件：

- Group A 完成
- Group B 已稳定

完成条件：

- 压力矩阵结构不再每轮重建
- 热膨胀源项不再每轮重复计算

### Group D：物性与能量方程层

包含任务：

- `5. _update_fluid_properties`
- `8. _step_energy_implicit`
- `9. _calc_enthalpy_time_derivative_explicit` 的完整版优化

本组目标：

- 清理对象写回路径
- 复用能量矩阵结构
- 减少大数组重复分配

本组禁止拆分：

- [ ] 不建议只改 `_step_energy_implicit`，但保持旧的对象同步方式
- [ ] 不建议只改 `_update_fluid_properties` 的回写策略，而不审查能量方程

本组建议内容：

- [ ] 将物性评估和对象回写拆开
- [ ] 使用切片赋值替代重新绑定
- [ ] 预构建焓方程稀疏结构
- [ ] 去掉 `T_old_vec.copy()` 和 `h_old_vec.copy()` 的不必要分配
- [ ] 完整向量化 `_calc_enthalpy_time_derivative_explicit`
- [ ] 统一由 `_sync_state_to_objects()` 负责最终对象同步

进入条件：

- Group B 完成
- Group C 完成

完成条件：

- 能量方程和物性更新都遵守统一同步规则
- 能量矩阵结构复用正常工作

### Group E：求解器层

包含任务：

- `10. _solve_linear_system`

本组目标：

- 只在前面结构和状态路径稳定后，再评估求解器替换或分解复用

本组禁止拆分：

- [ ] 不允许提前于 Group C / Group D 动求解器

本组建议内容：

- [ ] 评估 `factorized()` 或分解复用
- [ ] 评估是否值得切换迭代法
- [ ] 仅在矩阵性质确认稳定后再尝试

进入条件：

- Group C 完成
- Group D 完成

完成条件：

- 新求解器策略不破坏稳定性
- 确认收益明显大于维护成本

## 7. 不允许单独实施的任务

下面这些任务，单独做风险高：

- [ ] `3. _update_flow_rates` 单独实施
- [ ] `4. step_Picard` 单独实施
- [ ] `6. _calc_thermal_expansion_source` 单独实施
- [ ] `8. _step_energy_implicit` 在 Group B 之前单独实施
- [ ] `10. _solve_linear_system` 提前实施

## 8. 推荐实施顺序

1. Group A：基础缓存层
2. Group B：状态同步与液压迭代层
3. Group C：压力方程层
4. Group D：物性与能量方程层
5. Group E：求解器层

## 9. 每轮优化后必须检查

- [ ] `benchmark_open_loop_complex_300cv.py` 的 3 步 smoke test 是否通过
- [ ] 默认 2000 步 benchmark 是否通过
- [ ] `TEASAProfiler` 热点是否仍主要集中在 `HydraulicNetwork`
- [ ] 入口 / 出口流量是否仍量级合理
- [ ] 周期热源是否仍正常波动
- [ ] 出口温度和各支路温度是否仍在合理范围内

## 10. 备注

- 当前 TODO 已经不是单纯的函数清单，而是“可安全实施的分组清单”
- 后续执行时，优先以 Group 为单位推进，不建议再按单个函数零散修改
