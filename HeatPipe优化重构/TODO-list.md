# HeatPipe / RingHP TODO List

## 1. 范围

- 目标文件：
  - `Components/HPwithFin.py`
  - `Components/basicComponents/HeatPipe2D.py`
  - `Components/RingHP.py`
- 基线脚本：`HeatPipe优化重构/benchmark_ringhp_v5.py`
- 参考算例：`CoolantLoop/test_coolant_loop_v5.py`

## 2. 当前基线

- 基线配置：`dt = 0.05 s`，`n_steps = 20`，`inner_iter = 1`
- 基线输出：`HeatPipe优化重构/ringhp_v5_baseline_history.csv`
- 本轮 profiler 热点：

| 函数 | 调用次数 | 总耗时 |
| --- | ---: | ---: |
| `HeatPipe2D._update_properties` | 8769 | `2.111064 s` |
| `HeatPipe2D._compute_fluxes` | 8286 | `1.603796 s` |
| `HeatPipe2D._update_boundaries_state` | 8769 | `0.734677 s` |
| `RingHP.pre_step` | 40 | `0.692757 s` |
| `HPwithFin.pre_step` | 460 | `0.692228 s` |
| `HPwithFin._solve_fin_quasi_steady` | 920 | `0.666283 s` |

## 3. 静态分析发现

### P0

- [ ] 修复 `RingHP.py:215` 的闭包绑定错误。`constant_h_corr()` 在循环内定义但闭包捕获的是同一个 `proxy` 变量，循环结束后所有 `coupler_hp` 都会读到最后一个节点的 `proxy` 状态。建议改成工厂函数或 `functools.partial`，把当前节点的 `proxy` / `d_h` 显式绑定进去。
- [ ] 修复 `HPwithFin.py:297` 的翅片导热系数硬编码。当前 `k_fin_array = 348.9` 完全绕过了 `wall_mat` / `self.k_fin_mat`，一旦热管壁材不是这一个常数，翅片求解会直接失真。建议改为基于当前温度场调用材料导热率，并确认是否需要温度相关迭代。

### P1

- [ ] 修正 `RingHP.py:435` 和 `RingHP.py:443` 的总量统计口径。`get_total_heat_rejection()` 和 `get_total_external_heat_absorption()` 只累加代表热管单元，没有乘回 `hp_multipliers`，与 `SingleVolumeProxy.add_coupling_source_distribution()` 的总量折算口径不一致。建议统一为“真实节点总量”，或把当前接口改名为“representative-unit total”避免误用。
- [ ] 复核 `RingHP.py:194-207` 的阻塞/局阻模型。`A_proj_total` 和 `sigma` 按单根热管计算，但 `K_loss_val` 又通过 `N_eff` 部分引入 `N_hp`，导致“几何校核”和“阻力计算”不在同一物理口径上。建议先明确单排等效还是多排等效，再统一投影面积、孔隙率和阻力公式。
- [ ] 给 `HPwithFin.py:348-351` 的 Thomas 追赶过程增加数值保护。当前 `denom` 没有下限检查，若辐射线性化后主对角接近奇异，可能直接产生 `inf/nan` 并污染整个翅片等效热阻。建议增加最小阈值和回退策略。

### P2

- [ ] 减少 `HeatPipe2D.py:137-162` 的逐次分配。`_update_properties()` 每次都新建 `k_2d/rho_2d/cp_2d` 三个二维数组，再 flatten 回写，是目前最大的单点热点。建议引入类内 scratch buffer，并尽量原地更新。
- [ ] 减少 `HeatPipe2D.py:275-330` 的中间数组开销。`_compute_fluxes()` 每次都新建 `Q_net_2d` 并做多次切片累加，累计调用次数高，值得和 `_update_properties()` 一起做 buffer 化。
- [ ] 精简 `HPwithFin.py:323-356` 的矩阵临时量。当前每次翅片求解都重新申请 `a/b/c/d/c_prime/d_prime/T_new`，在 `920` 次调用下已形成稳定热点。建议把 Thomas 系数数组改成可复用缓存。
- [ ] 优化 `HPwithFin.py:399-405` 的双重求解。`pre_step()` 为了估算等效导热率，每个时间步对每根热管做两次 `_solve_fin_quasi_steady()`，代价几乎翻倍。建议优先评估解析线性化、一次求解后导出切线，或引入上一步导热率预测。
- [ ] 优化 `RingHP.py:49-68` 的单节点代理接口。`FluidSolidCouple.execute()` 高频读取 `temperature_vector / pressure_vector / density_vector / velocity_vector`，而 `SingleVolumeProxy` 每次都创建新的长度为 1 的数组。建议缓存 1 元素 buffer 并原位刷新，或扩展 `FluidSolidCouple` 对标量代理的支持。

### P3

- [ ] 评估 `HeatPipe2D.py:239-259` 和 `HeatPipe2D.py:315-327` 的切片索引缓存。当前每次边界更新/热流拼接都重复构造 `idx_eva/idx_aba/idx_con` 并切片，可以作为低风险的轻量优化收尾项。

## 4. 建议实施顺序

1. 先修复 `RingHP.py:215` 的闭包错误和 `HPwithFin.py:297` 的材料硬编码。
2. 再统一 `RingHP` 的总量统计口径与阻塞/局阻模型。
3. 最后进入 `HeatPipe2D` / `HPwithFin` 的 buffer 化和翅片求解降本。
