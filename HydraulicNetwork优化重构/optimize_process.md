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
