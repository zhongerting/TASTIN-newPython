# CoolantLoop - Codex 首读手册

> 基于当前工作区代码生成。每次让 Codex 处理 `CoolantLoop` 目录时，应先阅读本文件。
>
> 本目录是冷却回路用例层：包含当前集流环模型库、运行包装器、诊断测试、历史测试、独立分析脚本和既有运行产物。修改前仍应按任务定点回看源码。

## 1. Codex 首读指引

### 1.1 可直接依据本手册处理的任务

- 判断当前推荐入口，选择 `model_*.py`、`run_*.py` 或诊断脚本。
- 定位集流环、进口缓冲区、热腿、出口支路、manifold、断点文件和 profiler 报告。
- 解释两个生产模型的拓扑差别、默认工况、离散规模和续算链路。
- 为调试选择最小验证脚本、历史回归脚本或完整集流环测试。

### 1.2 必须回看源码的情况

遇到以下任务，不要只依赖本手册：

- 修改物理模型、几何参数、热管倍数、辐射边界或材料。
- 修改流体拓扑、连接顺序、`MacroFlowJunction` 倍数或节点映射。
- 修改自适应步长、Picard 内迭代、日志周期、CSV 字段或 profiler 行为。
- 修改 `build_model()` 返回字典、`run_case()` 参数或包装器续算链路。
- 修改上级 `Components`、`Solvers`、`Materials` 接口或状态保存字段。

## 2. 目录定位与文件地图

### 2.1 当前推荐入口

| 文件 | 类型 | 用途 | 当前入口建议 |
|---|---|---|---|
| [`model_collector_ring_6segment.py`](./model_collector_ring_6segment.py) | 生产模型库 | 构造 6 个拼接的 `1/6 RingHP` 扇区；支持运行、断点和 profiler | 优先用于分扇区定位、性能分析 |
| [`model_collector_ring_full_ringhp.py`](./model_collector_ring_full_ringhp.py) | 生产模型库 | 构造单一 360 度 `RingHP`；支持运行和断点 | 优先用于完整环模型 |
| [`run_collector_ring_6segment_160s.py`](./run_collector_ring_6segment_160s.py) | 运行包装器 | 从初始态运行 6segment 到 `160 s`，每 `40 s` 保存断点并输出 profiler | 可直接运行 |
| [`run_collector_ring_full_ringhp_200s_resume.py`](./run_collector_ring_full_ringhp_200s_resume.py) | 续算包装器 | 从 `60 s` 快照续算 full ringhp 到 `200 s` | 依赖指定输入快照 |
| [`run_collector_ring_full_ringhp_500s_resume.py`](./run_collector_ring_full_ringhp_500s_resume.py) | 续算包装器 | 从 `200 s` 快照续算 full ringhp 到 `500 s` | 可接续前一个包装器 |
| [`run_collector_ring_full_ringhp_1p3kgs_100s.py`](./run_collector_ring_full_ringhp_1p3kgs_100s.py) | 流量变体包装器 | 从初始态运行 full ringhp，总入口流量临时设为 `1.3 kg/s`，目标 `100 s` | 1.3 kg/s 短期测试入口；早期瞬态需小步长 |
| [`run_collector_ring_full_ringhp_1p3kgs_inlet_init_steady.py`](./run_collector_ring_full_ringhp_1p3kgs_inlet_init_steady.py) | 稳态搜索包装器 | 从 `t=0` 构造 full ringhp，总入口流量 `1.3 kg/s`，流体和集流环固体初温默认设为入口温度；热管初温默认同入口温度，也可用 `--hp-init-temp 800` 复现早期完整环口径；默认前 `0.1 s` 小步长启动 | 1.3 kg/s 真初态到稳态首选入口 |
| [`run_collector_ring_full_ringhp_geometry100hp_test.py`](./run_collector_ring_full_ringhp_geometry100hp_test.py) | 几何变体测试包装器 | 矩形集流环 `110 mm x 40 mm`、总热管数 `100`、热管 `0.1/0/0.5 m`、蒸汽半径 `7.5 mm`、吸液芯厚 `0.5 mm`、壁厚 `1 mm`、翅片 `20 mm x 0.4 mm` | 新结构错位检查和短计算测试 |
| [`run_collector_ring_full_ringhp_geometry100hp_potassium_test.py`](./run_collector_ring_full_ringhp_geometry100hp_potassium_test.py) | 钾热管几何变体包装器 | 复用 `geometry100hp` 的矩形集流环、100 根热管和几何/流量/发射率参数；临时重建 `PotassiumHP` 与钾热管 `WickMaterial` | 钾热管独立 A/B 对比入口；不要加载钠热管 restart |
| [`model_collector_ring_full_ringhp_geometry100hp_potassium_mixed.py`](./model_collector_ring_full_ringhp_geometry100hp_potassium_mixed.py) | 钾热管混合节点独立算例 | 固化 `geometry100hp`、100 根钾热管、总流量 `1.3 kg/s`、`eps_hp=eps_fin=0.6`、集流环壁面 `eps=0.05`、`T_space=200 K`；入口/出口均有小体积无热管混合节点 | 后续闭式回路和入口/出口混合节点策略的首选基线 |
| [`model_collector_ring_6segment_geometry100hp_potassium_mixed.py`](./model_collector_ring_6segment_geometry100hp_potassium_mixed.py) | 6 段钾热管混合节点独立算例 | 将 `geometry100hp` 拆为 6 个非闭合 `1/6 RingHP` 段，段间插入 3 个入口混合体和 3 个出口混合体；100 根钾热管按 `17/16/17/16/17/17` 分配，单段 3 个控制体 | 当前入口/出口连接混合体方案的首选实现；避免把混合体建成带热管的环节点 |
| [`run_collector_ring_6segment_geometry100hp_potassium_anisotropic_wick_test.py`](./run_collector_ring_6segment_geometry100hp_potassium_anisotropic_wick_test.py) | 吸液芯各向异性导热 A/B 测试包装器 | 从 6 段钾热管混合节点算例构建系统，运行前对所有代表热管调用 `set_wick_conductivity_mode(True)`；默认从 `500 s` restart 跑到 `501 s` 并统计 BDF warning | 用于验证“轴向赝热导、径向结构导热”是否降低 BDF warning 和改变稳态温差；不是生产默认模型 |
| [`run_collector_ring_full_ringhp_steady_debug.py`](./run_collector_ring_full_ringhp_steady_debug.py) | 稳态调试包装器 | 默认从最新 full ringhp 快照短续算，输出 history、restart 和 `T_out_avg` 斜率判据 | 集流环+热管稳态调试首选 |
| [`run_collector_ring_full_ringhp_energy_audit.py`](./run_collector_ring_full_ringhp_energy_audit.py) | 能量守恒审计包装器 | 从 full ringhp restart 续算并输出 `Q_loop - Q_rejection - dU/dt` 残差 | 集流环+热管整体能量守恒检查首选 |

`model_*.py` 可直接运行其内置 `50 s` 默认案例，也可被包装器 import。需要自定义参数时，优先调用 `run_case(...)`；需要检查或改造对象时，再单独调用 `build_model()`。

### 2.2 诊断测试与历史测试

| 文件 | 分类 | 用途 | 备注 |
|---|---|---|---|
| [`test_single_header_cell_one_hp.py`](./test_single_header_cell_one_hp.py) | 最小单元验证 | 单个 header cell 与单个带翅片热管，输出能量审计 CSV 和 PNG | 基础组件问题首选 |
| [`test_full_collector_ring.py`](./test_full_collector_ring.py) | 完整集流环验证 | 较早的完整集流环用例，检查流量分配、温度和散热 | 诊断参考，不是当前生产入口 |
| [`test_coolant_loop_v5.py`](./test_coolant_loop_v5.py) | 历史回归与诊断 | 双环支路、壁面辐射、热管散热分解、CSV 和断点 | 诊断价值高，不代表当前默认工况 |
| [`test_coolant_loop_v5_wall_radiation_only.py`](./test_coolant_loop_v5_wall_radiation_only.py) | 壁面辐射诊断 | 隔离壁面辐射路径，统计 `dU/dt` 和能量残差 | 排查辐射散热首选 |
| [`test_coolant_loop_v4.py`](./test_coolant_loop_v4.py) | 历史测试 | v4 演进版本 | 不要将参数当作生产默认值 |
| [`test_coolant_loop_v4_1.py`](./test_coolant_loop_v4_1.py) | 历史测试 | v4.1 演进版本 | 不要将 monkey-patch 当作当前接口 |
| [`test_coolant_loop_v4_2.py`](./test_coolant_loop_v4_2.py) | 历史回归与能量审计 | v4.2 分项储能、`dU/dt`、辐射分解和 CSV/PNG | 长尾能量问题的重要参考 |
| [`verify_dudt_long_tail_case.py`](./verify_dudt_long_tail_case.py) | 独立分析 | 用解析模型隔离验证小幅长期 `dU/dt` 尾部 | 不构造 TASTIN 回路 |

### 2.3 参考资料与生成产物

| 文件或模式 | 分类 | 用途 |
|---|---|---|
| [`bug_report.md`](./bug_report.md) | 历史问题记录 | 早期源码问题报告；部分判断已过时 |
| [`bug_report (1).md`](./bug_report%20%281%29.md) | 历史问题记录 | 更新版问题报告，包含对早期误诊的撤回说明 |
| [`collector_ring_6segment_buffered_160s_profiler_report.txt`](./collector_ring_6segment_buffered_160s_profiler_report.txt) | profiler 报告 | 已有性能分析结果，仅作参考 |
| `*_history.csv`、`*_diagnostics.csv`、`smoke_*.csv` | 运行产物 | 历史数据和冒烟测试输出，不是代码入口 |
| `*_restart.npz`、`*_restart_t????s.npz` | 断点产物 | `SystemManager` 全局状态快照 |
| `*.log` | 运行产物 | 包装器历史日志 |
| `*.png` | 报告图片 | 独立分析或诊断脚本生成图 |

不要把 `bug_report*.md`、profiler 报告、CSV、PNG、日志或 `.npz` 当作当前运行入口。它们用于回溯现象、对比结果和恢复计算。

## 3. 当前主干模型

### 3.1 共享默认工况

两个 `model_*.py` 使用同一组生产默认值：

| 参数 | 值 | 说明 |
|---|---:|---|
| `T_SPACE` | `3.0 K` | 辐射环境温度 |
| `T_INLET` | `843.0 K` | 进口边界温度 |
| `T_INIT` | `863.0 K` | 初始流体和固体温度 |
| `P_OUTLET` | `160000.0 Pa` | 出口定压边界 |
| `W_TOTAL` | `2.2 kg/s` | 进口总质量流量 |
| `W_BRANCH_TOTAL` | `W_TOTAL / 3` | 三条热腿的初始分支流量 |
| `DEFAULT_T_END` | `50.0 s` | 直接运行模型库时的默认终止时间 |

默认生产模型材料路径为：冷却剂 `SodiumPotassium78`，壁材 `SS316`，热管工质 `SodiumHP`，吸液芯等效材料 `WickMaterial`。钾热管独立包装器会在运行时临时替换为 `PotassiumHP` 和对应钾热管 `WickMaterial`，不修改默认模型库。

当前换热模型口径：

- 冷却剂与集流环内壁仍使用 `lyon_martinelli`，即 `Nu = 7 + 0.025 Pe^0.8`，适用于通道内液态金属对壁换热。
- 冷却剂与热管蒸发段外壁使用液态金属横掠单圆柱模型：`Nu_D = K_cal F_wake C sqrt(Pe_D)`，默认 `C=0.65`、`K_cal=1.0`、`F_wake=1.0`。
- 热管蒸发段特征长度取热管外径 `D_hp`，换热面积取 `pi D_hp L_eva`；局部速度用热管投影面积修正后的最小通流截面速度。
- 一个流体控制体内若代表多根热管，当前仍使用 `Q_total = N_hp * Q_single` 的集总源项放大，未解析控制体内部逐根热管串联升温。

### 3.2 共享几何与离散

| 区域 | 关键参数 | 离散 |
|---|---|---:|
| 进口缓冲区 | `L_INLET_BUFFER=0.20 m`，面积为三条热腿总面积 | `N_INLET_BUFFER=5` |
| 三条热腿 | 每条 `L_HOT_LEG=2.19632 m`，`R_IN_HOT_LEG=0.0138 m` | 每条 `N_HOT_LEG=28` |
| 每个 `1/6` 环段 | `L_SECTOR=0.793 m`，`AREA_RING=0.0016065 m^2` | `N_SECTOR=4` |
| 三条出口支路或 manifold | 每条 `L=0.40911 m`，内半径 `0.009 m` | 每条 `N=5` |
| 出口缓冲区 | `L_OUTLET_BUFFER=0.20 m`，面积为三条出口通道总面积 | `N_OUTLET_BUFFER=5` |
| 热管倍数 | 每个 `1/6` 环段为 `[6, 7, 6, 7]` | 每个环段共 `26` 根，整环共 `156` 根 |

每个环段壁面固体径向仅 `n_x=1`。6segment 为 6 个 `n_y=4` 的固体和 `RingHP`；full ringhp 为一个 `n_y=24` 的固体和 `RingHP`。热管冷凝段和翅片内部还有各自离散，修改时应进入 `Components` 源码核对。

### 3.3 `model_collector_ring_6segment.py`

该模型不是 6 条互不相连的并联扇区，而是 6 个首尾拼接的 `1/6 RingHP` 扇区。扇区端点形成 3 个进口节点和 3 个出口节点：

```text
InletBoundary
  -> InletBuffer
  -> 3 x HotLeg
  -> I1 / I2 / I3
  -> 6 个拼接环段: S1(I1->O1), S2(O1->I2), ... , S6(O3->I1)
  -> O1 / O2 / O3
  -> 3 x OutletBranch
  -> OutletBuffer
  -> OutletBoundary
```

跨热腿与环、环与出口支路的连接使用 `MacroFlowJunction(multiplier=2)`。诊断流量时要区分 `junc.W` 与 `junc.get_mass_flow_for(...)` 的口径，不要把代表侧和物理放大侧直接混用。当前对称建模口径是：程序只建一个物理集流环，动量方程和热管管束压降按单环流量计算；宏观端质量和能量按 2 倍代表两个对称集流环。

`build_model()` 返回字典的关键对象：

```text
inlet_boundary, outlet_boundary
inlet_buffer_channel, hot_legs, outlet_branches, outlet_buffer_channel
ring_nodes, sectors, solids, ring_hps, sector_specs
inlet_junction, inlet_buffer_to_hot_leg, hot_leg_to_ring
ring_to_outlet_branch, outlet_branch_to_outlet_buffer, outlet_junction
sector_link_junctions
all_vols, all_juncs, network, sys_mgr
```

定位建议：

- 修改单个 `1/6` 环段：看 `sector_specs`、`sectors`、`solids`、`ring_hps`。
- 修改环上接口：看 `ring_nodes` 和 `sector_link_junctions`。
- 修改支路流量口径：看 `hot_leg_to_ring` 与 `ring_to_outlet_branch`。
- 分析耗时：看该文件内 `PROFILER_KEY_FUNCTIONS` 和 profiler 输出。

### 3.4 `model_collector_ring_full_ringhp.py`

该模型把整环表示为一个 360 度流体通道、一个固体和一个 `RingHP`。环通道共 `N_RING=6*N_SECTOR=24` 个节点，并由 `ring_closure` 将末节点闭合回节点 `0`：

```text
InletBoundary
  -> InletBuffer
  -> 3 x HotLeg
  -> Ring nodes [0, 8, 16]
  -> 单一 360 度 RingHP 通道，末节点通过 ring_closure 回到节点 0
  -> Ring nodes [4, 12, 20]
  -> 3 x Manifold
  -> OutletBuffer
  -> OutletBoundary
```

`HP_MULTIPLIERS_RING = [6, 7, 6, 7] * 6`。跨热腿与环、环与 manifold 的连接同样使用 `MacroFlowJunction(multiplier=2)`。`ring_closure` 使用 `RingHP.outlet_dynamic_loss_params` 承接最后一个环节点的热管管束压降；不要只读取 `outlet_k_loss` 判断热管阻力。

`build_model()` 返回字典的关键对象：

```text
inlet_boundary, outlet_boundary
inlet_buffer_channel, hot_legs, manifolds, outlet_buffer_channel
ring_channel, ring_solid, ring_hp, ring_closure
inlet_junction, inlet_buffer_to_hot_leg, hot_leg_to_ring
ring_to_manifold, manifold_to_outlet_buffer, outlet_junction
all_vols, all_juncs, network, sys_mgr
```

### 3.5 两个模型如何选择

| 关注点 | `6segment` | `full_ringhp` |
|---|---|---|
| 环表示 | 6 个拼接的 `1/6 RingHP` | 1 个 360 度 `RingHP` |
| 环流体节点 | `6 x 4` | `24` |
| 环固体与宏观组件 | 6 组 | 1 组 |
| 出口通道命名 | `outlet_branches` | `manifolds` |
| 环闭合 | 6 个 `sector_link_junctions` | 1 个 `ring_closure` |
| profiler | 内置输出 | 无内置 profiler 输出 |
| 断点续算 | 支持 | 支持 |

## 4. 运行与断点续算

### 4.1 `run_case()` 接口

两个模型的 `run_case()` 都会自行调用 `build_model()`，最后返回 `(model, history)`。共同参数为：

```python
run_case(
    case_name=...,
    t_end=50.0,
    min_dt=1.0e-3,
    max_dt=0.5,
    safety_factor=1.0,
    inner_iter=2,
    print_every_time=1.0,
    csv_path=None,
    restart_from=None,
    restart_save_path=None,
    restart_save_every=10.0,
)
```

`model_collector_ring_6segment.py` 额外支持：

```python
profiler_summary_path=None
profiler_snapshot_path=None
profiler_report_path=None
```

运行循环通过 `sys_mgr.compute_adaptive_dt(min_dt=1e-3, max_dt=0.5, ...)` 计算步长，并在日志时刻、断点时刻和 `t_end` 前主动截短步长。`history` 与 CSV 按每个仿真步记录；`print_every_time` 控制终端日志周期，默认 `1 s`。

6segment 的 CSV 额外包含 CPU 和 wall time，并始终输出 profiler summary CSV、snapshot CSV 和文本报告。full ringhp 只输出历史 CSV。

### 4.2 `.npz` 断点 API

断点由 [`../Solvers/SystemManager.py`](../Solvers/SystemManager.py) 统一管理：

```python
sys_mgr.save_global_state("restart_state.npz")
sys_mgr.load_global_state("restart_state.npz")
```

不要使用旧手册中的 `save_snapshot()`、`load_snapshot()` 或 `.pkl` 模板。`save_global_state()` 会汇总系统时间、上一步长、流体网络、固体、宏观组件和可选点堆状态；`load_global_state()` 加载后还会同步固体时间、刷新边界缓存并重建耦合源项。

当 `restart_save_path` 非空时，运行结束会保存最终文件；当 `restart_save_every > 0` 时，还会生成带时间后缀的检查点，例如 `_restart_t0200s.npz`。恢复后，两个模型都会重新施加进口温度、进口压力和出口定压边界，再同步到水力网络。

### 4.3 三个包装器与续算链路

| 包装器 | 命令 | 输入 | 主要输出 |
|---|---|---|---|
| [`run_collector_ring_6segment_160s.py`](./run_collector_ring_6segment_160s.py) | `python run_collector_ring_6segment_160s.py` | 初始态 | `collector_ring_6segment_buffered_160s_history.csv`、每 `40 s` 检查点、最终 restart、三类 profiler 文件 |
| [`run_collector_ring_full_ringhp_200s_resume.py`](./run_collector_ring_full_ringhp_200s_resume.py) | `python run_collector_ring_full_ringhp_200s_resume.py` | `collector_ring_full_ringhp_buffered_200s_resume_from50s_restart_t0060s.npz` | 到 `200 s` 的历史 CSV、每 `20 s` 检查点、最终 restart |
| [`run_collector_ring_full_ringhp_500s_resume.py`](./run_collector_ring_full_ringhp_500s_resume.py) | `python run_collector_ring_full_ringhp_500s_resume.py` | 前一步生成的 `_restart_t0200s.npz` | 到 `500 s` 的历史 CSV、每 `100 s` 检查点、最终 restart |
| [`run_collector_ring_full_ringhp_1p3kgs_100s.py`](./run_collector_ring_full_ringhp_1p3kgs_100s.py) | `python run_collector_ring_full_ringhp_1p3kgs_100s.py` | 初始态；仅在 wrapper 内临时覆盖 `W_TOTAL=1.3 kg/s`；默认 `max_dt=5e-3 s`、`safety_factor=0.1` | `1p3kgs_100s` history、每 `20 s` 检查点、最终 restart；已验证可到 `0.05 s`，完整 `100 s` 属长时运行 |
| [`run_collector_ring_full_ringhp_1p3kgs_inlet_init_steady.py`](./run_collector_ring_full_ringhp_1p3kgs_inlet_init_steady.py) | `python run_collector_ring_full_ringhp_1p3kgs_inlet_init_steady.py` | 初始态；临时覆盖 `W_TOTAL=1.3 kg/s`、`T_INIT=T_INLET`，默认 `HP_INITIAL_TEMP=T_INLET`；若要参考 `test_full_collector_ring.py`，加 `--hp-init-temp 800`；默认 `startup_time=0.1 s`、`startup_max_dt=5e-3 s`、`startup_safety_factor=0.1`；支持 `--restart-from` 分段续算 | 到稳态或 `100 s` 上限；history、每 `20 s` 检查点、最终 restart；全入口温度已验证到 `0.1 s`；`--hp-init-temp 800` 已验证到 `0.2 s`，该口径启动瞬态需更长小步长 |
| [`run_collector_ring_full_ringhp_geometry100hp_test.py`](./run_collector_ring_full_ringhp_geometry100hp_test.py) | `python run_collector_ring_full_ringhp_geometry100hp_test.py` | 初始态；临时覆盖矩形通道、100 根热管和新热管/翅片几何；默认 `W_TOTAL=1.3 kg/s`、`T_INIT=T_INLET`、`HP_INITIAL_TEMP=800 K`、`t_end=0.05 s` | 新几何 history 和 restart；已验证可到 `0.05 s`，运行中可能出现不致停的 BDF 警告 |
| [`run_collector_ring_full_ringhp_geometry100hp_potassium_test.py`](./run_collector_ring_full_ringhp_geometry100hp_potassium_test.py) | `python run_collector_ring_full_ringhp_geometry100hp_potassium_test.py` | 初始态；复用 `geometry100hp` 覆盖参数，但热管工质临时替换为 `PotassiumHP(name="HP_Fluid_K")`，并重建钾热管 `WickMaterial_K` | 钾热管 history 和 restart；用于与钠热管几何变体 A/B 对比，restart 不能与钠热管算例混用 |
| [`run_collector_ring_full_ringhp_steady_debug.py`](./run_collector_ring_full_ringhp_steady_debug.py) | `python run_collector_ring_full_ringhp_steady_debug.py --duration 1` | 默认选择最新 full ringhp restart，优先 `500 s` 快照 | 短续算 history、最终 restart 和末段 `T_out_avg` 斜率判稳摘要 |
| [`run_collector_ring_full_ringhp_energy_audit.py`](./run_collector_ring_full_ringhp_energy_audit.py) | `python run_collector_ring_full_ringhp_energy_audit.py --duration 1` | 默认选择最新 full ringhp 稳态调试快照，回退到 `500 s` 快照 | 能量审计 CSV，最终 `Q_loop_flux`、辐射排热、储能率和残差摘要 |

当前目录中可见 `200 s` 和 `500 s` 链路产出的检查点与日志；`200 s` 包装器声明的 `60 s` 输入快照当前不在目录清单中。重新执行该包装器前，应先确认该输入文件存在。

## 5. 测试与调试导航

### 5.1 按目标选脚本

| 调试目标 | 优先阅读或运行 |
|---|---|
| 最小热管与 header cell 验证 | [`test_single_header_cell_one_hp.py`](./test_single_header_cell_one_hp.py) |
| 历史行为回归 | [`test_coolant_loop_v4.py`](./test_coolant_loop_v4.py)、[`test_coolant_loop_v4_1.py`](./test_coolant_loop_v4_1.py)、[`test_coolant_loop_v4_2.py`](./test_coolant_loop_v4_2.py)、[`test_coolant_loop_v5.py`](./test_coolant_loop_v5.py) |
| 壁面辐射隔离诊断 | [`test_coolant_loop_v5_wall_radiation_only.py`](./test_coolant_loop_v5_wall_radiation_only.py) |
| 完整集流环流量与温度分配 | [`test_full_collector_ring.py`](./test_full_collector_ring.py)，再回看两个生产模型 |
| 长时间小幅 `dU/dt` 尾部 | [`verify_dudt_long_tail_case.py`](./verify_dudt_long_tail_case.py)、[`test_coolant_loop_v4_2.py`](./test_coolant_loop_v4_2.py) |
| 性能瓶颈 | [`model_collector_ring_6segment.py`](./model_collector_ring_6segment.py)、[`../profiler.py`](../profiler.py)、既有 profiler 报告 |

### 5.2 常见问题定位

| 现象 | 首读文件 | 重点 |
|---|---|---|
| 能量守恒或 `dU/dt` 异常 | `run_collector_ring_full_ringhp_energy_audit.py`、`test_single_header_cell_one_hp.py`、`test_coolant_loop_v4_2.py`、`verify_dudt_long_tail_case.py` | full ringhp 先用审计包装器；底层问题再区分流体、header、热管和翅片储能；区分瞬态尾部与真实残差 |
| full ringhp 稳态判断 | `run_collector_ring_full_ringhp_steady_debug.py`、`model_collector_ring_full_ringhp.py` | 从最新 restart 短续算，先看 `T_out_avg` 末段斜率和流量闭合，再决定是否做长续算 |
| 辐射散热异常 | `test_coolant_loop_v5_wall_radiation_only.py`、`test_coolant_loop_v5.py` | 壁面辐射、热管裸壁辐射、翅片辐射和面积口径 |
| 流量分配异常 | 两个 `model_*.py`、`test_full_collector_ring.py` | `MacroFlowJunction` 倍数、环接口节点、代表侧和放大侧流量 |
| 断点恢复异常 | 两个 `model_*.py`、三个 `run_*.py`、`../Solvers/SystemManager.py` | `.npz` 路径、保存字段、加载后的边界重施加和耦合刷新 |
| 运行过慢 | `model_collector_ring_6segment.py`、`../profiler.py`、profiler 报告 | 自适应步长、水力、导热和耦合器耗时 |

`test_coolant_loop_v4*` 与 `test_coolant_loop_v5.py` 属于历史演进测试。它们保留了诊断价值，但其中工况、猴子补丁、输出字段和局部参数不得覆盖当前两个生产模型的默认值。

## 6. 上级依赖导航

本目录直接装配上级模块，不在这里复制其内部实现。修改组件内部行为时，先阅读对应上级手册，再进入源码。

| 上级模块 | 本目录直接使用的职责 | 首读文档 |
|---|---|---|
| `Components` | `RingHP` 组织代表性 `HPwithFin`，按 `hp_multipliers` 放大真实热管数量，并向 `SystemManager` 暴露内部固体和耦合器 | [`../Components/COMPONENTS_DETAILED_INTRO.md`](../Components/COMPONENTS_DETAILED_INTRO.md)、[`../Components/BASICCOMPONENTS_DETAILED_INTRO.md`](../Components/BASICCOMPONENTS_DETAILED_INTRO.md) |
| `Solvers/Hydrodynamics` | `IncompressibleFluidChannel`、边界控制体、`FlowJunction`、`MacroFlowJunction` 和 `HydraulicNetwork` 组成流体网络 | [`../Solvers/AI_AGENT_SOLVERS_ANALYSIS.md`](../Solvers/AI_AGENT_SOLVERS_ANALYSIS.md)、[`../Solvers/HYDRODYNAMICS_DETAILED_INTRO.md`](../Solvers/HYDRODYNAMICS_DETAILED_INTRO.md) |
| `Solvers/HeatConduction` | `Mesh2D`、`HeatConduction2D` 和动态辐射边界求解集流环壁面导热与散热 | [`../Solvers/HEATCONDUCTION_DETAILED_INTRO.md`](../Solvers/HEATCONDUCTION_DETAILED_INTRO.md) |
| `Solvers/SystemManager.py` | 初始化、全局步进、自适应步长、耦合刷新、回滚和 `.npz` 断点续算 | [`../Solvers/SYSTEMMANAGER_DETAILED_INTRO.md`](../Solvers/SYSTEMMANAGER_DETAILED_INTRO.md) |
| `Materials` | `SodiumPotassium78`、`SS316`、`SodiumHP`/`PotassiumHP`、`WickMaterial` 的温度相关物性 | [`../Materials/AI_AGENT_MATERIALS_ANALYSIS.md`](../Materials/AI_AGENT_MATERIALS_ANALYSIS.md) |

## 7. 维护约定

以下内容变化时，必须同步更新本手册：

1. 两个生产模型的默认常量、几何离散或热管倍数。
2. 环拓扑、节点映射、`MacroFlowJunction` 倍数或支路命名。
3. `build_model()` 返回字典。
4. `run_case()` 参数、默认步长范围、日志、CSV、profiler 或断点行为。
5. 三个 `run_*.py` 包装器的输入快照、终止时间、输出文件或续算链路。
6. 测试入口、独立分析脚本或上级首读文档。

更新后至少逐项核对两个 `model_*.py` 的常量、`build_model()` 返回值和 `run_case()` 签名，并核对三个包装器中的文件名与恢复链路。

## 8. 2026-06-04 geometry100hp emissivity scan note

- `Components/HPwithFin.py` now supports separate bare heat-pipe emissivity (`emissivity`) and fin emissivity (`fin_emissivity`). If `fin_emissivity` is omitted, it falls back to `emissivity`.
- `Components/RingHP.py` passes the optional `fin_emissivity` through to each representative `HPwithFin`.
- `CoolantLoop/model_collector_ring_full_ringhp.py` exposes `HP_EMISSIVITY` and `FIN_EMISSIVITY`; `run_collector_ring_full_ringhp_geometry100hp_test.py` exposes them as `--hp-emissivity` and `--fin-emissivity`.
- `CoolantLoop/scan_geometry100hp_emissivity.py` is an instantaneous diagnostic. It loads a restart, recomputes current fin/bare radiation for multiple emissivity pairs, and reports equivalent temperature drop. It does not advance time and must not be treated as a final steady-state result.
- Screening from `collector_ring_full_ringhp_geometry100hp_resume_100s_restart_t0040s.npz` showed that fin emissivity is slightly more sensitive than bare heat-pipe emissivity, but the `40-100 s` outlet-temperature evolution is still dominated by stored thermal energy. Dynamic continuations found:
  - `eps_hp=eps_fin=0.3`: `T_out_avg=759.838 K` at `70 s`, `DeltaT=83.162 K` (below the `85-95 K` target band).
  - `eps_hp=eps_fin=0.6`: `T_out_avg=752.123 K` at `70 s`, `DeltaT=90.877 K`; `T_out_avg=748.730 K` at `100 s`, `DeltaT=94.270 K` (inside but near the upper edge).
  - `eps_hp=eps_fin=0.85`: `T_out_avg=755.266 K` at `55 s`, `DeltaT=87.734 K`.
  - `eps_hp=eps_fin=1.0`: `T_out_avg=754.309 K` at `55 s`, `DeltaT=88.691 K`.
- For the current 100-heat-pipe geometry and `W_TOTAL=1.3 kg/s`, the best tested short-term candidate is `HP_EMISSIVITY=0.6`, `FIN_EMISSIVITY=0.6`. A formal from-zero long run should still test the `0.55-0.60` range because the `100 s` result is close to the `95 K` upper bound and was continued from the existing `40 s` transient restart.
- Follow-up 200 s continuations were run in parallel with equal bare heat-pipe and fin emissivity (`eps_hp=eps_fin=eps`) from the existing transient restarts. At `200 s`, the outlet temperature and inlet-outlet drop were:
  - `eps=0.40`: `T_out_avg=762.592 K`, `DeltaT=80.408 K` (below target).
  - `eps=0.45`: `T_out_avg=757.854 K`, `DeltaT=85.146 K` (lower edge).
  - `eps=0.50`: `T_out_avg=753.894 K`, `DeltaT=89.106 K` (center of target band).
  - `eps=0.55`: `T_out_avg=750.740 K`, `DeltaT=92.260 K` (inside target band).
  - `eps=0.60`: `T_out_avg=748.199 K`, `DeltaT=94.801 K` (upper edge).
  The `160-200 s` outlet-temperature changes were small for `eps=0.50-0.60`; `eps=0.50` is the most balanced tested value for the requested `85-95 K` drop. Radiation diagnostics at `200 s` gave `Q_bare/Q_fin` pairs of approximately `17.50/23.33 kW` for `eps=0.50`, `18.47/24.52 kW` for `eps=0.55`, and `19.31/25.54 kW` for `eps=0.60`; fin radiation remains about `1.32x` bare heat-pipe radiation when both emissivities are adjusted together.

## 9. 2026-06-06 tube-bank pressure-drop note

- `RingHP` no longer uses the old fixed heat-pipe obstruction estimate `C_D * A_proj/A_flow * N_eff` for nodes with heat pipes. It configures `dynamic_loss_params["model"] == "inline_tube_bank_euler"` on the corresponding ring junctions.
- The model uses the inline tube-bank Euler-number relation `Eu = 0.67 * (S_T/d - 1)^(-0.5) * Re^(-0.15)` and converts it to the solver's `K_eq` with `K_eq = Eu * N_hp * (A_flow/A_min)^2`.
- For `geometry100hp`, use `A_flow = 0.110 * 0.040 = 0.0044 m2`, `d = 0.018 m`, `L_eva = 0.100 m`, so `S_T = A_flow/L_eva = 0.044 m` and `S_T/d = 2.44`. This value is tied to the same minimum-area definition used for `v_max`; do not substitute `110/18` or `40/18`.
- `HP_MULTIPLIERS_RING` is interpreted as the node's streamwise tube-row count for pressure drop and as the node heat-pipe multiplier for heat transfer. With the current `geometry100hp` distribution, `N_hp` is `4` or `5` by node.

## 10. 2026-06-06 geometry100hp potassium mixing-node case

- `model_collector_ring_full_ringhp_geometry100hp_potassium_mixed.py` is a standalone case for the 100-potassium-heat-pipe geometry with inlet/outlet mixing nodes. It should be preferred over temporary runtime overrides when reproducing the `eps_hp=eps_fin=0.6`, ring-wall `eps=0.05`, `T_space=200 K` condition.
- Topology:
  ```text
  inlet_boundary -> inlet_buffer -> MacroFlowJunction x3
    -> single-ring hot_leg x3 -> inlet_mix_node x3 -> ring nodes [0, 8, 16]
  ring nodes [4, 12, 20] -> outlet_mix_node x3 -> single-ring manifold x3
    -> MacroFlowJunction x3 -> outlet_buffer -> outlet_boundary
  ```
- The `MacroFlowJunction(multiplier=2)` is now placed between the macro header/buffer and the single-ring inlet/outlet pipe. Thus the inlet/outlet buffers see total physical flow, while each hot leg, mixing node, manifold, and ring channel sees the single modeled ring flow.
- The mixing nodes use `MIX_NODE_VOLUME_FACTOR = 0.10` times one ring-node volume, have no heat pipes, and initially add no local loss. They are intended to separate inlet/outlet injection and extraction from the heat-pipe-bearing ring control volumes.

## 11. 2026-06-06 6segment geometry100hp potassium mixed-node case

- `model_collector_ring_6segment_geometry100hp_potassium_mixed.py` is the segmented implementation of the revised inlet/outlet connection strategy. It reuses the `geometry100hp` potassium heat-pipe geometry, `W_TOTAL=1.3 kg/s`, `eps_hp=eps_fin=0.6`, ring-wall `eps=0.05`, and `T_space=200 K`.
- The hydraulic ring is built from six non-closed `1/6 RingHP` segments instead of one closed ring. The segment connection order is:
  ```text
  I1 -> A1 -> O1 -> A2 -> I2 -> A3 -> O2
     -> A4 -> I3 -> A5 -> O3 -> A6 -> I1
  ```
  where `I1/I2/I3` are inlet mixing volumes and `O1/O2/O3` are outlet mixing volumes.
- Each segment has three coolant control volumes. Heat-pipe multipliers are `A1=[5,6,6]`, `A2=[5,5,6]`, `A3=[5,6,6]`, `A4=[5,5,6]`, `A5=[5,6,6]`, `A6=[5,6,6]`, giving segment totals `17/16/17/16/17/17` and total heat-pipe count `100`.
- The mixing volumes are not heat-pipe-bearing ring nodes. Inlet mixing volumes use `AREA_HOT_LEG * (2*DH_HOT_LEG)` and outlet mixing volumes use `AREA_MANIFOLD * (2*DH_MANIFOLD)`. Links from mixing volumes into/out of RingHP segments use `AREA_RING`; no additional local loss is added initially.
- `MacroFlowJunction(multiplier=2)` is only used between macro buffers and single-ring branches. A smoke run to `0.005 s` verified `W_in_total=1.300000 kg/s`, `W_ring_in_total=0.650000 kg/s`, and near-even inlet branch flows of about `0.2167 kg/s` each.

## 12. 2026-06-07 anisotropic wick conductivity test wrapper

- `run_collector_ring_6segment_geometry100hp_potassium_anisotropic_wick_test.py` is a diagnostic wrapper for the 6-segment potassium case. It monkey-patches only the constructed model instance: after `build_model()`, each representative heat-pipe solid receives `hp.set_wick_conductivity_mode(True)`.
- The anisotropic wick mode keeps axial conductivity as `WickMaterial.conductivity_axial()` including pseudothermal conductivity, while radial conductivity uses `WickMaterial.conductivity_radial()`, currently structural conductivity only. The base production case remains isotropic by default.
- Default command:
  ```powershell
  & "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" CoolantLoop\run_collector_ring_6segment_geometry100hp_potassium_anisotropic_wick_test.py
  ```
  This loads `collector_ring_6segment_geometry100hp_potassium_mixed_500s_resume_restart.npz`, runs to `501 s`, captures BDF RuntimeWarnings, and prints final `T_out_avg`, inlet-outlet temperature drop, flow closure, and top wick-conductivity diagnostics.
- For A/B control, pass `--isotropic-control`. For long re-equilibration, pass e.g. `--t-end 650 --print-every-time 10 --restart-save-every 25`. Do not treat anisotropic results as production defaults until the warning rate, final temperature drop, and representative heat-pipe axial/radial temperature distributions have been compared.
- A 2026-06-07 `500->501 s` short A/B run confirmed the wrapper path and diagnostics. In anisotropic mode, representative high-temperature heat pipes used `k_axial≈3.8e7 W/(m*K)` but `k_radial_model≈88 W/(m*K)`; in isotropic control, `k_radial_model≈k_axial`. The BDF warning is intermittent over 1 s windows, so absence or presence in one short run is not by itself a production decision criterion.
