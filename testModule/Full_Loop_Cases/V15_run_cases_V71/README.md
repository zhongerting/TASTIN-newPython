# V15_run_cases_V71 运行算例说明

本目录是 `V15_run_cases` 的 V71 副本。拓扑、水力 smoke 设置、冷却剂、TEC 开关和管翅式辐射器参数沿用 V15；主要物理差异包括堆芯轴向核功率分布改为中心 `0.30 m` 均匀加热。V71 另将总泵压头设为 `20200 Pa`（两台串联泵各 `10100 Pa`），并将辐射管进、出口局部阻力系数均设为 `K=38`，用于把当前热态总流量标定到约 `1.18 kg/s`。

## 与 V15 的差异

- 构建入口：`testModule.Full_Loop_Cases.build_v15_v71_case_a_system(...)`
- case_version：`v15_v71_center0p30_uniform_pipefin_full_loop`
- 总泵压头：`20200 Pa`，两台串联泵各 `10100 Pa`；辐射管进、出口局部阻力系数均为 `K=38`。`FullLoopFlowConfig.total_flow_kg_s=1.3 kg/s` 仅用于设计/初始流量，实际总流量由水力求解得到约 `1.18 kg/s`。
- 轴向功率分布：`center_0p30m_uniform`
- profile 生成规则：在 TFE 全长 `0.065 + 0.377 + 0.065 = 0.507 m` 中居中取 `0.30 m`，按轴向网格单元与该区间的重叠长度分配功率，并归一化到总和 `1.0`。

环功率份额、代表 TFE 倍率、三条冷回流支路、78 根辐射管和闭式流网拓扑不变。

## 水力-only smoke

结果文件：`v15_v71_flow_path_smoke_result.json`

生成命令：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m unittest testModule.Full_Loop_Cases.V15_run_cases_V71.test_v15_v71_flow_path_smoke
```

测试范围与 V15 smoke 相同：只验证闭式流网连通性、有限压力/流量、唯一压力参考和 V71 profile 元数据；不调用 `SystemManager.initialize_system()`，不执行热耦合器，也不推进固体导热模型。
## 723 K 冷启动与单次 TEC 求解

`run_v15_v71_cold_start.py` 提供两个彼此隔离的阶段：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m testModule.Full_Loop_Cases.V15_run_cases_V71.run_v15_v71_cold_start thermal
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -m testModule.Full_Loop_Cases.V15_run_cases_V71.run_v15_v71_cold_start tec-once --restart testModule\Full_Loop_Cases\V15_run_cases_V71\runs\cold_start_723k_1000s\thermal_1000s_restart.npz
```

固定工况为全系统 `723 K`、空间背景 `4 K`、堆芯总功率 `106 kW`、TEC 在热工阶段关闭、固体导热 `implicit_euler`。热工阶段推进到 `1000 s` 后保存 restart；`tec-once` 重载该状态，启用查表并只执行一次主串联电路 `fixed_u=27.2 V` 求解，不继续推进时间。
