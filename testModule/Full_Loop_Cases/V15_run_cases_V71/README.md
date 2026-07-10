# V15_run_cases_V71 运行算例说明

本目录是 `V15_run_cases` 的 V71 副本。拓扑、水力 smoke 设置、泵压头、冷却剂、TEC 开关和管翅式辐射器参数沿用 V15；唯一物理差异是堆芯轴向核功率分布改为中心 `0.30 m` 均匀加热。

## 与 V15 的差异

- 构建入口：`testModule.Full_Loop_Cases.build_v15_v71_case_a_system(...)`
- case_version：`v15_v71_center0p30_uniform_pipefin_full_loop`
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
