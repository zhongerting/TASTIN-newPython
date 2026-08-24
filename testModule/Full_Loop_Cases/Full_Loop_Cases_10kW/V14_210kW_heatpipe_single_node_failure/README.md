# V14 210 kW 上下对应单节点热管完全失效

本算例从 `V14_210kW_fixed_power_external_heat_2orbits/runs/two_orbits_from13864_20260720/checkpoint_t019865s.npz` 重启动。

- 失效位置：A5（`sector_index=4`）的 `local_node=2`（`theta_014`）。
- 失效定义：上、下集流环对应节点的有效传热能力均为 0。
- 集流环额定热管数量、局部阻力、流量网络、热管固体和外热流边界均保留。
- TEC 正常计算；遮热罩关闭；冷却剂温度只记录，不作为 1058 K 停堆条件。
- 初始固定功率 210 kW；五项固体温限任一达到后切换点堆并施加 -2$。
- `history.csv`、`history_coolant.csv`、`history_solids.csv`、`history_electrical.csv`、`history_reactivity.csv` 每 1 s 记录。
- 周期重启动每 100 s；另存 `accident_start_restart.npz`、`scram_restart.npz` 和 `final_restart.npz`。

运行：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" `
  testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_heatpipe_single_node_failure\run_v14_heatpipe_single_node_failure.py
```

默认运行一个外热流周期；若触发停堆，至少延长计算半个周期。可用 `--duration`、`--dt`、`--output-dir` 覆盖默认参数。
