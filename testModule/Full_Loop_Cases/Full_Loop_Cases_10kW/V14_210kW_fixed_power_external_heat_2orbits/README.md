# V14 固定 210 kW + N18 外热流双周期算例

该算例从当前 `checkpoint_t013864s.npz` 继续计算，但保持堆芯总功率固定为
`210000 W`，不启用点堆功率演化。加载 restart 的时刻作为 N18 外热历史的
相位 `0 s`。

- 单周期：`5668.144369 s`
- 总计算时长：`11336.288738 s`
- 上、下集流环分别复用 N18 的 0～17 列
- TEC 设置沿用 restart 邻近的 `run_config.json`
- 默认时间步：`0.05 s`
- 默认每 `10 s` 记录一次，每 `600 s` 保存 checkpoint

运行命令：

```powershell
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_fixed_power_external_heat_2orbits\run_v14_210kw_fixed_power_external_heat_2orbits.py
```
