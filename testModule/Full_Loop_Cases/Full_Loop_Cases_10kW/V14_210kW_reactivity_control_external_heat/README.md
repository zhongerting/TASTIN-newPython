# V14 210 kW 反应性控制与外热流算例

该算例复用 `V14_210kW_reactivity_control`，并显式启用
`Components/ExternalHeatSources/is58p5_w0_8p12_N18_sum.csv`。

- 360 个角度采样点映射到完整的 `5668.144369 s` 轨道周期。
- 上环和下环分别复用 N18 的第 0～17 列。
- 每个 60° 扇区的 3 个代表热管节点依次使用 `3n`、`3n+1`、`3n+2` 列。
- 加载的 restart 绝对时间保持不变，但被定义为本次外热历史的相位 `0 s`。
- 管壁按半圆投影吸热；翅片按单侧投影面积吸热；吸收功率乘以 `0.992` 和表面发射率。

默认入口：

```powershell
& 'E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe' testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_reactivity_control_external_heat\run_v14_210kw_reactivity_control_external_heat.py
```
