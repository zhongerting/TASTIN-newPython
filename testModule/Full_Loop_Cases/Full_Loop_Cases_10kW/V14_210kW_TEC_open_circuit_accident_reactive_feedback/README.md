# V14 210 kW TEC 开路事故：温度反应性反馈

初态、TEC 开路模型、轨道外热、温限、checkpoint 和 history 记录计划与固定功率算例相同。

本算例在事故开始时以 210 kW 和 restart 温度初始化点堆，将初始温度反馈校零；随后功率由温度反馈演化。首次达到任一温限时，在已有温度反馈基础上持续叠加 -2$ 外加反应性。冷却剂流体温度只记录，不触发停堆。

默认正式输出目录为 `runs/final_from_t019865s/`，记录计划、详细历史和重启动产物与固定功率算例相同。

入口：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -u testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_TEC_open_circuit_accident_reactive_feedback\run_v14_tec_open_circuit_reactive_feedback.py
```
## 正式计算结果

runs/final_from_t019865s/ 已完成一个轨道周期 5668.144369 s，全过程未触发停堆且水力收敛。全程峰值为通道壁 875.492 K、燃料 2459.505 K、接收极 898.084 K、慢化剂 850.731 K、反射层 797.831 K，均低于限值；末态功率为 8.504 kW。最终重启动 SHA256 为 20137C935CB03D85DAF42699FA9ACC22546D615EA3D0FE04CBAEB5DB39A41212。

本次正式历史具有与固定功率算例相同的 100 s 浮点边界额外记录点；运行器中的节拍判定已修正。
