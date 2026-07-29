# V14 210 kW TEC 开路事故：固定功率控制

初态为 `checkpoint_t019865s.npz`。事故时全部 58 根 TFE 的主 TEC 电路永久开路，电流、电功率、电子热流和焦耳热清零，保留 TEC 间隙被动导热。轨道外热保持开启并沿用源 restart 的周期和相位。

越限前控制系统将堆芯功率固定在 210 kW。首次达到通道壁 1058 K、芯块 2700 K、接收极 1023 K、慢化剂 930 K 或反射层 1000 K 时，以该时刻温度和 210 kW 初始化点堆并校零反馈，随后持续施加 -2$ 外加反应性。冷却剂流体温度只记录，不触发停堆。

默认步长 0.05 s，checkpoint 每 50 s。history 在事故后 0-20 s 每 0.1 s、20-100 s 每 1 s、之后每 10 s；停堆后重新执行同一记录计划。总时长同时满足事故后至少一个轨道周期和停堆后至少半个轨道周期。

默认正式输出目录为 `runs/final_from_t019865s/`，包括紧凑历史、冷却剂/固体/电气/反应性详细历史、50 s 重启动、停堆事件与最终重启动。

入口：

```powershell
& "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -u testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_TEC_open_circuit_accident_fixed_power\run_v14_tec_open_circuit_fixed_power.py
```
## 正式计算结果

runs/final_from_t019865s/ 已完成一个轨道周期 5668.144369 s，全过程未触发停堆且水力收敛。全程峰值为通道壁 893.362 K、燃料 2517.546 K、接收极 917.759 K、慢化剂 867.203 K、反射层 808.874 K，均低于限值。最终重启动 SHA256 为 DFDDBB88094BB50124C70399B30526B83273D830E3F47439931A14B937FD8BDA。

正式历史受累计浮点误差影响，在 100 s 边界后额外记录了 101 s 点，随后保持 10 s 间隔；物理推进、温限判断、50 s checkpoint 和最终 restart 不受影响。运行器已修正该边界判定，后续运行从 100 s 直接进入 10 s 节拍。
