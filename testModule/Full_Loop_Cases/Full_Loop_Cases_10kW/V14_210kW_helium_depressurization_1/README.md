# V14 210 kW 氦气瞬时失压事故（轨道外热初态）

## 算例定义

- 初态：V14_210kW_fixed_power_external_heat_2orbits/runs/two_orbits_from13864_20260720/checkpoint_t019265s.npz。
- 事故时刻：5个代表性TFE的接收极—内套管氦气隙同时令 k_gas=0，即 h_He 从 5678 降为 0 W/(m2*K)；间隙辐射保留。
- 保留：NaK水力、流固换热、固体导热、热管、辐射器、轨道外热、TEC和温度反应性反馈。
- 功率：从210 kW固定功率checkpoint接管为点堆动力学，不再施加固定功率源。
- 两种工况：自然发展；事故后5 s施加 -2 dollar。停堆后主TEC电流不高于0.01 A时转为永久开路并停止TEC计算。
- 最长事故时间：2000 s；温限先到则终止。

温限：接收极1500 K、发射极3000 K、冷却剂1058 K、慢化剂930 K、反射层1000 K。

记录频率与LOCA一致：0–20 s每0.5 s，20–100 s每2 s，100–400 s每5 s，400–600 s每10 s，600 s以后每20 s。输出包括history.csv以及冷却剂、固体、电学、反应性四个独立history文件。

## 正式运行

自然发展：

    & "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -u testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization_1\run_v14_helium_depressurization_1.py --duration 2000 --output-dir testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization_1\runs\feedback_noscram_2000s_staged

5 s停堆：

    & "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" -u testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization_1\run_v14_helium_depressurization_1.py --duration 2000 --scram-time 5 --scram-reactivity-dollars -2 --output-dir testModule\Full_Loop_Cases\Full_Loop_Cases_10kW\V14_210kW_helium_depressurization_1\runs\feedback_scram5s_minus2dollar_2000s_staged_tecopen001A
