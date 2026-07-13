# ThermoCalc 串联 TEC 固定电流设计

## 目标

为现有串联 TEC 电路增加 `fixed_i` 模式。用户给定串联总电流，模型计算端电压和各 TEC 内部电势、焦耳热及热流。现有 `fixed_u`、`fixed_r` 和并联模式保持不变。

## 实现边界

- 在 `CalculationMode` 中增加 `FixedCurrent`。
- 在 `circuitTECs` 中增加 `isFixedI` 和 `iFixedCircuitCalc()`。
- `iFixedCircuitCalc()` 复用现有 `circuitCalc(I)`，不复制串联电路算法。
- `ThermoCalcWrapper.setup_circuit_mode("fixed_i", target_current)` 映射到新模式。
- 不修改上层 ReactorCore 电路拓扑或默认模式。

## 计算与失败语义

给定 `Itarget >= 0` 后执行 `circuitCalc(Itarget)`。

当内部求解收敛、端电压有限且 `Uout > 0` 时：

```text
Iout = Itarget
Uout = circuitCalc(Itarget) 的端电压
converged = true
iteration_count = 1
```

当目标电流求解失败、端电压非有限或 `Uout <= 0` 时，不允许进入外加电流耗电状态，随后执行 `circuitCalc(0)` 恢复开路物理状态：

```text
Iout = 0
Uout = 开路电压
converged = false
iteration_count = 1
```

如果 `circuitCalc(0)` 也失败或得到非有限端电压，则安全返回：

```text
Iout = 0
Uout = 0
converged = false
iteration_count = 1
```

同时输出一次明确警告，避免非有限电输出进入系统计算。

## 验证顺序

1. 先直接验证 `circuitCalc(0)`：覆盖单根、多根串联、解析发射和查表发射，检查端电压、电势、焦耳热和热流全部有限。
2. 验证可实现固定电流：`Iout` 等于目标值、`Uout > 0`、内部结果有限。
3. 验证不可实现固定电流：回退到零电流开路状态并标记未收敛。
4. 验证双重失败安全路径：输出为有限零值并发出警告。
5. 回归 `fixed_u`、`fixed_r` 和三种并联模式。

## 构建隔离

修改后只编译到新的独立构建目录，并通过 `THERMOCALC_PYD_DIR` 加载测试扩展。不得覆盖、重命名或删除 `ThermoCalc/te_solver.cp312-win_amd64.pyd`，因为该生产扩展正在被其他算例使用。