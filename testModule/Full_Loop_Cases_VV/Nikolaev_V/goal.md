# Nikolaev_V 验证算例目标

最终目标是在 `testModule/Full_Loop_Cases_VV/Nikolaev_V` 下完成 Nikolaev 1995 SPACE-R 单电池 TFE mock-up 文献的验证算例、参数调整记录和验证报告。

## 当前阶段

首版目标是完成表格级验证：

- 固化 Table 1-4 文献数据。
- 建立 TOPAZ-II 相似的单 TFE mock-up 积分模型。
- 允许对文献未说明的闭合参数做全局调整。
- 输出 CSV、JSON 和 Markdown 验证报告。

## 约束

- 不修改主程序 Python/C++ 文件。
- 不用模型输出覆盖文献原始表格数据。
- 所有参数调整都必须记录在 `PARAMETER_ADJUSTMENT_GUIDE.md` 或后续 adjustment log 中。
- 如果后续 Figure 4 数字化后发现当前模型无法闭合，应如实报告并说明缺失边界。

## 完成判据

- `run_nikolaev_validation.py` 能生成一次完整运行目录。
- `validation_report.md` 包含源文献信息、输入参数、Table 2 对比和误差指标。
- 单元测试和 `py_compile` 通过。
