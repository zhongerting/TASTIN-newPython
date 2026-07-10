# TOPAZ-II V&V 图像包

本图像包为论文 V&V 调研配套图片，主要用于 Origin/WebPlotDigitizer 数字化曲线。图片均来自公开文献 PDF 的图页裁剪/渲染，不包含全文 PDF。

## 目录

- `Paramonov_ElGenk_1994_V71_system/`：TOPAZ-II V-71/TSET 系统级电加热试验图线和回路几何表，适用于系统热工水力、压力、流量、局部温度/欠热度和系统级电输出验证。
- `Venable_1995_single_cell_TFE_main/`：NPS 单电池 TFE 试验台结构图、I-V 曲线、功率-电流曲线、效率表，适用于单 TFE 电输出与热阻模型验证。
- `Venable_1995_single_cell_TFE_appendixB_IV/`：Venable 附录 B 中按铯压排列的更多 I-V 上扫曲线，适用于扩展电特性验证。
- `00_manifest.csv`：所有图片文件、像素尺寸和来源说明。

## 数字化建议

1. 优先数字化 V-71 的堆芯入口/出口温度、系统压力、106 kW 局部温度/欠热度。
2. TFE 电特性优先数字化 Venable Fig. 6-1、Fig. 6-2、Fig. 6-5 至 Fig. 6-8。
3. 每条曲线至少保存：文献、图号、轴范围、单位、数字化日期、工具版本和估计读取误差。
4. 多曲线图建议先手动拾取实验符号点，再视图像质量决定是否自动识别。

## 图片数量

共 `50` 张 PNG。
