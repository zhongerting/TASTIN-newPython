# Nikolaev 1995 源数据提取记录

## PDF

- 文件：`e:\文献阅读\nikolaev1995.pdf`
- 页数：7 页
- 题名：**A single-cell TFE mock-up of the thermionic nuclear power system "Space-R"**
- 出处：AIP Conference Proceedings 324, 815 (1995)
- DOI：`10.1063/1.47120`

## 可定量验证数据

| 文献位置 | 内容 | 当前处理 |
| --- | --- | --- |
| Table 1 | TOPAZ-II experimental 与 SPACE-R prototype calculation 主要参数 | 已录入 `nikolaev_source_data.py` |
| Table 2 | 平均 TFE 在 0.7/0.8/0.9 V 下的 I、Q、Te、效率 | 已录入并作为主要验证表 |
| Table 3 | 燃料最高温度，`Vfree=20/30/40%`，`Kr=1/1.15` | 已录入并用二维插值面闭合 |
| Table 4 | 最大 capillary diameter，`Vfree=20/30/40%`，`Kr=1/1.15` | 已录入并用二维插值面闭合 |
| Figure 4 | Q=2.5-5.0 kW 的 calculated VAC 曲线 | 未数字化，当前只记录热功率范围 |

## OCR/原文歧义

- Table 1 中 `H~n (era)` 由 OCR 得到，当前按 effective height `46/70 cm` 记录，但需要原文图像复核。
- Table 1 中 `Emitter diameter (mm) 1.15/2.3` 与正文“emitter cladding thickening up to 2.3 mm”更一致，当前按 emitter cladding thickness 解释。
- 第 5 页关于 mock-up collector 的句子 OCR 为 `collector with the thickness 17.3 mm`，物理上可能是 collector diameter 或抽取错误；当前未作为强约束。

## 建模结论

这篇文献不是完整实验数据报告，而是 SPACE-R 单电池 TFE mock-up 的设计、计算和初步试验说明。最适合当前阶段验证的是 Table 2 的积分电热工作点，以及 Table 3/4 的设计计算表格。VAC 曲线需要后续数字化后才能做真正的曲线级验证。
