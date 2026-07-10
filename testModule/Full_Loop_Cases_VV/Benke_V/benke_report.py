from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Sequence

CASE_DIR = Path(__file__).resolve().parent


def _fmt(value: Any, digits: int = 6) -> str:
    if isinstance(value, float):
        return f"{value:.{digits}g}"
    if isinstance(value, int):
        return str(value)
    if value is None:
        return "N/A"
    return str(value)


def _range_rows(range_checks: dict[str, Any]) -> list[str]:
    rows = ["| 项目 | 计算值 | 范围/容差 | 状态 |", "| --- | ---:| ---:| --- |"]
    for name, check in range_checks.items():
        value = _fmt(check.get("value"))
        if "min" in check and "max" in check:
            bound = f"{_fmt(check['min'])} - {_fmt(check['max'])}"
        else:
            bound = f"abs <= {_fmt(check.get('tolerance_abs_w'))}"
        status = "PASS" if check.get("passed") else "FAIL"
        rows.append(f"| `{name}` | {value} | {bound} | {status} |")
    return rows


def _comparison_section(title: str, comparison: dict[str, Any]) -> list[str]:
    status = comparison.get("status", "unknown")
    lines = [f"## {title}", "", f"状态：`{status}`", ""]
    if status == "missing":
        lines.append(f"缺失文件：`{comparison.get('expected_file')}`")
        lines.append("")
        lines.append(comparison.get("message", "No digitized data available."))
        return lines
    if status != "compared":
        lines.append(json.dumps(comparison, indent=2, ensure_ascii=False))
        return lines

    source_file = comparison.get("source_file")
    if source_file:
        lines.append(f"数据文件：`{source_file}`")
        lines.append("")
    metric_rows = []
    for key, value in comparison.items():
        if key in {"status", "source_file"}:
            continue
        metric_rows.append((key, value))
    if metric_rows:
        lines.extend(["| 指标 | 数值 |", "| --- | ---:|"])
        for key, value in metric_rows:
            lines.append(f"| `{key}` | {_fmt(value)} |")
    return lines


def build_markdown_report(summary: dict[str, Any], run_dir: Path | None = None) -> str:
    validation = summary.get("validation", {})
    case = summary.get("case", {})
    config = summary.get("config", {})
    title = "Benke 热工水力验证结果报告"
    lines = [f"# {title}", ""]
    if run_dir is not None:
        lines.extend([f"运行目录：`{run_dir}`", ""])
    lines.extend(
        [
            "## 验证状态",
            "",
            f"状态：`{validation.get('status', 'unknown')}`",
            "",
            f"范围校核：`{validation.get('range_check_status', 'unknown')}`",
            "",
        ]
    )
    note = validation.get("missing_data_note")
    if note:
        lines.extend([f"说明：{note}", ""])

    lines.extend(
        [
            "## 输入参数",
            "",
            "| 参数 | 数值 |",
            "| --- | ---:|",
            f"| TISA input power W | {_fmt(case.get('tisa_power_w'))} |",
            f"| regulated He pressure torr | {_fmt(case.get('regulated_he_pressure_torr'))} |",
            f"| active length m | {_fmt(config.get('active_length_m'))} |",
            f"| TISA heated length m | {_fmt(config.get('tisa_heated_length_m'))} |",
            f"| coolant heat fraction | {_fmt(config.get('coolant_heat_fraction'))} |",
            f"| water inlet K | {_fmt(config.get('water_inlet_temperature_k'))} |",
            f"| water mass flow kg/s | {_fmt(config.get('water_mass_flow_kg_s'))} |",
            f"| water h W/(m2 K) | {_fmt(config.get('water_h_w_m2_k'))} |",
            f"| regulated He effective k W/(m K) | {_fmt(config.get('regulated_he_effective_k_w_m_k'))} |",
            f"| extra resistance K/W | {_fmt(config.get('extra_resistance_k_per_w'))} |",
            "",
            "## 主要输出",
            "",
            "| 输出 | 数值 |",
            "| --- | ---:|",
            f"| active-zone power W | {_fmt(summary.get('active_zone_power_w'))} |",
            f"| water outlet K | {_fmt(summary.get('water_bulk_outlet_k'))} |",
            f"| water delta-T K | {_fmt(summary.get('water_delta_t_k'))} |",
            f"| energy balance error W | {_fmt(summary.get('energy_balance_error_w'))} |",
            f"| collector inner mean/max K | {_fmt(summary.get('collector_inner_mean_k'))} / {_fmt(summary.get('collector_inner_max_k'))} |",
            f"| sleeve outer mean/max K | {_fmt(summary.get('sleeve_outer_mean_k'))} / {_fmt(summary.get('sleeve_outer_max_k'))} |",
            "",
            "## 文献范围校核",
            "",
        ]
    )
    lines.extend(_range_rows(validation.get("range_checks", {})))
    lines.append("")
    lines.extend(_comparison_section("套筒 12 点热电偶对比", validation.get("sleeve_thermocouple_comparison", {})))
    lines.append("")
    lines.extend(_comparison_section("水侧热平衡对比", validation.get("water_balance_comparison", {})))
    lines.extend(
        [
            "",
            "## 结论",
            "",
        ]
    )
    status = validation.get("status")
    if status == "complete_with_digitized_data":
        lines.append("套筒温度和水侧热平衡均已接入真实数字化数据，可作为 Benke 热工水力量化验证结果。")
    elif status == "quantitative_partial_with_digitized_data":
        lines.append("已有部分真实数字化数据接入，但 Benke 热工水力验证仍不完整；需要补齐缺失的套筒或水侧数据。")
    else:
        lines.append("当前仅完成文献范围校核和模型能量闭合，尚未完成曲线/测点级 Benke 实验验证。")
    lines.append("")
    return "\n".join(lines)


def generate_markdown_report(summary_path: Path | str, output_path: Path | str | None = None) -> Path:
    summary_path = Path(summary_path)
    with summary_path.open("r", encoding="utf-8") as handle:
        summary = json.load(handle)
    output = summary_path.with_name("validation_report.md") if output_path is None else Path(output_path)
    output.write_text(build_markdown_report(summary, run_dir=summary_path.parent), encoding="utf-8")
    return output


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate a Markdown report from a Benke run_summary.json file.")
    parser.add_argument("summary_path", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args(argv)
    output = generate_markdown_report(args.summary_path, args.output)
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
