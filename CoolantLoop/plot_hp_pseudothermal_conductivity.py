from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Materials.Solids.KHP import PotassiumHP
from Materials.Solids.NaHP import SodiumHP
from Materials.Solids.WallMaterial import SS316
from Materials.Solids.WickMaterial import WickMaterial


R_VAPOR = 7.5e-3
R_IN_WALL = 8.0e-3
POROSITY = 0.6


def make_wick(fluid):
    return WickMaterial(
        name=f"{fluid.name}_wick_check",
        solid_mat=SS316(name="SS316_wick_skeleton"),
        fluid_mat=fluid,
        porosity=POROSITY,
        r_vapor=R_VAPOR,
        r_in_wall=R_IN_WALL,
    )


def _polyline(points):
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def _nice_log_ticks(y_min, y_max):
    lo = int(np.floor(np.log10(y_min)))
    hi = int(np.ceil(np.log10(y_max)))
    return [10.0 ** i for i in range(lo, hi + 1)]


def _draw_panel(x, y_na, y_k, x_min, x_max, y_min, y_max, ox, oy, width, height, title, shade=None):
    def sx(v):
        return ox + (v - x_min) / (x_max - x_min) * width

    def sy(v):
        return oy + height - (np.log10(v) - np.log10(y_min)) / (np.log10(y_max) - np.log10(y_min)) * height

    parts = []
    parts.append(f'<text x="{ox + width / 2:.1f}" y="{oy - 26:.1f}" text-anchor="middle" class="title">{title}</text>')

    if shade is not None:
        x0, x1, label = shade
        parts.append(
            f'<rect x="{sx(x0):.2f}" y="{oy:.2f}" width="{sx(x1)-sx(x0):.2f}" '
            f'height="{height:.2f}" fill="#d62728" opacity="0.08"/>'
        )
        parts.append(f'<text x="{sx((x0+x1)/2):.1f}" y="{oy + 18:.1f}" text-anchor="middle" class="note">{label}</text>')

    for tick in np.linspace(x_min, x_max, 7):
        xt = sx(tick)
        parts.append(f'<line x1="{xt:.2f}" y1="{oy:.2f}" x2="{xt:.2f}" y2="{oy+height:.2f}" class="grid"/>')
        parts.append(f'<text x="{xt:.1f}" y="{oy+height+20:.1f}" text-anchor="middle" class="tick">{tick:.0f}</text>')

    for tick in _nice_log_ticks(y_min, y_max):
        yt = sy(tick)
        parts.append(f'<line x1="{ox:.2f}" y1="{yt:.2f}" x2="{ox+width:.2f}" y2="{yt:.2f}" class="grid"/>')
        parts.append(f'<text x="{ox-8:.1f}" y="{yt+4:.1f}" text-anchor="end" class="tick">1e{int(np.log10(tick))}</text>')

    parts.append(f'<rect x="{ox:.2f}" y="{oy:.2f}" width="{width:.2f}" height="{height:.2f}" fill="none" class="axis"/>')

    mask_na = np.isfinite(y_na) & (y_na > 0.0)
    mask_k = np.isfinite(y_k) & (y_k > 0.0)
    pts_na = [(sx(a), sy(b)) for a, b in zip(x[mask_na], y_na[mask_na])]
    pts_k = [(sx(a), sy(b)) for a, b in zip(x[mask_k], y_k[mask_k])]
    parts.append(f'<polyline points="{_polyline(pts_na)}" class="line na"/>')
    parts.append(f'<polyline points="{_polyline(pts_k)}" class="line k"/>')

    parts.append(f'<text x="{ox + width / 2:.1f}" y="{oy + height + 48:.1f}" text-anchor="middle" class="label">Temperature [K]</text>')
    parts.append(
        f'<text x="{ox - 58:.1f}" y="{oy + height / 2:.1f}" text-anchor="middle" '
        f'transform="rotate(-90 {ox - 58:.1f} {oy + height / 2:.1f})" class="label">'
        f'k_pse [W/(m*K)]</text>'
    )
    return "\n".join(parts)


def write_svg(out, t_all, na_all, k_all, t_valid, na_valid, k_valid):
    vals = np.concatenate([na_all, k_all, na_valid, k_valid])
    vals = vals[np.isfinite(vals) & (vals > 0.0)]
    y_min = 10.0 ** np.floor(np.log10(vals.min()))
    y_max = 10.0 ** np.ceil(np.log10(vals.max()))

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1380" height="560" viewBox="0 0 1380 560">
<style>
  text {{ font-family: Arial, sans-serif; fill: #222; }}
  .title {{ font-size: 20px; font-weight: 700; }}
  .label {{ font-size: 15px; }}
  .tick {{ font-size: 12px; fill: #444; }}
  .note {{ font-size: 12px; fill: #7a1c1c; }}
  .grid {{ stroke: #c8c8c8; stroke-width: 1; opacity: 0.55; }}
  .axis {{ stroke: #333; stroke-width: 1.3; }}
  .line {{ fill: none; stroke-width: 3; }}
  .na {{ stroke: #1f77b4; }}
  .k {{ stroke: #d62728; }}
</style>
<rect x="0" y="0" width="1380" height="560" fill="#ffffff"/>
<text x="690" y="34" text-anchor="middle" class="title">Na vs K heat-pipe pseudothermal conductivity</text>
{_draw_panel(t_all, na_all, k_all, 600.0, 1200.0, y_min, y_max, 95, 92, 540, 340, "Full view: 600-1200 K", shade=(700.0, 1033.0, "K main valid range"))}
{_draw_panel(t_valid, na_valid, k_valid, 700.0, 1033.0, y_min, y_max, 765, 92, 540, 340, "Overlap view: 700-1033 K")}
<line x1="536" y1="500" x2="576" y2="500" class="line na"/>
<text x="586" y="505" class="label">Na heat-pipe fluid</text>
<line x1="750" y1="500" x2="790" y2="500" class="line k"/>
<text x="800" y="505" class="label">K heat-pipe fluid</text>
<text x="690" y="536" text-anchor="middle" class="note">Geometry: r_vapor=7.5 mm, r_in_wall=8.0 mm; plotted with WickMaterial._conductivity_pseudothermal_direct()</text>
</svg>
'''
    out.write_text(svg, encoding="utf-8")


def write_focus_svg(out, x, y_na, y_k, x_min, x_max):
    vals = np.concatenate([y_na, y_k])
    vals = vals[np.isfinite(vals) & (vals > 0.0)]
    y_min = 10.0 ** np.floor(np.log10(vals.min()))
    y_max = 10.0 ** np.ceil(np.log10(vals.max()))

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="760" height="560" viewBox="0 0 760 560">
<style>
  text {{ font-family: Arial, sans-serif; fill: #222; }}
  .title {{ font-size: 20px; font-weight: 700; }}
  .label {{ font-size: 15px; }}
  .tick {{ font-size: 12px; fill: #444; }}
  .note {{ font-size: 12px; fill: #555; }}
  .grid {{ stroke: #c8c8c8; stroke-width: 1; opacity: 0.55; }}
  .axis {{ stroke: #333; stroke-width: 1.3; }}
  .line {{ fill: none; stroke-width: 3; }}
  .na {{ stroke: #1f77b4; }}
  .k {{ stroke: #d62728; }}
</style>
<rect x="0" y="0" width="760" height="560" fill="#ffffff"/>
<text x="380" y="34" text-anchor="middle" class="title">Na vs K pseudothermal conductivity: 550-850 K</text>
{_draw_panel(x, y_na, y_k, x_min, x_max, y_min, y_max, 115, 92, 560, 340, "Focused view")}
<line x1="220" y1="500" x2="260" y2="500" class="line na"/>
<text x="270" y="505" class="label">Na heat-pipe fluid</text>
<line x1="430" y1="500" x2="470" y2="500" class="line k"/>
<text x="480" y="505" class="label">K heat-pipe fluid</text>
<text x="380" y="536" text-anchor="middle" class="note">Geometry: r_vapor=7.5 mm, r_in_wall=8.0 mm; raw pseudothermal conductivity</text>
</svg>
'''
    out.write_text(svg, encoding="utf-8")


def main():
    t_all = np.linspace(600.0, 1200.0, 601)
    t_valid_k = np.linspace(700.0, 1033.0, 334)
    t_focus = np.linspace(550.0, 850.0, 301)

    na = SodiumHP(name="HP_Fluid_Na")
    k = PotassiumHP(name="HP_Fluid_K")
    wick_na = make_wick(na)
    wick_k = make_wick(k)

    k_pse_na = wick_na._conductivity_pseudothermal_direct(t_all)
    k_pse_k = wick_k._conductivity_pseudothermal_direct(t_all)
    k_pse_na_valid = wick_na._conductivity_pseudothermal_direct(t_valid_k)
    k_pse_k_valid = wick_k._conductivity_pseudothermal_direct(t_valid_k)
    k_pse_na_focus = wick_na._conductivity_pseudothermal_direct(t_focus)
    k_pse_k_focus = wick_k._conductivity_pseudothermal_direct(t_focus)

    out = Path(__file__).with_name("hp_pseudothermal_conductivity_na_vs_k.svg")
    write_svg(out, t_all, k_pse_na, k_pse_k, t_valid_k, k_pse_na_valid, k_pse_k_valid)
    focus_out = Path(__file__).with_name("hp_pseudothermal_conductivity_na_vs_k_550_850K.svg")
    write_focus_svg(focus_out, t_focus, k_pse_na_focus, k_pse_k_focus, 550.0, 850.0)

    sample_t = np.array([700.0, 800.0, 900.0, 1000.0, 1033.0])
    sample_na = wick_na._conductivity_pseudothermal_direct(sample_t)
    sample_k = wick_k._conductivity_pseudothermal_direct(sample_t)
    print(f"saved={out}")
    print(f"saved_focus={focus_out}")
    print("T[K],k_pse_Na[W/m/K],k_pse_K[W/m/K],K/Na")
    for temp, val_na, val_k in zip(sample_t, sample_na, sample_k):
        print(f"{temp:.1f},{val_na:.9e},{val_k:.9e},{val_k / val_na:.9e}")


if __name__ == "__main__":
    main()
