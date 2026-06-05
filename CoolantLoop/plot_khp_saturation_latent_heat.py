from pathlib import Path
import sys

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from Materials.Solids.KHP import PotassiumHP


def _polyline(points):
    return " ".join(f"{x:.2f},{y:.2f}" for x, y in points)


def _ticks_linear(v_min, v_max, n=7):
    return np.linspace(v_min, v_max, n)


def _ticks_log(v_min, v_max):
    lo = int(np.floor(np.log10(v_min)))
    hi = int(np.ceil(np.log10(v_max)))
    return [10.0 ** i for i in range(lo, hi + 1)]


def _draw_panel(x, series, x_min, x_max, y_min, y_max, ox, oy, width, height, title, y_label, log_y=False):
    def sx(v):
        return ox + (v - x_min) / (x_max - x_min) * width

    if log_y:
        ly_min = np.log10(y_min)
        ly_max = np.log10(y_max)

        def sy(v):
            return oy + height - (np.log10(v) - ly_min) / (ly_max - ly_min) * height

        y_ticks = _ticks_log(y_min, y_max)
        y_label_fmt = lambda v: f"1e{int(np.log10(v))}"
    else:
        def sy(v):
            return oy + height - (v - y_min) / (y_max - y_min) * height

        y_ticks = _ticks_linear(y_min, y_max)
        y_label_fmt = lambda v: f"{v/1e6:.2f}"

    parts = [f'<text x="{ox + width / 2:.1f}" y="{oy - 28:.1f}" text-anchor="middle" class="title">{title}</text>']
    for x_tick in _ticks_linear(x_min, x_max):
        xt = sx(x_tick)
        parts.append(f'<line x1="{xt:.2f}" y1="{oy:.2f}" x2="{xt:.2f}" y2="{oy+height:.2f}" class="grid"/>')
        parts.append(f'<text x="{xt:.1f}" y="{oy+height+20:.1f}" text-anchor="middle" class="tick">{x_tick:.0f}</text>')

    for y_tick in y_ticks:
        yt = sy(y_tick)
        parts.append(f'<line x1="{ox:.2f}" y1="{yt:.2f}" x2="{ox+width:.2f}" y2="{yt:.2f}" class="grid"/>')
        parts.append(f'<text x="{ox-8:.1f}" y="{yt+4:.1f}" text-anchor="end" class="tick">{y_label_fmt(y_tick)}</text>')

    parts.append(f'<rect x="{ox:.2f}" y="{oy:.2f}" width="{width:.2f}" height="{height:.2f}" fill="none" class="axis"/>')

    for cls, y in series:
        mask = np.isfinite(y)
        if log_y:
            mask &= y > 0.0
        pts = [(sx(a), sy(b)) for a, b in zip(x[mask], y[mask])]
        parts.append(f'<polyline points="{_polyline(pts)}" class="{cls}"/>')

    parts.append(f'<text x="{ox + width / 2:.1f}" y="{oy + height + 48:.1f}" text-anchor="middle" class="label">Temperature [K]</text>')
    parts.append(
        f'<text x="{ox - 66:.1f}" y="{oy + height / 2:.1f}" text-anchor="middle" '
        f'transform="rotate(-90 {ox - 66:.1f} {oy + height / 2:.1f})" class="label">{y_label}</text>'
    )
    return "\n".join(parts)


def write_svg(out, t, p_current, h_current):
    p_vals = p_current
    p_vals = p_vals[np.isfinite(p_vals) & (p_vals > 0.0)]
    p_min = 10.0 ** np.floor(np.log10(p_vals.min()))
    p_max = 10.0 ** np.ceil(np.log10(p_vals.max()))

    h_vals = h_current
    h_min = np.floor(h_vals.min() / 5.0e4) * 5.0e4
    h_max = np.ceil(h_vals.max() / 5.0e4) * 5.0e4

    svg = f'''<svg xmlns="http://www.w3.org/2000/svg" width="1380" height="620" viewBox="0 0 1380 620">
<style>
  text {{ font-family: Arial, sans-serif; fill: #222; }}
  .title {{ font-size: 20px; font-weight: 700; }}
  .label {{ font-size: 15px; }}
  .tick {{ font-size: 12px; fill: #444; }}
  .note {{ font-size: 12px; fill: #555; }}
  .grid {{ stroke: #c8c8c8; stroke-width: 1; opacity: 0.55; }}
  .axis {{ stroke: #333; stroke-width: 1.3; }}
  .current {{ fill: none; stroke: #111; stroke-width: 3.5; }}
</style>
<rect x="0" y="0" width="1380" height="620" fill="#ffffff"/>
<text x="690" y="34" text-anchor="middle" class="title">Potassium saturation pressure and latent heat</text>
<line x1="600" y1="76" x2="640" y2="76" class="current"/><text x="650" y="81" class="note">current KHP.py production formula</text>
{_draw_panel(t, [("current", p_current)], 350.0, 1000.0, p_min, p_max, 95, 126, 540, 340, "Saturation pressure", "p_sat [Pa]", log_y=True)}
{_draw_panel(t, [("current", h_current)], 350.0, 1000.0, h_min, h_max, 765, 126, 540, 340, "Latent heat of vaporization", "h_fg [MJ/kg]", log_y=False)}
<text x="690" y="588" text-anchor="middle" class="note">Pressure: ln(p_sat)=25.109-10488/T-0.448 ln(T), 350-1000 K. Latent heat: quadratic fit, 350-900 K.</text>
</svg>
'''
    out.write_text(svg, encoding="utf-8")


def main():
    khp = PotassiumHP()
    t = np.linspace(350.0, 1000.0, 651)
    p_current = khp.saturation_pressure(t)
    h_current = khp.latent_heat(t)

    out = Path(__file__).with_name("khp_saturation_pressure_latent_heat.svg")
    write_svg(out, t, p_current, h_current)

    sample_t = np.array([350.0, 550.0, 600.0, 750.0, 900.0, 1000.0])
    print(f"saved={out}")
    print("T[K],p_sat[Pa],h_fg[J/kg]")
    for temp, p_val, h_val in zip(sample_t, khp.saturation_pressure(sample_t), khp.latent_heat(sample_t)):
        print(f"{temp:.1f},{float(p_val):.9e},{float(h_val):.9e}")


if __name__ == "__main__":
    main()
