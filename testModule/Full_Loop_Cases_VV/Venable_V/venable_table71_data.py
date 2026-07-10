from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Tuple


@dataclass(frozen=True)
class VenableTable71Case:
    case_id: str
    q_az_w: float
    p_out_exp_w: float
    eta_exp_percent: float
    pcs_torr: float


def pcs_torr_for_qaz(q_az_w: float) -> float:
    """Return the optimal cesium pressure segment reported for Table 7-1."""
    if 892.0 <= q_az_w <= 1405.0:
        return 0.4
    if 1580.0 <= q_az_w <= 2112.0:
        return 0.5
    if 2281.0 <= q_az_w <= 2637.0:
        return 0.8
    if 2813.0 <= q_az_w <= 3162.0:
        return 1.0
    raise ValueError(f"No Venable Table 7-1 optimal Cs pressure for Q_az={q_az_w:g} W")


_TABLE71_POWER_ROWS: Tuple[Tuple[float, float, float], ...] = (
    (892.0, 10.23, 1.15),
    (1062.0, 17.80, 1.68),
    (1237.0, 30.13, 2.44),
    (1405.0, 45.00, 3.20),
    (1580.0, 63.25, 4.01),
    (1755.0, 77.28, 4.40),
    (1933.0, 86.26, 4.46),
    (2112.0, 103.97, 4.92),
    (2281.0, 115.44, 5.06),
    (2474.0, 129.87, 5.25),
    (2637.0, 146.75, 5.57),
    (2813.0, 167.06, 5.94),
    (2999.0, 178.16, 5.94),
    (3162.0, 192.46, 6.09),
)


def _case_id(q_az_w: float) -> str:
    return f"venable_table71_qaz_{int(round(q_az_w)):04d}w"


def iter_table71_cases() -> Iterable[VenableTable71Case]:
    """Yield Venable Table 7-1 maximum-power points.

    The table's active-zone power is already the active-zone heat input Q_az.
    It must not be multiplied by Benke's 0.88 active-zone correction again.
    """
    for q_az_w, p_out_exp_w, eta_exp_percent in _TABLE71_POWER_ROWS:
        yield VenableTable71Case(
            case_id=_case_id(q_az_w),
            q_az_w=q_az_w,
            p_out_exp_w=p_out_exp_w,
            eta_exp_percent=eta_exp_percent,
            pcs_torr=pcs_torr_for_qaz(q_az_w),
        )


TABLE71_CASES: Tuple[VenableTable71Case, ...] = tuple(iter_table71_cases())
