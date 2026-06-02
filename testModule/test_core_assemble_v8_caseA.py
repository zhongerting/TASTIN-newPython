import os
import sys
from typing import Any, Dict, Optional, Sequence

import numpy as np

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from test_core_assemble_v5 import COREINPUT_RING_SHARES_RAW
from test_core_assemble_v7_caseA import (
    CASE_A_DESIGN_TOTAL_FLOW_KG_S,
    _case_a_electric_diagnostics,
    _case_a_flow_diagnostics,
    _case_a_reset_design_flows_after_restart,
    build_v7_case_a_system,
)


V8_REPRESENTATIVE_NAMES = ("Center", "Ring1", "Ring2", "Ring3_TEC", "Ring3_Open")
V8_RING_MULTIPLIERS = (1, 6, 12, 15, 3)
V8_TEC_RING_MULTIPLIERS = (1, 6, 12, 15, 0)
V8_RING_MAPPING = {
    "Center": 0,
    "Ring1": 1,
    "Ring2": 2,
    "Ring3_TEC": 3,
    "Ring3_Open": 3,
}
V8_PHYSICAL_RING_COUNT = 4


def build_v8_power_factors(
    representative_names: Sequence[str],
    multipliers: Sequence[int],
    ring_mapping: Dict[str, int],
) -> Dict[str, float]:
    """Assign one per-real-TFE power factor to every representative in a physical ring."""
    shares = np.asarray(COREINPUT_RING_SHARES_RAW, dtype=float)
    shares = shares / float(np.sum(shares))
    ring_totals = {idx: 0 for idx in range(V8_PHYSICAL_RING_COUNT)}
    for name, multiplier in zip(representative_names, multipliers):
        ring_totals[int(ring_mapping[name])] += int(multiplier)

    return {
        name: float(shares[int(ring_mapping[name])] / ring_totals[int(ring_mapping[name])])
        for name in representative_names
    }


def build_v8_case_a_system(
    ring_multipliers: Optional[Sequence[int]] = None,
    tec_ring_multipliers: Optional[Sequence[int]] = None,
    **kwargs: Any,
) -> Dict[str, Any]:
    """Build V8 CaseA with the outer physical ring split into TEC and passive representatives."""
    multipliers = list(V8_RING_MULTIPLIERS if ring_multipliers is None else ring_multipliers)
    tec_multipliers = list(
        V8_TEC_RING_MULTIPLIERS if tec_ring_multipliers is None else tec_ring_multipliers
    )
    if len(multipliers) != len(V8_REPRESENTATIVE_NAMES):
        raise ValueError("V8 ring_multipliers must contain five values.")
    if len(tec_multipliers) != len(V8_REPRESENTATIVE_NAMES):
        raise ValueError("V8 tec_ring_multipliers must contain five values.")

    build = build_v7_case_a_system(
        ring_multipliers=multipliers,
        tec_ring_multipliers=tec_multipliers,
        representative_names=V8_REPRESENTATIVE_NAMES,
        representative_ring_mapping=V8_RING_MAPPING,
        representative_power_factors=build_v8_power_factors(
            V8_REPRESENTATIVE_NAMES,
            multipliers,
            V8_RING_MAPPING,
        ),
        physical_ring_count=V8_PHYSICAL_RING_COUNT,
        core_name="TASTIN_Core_V8_CaseA",
        **kwargs,
    )
    build["case_version"] = "v8"
    build["physical_ring_count"] = V8_PHYSICAL_RING_COUNT
    build["passive_tfe_names"] = ["Ring3_Open"]
    return build


__all__ = [
    "CASE_A_DESIGN_TOTAL_FLOW_KG_S",
    "V8_PHYSICAL_RING_COUNT",
    "V8_REPRESENTATIVE_NAMES",
    "V8_RING_MAPPING",
    "V8_RING_MULTIPLIERS",
    "V8_TEC_RING_MULTIPLIERS",
    "_case_a_electric_diagnostics",
    "_case_a_flow_diagnostics",
    "_case_a_reset_design_flows_after_restart",
    "build_v8_case_a_system",
    "build_v8_power_factors",
]
