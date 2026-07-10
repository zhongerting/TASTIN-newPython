"""V15 V71 variant with center-0.30 m uniform core heating."""

from typing import Any, Dict, Optional, Sequence

import numpy as np

from .common_config import FullLoopCoreConfig, FullLoopFlowConfig, FullLoopPumpConfig
from .v15_case import build_v15_case_a_system, v15_basic_diagnostics
from .v15_pipefin_radiator import V15PipeFinRadiatorConfig


V15_V71_PROFILE_NAME = "center_0p30m_uniform"


def build_center_uniform_axial_power_profile(
    n_lower: int,
    n_active: int,
    n_upper: int,
    heater_length_m: float,
    axial_lengths_m: Sequence[float] = (0.065, 0.377, 0.065),
) -> np.ndarray:
    lengths = tuple(float(value) for value in axial_lengths_m)
    nodes = (int(n_lower), int(n_active), int(n_upper))
    if len(lengths) != 3 or any(value <= 0.0 for value in lengths):
        raise ValueError("axial_lengths_m must contain three positive lengths.")
    if any(value <= 0 for value in nodes):
        raise ValueError("n_lower, n_active, and n_upper must be positive.")
    heater_length = float(heater_length_m)
    full_length = float(sum(lengths))
    if heater_length <= 0.0 or heater_length > full_length:
        raise ValueError("heater_length_m must be positive and no larger than the full axial length.")

    node_lengths = np.array(
        [lengths[0] / nodes[0]] * nodes[0]
        + [lengths[1] / nodes[1]] * nodes[1]
        + [lengths[2] / nodes[2]] * nodes[2],
        dtype=float,
    )
    faces = np.insert(np.cumsum(node_lengths), 0, 0.0)
    heater_start = 0.5 * (full_length - heater_length)
    heater_end = heater_start + heater_length
    overlap = np.maximum(0.0, np.minimum(faces[1:], heater_end) - np.maximum(faces[:-1], heater_start))
    total = float(np.sum(overlap))
    if total <= 0.0:
        raise ValueError("center uniform axial power profile has zero heated length.")
    return overlap / total


def build_v15_v71_case_a_system(
    core_config: Optional[FullLoopCoreConfig] = None,
    flow_config: Optional[FullLoopFlowConfig] = None,
    pump_config: Optional[FullLoopPumpConfig] = None,
    radiator_config: Optional[V15PipeFinRadiatorConfig] = None,
) -> Dict[str, Any]:
    build = build_v15_case_a_system(
        core_config=core_config,
        flow_config=flow_config,
        pump_config=pump_config,
        radiator_config=radiator_config,
    )
    profile = build_center_uniform_axial_power_profile(6, 25, 6, 0.30)
    for tfe in build["tfes"].values():
        tfe.axial_power_profile = profile.copy()
        tfe.solids["pellet"].set_axial_power_profile(profile)

    build["case_version"] = "v15_v71_center0p30_uniform_pipefin_full_loop"
    build["axial_power_profile_name"] = V15_V71_PROFILE_NAME
    build["axial_power_profile"] = profile
    return build


def v15_v71_basic_diagnostics(build: Dict[str, Any]) -> Dict[str, float]:
    diagnostics = v15_basic_diagnostics(build)
    diagnostics["axial_power_profile_sum"] = float(np.sum(build["axial_power_profile"]))
    diagnostics["axial_power_profile_heated_node_count"] = int(np.count_nonzero(build["axial_power_profile"] > 0.0))
    return diagnostics
