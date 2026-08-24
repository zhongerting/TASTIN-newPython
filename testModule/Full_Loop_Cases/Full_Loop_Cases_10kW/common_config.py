from dataclasses import dataclass, field
from typing import Dict, Optional, Sequence, Tuple


REPRESENTATIVE_NAMES = ("Center", "Ring1", "Ring2", "Ring3", "Ring4")
DEFAULT_RING_MULTIPLIERS = (1, 6, 9, 18, 24)
DEFAULT_TEC_RING_MULTIPLIERS = (1, 6, 9, 18, 24)
DEFAULT_RING_MAPPING = {
    "Center": 0,
    "Ring1": 1,
    "Ring2": 2,
    "Ring3": 3,
    "Ring4": 4,
}


@dataclass(frozen=True)
class ReservedParallelTecConfig:
    enabled: bool = False
    mode: str = "fixed_u"
    target_value: float = 0.8
    current_guess_a: float = 6000.0
    multipliers: Optional[Dict[str, int]] = None


@dataclass(frozen=True)
class FullLoopCoreConfig:
    representative_names: Sequence[str] = REPRESENTATIVE_NAMES
    ring_multipliers: Sequence[int] = DEFAULT_RING_MULTIPLIERS
    tec_ring_multipliers: Sequence[int] = DEFAULT_TEC_RING_MULTIPLIERS
    representative_ring_mapping: Dict[str, int] = field(default_factory=lambda: dict(DEFAULT_RING_MAPPING))
    physical_ring_count: int = 5
    coolant_material: str = "SodiumPotassium78"
    inlet_temperature_k: float = 727.0
    reference_pressure_pa: float = 205000.0
    core_name: str = "TOPAZ2_FullLoop_Core"
    solid_heat_capacity_scale: float = 1.0
    solid_heat_capacity_scale_scope: str = "global_outer"
    tec_gap_h_eq_w_m2_k: float = 29.0
    tec_gap_gas: str = "Cesium"
    global_outer_heat_capacity_scale: Optional[float] = None
    main_tec_enabled: bool = True
    main_tec_mode: str = "fixed_u"
    main_tec_target_value: float = 50.5
    main_tec_current_guess_a: float = 150.0
    main_tec_topology: str = "series"
    tec_lookup_enabled: Optional[bool] = None
    tec_lookup_db: Optional[str] = None
    tec_lookup_regions: Optional[Tuple[str, ...]] = None
    reserved_parallel_tec: ReservedParallelTecConfig = field(default_factory=ReservedParallelTecConfig)

@dataclass(frozen=True)
class FullLoopFlowConfig:
    total_flow_kg_s: float = 1.3
    core_connector_volume_m3: float = 1.0e-5
    core_connector_length_m: float = 0.02
    radiator_header_area_m2: float = 3.8e-4
    radiator_header_dh_m: float = 0.014
    core_inlet_segment_length_m: float = 0.13
    core_inlet_segment_area_m2: float = 3.8e-4
    core_inlet_segment_dh_m: float = 0.014
    core_inlet_segment_n_nodes: int = 1
    radiator_inlet_header_length_m: float = 0.13
    radiator_outlet_header_length_m: float = 0.13
    radiator_header_n_nodes: int = 1
    placeholder_bridge_length_m: float = 0.02
    connector_k_loss: float = 0.0
    radiator_bridge_k_loss: float = 0.0


@dataclass(frozen=True)
class FullLoopPumpConfig:
    pump_total_head_pa: float = 7900.0
    pump_count: int = 2
    pump_flow_control: bool = False
    target_flow_kg_s: Optional[float] = None
    pump_node_volume_m3: float = 1.0e-5
    pump_node_length_m: float = 0.02
    pump_area_m2: Optional[float] = None
    pump_dh_m: Optional[float] = None


def validate_sequence_lengths(core_config: FullLoopCoreConfig) -> Tuple[Tuple[str, ...], Tuple[int, ...], Tuple[int, ...]]:
    names = tuple(str(name) for name in core_config.representative_names)
    multipliers = tuple(int(value) for value in core_config.ring_multipliers)
    tec_multipliers = tuple(int(value) for value in core_config.tec_ring_multipliers)
    if not names or len(set(names)) != len(names):
        raise ValueError("representative_names must contain unique names.")
    if len(multipliers) != len(names):
        raise ValueError("ring_multipliers must match representative_names.")
    if len(tec_multipliers) != len(names):
        raise ValueError("tec_ring_multipliers must match representative_names.")
    if any(value <= 0 for value in multipliers):
        raise ValueError("ring_multipliers must be positive.")
    if any(value < 0 for value in tec_multipliers):
        raise ValueError("tec_ring_multipliers must be non-negative.")
    for mult, tec_mult in zip(multipliers, tec_multipliers):
        if tec_mult > mult:
            raise ValueError("tec_ring_multipliers cannot exceed ring_multipliers.")
    if set(core_config.representative_ring_mapping) != set(names):
        raise ValueError("representative_ring_mapping keys must match representative_names.")
    return names, multipliers, tec_multipliers
