from typing import Any, Dict


def full_loop_common_diagnostics(build: Dict[str, Any]) -> Dict[str, Any]:
    net = build["system"].fluid_solver
    fixed_pressure = [
        getattr(vol, "name", "")
        for vol in net.volumes_obj
        if bool(getattr(vol, "is_pressure_boundary", False))
    ]
    pressure_reference = [
        getattr(vol, "name", "")
        for vol in net.volumes_obj
        if bool(getattr(vol, "is_pressure_reference", False))
    ]
    pump_a = build["pump_a"]
    pump_b = build["pump_b"]
    return {
        "case_version": build["case_version"],
        "volume_count": len(net.volumes_obj),
        "junction_count": len(net.junctions_obj),
        "fixed_pressure_boundaries": fixed_pressure,
        "fixed_pressure_boundary_count": len(fixed_pressure),
        "pressure_reference": pressure_reference,
        "pressure_reference_count": len(pressure_reference),
        "pump_total_head_pa": float(pump_a.delta_p + pump_b.delta_p),
        "pump_a_head_pa": float(pump_a.delta_p),
        "pump_b_head_pa": float(pump_b.delta_p),
        "total_flow_design_kg_s": float(build["total_flow_design_kg_s"]),
        "single_tfe_flow_design_kg_s": float(build["single_tfe_flow_design_kg_s"]),
        "ring_multipliers": dict(build["ring_multipliers"]),
        "tec_ring_multipliers": dict(build["tec_ring_multipliers"]),
        "radiator_inlet_header_t_k": float(build["radiator_inlet_header"].volumes[0].T),
        "radiator_outlet_header_t_k": float(build["radiator_outlet_header"].volumes[-1].T),
    }
