import logging
import os
import sys
import numpy as np

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)

current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Solvers.SystemManager import SystemManager
from Solvers.Hydrodynamics.HydraulicNetwork import HydraulicNetwork
from Solvers.Hydrodynamics.Components import (
    IncompressibleFluidVolume,
    IncompressibleFluidChannel,
    FlowJunction,
    MacroFlowJunction,
    PumpJunction,
    PressurizerVolume,
)
from Solvers.HeatConduction.HeatConduction import HeatConduction2D
from Solvers.HeatConduction.Mesh import Mesh2D
from Materials.Solids.WallMaterial import SS316

import CoolantLoop.model_collector_ring_full_ringhp_geometry100hp_potassium_mixed as cfg
import CoolantLoop.model_collector_ring_6segment_geometry100hp_potassium_mixed as open_case


logging.basicConfig(level=logging.WARNING)


# ---------------------------------------------------------
# Closed-loop operating parameters
# ---------------------------------------------------------
T_SPACE = open_case.T_SPACE
T_INLET_TARGET = open_case.T_INLET
P_REF = open_case.P_OUTLET
T_INIT = open_case.T_INIT
W_TOTAL_TARGET = open_case.W_TOTAL
OPEN_LOOP_STEADY_OUTLET_T = 756.5365807303224

DEFAULT_PUMP_DELTA_P = 2070.0
DEFAULT_HEATER_POWER = 98.0e3
DEFAULT_T_END = 10.0
DEFAULT_PRINT_EVERY_TIME = 1.0
DEFAULT_RESTART_SAVE_EVERY = 10.0
DEFAULT_THERMAL_INIT_RESTART = os.path.join(
    current_dir,
    "collector_ring_6segment_geometry100hp_potassium_mixed_local_loss_630s_from510s_restart.npz",
)

L_HEATER = 1.0
N_HEATER = 10
AREA_HEATER = open_case.AREA_INLET_BUFFER
DH_HEATER = open_case.DH_INLET_BUFFER

L_RETURN = 1.0
N_RETURN = 10
AREA_RETURN = open_case.AREA_OUTLET_BUFFER
DH_RETURN = open_case.DH_OUTLET_BUFFER

PRESSURIZER_LENGTH = 0.20
PRESSURIZER_AREA = AREA_RETURN
PRESSURIZER_DH = DH_RETURN
PRESSURIZER_VOLUME = PRESSURIZER_AREA * PRESSURIZER_LENGTH

K_PUMP_TO_HEATER = 0.0
K_HEATER_TO_INLET_BUFFER = 0.2
K_OUTLET_BUFFER_TO_RETURN = 0.2
K_RETURN_TO_PRESSURIZER = 0.2


nak = open_case.nak


def set_volume_temperature(vol, temperature):
    temperature = float(temperature)
    pressure = float(getattr(vol, "P", P_REF))
    vol.T = temperature
    vol.h = vol.material.enthalpy(temperature, pressure)
    vol.rho = vol.material.density(temperature, pressure)
    vol.mu = vol.material.viscosity(temperature, pressure)


def set_channel_temperature_profile(channel, t_start, t_end):
    if len(channel.volumes) == 1:
        set_volume_temperature(channel.volumes[0], 0.5 * (t_start + t_end))
        return
    for idx, vol in enumerate(channel.volumes):
        frac = idx / (len(channel.volumes) - 1)
        set_volume_temperature(vol, t_start + frac * (t_end - t_start))


def build_sector_solid(name, initial_temp):
    mesh = Mesh2D(
        x_dim=open_case.R_OUT_RING - open_case.R_IN_RING,
        n_x=1,
        y_dim=open_case.L_SECTOR,
        n_y=open_case.N_SECTOR,
        geometry_type="cylindrical",
        inner_radius=open_case.R_IN_RING,
    )
    solid = HeatConduction2D(
        mesh=mesh,
        material=SS316(),
        name=name,
        initial_temp=initial_temp,
    )
    solid.set_ode_method(open_case.RING_WALL_ODE_METHOD)
    bare_area_array = solid.boundaries["right"].area / 2.0
    solid.boundaries["right"].add_dynamic_radiation_condition(
        emissivity=open_case.RING_EMISSIVITY,
        bare_area_array=bare_area_array,
        T_env=T_SPACE,
    )
    return solid


def build_mix_node(name, kind, initial_temp):
    if kind == "inlet":
        area = open_case.AREA_HOT_LEG
        dh = open_case.DH_HOT_LEG
    elif kind == "outlet":
        area = open_case.AREA_MANIFOLD
        dh = open_case.DH_MANIFOLD
    else:
        raise ValueError(f"Unknown mix node kind: {kind}")

    length = 2.0 * dh
    return IncompressibleFluidVolume(
        name=name,
        volume=area * length,
        length=length,
        flow_area=area,
        hydraulic_diam=dh,
        initial_P=P_REF,
        initial_T=initial_temp,
        material=nak,
    )


def _set_channel_flow_guess(channel, mass_flow):
    for junction in channel.internal_junctions:
        junction.W = mass_flow


def _total_fluid_inventory_energy(volumes):
    total = 0.0
    for vol in volumes:
        total += float(vol.rho) * float(vol.vol) * float(vol.h)
    return total


def load_solid_thermal_state(sys_mgr, restart_path):
    if restart_path is None:
        return 0
    if not os.path.exists(restart_path):
        print(f"Thermal initialization restart not found, using analytic initial profile: {restart_path}")
        return 0

    loaded = 0
    with np.load(restart_path, allow_pickle=False) as data:
        data_dict = dict(data)
        for name, solid in sys_mgr.solid_components.items():
            prefix = f"Solid_{name}"
            if f"{prefix}/T" not in data_dict:
                continue
            solid.load_state_dict(data_dict, prefix=prefix)
            loaded += 1
    sys_mgr._sync_solid_times_to_global()
    return loaded


def build_model(
    pump_delta_p=DEFAULT_PUMP_DELTA_P,
    heater_power=DEFAULT_HEATER_POWER,
    initial_temp=T_INIT,
    initialize_from_open_loop_profile=True,
    initial_inlet_temp=T_INLET_TARGET,
    initial_outlet_temp=OPEN_LOOP_STEADY_OUTLET_T,
    thermal_restart_from=DEFAULT_THERMAL_INIT_RESTART,
):
    pressurizer = PressurizerVolume(
        name="Pressurizer",
        volume=PRESSURIZER_VOLUME,
        length=PRESSURIZER_LENGTH,
        flow_area=PRESSURIZER_AREA,
        hydraulic_diam=PRESSURIZER_DH,
        initial_P=P_REF,
        initial_T=initial_temp,
        material=nak,
    )
    pressurizer.set_pressure(P_REF)

    heater_channel = IncompressibleFluidChannel(
        name="HeaterChannel",
        n_nodes=N_HEATER,
        total_length=L_HEATER,
        flow_area=AREA_HEATER,
        hydraulic_diam=DH_HEATER,
        initial_P=P_REF,
        initial_T=initial_temp,
        material=nak,
    )
    inlet_buffer_channel = IncompressibleFluidChannel(
        name="InletBuffer",
        n_nodes=open_case.N_INLET_BUFFER,
        total_length=open_case.L_INLET_BUFFER,
        flow_area=open_case.AREA_INLET_BUFFER,
        hydraulic_diam=open_case.DH_INLET_BUFFER,
        initial_P=P_REF,
        initial_T=initial_temp,
        material=nak,
    )
    outlet_buffer_channel = IncompressibleFluidChannel(
        name="OutletBuffer",
        n_nodes=open_case.N_OUTLET_BUFFER,
        total_length=open_case.L_OUTLET_BUFFER,
        flow_area=open_case.AREA_OUTLET_BUFFER,
        hydraulic_diam=open_case.DH_OUTLET_BUFFER,
        initial_P=P_REF,
        initial_T=initial_temp,
        material=nak,
    )
    return_channel = IncompressibleFluidChannel(
        name="ReturnChannel",
        n_nodes=N_RETURN,
        total_length=L_RETURN,
        flow_area=AREA_RETURN,
        hydraulic_diam=DH_RETURN,
        initial_P=P_REF,
        initial_T=initial_temp,
        material=nak,
    )

    hot_legs = [
        IncompressibleFluidChannel(
            name=f"HotLeg_{idx}",
            n_nodes=open_case.N_HOT_LEG,
            total_length=open_case.L_HOT_LEG,
            flow_area=open_case.AREA_HOT_LEG,
            hydraulic_diam=open_case.DH_HOT_LEG,
            initial_P=P_REF,
            initial_T=initial_temp,
            material=nak,
        )
        for idx in range(1, 4)
    ]
    manifolds = [
        IncompressibleFluidChannel(
            name=f"Manifold_{idx}",
            n_nodes=open_case.N_MANIFOLD,
            total_length=open_case.L_MANIFOLD,
            flow_area=open_case.AREA_MANIFOLD,
            hydraulic_diam=open_case.DH_MANIFOLD,
            initial_P=P_REF,
            initial_T=initial_temp,
            material=nak,
        )
        for idx in range(1, 4)
    ]

    inlet_mix_nodes = {
        key: build_mix_node(
            f"InletMix_{key}",
            "inlet",
            initial_inlet_temp if initialize_from_open_loop_profile else initial_temp,
        )
        for key in open_case.INLET_MIX_KEYS
    }
    outlet_mix_nodes = {
        key: build_mix_node(
            f"OutletMix_{key}",
            "outlet",
            initial_outlet_temp if initialize_from_open_loop_profile else initial_temp,
        )
        for key in open_case.OUTLET_MIX_KEYS
    }
    mix_nodes = {**inlet_mix_nodes, **outlet_mix_nodes}

    sectors = []
    solids = []
    ring_hps = []
    segment_links = []
    segment_entry_links = []
    segment_exit_links = []

    for sector_name, start_key, end_key, multipliers in open_case.SEGMENT_SPECS:
        if initialize_from_open_loop_profile:
            t_sector_start = (
                initial_inlet_temp if start_key in open_case.INLET_MIX_KEYS else initial_outlet_temp
            )
            t_sector_end = (
                initial_inlet_temp if end_key in open_case.INLET_MIX_KEYS else initial_outlet_temp
            )
            t_sector_solid = 0.5 * (t_sector_start + t_sector_end)
        else:
            t_sector_start = initial_temp
            t_sector_end = initial_temp
            t_sector_solid = initial_temp

        channel = IncompressibleFluidChannel(
            name=f"{sector_name}_Channel",
            n_nodes=open_case.N_SECTOR,
            total_length=open_case.L_SECTOR,
            flow_area=open_case.AREA_RING,
            hydraulic_diam=open_case.DH_RING,
            initial_P=P_REF,
            initial_T=t_sector_solid,
            material=nak,
        )
        if initialize_from_open_loop_profile:
            set_channel_temperature_profile(channel, t_sector_start, t_sector_end)
        solid = build_sector_solid(f"{sector_name}_Solid", t_sector_solid)
        ring_hp = cfg.build_ring_hp(
            name=f"{sector_name}_RingHP",
            fluid_channel=channel,
            solid_header=solid,
            hp_multipliers=multipliers,
        )
        open_case.configure_ring_hp_heat_pipe_solver(ring_hp)
        sectors.append(channel)
        solids.append(solid)
        ring_hps.append(ring_hp)

        entry_link = FlowJunction(
            name=f"J_{start_key}_to_{sector_name}",
            from_vol=mix_nodes[start_key],
            to_vol=channel.volumes[0],
            flow_area=open_case.AREA_RING,
            k_loss=open_case.K_INLET_MIX_TO_RING_SEGMENT,
        )
        exit_link = FlowJunction(
            name=f"J_{sector_name}_to_{end_key}",
            from_vol=channel.volumes[-1],
            to_vol=mix_nodes[end_key],
            flow_area=open_case.AREA_RING,
            k_loss=ring_hp.outlet_k_loss + open_case.K_RING_SEGMENT_TO_OUTLET_MIX,
            dynamic_loss_params=ring_hp.outlet_dynamic_loss_params,
        )
        segment_entry_links.append(entry_link)
        segment_exit_links.append(exit_link)
        segment_links.extend([entry_link, exit_link])

    pump_junction = PumpJunction(
        name="J_Pressurizer_Pump_Heater",
        from_vol=pressurizer,
        to_vol=heater_channel.volumes[0],
        flow_area=AREA_HEATER,
        k_loss=K_PUMP_TO_HEATER,
        delta_p=pump_delta_p,
    )
    heater_to_inlet_buffer = FlowJunction(
        name="J_Heater_InletBuffer",
        from_vol=heater_channel.volumes[-1],
        to_vol=inlet_buffer_channel.volumes[0],
        flow_area=AREA_HEATER,
        k_loss=K_HEATER_TO_INLET_BUFFER,
    )
    outlet_buffer_to_return = FlowJunction(
        name="J_OutletBuffer_Return",
        from_vol=outlet_buffer_channel.volumes[-1],
        to_vol=return_channel.volumes[0],
        flow_area=AREA_RETURN,
        k_loss=K_OUTLET_BUFFER_TO_RETURN,
    )
    return_to_pressurizer = FlowJunction(
        name="J_Return_Pressurizer",
        from_vol=return_channel.volumes[-1],
        to_vol=pressurizer,
        flow_area=AREA_RETURN,
        k_loss=K_RETURN_TO_PRESSURIZER,
    )

    inlet_buffer_to_hot_leg = []
    hot_leg_to_inlet_mix = []
    outlet_mix_to_manifold = []
    manifold_to_outlet_buffer = []

    for idx, key in enumerate(open_case.INLET_MIX_KEYS):
        inlet_buffer_to_hot_leg.append(
            MacroFlowJunction(
                name=f"J_InletBuffer_HotLeg_{idx + 1}",
                from_vol=inlet_buffer_channel.volumes[-1],
                to_vol=hot_legs[idx].volumes[0],
                macro_vol=inlet_buffer_channel.volumes[-1],
                multiplier=2,
                flow_area=open_case.AREA_HOT_LEG,
                k_loss=open_case.K_INLET_HEADER_TO_HOT_LEG,
            )
        )
        hot_leg_to_inlet_mix.append(
            FlowJunction(
                name=f"J_HotLeg_{idx + 1}_InletMix_{key}",
                from_vol=hot_legs[idx].volumes[-1],
                to_vol=inlet_mix_nodes[key],
                flow_area=open_case.AREA_HOT_LEG,
                k_loss=open_case.K_HOT_LEG_TO_INLET_MIX,
            )
        )

    for idx, key in enumerate(open_case.OUTLET_MIX_KEYS):
        outlet_mix_to_manifold.append(
            FlowJunction(
                name=f"J_OutletMix_{key}_Manifold_{idx + 1}",
                from_vol=outlet_mix_nodes[key],
                to_vol=manifolds[idx].volumes[0],
                flow_area=open_case.AREA_MANIFOLD,
                k_loss=open_case.K_OUTLET_MIX_TO_MANIFOLD,
            )
        )
        manifold_to_outlet_buffer.append(
            MacroFlowJunction(
                name=f"J_Manifold_{idx + 1}_OutletBuffer",
                from_vol=manifolds[idx].volumes[-1],
                to_vol=outlet_buffer_channel.volumes[0],
                macro_vol=outlet_buffer_channel.volumes[0],
                multiplier=2,
                flow_area=open_case.AREA_MANIFOLD,
                k_loss=open_case.K_MANIFOLD_TO_OUTLET_HEADER,
            )
        )

    all_vols = [pressurizer]
    all_vols.extend(heater_channel.volumes)
    all_vols.extend(inlet_buffer_channel.volumes)
    for channel in hot_legs:
        all_vols.extend(channel.volumes)
    all_vols.extend(inlet_mix_nodes.values())
    all_vols.extend(outlet_mix_nodes.values())
    for channel in sectors:
        all_vols.extend(channel.volumes)
    for channel in manifolds:
        all_vols.extend(channel.volumes)
    all_vols.extend(outlet_buffer_channel.volumes)
    all_vols.extend(return_channel.volumes)

    all_juncs = [pump_junction]
    all_juncs.extend(heater_channel.internal_junctions)
    all_juncs.append(heater_to_inlet_buffer)
    all_juncs.extend(inlet_buffer_channel.internal_junctions)
    all_juncs.extend(inlet_buffer_to_hot_leg)
    all_juncs.extend(hot_leg_to_inlet_mix)
    all_juncs.extend(segment_links)
    all_juncs.extend(outlet_mix_to_manifold)
    all_juncs.extend(manifold_to_outlet_buffer)
    all_juncs.extend(outlet_buffer_channel.internal_junctions)
    all_juncs.append(outlet_buffer_to_return)
    all_juncs.extend(return_channel.internal_junctions)
    all_juncs.append(return_to_pressurizer)
    for channel in hot_legs:
        all_juncs.extend(channel.internal_junctions)
    for channel in sectors:
        all_juncs.extend(channel.internal_junctions)
    for channel in manifolds:
        all_juncs.extend(channel.internal_junctions)

    if initialize_from_open_loop_profile:
        set_volume_temperature(pressurizer, initial_outlet_temp)
        set_channel_temperature_profile(heater_channel, initial_outlet_temp, initial_inlet_temp)
        set_channel_temperature_profile(inlet_buffer_channel, initial_inlet_temp, initial_inlet_temp)
        set_channel_temperature_profile(outlet_buffer_channel, initial_outlet_temp, initial_outlet_temp)
        set_channel_temperature_profile(return_channel, initial_outlet_temp, initial_outlet_temp)
        for channel in hot_legs:
            set_channel_temperature_profile(channel, initial_inlet_temp, initial_inlet_temp)
        for channel in manifolds:
            set_channel_temperature_profile(channel, initial_outlet_temp, initial_outlet_temp)

    w_hot_leg = W_TOTAL_TARGET / 6.0
    pump_junction.W = W_TOTAL_TARGET
    _set_channel_flow_guess(heater_channel, W_TOTAL_TARGET)
    heater_to_inlet_buffer.W = W_TOTAL_TARGET
    _set_channel_flow_guess(inlet_buffer_channel, W_TOTAL_TARGET)
    outlet_buffer_to_return.W = W_TOTAL_TARGET
    _set_channel_flow_guess(return_channel, W_TOTAL_TARGET)
    return_to_pressurizer.W = W_TOTAL_TARGET
    for junction in inlet_buffer_to_hot_leg:
        junction.W = w_hot_leg
    for junction in hot_leg_to_inlet_mix:
        junction.W = w_hot_leg
    for channel in hot_legs:
        _set_channel_flow_guess(channel, w_hot_leg)
    for junction in outlet_mix_to_manifold:
        junction.W = w_hot_leg
    for junction in manifold_to_outlet_buffer:
        junction.W = w_hot_leg
    for channel in manifolds:
        _set_channel_flow_guess(channel, w_hot_leg)

    network = HydraulicNetwork(all_vols, all_juncs, gravity_vector=0.0)
    sys_mgr = SystemManager(fluid_network=network)
    for ring_hp in ring_hps:
        sys_mgr.add_component(ring_hp)
    loaded_thermal_solids = load_solid_thermal_state(sys_mgr, thermal_restart_from)

    heater_power_ref = {"value": float(heater_power)}

    def apply_heater_source(_sys_mgr):
        q_each = heater_power_ref["value"] / len(heater_channel.volumes)
        for vol in heater_channel.volumes:
            vol.Q_vol += q_each

    sys_mgr.add_persistent_fluid_source(apply_heater_source)

    return {
        "pressurizer": pressurizer,
        "heater_channel": heater_channel,
        "inlet_buffer_channel": inlet_buffer_channel,
        "outlet_buffer_channel": outlet_buffer_channel,
        "return_channel": return_channel,
        "hot_legs": hot_legs,
        "manifolds": manifolds,
        "inlet_mix_nodes": inlet_mix_nodes,
        "outlet_mix_nodes": outlet_mix_nodes,
        "sectors": sectors,
        "solids": solids,
        "ring_hps": ring_hps,
        "pump_junction": pump_junction,
        "heater_to_inlet_buffer": heater_to_inlet_buffer,
        "outlet_buffer_to_return": outlet_buffer_to_return,
        "return_to_pressurizer": return_to_pressurizer,
        "inlet_buffer_to_hot_leg": inlet_buffer_to_hot_leg,
        "hot_leg_to_inlet_mix": hot_leg_to_inlet_mix,
        "segment_links": segment_links,
        "segment_entry_links": segment_entry_links,
        "segment_exit_links": segment_exit_links,
        "outlet_mix_to_manifold": outlet_mix_to_manifold,
        "manifold_to_outlet_buffer": manifold_to_outlet_buffer,
        "all_vols": all_vols,
        "all_juncs": all_juncs,
        "network": network,
        "sys_mgr": sys_mgr,
        "heater_power_ref": heater_power_ref,
        "last_fluid_energy": None,
        "thermal_restart_from": thermal_restart_from,
        "loaded_thermal_solids": loaded_thermal_solids,
    }


def get_model_statistics(model):
    solid_nodes = int(sum(np.size(solid.T) for solid in model["solids"]))
    fluid_nodes = int(len(model["all_vols"]))
    flow_junctions = int(len(model["all_juncs"]))
    return solid_nodes, fluid_nodes, flow_junctions


def print_model_summary(model):
    print("Model: closed-loop 6 RingHP segments + heater + pump + pressurizer")
    print(f"Volumes: {len(model['all_vols'])}")
    print(f"Junctions: {len(model['all_juncs'])}")
    print(f"Pressure reference  : Pressurizer at {P_REF:.1f} Pa")
    print(f"Pump delta-p        : {model['pump_junction'].delta_p:.3f} Pa")
    print(f"Heater power        : {model['heater_power_ref']['value']:.3f} W")
    print(f"Thermal restart     : {model['thermal_restart_from']}")
    print(f"Loaded thermal sols : {model['loaded_thermal_solids']}")
    print(f"Target total flow   : {W_TOTAL_TARGET:.6f} kg/s")
    print(f"RingHP segment nodes: {open_case.N_SECTOR} per segment, {6 * open_case.N_SECTOR} total")
    print(f"HP segment totals   : {[sum(spec[3]) for spec in open_case.SEGMENT_SPECS]}")
    print(f"HP/fin emissivity   : {open_case.HP_EMISSIVITY:.3f} / {open_case.FIN_EMISSIVITY:.3f}")
    print(f"Ring emissivity     : {open_case.RING_EMISSIVITY:.3f}, T_space = {T_SPACE:.1f} K")
    print(f"Ring wall ODE method: {open_case.RING_WALL_ODE_METHOD}")
    print(
        "HP solver           : "
        f"wick_anisotropic={open_case.HP_WICK_ANISOTROPIC}, "
        f"k_cap={open_case.HP_WICK_CONDUCTIVITY_CAP}, "
        f"integrator={open_case.HP_TIME_INTEGRATOR}, "
        f"theta={open_case.HP_THETA_IMPLICIT_VALUE:.3f}, "
        f"implicit_boundary={open_case.HP_IMPLICIT_BOUNDARY_LINEARIZATION}"
    )


def print_pre_run_summary(model, case_name):
    solid_nodes, fluid_nodes, flow_junctions = get_model_statistics(model)
    print("=" * 78)
    print(f"Case: {case_name}")
    print("=" * 78)
    print_model_summary(model)
    print("-" * 78)
    print("Pre-run Summary")
    print(f"  Solid nodes    : {solid_nodes}")
    print(f"  Fluid nodes    : {fluid_nodes}")
    print(f"  Flow junctions : {flow_junctions}")
    print("=" * 78)


def _record_history_row(model, current_t, dt):
    pressurizer = model["pressurizer"]
    heater_channel = model["heater_channel"]
    inlet_buffer_channel = model["inlet_buffer_channel"]
    outlet_buffer_channel = model["outlet_buffer_channel"]
    return_channel = model["return_channel"]
    pump_junction = model["pump_junction"]
    heater_to_inlet_buffer = model["heater_to_inlet_buffer"]
    outlet_buffer_to_return = model["outlet_buffer_to_return"]
    return_to_pressurizer = model["return_to_pressurizer"]
    hot_leg_to_inlet_mix = model["hot_leg_to_inlet_mix"]
    outlet_mix_to_manifold = model["outlet_mix_to_manifold"]
    inlet_buffer_to_hot_leg = model["inlet_buffer_to_hot_leg"]
    manifold_to_outlet_buffer = model["manifold_to_outlet_buffer"]
    inlet_mix_nodes = model["inlet_mix_nodes"]
    outlet_mix_nodes = model["outlet_mix_nodes"]
    hot_legs = model["hot_legs"]
    manifolds = model["manifolds"]
    sectors = model["sectors"]
    segment_entry_links = model["segment_entry_links"]
    segment_exit_links = model["segment_exit_links"]

    w_ring_in_total = float(sum(j.W for j in hot_leg_to_inlet_mix))
    w_ring_out_total = float(sum(j.W for j in outlet_mix_to_manifold))
    t_out_avg = float(np.mean([channel.volumes[-1].T for channel in manifolds]))
    h_in = float(inlet_buffer_channel.volumes[-1].h)
    h_out = float(outlet_buffer_channel.volumes[-1].h)
    q_reject_est = float(pump_junction.W) * (h_in - h_out)

    fluid_energy = _total_fluid_inventory_energy(model["all_vols"])
    if model["last_fluid_energy"] is None or dt <= 0.0:
        dudt_fluid = 0.0
    else:
        dudt_fluid = (fluid_energy - model["last_fluid_energy"]) / dt
    model["last_fluid_energy"] = fluid_energy

    fluid_temps = [float(vol.T) for vol in model["all_vols"]]
    q_heater = float(model["heater_power_ref"]["value"])
    row = {
        "time": current_t,
        "dt": dt,
        "W_total_loop": float(pump_junction.W),
        "W_ring_in_total": w_ring_in_total,
        "W_ring_out_total": w_ring_out_total,
        "P_pressurizer": float(pressurizer.P),
        "P_pump_in": float(pump_junction.from_vol.P),
        "P_pump_out": float(pump_junction.to_vol.P),
        "dP_pump": float(pump_junction.delta_p),
        "T_heater_in": float(heater_channel.volumes[0].T),
        "T_heater_out": float(heater_channel.volumes[-1].T),
        "T_inlet_buffer_out": float(inlet_buffer_channel.volumes[-1].T),
        "T_out_avg": t_out_avg,
        "T_return_out": float(return_channel.volumes[-1].T),
        "Q_heater": q_heater,
        "Q_reject_est": q_reject_est,
        "dUdt_fluid_est": dudt_fluid,
        "Q_balance_residual": q_heater - q_reject_est - dudt_fluid,
        "T_loop_min": float(min(fluid_temps)),
        "T_loop_max": float(max(fluid_temps)),
        "P_inlet_buffer_out": float(inlet_buffer_channel.volumes[-1].P),
        "P_outlet_buffer_in": float(outlet_buffer_channel.volumes[0].P),
        "P_outlet_buffer_out": float(outlet_buffer_channel.volumes[-1].P),
        "dP_heater_to_inlet_buffer": open_case.pressure_drop_along_flow(heater_to_inlet_buffer),
        "dP_outlet_buffer_to_return": open_case.pressure_drop_along_flow(outlet_buffer_to_return),
        "dP_return_to_pressurizer": open_case.pressure_drop_along_flow(return_to_pressurizer),
    }

    for idx, key in enumerate(open_case.INLET_MIX_KEYS, start=1):
        node = inlet_mix_nodes[key]
        row[f"T_inlet_mix_{key}"] = float(node.T)
        row[f"P_inlet_mix_{key}"] = float(node.P)
        row[f"W_macro_inlet_to_hotleg_{idx}"] = float(
            inlet_buffer_to_hot_leg[idx - 1].get_mass_flow_for(inlet_buffer_channel.volumes[-1])
        )
        row[f"W_hotleg_to_inlet_mix_{key}"] = float(hot_leg_to_inlet_mix[idx - 1].W)
        row[f"dP_macro_inlet_to_hotleg_{idx}"] = open_case.pressure_drop_along_flow(
            inlet_buffer_to_hot_leg[idx - 1]
        )
        row[f"dP_hotleg_to_inlet_mix_{key}"] = open_case.pressure_drop_along_flow(
            hot_leg_to_inlet_mix[idx - 1]
        )
    for idx, key in enumerate(open_case.OUTLET_MIX_KEYS, start=1):
        node = outlet_mix_nodes[key]
        row[f"T_outlet_mix_{key}"] = float(node.T)
        row[f"P_outlet_mix_{key}"] = float(node.P)
        row[f"W_outlet_mix_to_manifold_{key}"] = float(outlet_mix_to_manifold[idx - 1].W)
        row[f"W_macro_manifold_to_outlet_{idx}"] = float(
            manifold_to_outlet_buffer[idx - 1].get_mass_flow_for(outlet_buffer_channel.volumes[0])
        )
        row[f"dP_outlet_mix_to_manifold_{key}"] = open_case.pressure_drop_along_flow(
            outlet_mix_to_manifold[idx - 1]
        )
        row[f"dP_macro_manifold_to_outlet_{idx}"] = open_case.pressure_drop_along_flow(
            manifold_to_outlet_buffer[idx - 1]
        )
    for idx, channel in enumerate(hot_legs, start=1):
        row[f"T_hotleg_{idx}_out"] = float(channel.volumes[-1].T)
        row[f"P_hotleg_{idx}_in"] = float(channel.volumes[0].P)
        row[f"P_hotleg_{idx}_out"] = float(channel.volumes[-1].P)
    for idx, channel in enumerate(manifolds, start=1):
        row[f"T_manifold_{idx}_out"] = float(channel.volumes[-1].T)
        row[f"P_manifold_{idx}_in"] = float(channel.volumes[0].P)
        row[f"P_manifold_{idx}_out"] = float(channel.volumes[-1].P)
    for idx, channel in enumerate(sectors, start=1):
        row[f"T_A{idx}_in"] = float(channel.volumes[0].T)
        row[f"T_A{idx}_out"] = float(channel.volumes[-1].T)
        row[f"P_A{idx}_in"] = float(channel.volumes[0].P)
        row[f"P_A{idx}_out"] = float(channel.volumes[-1].P)
        row[f"W_A{idx}_entry"] = float(segment_entry_links[idx - 1].W)
        row[f"W_A{idx}_exit"] = float(segment_exit_links[idx - 1].W)
        row[f"dP_A{idx}_entry"] = open_case.pressure_drop_along_flow(segment_entry_links[idx - 1])
        row[f"dP_A{idx}_exit"] = open_case.pressure_drop_along_flow(segment_exit_links[idx - 1])

    return row


def run_case(
    case_name="collector_ring_6segment_geometry100hp_potassium_closed_loop_10s",
    t_end=DEFAULT_T_END,
    min_dt=1.0e-3,
    max_dt=0.2,
    safety_factor=1.0,
    inner_iter=2,
    print_every_time=DEFAULT_PRINT_EVERY_TIME,
    csv_path=None,
    restart_from=None,
    restart_save_path=None,
    restart_save_every=DEFAULT_RESTART_SAVE_EVERY,
    pump_delta_p=DEFAULT_PUMP_DELTA_P,
    heater_power=DEFAULT_HEATER_POWER,
    initial_temp=T_INIT,
    initialize_from_open_loop_profile=True,
    thermal_restart_from=DEFAULT_THERMAL_INIT_RESTART,
):
    model = build_model(
        pump_delta_p=pump_delta_p,
        heater_power=heater_power,
        initial_temp=initial_temp,
        initialize_from_open_loop_profile=initialize_from_open_loop_profile,
        thermal_restart_from=thermal_restart_from,
    )
    sys_mgr = model["sys_mgr"]

    print_pre_run_summary(model, case_name)

    if restart_from is not None:
        sys_mgr.load_global_state(restart_from)
        print(f"Restart loaded from: {restart_from}")
        print(f"Restart time: {sys_mgr.global_time:.6f} s")
    else:
        sys_mgr.initialize_system()
        print("System initialized from closed-loop initial condition.")

    history = []
    next_print_time = cfg.next_event_time(sys_mgr.global_time, print_every_time)
    next_restart_save_time = cfg.next_event_time(sys_mgr.global_time, restart_save_every)

    while sys_mgr.global_time < t_end:
        dt = sys_mgr.compute_adaptive_dt(
            min_dt=min_dt,
            max_dt=max_dt,
            safety_factor=safety_factor,
        )
        if next_restart_save_time is not None:
            if sys_mgr.global_time < next_restart_save_time < sys_mgr.global_time + dt:
                dt = next_restart_save_time - sys_mgr.global_time
        if next_print_time is not None:
            if sys_mgr.global_time < next_print_time < sys_mgr.global_time + dt:
                dt = min(dt, next_print_time - sys_mgr.global_time)
        dt = min(dt, t_end - sys_mgr.global_time)
        if dt <= 1.0e-12:
            break

        sys_mgr.step(dt=dt, inner_iter=inner_iter)
        current_t = sys_mgr.global_time
        row = _record_history_row(model, current_t, dt)
        history.append(row)

        should_print = current_t >= t_end
        if next_print_time is not None and current_t >= next_print_time - 1.0e-12:
            should_print = True
        if should_print:
            print(
                f"t = {current_t:8.3f} s | "
                f"T_in = {row['T_inlet_buffer_out']:.3f} K | "
                f"T_out_avg = {row['T_out_avg']:.3f} K | "
                f"W_loop = {row['W_total_loop']:.4f} kg/s | "
                f"Q_res = {row['Q_balance_residual'] / 1000.0:.3f} kW"
            )
            while next_print_time is not None and current_t >= next_print_time - 1.0e-12:
                next_print_time += print_every_time

        if (
            restart_save_path is not None
            and restart_save_every > 0.0
            and next_restart_save_time is not None
            and current_t >= next_restart_save_time - 1.0e-12
        ):
            checkpoint_path = cfg.restart_checkpoint_path(restart_save_path, next_restart_save_time)
            sys_mgr.save_global_state(checkpoint_path)
            print(f"Restart saved at t={current_t:.3f} s: {checkpoint_path}")
            while next_restart_save_time is not None and current_t >= next_restart_save_time - 1.0e-12:
                next_restart_save_time += restart_save_every

    if csv_path is None:
        csv_path = os.path.join(current_dir, f"{case_name}_history.csv")
    cfg.write_history_csv(csv_path, history)

    if restart_save_path is not None:
        sys_mgr.save_global_state(restart_save_path)
        print(f"Final restart saved: {restart_save_path}")

    print("=" * 70)
    print(f"Case completed: {case_name}")
    print("=" * 70)
    return model, history


if __name__ == "__main__":
    run_case(
        case_name="collector_ring_6segment_geometry100hp_potassium_closed_loop_10s",
        t_end=10.0,
        max_dt=0.2,
        print_every_time=1.0,
        csv_path=os.path.join(current_dir, "collector_ring_6segment_geometry100hp_potassium_closed_loop_10s_history.csv"),
        restart_save_path=os.path.join(current_dir, "collector_ring_6segment_geometry100hp_potassium_closed_loop_10s_restart.npz"),
        restart_save_every=5.0,
    )
