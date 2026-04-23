import json
import os
from typing import Any, Dict

import numpy as np


def _to_list(arr: np.ndarray):
    return np.asarray(arr, dtype=float).tolist()


def extract_center_channel_temperature_field(
    restart_file: str = None,
    output_file: str = None,
) -> str:
    workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..'))

    if restart_file is None:
        restart_file = os.path.join(workspace_root, 'test_core_assemble_v5_restart_t5000.npz')
    if output_file is None:
        output_file = os.path.join(workspace_root, 'center_channel_temperature_field_v5_t5000.json')

    data = np.load(restart_file, allow_pickle=False)

    geom = {
        'r_pellet_inner_m': 4.0e-3,
        'r_pellet_outer_m': 8.5e-3,
        'r_fission_gas_outer_m': 8.65e-3,
        'r_emitter_outer_m': 9.8e-3,
        'r_collector_inner_m': 10.3e-3,
        'r_collector_outer_m': 11.85e-3,
        'r_inner_clad_inner_m': 11.90e-3,
        'r_inner_clad_outer_m': 12.25e-3,
        'r_coolant_inner_m': 12.25e-3,
        'r_coolant_outer_m': 12.95e-3,
        'r_outer_clad_outer_m': 13.30e-3,
        'r_moderator_inner_m': 13.52e-3,
        'r_moderator_outer_m': 16.27e-3,
        'height_m': 0.507,
    }

    axial_lengths = np.array([0.065 / 6] * 6 + [0.377 / 25] * 25 + [0.065 / 6] * 6, dtype=float)
    axial_centers_m = np.insert(np.cumsum(axial_lengths), 0, 0.0)[:-1] + 0.5 * axial_lengths

    pellet_edges = np.linspace(geom['r_pellet_inner_m'], geom['r_pellet_outer_m'], 6)
    pellet_centers_m = 0.5 * (pellet_edges[:-1] + pellet_edges[1:])
    moderator_edges = np.linspace(geom['r_moderator_inner_m'], geom['r_moderator_outer_m'], 4)
    moderator_centers_m = 0.5 * (moderator_edges[:-1] + moderator_edges[1:])

    single_centers_m = {
        'fission_gas_gap': 0.5 * (geom['r_pellet_outer_m'] + geom['r_fission_gas_outer_m']),
        'emitter': 0.5 * (geom['r_fission_gas_outer_m'] + geom['r_emitter_outer_m']),
        'cs_gap': 0.5 * (geom['r_emitter_outer_m'] + geom['r_collector_inner_m']),
        'collector': 0.5 * (geom['r_collector_inner_m'] + geom['r_collector_outer_m']),
        'he_gap': 0.5 * (geom['r_collector_outer_m'] + geom['r_inner_clad_inner_m']),
        'inner_clad': 0.5 * (geom['r_inner_clad_inner_m'] + geom['r_inner_clad_outer_m']),
        'coolant': 0.5 * (geom['r_coolant_inner_m'] + geom['r_coolant_outer_m']),
        'outer_clad': 0.5 * (geom['r_coolant_outer_m'] + geom['r_outer_clad_outer_m']),
        'co2_gap': 0.5 * (geom['r_outer_clad_outer_m'] + geom['r_moderator_inner_m']),
    }

    def arr(key: str) -> np.ndarray:
        return np.asarray(data[key], dtype=float)

    pellet = arr('Solid_Center_Pellet/T').reshape((5, 37))
    emitter = arr('Solid_Center_Emitter/T').reshape((1, 37))[0]
    collector = arr('Solid_Center_Collector/T').reshape((1, 37))[0]
    inner_clad = arr('Solid_Center_InnerClad/T').reshape((1, 37))[0]
    outer_clad = arr('Solid_Center_OuterClad/T').reshape((1, 37))[0]
    moderator = arr('Solid_Center_Moderator/T').reshape((3, 37))
    coolant = arr('Fluid/T_vec')[2:39]

    pellet_right = arr('Solid_Center_Pellet/boundaries/right/T_surface')
    emitter_left = arr('Solid_Center_Emitter/boundaries/left/T_surface')
    emitter_right = arr('Solid_Center_Emitter/boundaries/right/T_surface')
    collector_left = arr('Solid_Center_Collector/boundaries/left/T_surface')
    collector_right = arr('Solid_Center_Collector/boundaries/right/T_surface')
    iclad_left = arr('Solid_Center_InnerClad/boundaries/left/T_surface')
    oclad_right = arr('Solid_Center_OuterClad/boundaries/right/T_surface')
    moderator_left = arr('Solid_Center_Moderator/boundaries/left/T_surface')

    fission_gas_gap = 0.5 * (pellet_right + emitter_left)
    cs_gap = 0.5 * (emitter_right + collector_left)
    he_gap = 0.5 * (collector_right + iclad_left)
    co2_gap = 0.5 * (oclad_right + moderator_left)

    exported: Dict[str, Any] = {
        'source_restart_file': os.path.abspath(restart_file),
        'steady_time_s': float(data['System/global_time'][0]),
        'component': 'Center',
        'notes': [
            'This file contains the full center-channel temperature field extracted from test_core_assemble_v5.',
            'Gap regions are simplified in v5, so gap temperatures are approximated by averaging the adjacent boundary surface temperatures at each axial node.',
        ],
        'geometry': geom,
        'axial_node_lengths_m': _to_list(axial_lengths),
        'axial_node_centers_m': _to_list(axial_centers_m),
        'regions': {
            'pellet': {
                'radial_node_centers_m': _to_list(pellet_centers_m),
                'temperature_K': _to_list(pellet),
            },
            'fission_gas_gap': {
                'radial_center_m': float(single_centers_m['fission_gas_gap']),
                'temperature_K_approx': _to_list(fission_gas_gap),
            },
            'emitter': {
                'radial_center_m': float(single_centers_m['emitter']),
                'temperature_K': _to_list(emitter),
            },
            'cs_gap': {
                'radial_center_m': float(single_centers_m['cs_gap']),
                'temperature_K_approx': _to_list(cs_gap),
            },
            'collector': {
                'radial_center_m': float(single_centers_m['collector']),
                'temperature_K': _to_list(collector),
            },
            'he_gap': {
                'radial_center_m': float(single_centers_m['he_gap']),
                'temperature_K_approx': _to_list(he_gap),
            },
            'inner_clad': {
                'radial_center_m': float(single_centers_m['inner_clad']),
                'temperature_K': _to_list(inner_clad),
            },
            'coolant': {
                'radial_center_m': float(single_centers_m['coolant']),
                'temperature_K': _to_list(coolant),
            },
            'outer_clad': {
                'radial_center_m': float(single_centers_m['outer_clad']),
                'temperature_K': _to_list(outer_clad),
            },
            'co2_gap': {
                'radial_center_m': float(single_centers_m['co2_gap']),
                'temperature_K_approx': _to_list(co2_gap),
            },
            'moderator': {
                'radial_node_centers_m': _to_list(moderator_centers_m),
                'temperature_K': _to_list(moderator),
            },
        },
    }

    with open(output_file, 'w', encoding='utf-8') as fp:
        json.dump(exported, fp, ensure_ascii=False, indent=2)

    return output_file


if __name__ == '__main__':
    out = extract_center_channel_temperature_field()
    print(out)
