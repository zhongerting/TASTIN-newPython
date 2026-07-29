import numpy as np

from .common_config import FullLoopCoreConfig
from .v14_case import build_v14_case_a_system
from .v14_heatpipe_radiator import (
    LOWER_HP_MULTIPLIERS,
    UPPER_HP_MULTIPLIERS,
    V14HeatPipeRadiatorConfig,
)


def test_v14_thermal_shield_uses_two_ring_sector_mapping_and_hp_multipliers():
    build = build_v14_case_a_system(
        core_config=FullLoopCoreConfig(main_tec_enabled=False),
        radiator_config=V14HeatPipeRadiatorConfig(
            hot_branch_n_nodes=1,
            manifold_node_counts=(1, 1, 1),
            hp_n_con=2,
            n_fin_height=3,
            external_heat_enabled=True,
            thermal_shield_enabled=True,
            thermal_shield_initially_active=True,
        ),
    )

    shield = build['radiator_thermal_shield']
    units = build['radiator_units']
    expected_multipliers = np.asarray(
        [value for sector in UPPER_HP_MULTIPLIERS + LOWER_HP_MULTIPLIERS for value in sector],
        dtype=float,
    )

    assert shield is not None
    assert shield.external_heat_source.table_library.get_table(0).time[-1] == 5668.14
    assert len(units) == 36
    np.testing.assert_allclose(
        [unit.radiation_area_multiplier for unit in units],
        expected_multipliers,
    )
    np.testing.assert_array_equal(
        shield._shield2_sector_indices(len(units)),
        np.repeat(np.arange(12), 3),
    )
    assert sum(unit.radiation_area_multiplier for unit in units) == 340.0
    assert build['system'].components.index(shield) < build['system'].components.index(build['ring_hps'][0])

    temperatures = (700.0, 800.0, 900.0)
    for unit, temperature in zip(units[:3], temperatures):
        unit.last_fin_temperature = np.full((2, 2), temperature)
        unit._has_valid_fin_temperature = True
    tc4_rad, _ = shield._shield2_tc4_rad()
    expected_t4 = np.average(np.asarray(temperatures) ** 4, weights=expected_multipliers[:3])
    np.testing.assert_allclose(tc4_rad[0], expected_t4, rtol=1.0e-12)

    current_time = 123.0
    raw_n6 = shield.external_heat_source.get_heat_flux(current_time)
    shield.pre_step(dt=0.1, current_time=current_time)
    assert shield.last_active is True
    np.testing.assert_allclose(shield.qsss_w_m2[:6], 0.992 * raw_n6)
    np.testing.assert_array_equal(shield.qsss_w_m2[6:], 0.0)
    assert all(source.scale_factor == 0.0 for source in shield.direct_external_heat_sources)
    assert all(np.isfinite(unit.T_space) and unit.T_space > 200.0 for unit in units)

    shield.set_active(False)
    shield.pre_step(dt=0.1, current_time=current_time)
    np.testing.assert_array_equal(shield.qsss_w_m2, 0.0)
    assert all(source.scale_factor == 1.0 for source in shield.direct_external_heat_sources)

    shield.external_heat_source.time_origin_s = 42.0
    state = shield.get_state_dict("Macro_TestShield")
    shield.set_active(True)
    shield.external_heat_source.time_origin_s = 0.0
    shield.load_state_dict(state, "Macro_TestShield")
    assert shield.active_override is False
    assert shield.external_heat_source.time_origin_s == 42.0
    assert all(source.time_origin_s == 42.0 for source in shield.direct_external_heat_sources)


if __name__ == '__main__':
    test_v14_thermal_shield_uses_two_ring_sector_mapping_and_hp_multipliers()
    print('V14 thermal-shield coupling check passed.')
