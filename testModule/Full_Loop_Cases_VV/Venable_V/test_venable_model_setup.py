from __future__ import annotations

import csv
import json
import sys
import unittest
from pathlib import Path

import numpy as np

# 获取当前文件所在目录并添加到系统路径
CASE_DIR = Path(__file__).resolve().parent
if str(CASE_DIR) not in sys.path:
    sys.path.insert(0, str(CASE_DIR))


class VenableTable71DataTests(unittest.TestCase):
    def test_table71_rows_efficiency_and_cs_pressure_segments(self) -> None:
        # 测试Table 7.1数据的行数、效率和铯压力分段
        from venable_table71_data import iter_table71_cases, pcs_torr_for_qaz

        cases = list(iter_table71_cases())

        self.assertEqual(len(cases), 14)
        self.assertEqual(cases[0].q_az_w, 892.0)
        self.assertEqual(cases[-1].q_az_w, 3162.0)

        for case in cases:
            eta_percent = 100.0 * case.p_out_exp_w / case.q_az_w
            self.assertAlmostEqual(case.eta_exp_percent, eta_percent, delta=0.015)
            self.assertEqual(case.pcs_torr, pcs_torr_for_qaz(case.q_az_w))

        self.assertEqual(pcs_torr_for_qaz(1405.0), 0.4)
        self.assertEqual(pcs_torr_for_qaz(1580.0), 0.5)
        self.assertEqual(pcs_torr_for_qaz(2637.0), 0.8)
        self.assertEqual(pcs_torr_for_qaz(2813.0), 1.0)

    def test_qaz_is_not_reduced_by_benke_active_zone_factor(self) -> None:
        # 测试Q_az是否未被Benke有效区域因子减小
        from venable_table71_data import iter_table71_cases

        cases = list(iter_table71_cases())

        self.assertEqual(cases[0].q_az_w, 892.0)
        self.assertNotAlmostEqual(cases[0].q_az_w, 892.0 * 0.88)


class VenableSingleTfeModelTests(unittest.TestCase):
    def test_thermal_network_heat_source_is_centered_and_conserves_qaz(self) -> None:
        # 测试热网络热源是否居中并守恒Q_az
        from venable_thermal_network import build_tisa_heat_source_w

        heat = build_tisa_heat_source_w(
            q_total_w=3000.0,
            n_nodes=50,
            active_length_m=0.375,
            heated_length_m=0.300,
        )

        self.assertEqual(heat.shape, (50,))
        self.assertAlmostEqual(float(heat.sum()), 3000.0, places=9)
        self.assertLess(float(heat[0]), float(heat[len(heat) // 2]))
        self.assertLess(float(heat[-1]), float(heat[len(heat) // 2]))
        self.assertAlmostEqual(float(heat[0]), float(heat[-1]), places=12)

    def test_thermal_network_energy_balance_and_flow_sensitivity(self) -> None:
        # 测试热网络能量平衡和流量敏感性
        from venable_single_tfe_model import DEFAULT_GEOMETRY
        from venable_table71_data import TABLE71_CASES
        from venable_thermal_network import VenableThermalNetworkConfig, solve_thermal_network

        case = TABLE71_CASES[-1]
        base = solve_thermal_network(case, DEFAULT_GEOMETRY, VenableThermalNetworkConfig())
        high_flow = solve_thermal_network(
            case,
            DEFAULT_GEOMETRY,
            VenableThermalNetworkConfig(cooling_water_mass_flow_kg_s=0.09),
        )

        expected_delta_k = case.q_az_w / (
            base.config.cooling_water_mass_flow_kg_s * base.config.cooling_water_cp_j_kg_k
        )
        self.assertAlmostEqual(
            float(base.water_bulk_outlet_k - base.config.cooling_water_inlet_temperature_k),
            expected_delta_k,
            places=9,
        )
        self.assertLess(abs(base.energy_balance_error_w), 1.0e-9)
        self.assertTrue(np.all(base.emitter_temperature_k > base.collector_temperature_k))
        self.assertTrue(np.all(base.collector_temperature_k > base.water_bulk_temperature_k))
        self.assertLess(float(high_flow.emitter_temperature_k.mean()), float(base.emitter_temperature_k.mean()))
        self.assertLess(float(high_flow.collector_temperature_k.mean()), float(base.collector_temperature_k.mean()))

    def test_build_case_model_supports_thermal_network_mode(self) -> None:
        # 测试build_case_model是否支持热网络模式
        from venable_single_tfe_model import (
            THERMAL_MODEL_THERMAL_NETWORK_V1,
            VenableThermalClosure,
            build_case_model,
        )
        from venable_table71_data import TABLE71_CASES

        model = build_case_model(
            TABLE71_CASES[-1],
            thermal_closure=VenableThermalClosure(
                thermal_model_mode=THERMAL_MODEL_THERMAL_NETWORK_V1,
            ),
        )

        self.assertIsNotNone(model.thermal_network_result)
        self.assertEqual(model.arrays.temitter_k.shape, (1, model.geometry.n_nodes))
        self.assertTrue(np.all(np.isfinite(model.arrays.temitter_k)))
        self.assertTrue(np.all(model.arrays.temitter_k > model.arrays.tcollector_k))

    def test_cesium_pressure_formula_matches_documented_thermocalc_range(self) -> None:
        # 测试铯压力公式是否与文档化的ThermoCalc范围匹配
        from venable_single_tfe_model import (
            cesium_pressure_from_tcs,
            tcs_from_cesium_pressure,
        )

        self.assertAlmostEqual(tcs_from_cesium_pressure(0.02), 441.44, delta=0.02)
        self.assertAlmostEqual(tcs_from_cesium_pressure(5.0), 614.62, delta=0.02)
        self.assertAlmostEqual(cesium_pressure_from_tcs(441.44), 0.02, delta=5.0e-5)
        self.assertAlmostEqual(cesium_pressure_from_tcs(614.62), 5.0, delta=5.0e-3)

    def test_tisa_axial_profile_uses_centered_300mm_heated_length_and_preserves_mean(self) -> None:
        # 测试TISA轴向分布是否使用居中的300mm加热长度并保持平均值
        from venable_single_tfe_model import VenableThermalClosure, build_case_model
        from venable_table71_data import TABLE71_CASES

        closure = VenableThermalClosure(
            emitter_mean_min_k=1500.0,
            emitter_mean_max_k=1500.0,
            collector_mean_min_k=750.0,
            collector_mean_max_k=750.0,
            axial_shape_amplitude=0.2,
            axial_profile_mode="tisa_300mm",
        )
        model = build_case_model(TABLE71_CASES[0], thermal_closure=closure)
        emitter = model.arrays.temitter_k[0]

        self.assertAlmostEqual(float(emitter.mean()), 1500.0, places=9)
        self.assertLess(float(emitter[0]), 1500.0)
        self.assertLess(float(emitter[-1]), 1500.0)
        self.assertGreater(float(emitter[len(emitter) // 2]), 1500.0)

    def test_benke_water_jacket_collector_boundary_responds_to_flow_and_inlet_temperature(self) -> None:
        # 测试Benke水套集电极边界是否响应流量和入口温度
        from venable_single_tfe_model import VenableThermalClosure, build_case_model
        from venable_table71_data import TABLE71_CASES

        base = VenableThermalClosure(
            emitter_mean_min_k=1500.0,
            emitter_mean_max_k=1500.0,
            collector_boundary_mode="benke_water_jacket",
            cooling_water_inlet_temperature_k=310.0,
            cooling_water_mass_flow_kg_s=0.03,
            water_heat_transfer_coefficient_w_m2_k=800.0,
            coolant_heat_pickup_fraction=0.4,
            regulated_he_gap_effective_k_w_m_k=0.08,
        )
        high_flow = VenableThermalClosure(
            emitter_mean_min_k=1500.0,
            emitter_mean_max_k=1500.0,
            collector_boundary_mode="benke_water_jacket",
            cooling_water_inlet_temperature_k=310.0,
            cooling_water_mass_flow_kg_s=0.09,
            water_heat_transfer_coefficient_w_m2_k=800.0,
            coolant_heat_pickup_fraction=0.4,
            regulated_he_gap_effective_k_w_m_k=0.08,
        )
        hot_inlet = VenableThermalClosure(
            emitter_mean_min_k=1500.0,
            emitter_mean_max_k=1500.0,
            collector_boundary_mode="benke_water_jacket",
            cooling_water_inlet_temperature_k=330.0,
            cooling_water_mass_flow_kg_s=0.03,
            water_heat_transfer_coefficient_w_m2_k=800.0,
            coolant_heat_pickup_fraction=0.4,
            regulated_he_gap_effective_k_w_m_k=0.08,
        )

        case = TABLE71_CASES[-1]
        base_model = build_case_model(case, thermal_closure=base)
        high_flow_model = build_case_model(case, thermal_closure=high_flow)
        hot_inlet_model = build_case_model(case, thermal_closure=hot_inlet)

        self.assertLess(
            float(high_flow_model.arrays.tcollector_k.mean()),
            float(base_model.arrays.tcollector_k.mean()),
        )
        self.assertGreater(
            float(hot_inlet_model.arrays.tcollector_k.mean()),
            float(base_model.arrays.tcollector_k.mean()),
        )


if __name__ == "__main__":
    unittest.main()
