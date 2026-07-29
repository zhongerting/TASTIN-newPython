import unittest

try:
    from . import analyze_candidates
except ImportError:
    analyze_candidates = None


class CandidateAnalysisTests(unittest.TestCase):
    def test_stable_target_passes_and_drift_or_band_miss_fails(self):
        self.assertIsNotNone(analyze_candidates)
        stable = [
            {
                "elapsed_s": float(t),
                "tec_main_electric_power_W": 2100.0 + 5.0 * ((t // 10) % 2),
                "tec_main_converged": "True",
                "fluid_converged": "True",
                "core_outlet_T_K": 760.0 + 0.01 * t,
                "hp_condenser_temperature_mean_K": 690.0 + 0.005 * t,
            }
            for t in range(0, 301, 10)
        ]
        accepted = analyze_candidates.analyze_rows(stable, window_s=300.0)
        self.assertTrue(accepted["accepted"])
        self.assertFalse(
            analyze_candidates.analyze_rows(stable[:10], window_s=300.0)["accepted"]
        )

        drifting = [dict(row) for row in stable]
        for row in drifting:
            row["tec_main_electric_power_W"] += 0.2 * row["elapsed_s"]
        self.assertFalse(
            analyze_candidates.analyze_rows(drifting, window_s=300.0)["accepted"]
        )

        out_of_band = [dict(row, tec_main_electric_power_W=2300.0) for row in stable]
        self.assertFalse(
            analyze_candidates.analyze_rows(out_of_band, window_s=300.0)["accepted"]
        )


if __name__ == "__main__":
    unittest.main()
