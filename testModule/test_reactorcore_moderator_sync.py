import unittest
from unittest.mock import patch

import numpy as np

from Solvers.HeatConduction.Boundary import BoundaryRegion, ResistanceBC
from testModule.run_v13_caseA_closed_loop import build_case, parse_args


class ReactorCoreModeratorSyncTests(unittest.TestCase):
    def test_v13_build_does_not_compute_construction_placeholder_gap_bc(self):
        with patch(
            "sys.argv",
            [
                "run_v13_caseA_closed_loop.py",
                "--duration",
                "1",
                "--disable-tec-coupled",
                "--solid-ode-method",
                "implicit_euler",
            ],
        ):
            args = parse_args()

        placeholder_flux_events = []
        original_compute = BoundaryRegion.compute_net_flux_for_solver

        def monitored_compute(boundary):
            for condition in getattr(boundary, "conditions", []):
                if not isinstance(condition, ResistanceBC):
                    continue
                r_total = (
                    np.asarray(condition.R_ext, dtype=float)
                    + np.asarray(condition.R_add, dtype=float)
                )
                t_ext = np.asarray(condition.T_ext, dtype=float)
                placeholder = np.isfinite(r_total) & (r_total == 0.0) & (t_ext == 300.0)
                if np.any(placeholder):
                    placeholder_flux_events.append(tuple(getattr(boundary, "shape", ())))
            return original_compute(boundary)

        with patch.object(BoundaryRegion, "compute_net_flux_for_solver", monitored_compute):
            build_case(args)

        self.assertEqual(
            placeholder_flux_events,
            [],
            "V13 build must not compute gap-coupler construction placeholders "
            "with R_ext=0 and T_ext=300 K.",
        )

    def test_pre_step_syncs_tfe_gap_couplers_before_virtual_moderator_flux(self):
        with patch(
            "sys.argv",
            [
                "run_v13_caseA_closed_loop.py",
                "--duration",
                "1",
                "--disable-tec-coupled",
                "--solid-ode-method",
                "implicit_euler",
            ],
        ):
            args = parse_args()
        build = build_case(args)
        core = build["core"]

        zero_resistance_flux_events = []
        original_compute = BoundaryRegion.compute_net_flux_for_solver

        def monitored_compute(boundary):
            for condition in getattr(boundary, "conditions", []):
                if not isinstance(condition, ResistanceBC):
                    continue
                r_total = (
                    np.asarray(condition.R_ext, dtype=float)
                    + np.asarray(condition.R_add, dtype=float)
                )
                if np.any(np.isfinite(r_total) & (r_total == 0.0)):
                    zero_resistance_flux_events.append(tuple(getattr(boundary, "shape", ())))
            return original_compute(boundary)

        with patch.object(BoundaryRegion, "compute_net_flux_for_solver", monitored_compute):
            core.pre_step(dt=0.1, current_time=0.0)

        self.assertEqual(
            zero_resistance_flux_events,
            [],
            "ReactorCore.pre_step() must not compute virtual moderator flux "
            "before TFE gap/solid couplers replace their zero-resistance "
            "placeholder boundary conditions.",
        )


if __name__ == "__main__":
    unittest.main()
