import unittest

import numpy as np

from Components.ExternalHeatSources import OrbitalTableHeatSource
from Components.ExternalHeatSources.embedded_flux_tables import (
    EmbeddedFluxTable,
    EmbeddedFluxTableLibrary,
)
from testModule.Full_Loop_Cases.Full_Loop_Cases_10kW.v14_heatpipe_radiator import (
    V14HeatPipeRadiatorConfig,
    _external_heat_config_for_sector,
)


class ExternalHeatCaseTests(unittest.TestCase):
    def test_restart_time_is_external_heat_phase_zero(self):
        library = EmbeddedFluxTableLibrary({
            0: EmbeddedFluxTable(
                table_id=0,
                name="test",
                time=np.array([0.0, 10.0]),
                values=np.array([100.0, 200.0]),
                periodic=True,
            ),
        })
        source = OrbitalTableHeatSource(
            shape=(1,),
            table_ids=0,
            table_library=library,
            time_origin_s=13864.2,
        )

        self.assertAlmostEqual(source.get_heat_flux(13864.2)[0], 100.0)
        self.assertAlmostEqual(source.get_heat_flux(13869.2)[0], 150.0)

    def test_upper_and_lower_sectors_reuse_the_same_n18_columns(self):
        config = V14HeatPipeRadiatorConfig(
            external_heat_enabled=True,
            external_heat_period_s=5668.144369,
            external_heat_time_origin_s=13864.2,
        )
        marker_library = object()

        upper = _external_heat_config_for_sector(marker_library, 4, config)
        lower = _external_heat_config_for_sector(marker_library, 4, config)

        self.assertEqual(upper["table_ids_by_node"], [12, 13, 14])
        self.assertEqual(lower["table_ids_by_node"], [12, 13, 14])
        self.assertEqual(upper["time_origin_s"], 13864.2)
        self.assertEqual(config.external_heat_period_s, 5668.144369)


if __name__ == "__main__":
    unittest.main()
