from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


RUNNER = (
    Path(__file__).resolve().parents[1]
    / "validation"
    / "run_halpha_three_component_calibration.py"
)
SPEC = importlib.util.spec_from_file_location("three_component_calibration", RUNNER)
assert SPEC is not None and SPEC.loader is not None
calibration = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(calibration)


class ThreeComponentCalibrationTests(unittest.TestCase):
    def test_matrix_definitions_are_unique_and_contain_three_components(self) -> None:
        cases = calibration.definitions(2)
        self.assertEqual(len(cases), 18)
        self.assertEqual(len({case["name"] for case in cases}), 18)
        self.assertTrue(all(len(case["components"]) == 3 for case in cases))
        self.assertTrue(
            all(
                [item.velocity_kms for item in case["components"]]
                == [-case["separation_kms"], 0.0, case["separation_kms"]]
                for case in cases
            )
        )

    def test_profiles_allow_three_components_and_differ_in_rigor(self) -> None:
        pilot = calibration.fit_config("pilot")
        tight = calibration.fit_config("tight")
        self.assertEqual(pilot["kinematics"]["max_components"], 3)
        self.assertEqual(tight["kinematics"]["max_components"], 3)
        self.assertGreater(
            tight["sampling"]["min_num_live_points"],
            pilot["sampling"]["min_num_live_points"],
        )
        self.assertLess(tight["sampling"]["dlogz"], pilot["sampling"]["dlogz"])


if __name__ == "__main__":
    unittest.main()
