from __future__ import annotations

import unittest

from validation.run_survey_1d_regression import (
    N_MALFORMED,
    N_VALID,
    config,
    identifiers,
)


class SurveyRegressionTests(unittest.TestCase):
    def test_large_matrix_counts_and_ids_are_deterministic(self) -> None:
        valid, malformed = identifiers()
        self.assertEqual(len(valid), N_VALID)
        self.assertEqual(len(malformed), N_MALFORMED)
        self.assertEqual(len(set(valid) | set(malformed)), N_VALID + N_MALFORMED)
        self.assertEqual(identifiers(), (valid, malformed))

    def test_workflow_and_science_configs_use_production_survey_input(self) -> None:
        from pathlib import Path

        workflow = config(Path("input.fits"), Path("output"), 0, 1)
        science = config(Path("input.fits"), Path("output"), 2, 1)
        self.assertEqual(workflow["input"]["kind"], "survey_table")
        self.assertEqual(workflow["fit"]["kinematics"]["max_components"], 0)
        self.assertEqual(science["fit"]["kinematics"]["max_components"], 2)
        self.assertEqual(science["fit"]["sampling"]["seed"], 1)
        self.assertEqual(science["fit"]["lines"][2]["ratio_to"], "oiii5007")


if __name__ == "__main__":
    unittest.main()
