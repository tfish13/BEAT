from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


RUNNER = Path(__file__).resolve().parents[1] / "validation" / "run_nirspec_injection_recovery.py"
SPEC = importlib.util.spec_from_file_location("nirspec_validation", RUNNER)
validation = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(validation)


class NIRSpecValidationTests(unittest.TestCase):
    def test_standard_profile_is_tighter_and_contains_all_count_classes(self) -> None:
        specification = validation.DATASETS["ngc4151_sivi"]
        standard = validation.validation_fit(specification, "standard")
        pilot = validation.validation_fit(specification, "pilot")
        self.assertEqual(standard["sampling"]["nsteps"], 20)
        self.assertGreater(
            standard["sampling"]["min_num_live_points"],
            pilot["sampling"]["min_num_live_points"],
        )
        cases = validation.case_definitions("ngc4151_sivi", "standard")
        self.assertEqual(len(cases), 8)
        self.assertEqual(
            {len(case["components"]) for case in cases}, {0, 1, 2, 3}
        )

    def test_boundary_matrix_covers_requested_stressors(self) -> None:
        cases = validation.case_definitions(
            "ngc4151_sivi", "expanded", "boundary"
        )
        self.assertEqual(len(cases), 10)
        names = {case["name"] for case in cases}
        self.assertTrue(any("snr5" in name for name in names))
        self.assertTrue(any("sigma250" in name for name in names))
        self.assertTrue(any("sep150" in name for name in names))
        self.assertTrue(any("ratio025" in name for name in names))
        self.assertTrue(any("triple_adj300" in name for name in names))
        self.assertEqual(
            sum(case.get("donor") == "alternate" for case in cases), 4
        )

    def test_donor_replicates_balance_donors_and_truth_classes(self) -> None:
        cases = validation.case_definitions(
            "ngc4151_sivi", "expanded", "donor_replicates"
        )
        self.assertEqual(len(cases), 20)
        for donor in ("primary", "alternate"):
            subset = [case for case in cases if case["donor"] == donor]
            self.assertEqual(len(subset), 10)
            self.assertEqual(
                sum(len(case["components"]) == 0 for case in subset), 5
            )
            self.assertEqual(
                sum(len(case["components"]) == 1 for case in subset), 5
            )

    def test_powered_partial_matrix_passes_only_its_scope(self) -> None:
        scores = []
        for count in (0, 1):
            for replicate in range(20):
                scores.append(
                    {
                        "spectrum_id": f"case-{count}-{replicate}",
                        "true_components": count,
                        "recovered_components": count,
                        "component_count_correct": True,
                        "selection_status": "accepted",
                        "runtime_seconds": 1.0,
                        "parameters": []
                        if count == 0
                        else [
                            {
                                "velocity_error_kms": 2.0,
                                "sigma_fractional_error": 0.05,
                                "flux_fractional_error": 0.04,
                            }
                        ],
                    }
                )
        summary = validation.summarize(scores)
        self.assertEqual(summary["provisional_gate_status"], "pass_scope")
        self.assertTrue(summary["statistically_powered"])
        self.assertFalse(summary["full_component_count_coverage"])
        self.assertIsNone(summary["provisional_gates"]["two_component_accuracy"])

    def test_numerical_pass_is_reported_as_underpowered(self) -> None:
        scores = []
        for count in range(4):
            scores.append(
                {
                    "spectrum_id": f"case-{count}",
                    "true_components": count,
                    "recovered_components": count,
                    "component_count_correct": True,
                    "selection_status": "accepted",
                    "runtime_seconds": 1.0,
                    "parameters": []
                    if count == 0
                    else [
                        {
                            "velocity_error_kms": 2.0,
                            "sigma_fractional_error": 0.05,
                            "flux_fractional_error": 0.04,
                        }
                    ],
                }
            )
        summary = validation.summarize(scores)
        self.assertEqual(summary["provisional_gate_status"], "pass_underpowered")
        self.assertFalse(summary["statistically_powered"])
        self.assertTrue(all(summary["provisional_gates"].values()))

    def test_ambiguous_selection_fails_reliability_gate(self) -> None:
        scores = []
        for count in range(4):
            scores.append(
                {
                    "spectrum_id": f"case-{count}",
                    "true_components": count,
                    "recovered_components": count,
                    "component_count_correct": True,
                    "selection_status": "ambiguous" if count == 0 else "accepted",
                    "runtime_seconds": 1.0,
                    "parameters": []
                    if count == 0
                    else [
                        {
                            "velocity_error_kms": 2.0,
                            "sigma_fractional_error": 0.05,
                            "flux_fractional_error": 0.04,
                        }
                    ],
                }
            )
        summary = validation.summarize(scores)
        self.assertEqual(summary["provisional_gate_status"], "fail")
        self.assertFalse(summary["provisional_gates"]["selection_reliability"])


if __name__ == "__main__":
    unittest.main()
