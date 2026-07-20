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

    def test_powered_components_balance_classes_donors_and_stressors(self) -> None:
        cases = validation.case_definitions(
            "ngc4151_sivi", "expanded", "powered_components"
        )
        self.assertEqual(len(cases), 20)
        self.assertEqual(
            sum(len(case["components"]) == 2 for case in cases), 10
        )
        self.assertEqual(
            sum(len(case["components"]) == 3 for case in cases), 10
        )
        for donor in ("primary", "alternate"):
            subset = [case for case in cases if case["donor"] == donor]
            self.assertEqual(len(subset), 10)
            self.assertEqual(
                sum(len(case["components"]) == 2 for case in subset), 5
            )
            self.assertEqual(
                sum(len(case["components"]) == 3 for case in subset), 5
            )
        names = {case["name"] for case in cases}
        self.assertTrue(any("sep300_equal" in name for name in names))
        self.assertTrue(any("sep400_ratio025" in name for name in names))
        self.assertTrue(any("adj300_equal" in name for name in names))
        self.assertTrue(any("adj500_weak" in name for name in names))

    def test_g235h_weak_boundary_is_balanced_and_g395h_is_excluded(self) -> None:
        cases = validation.case_definitions(
            "ngc4151_sivi", "expanded", "g235h_weak_boundary"
        )
        self.assertEqual(len(cases), 32)
        self.assertEqual(
            validation.case_definitions(
                "ic5063_bralpha", "expanded", "g235h_weak_boundary"
            ),
            [],
        )
        self.assertEqual(
            {case["effective_weak_component_snr"] for case in cases},
            {3.75, 5.0, 6.25, 10.0},
        )
        for donor in ("primary", "alternate"):
            subset = [case for case in cases if case["donor"] == donor]
            self.assertEqual(len(subset), 16)
            self.assertEqual(
                sum(len(case["components"]) == 2 for case in subset), 8
            )
            self.assertEqual(
                sum(len(case["components"]) == 3 for case in subset), 8
            )
        for donor in ("primary", "alternate"):
            for family, separation in (
                ("double", "sep300"),
                ("double", "sep400"),
                ("triple", "adj400"),
                ("triple", "adj500"),
            ):
                paired = [
                    case
                    for case in cases
                    if case["donor"] == donor
                    and family in case["name"]
                    and separation in case["name"]
                ]
                self.assertEqual(len(paired), 4)
                self.assertEqual(len({case["seed"] for case in paired}), 1)

    def test_automatic_rerun_profile_preserves_standard_rigor(self) -> None:
        specification = validation.DATASETS["ngc4151_sivi"]
        fit = validation.validation_fit(
            specification,
            "standard",
            correlated_noise=True,
            automatic_rerun=True,
        )
        self.assertEqual(fit["selection"]["audit"]["mode"], "rerun")
        audit = fit["selection"]["audit"]["sampling"]
        self.assertGreaterEqual(
            audit["min_num_live_points"],
            fit["sampling"]["min_num_live_points"],
        )
        self.assertLessEqual(audit["dlogz"], fit["sampling"]["dlogz"])

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

    def test_manual_audit_adjudication_matches_rerun_status_semantics(self) -> None:
        accepted = validation.adjudicate_manual_audit(
            {
                "spectrum_id": "accepted",
                "selection_status": "accepted",
                "selection_reasons": [],
            }
        )
        moderate_max = validation.adjudicate_manual_audit(
            {
                "spectrum_id": "moderate-max",
                "selection_status": "convergence_unverified",
                "selection_reasons": [validation.MODERATE_MAX_REASON],
            }
        )
        ambiguous = validation.adjudicate_manual_audit(
            {
                "spectrum_id": "ambiguous",
                "selection_status": "ambiguous",
                "selection_reasons": [
                    "evidence uncertainty overlaps the selection threshold"
                ],
            }
        )
        self.assertEqual(accepted["selection_status"], "accepted_after_audit")
        self.assertEqual(
            moderate_max["selection_status"], "accepted_after_audit"
        )
        self.assertEqual(ambiguous["selection_status"], "ambiguous")

    def test_accepted_after_audit_passes_reliability_gate(self) -> None:
        scores = []
        for replicate in range(20):
            scores.append(
                {
                    "spectrum_id": f"double-{replicate}",
                    "true_components": 2,
                    "recovered_components": 2,
                    "component_count_correct": True,
                    "selection_status": "accepted_after_audit",
                    "runtime_seconds": 1.0,
                    "parameters": [
                        {
                            "velocity_error_kms": 2.0,
                            "sigma_fractional_error": 0.05,
                            "flux_fractional_error": 0.04,
                        }
                    ],
                }
            )
        summary = validation.summarize(scores)
        self.assertTrue(summary["provisional_gates"]["selection_reliability"])
        self.assertEqual(summary["provisional_gate_status"], "pass_scope")


if __name__ == "__main__":
    unittest.main()
