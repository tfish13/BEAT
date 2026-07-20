from __future__ import annotations

import unittest
from collections import Counter
from pathlib import Path

from validation.run_miri_bounded_validation import (
    REPRESENTATIVE_SEGMENTS,
    SEGMENTS,
    THRESHOLDS,
    definitions,
    validation_fit,
)


class MiriBoundedValidationTests(unittest.TestCase):
    def test_all_segments_and_bounded_case_balance(self) -> None:
        self.assertEqual(len(SEGMENTS), 12)
        self.assertEqual(set(REPRESENTATIVE_SEGMENTS), {"2A", "3B", "4C"})
        cases = definitions()
        self.assertEqual(len(cases), 26)
        reference = [case for case in cases if not case["profile_mismatch"]]
        self.assertEqual(len(reference), 24)
        counts = Counter(case["kind"] for case in reference)
        self.assertEqual(counts, Counter({
            "blank": 6,
            "single_narrow_pixel": 6,
            "single_resolved_half": 6,
            "double_resolved": 6,
        }))
        self.assertEqual(sum(case["profile_mismatch"] for case in cases), 2)

    def test_fit_uses_instrument_lsf_and_two_component_ceiling(self) -> None:
        fit = validation_fit(100000.0, 1000.0, "screening")
        self.assertEqual(fit["lsf"], {"model": "instrument"})
        self.assertEqual(fit["kinematics"]["max_components"], 2)
        self.assertEqual(fit["selection"]["delta_logz"], 5.0)
        self.assertEqual(fit["noise"]["model"], "independent")

    def test_thresholds_are_predeclared_for_bounded_scope(self) -> None:
        self.assertEqual(THRESHOLDS["segment_selection_correct_min"], 24)
        self.assertEqual(THRESHOLDS["single_correct_min"], 12)
        self.assertEqual(THRESHOLDS["double_correct_min"], 5)
        self.assertEqual(THRESHOLDS["reference_evidence_flags_max"], 0)

    def test_runner_freezes_the_complete_passing_matrix(self) -> None:
        source = Path("validation/run_miri_bounded_validation.py").read_text(
            encoding="utf-8"
        )
        self.assertIn('summary["all_bounded_gates_pass"]', source)
        self.assertIn("len(scores) == len(cases)", source)
        self.assertIn('"frozen_gate_manifest.json"', source)


if __name__ == "__main__":
    unittest.main()
