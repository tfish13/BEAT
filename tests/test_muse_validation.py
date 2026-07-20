from __future__ import annotations

import unittest
from collections import Counter

from beat.injection import muse_hbeta_oiii_fit
from validation.run_muse_powered_validation import (
    THRESHOLDS,
    definitions,
    fit_config,
    in_supported_red_domain,
)


class MusePoweredValidationTests(unittest.TestCase):
    def test_red_matrix_meets_powered_counts_and_balances_donors(self) -> None:
        cases = definitions()
        counts = Counter(
            (case["band"], case["component_class"]) for case in cases
        )
        self.assertEqual(counts[("red", "blank")], 60)
        self.assertEqual(counts[("red", "single")], 20)
        self.assertEqual(counts[("red", "double")], 24)
        self.assertEqual(counts[("red", "triple")], 20)
        donor_counts = Counter(
            (case["target"], case["band"], case["component_class"])
            for case in cases
        )
        for target in ("ngc2992", "ngc3393"):
            self.assertEqual(donor_counts[(target, "red", "blank")], 30)
            self.assertEqual(donor_counts[(target, "red", "single")], 10)
            self.assertEqual(donor_counts[(target, "red", "double")], 12)
            self.assertEqual(donor_counts[(target, "red", "triple")], 10)

    def test_red_matrix_covers_declared_triple_domain_and_parameter_range(self) -> None:
        triples = [
            case for case in definitions()
            if case["band"] == "red" and case["component_class"] == "triple"
        ]
        spacings = {
            case["components"][1].velocity_kms
            - case["components"][0].velocity_kms
            for case in triples
        }
        widths = {
            component.sigma_kms for case in triples for component in case["components"]
        }
        weak_snrs = {
            min(component.peak_snr * component.flux_scale for component in case["components"])
            for case in triples
        }
        self.assertEqual(spacings, {400.0, 500.0})
        self.assertTrue({60.0, 80.0, 100.0, 120.0, 160.0}.issubset(widths))
        self.assertGreaterEqual(min(weak_snrs), 10.0)

    def test_frozen_supported_red_definition_contains_120_cases(self) -> None:
        supported = [case for case in definitions() if in_supported_red_domain(case)]
        self.assertEqual(len(supported), 120)
        counts = Counter(case["component_class"] for case in supported)
        self.assertEqual(
            counts,
            Counter({"blank": 60, "single": 20, "double": 20, "triple": 20}),
        )

    def test_blue_check_contains_hbeta_oiii_and_wavelength_dependent_lsf(self) -> None:
        fit = muse_hbeta_oiii_fit(max_components=3)
        self.assertEqual(
            [line["name"] for line in fit["lines"]],
            ["hbeta", "oiii5007", "oiii4959"],
        )
        powered = fit_config("blue", "screening")
        self.assertEqual(powered["noise"]["model"], "ar1")
        self.assertEqual(powered["lsf"]["model"], "polynomial_fwhm_angstrom")
        self.assertEqual(THRESHOLDS["minimum_red_cases_per_class"], 20)


if __name__ == "__main__":
    unittest.main()
