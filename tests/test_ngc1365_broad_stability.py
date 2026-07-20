from __future__ import annotations

import unittest
from collections import Counter

from validation.run_ngc1365_broad_stability import (
    AUDIT_STOPPING_RULE,
    CONTINUA,
    CRITERIA,
    MODELS,
    case_definitions,
    fit_config,
)


class NGC1365BroadStabilityTests(unittest.TestCase):
    def test_matrix_has_six_models_for_each_of_nine_spaxels(self) -> None:
        cases = case_definitions()
        self.assertEqual(len(cases), 54)
        counts = Counter((case["x"], case["y"]) for case in cases)
        self.assertEqual(set(counts.values()), {6})
        self.assertEqual(set(CONTINUA), {"linear", "quadratic"})
        self.assertEqual(set(MODELS), {"narrow", "broad1", "broad2"})

    def test_broad_priors_are_wider_and_labels_cannot_exchange(self) -> None:
        one = fit_config("linear", "broad1", "screening")
        two = fit_config("linear", "broad2", "screening")
        self.assertEqual(one["broad_components"][0]["sigma_kms"], [400.0, 6000.0])
        self.assertEqual(two["broad_components"][0]["sigma_kms"], [400.0, 1800.0])
        self.assertEqual(two["broad_components"][1]["sigma_kms"], [1800.0, 6000.0])

    def test_continuum_and_noise_stress_are_scoped_consistently(self) -> None:
        linear = fit_config("linear", "broad2", "screening")
        quadratic = fit_config("quadratic", "broad2", "screening")
        self.assertEqual(linear["continuum"]["degree"], 1)
        self.assertEqual(quadratic["continuum"]["degree"], 2)
        self.assertEqual(linear["continuum"]["windows"], quadratic["continuum"]["windows"])
        self.assertEqual(linear["noise"], {"model": "ar1", "rho": "auto", "marginal_scale": "auto"})
        self.assertEqual(linear["kinematics"]["max_components"], 2)

    def test_interpretation_thresholds_are_predeclared(self) -> None:
        self.assertEqual(CRITERIA["minimum_spaxels"], 7)
        self.assertEqual(CRITERIA["continuum_sign_agreement_min"], 8)
        self.assertEqual(CRITERIA["delta_logz_threshold"], 5.0)
        self.assertEqual(CRITERIA["prior_boundary_fraction"], 0.05)
        self.assertIn("screening matrix is terminal", AUDIT_STOPPING_RULE)


if __name__ == "__main__":
    unittest.main()
