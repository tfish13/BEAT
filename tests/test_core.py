from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np

from beat.config import ConfigError, load_config, validate_config
from beat.data import linear_wavelength
from beat.fitting import fit_spectrum
from beat.injection import (
    InjectedComponent,
    analytic_injection_spectrum,
    muse_oiii_fit,
    real_residual_baseline,
    score_recovery,
)
from beat.lsf import FWHM_TO_SIGMA, lsf_sigma_angstrom
from beat.model import ModelDefinition, convolved_sigma, prepare_spectrum, relativistic_doppler
from beat.pipeline import run_pipeline
from beat.spectrum import Spectrum


def fit_config(max_components: int = 2) -> dict:
    return {
        "frame": "rest",
        "window": [4800.0, 5100.0],
        "minimum_valid_pixels": 20,
        "continuum": {"degree": 0, "windows": [[4800, 4830], [5070, 5100]]},
        "kinematics": {
            "max_components": max_components,
            "velocity_kms": [-500.0, 500.0],
            "sigma_kms": [30.0, 500.0],
        },
        "lines": [
            {"name": "hbeta", "wavelength": 4861.33},
            {"name": "oiii5007", "wavelength": 5006.84},
            {
                "name": "oiii4959",
                "wavelength": 4958.92,
                "ratio_to": "oiii5007",
                "ratio": 1.0 / 3.0,
            },
        ],
        "flux_prior": {"min_snr": 0.1, "max_signal_factor": 20},
        "selection": {"delta_logz": 5.0, "stop_when_not_preferred": True},
        "sampling": {"min_num_live_points": 40, "dlogz": 1.0},
    }


def blank_spectrum() -> Spectrum:
    wavelength = np.linspace(4750.0, 5150.0, 4001)
    return Spectrum(
        spectrum_id="synthetic",
        wavelength=wavelength,
        flux=np.zeros_like(wavelength),
        uncertainty=np.ones_like(wavelength),
        redshift=0.0,
    )


class SpectrumTests(unittest.TestCase):
    def test_all_invalid_pixels_raise_a_clear_validation_error(self) -> None:
        spectrum = Spectrum(
            "invalid", np.arange(5.0), np.full(5, np.nan), np.ones(5), 0.1
        )
        with self.assertRaisesRegex(ValueError, "fewer than two valid"):
            spectrum.prepared()

    def test_prepared_sorts_masks_and_deduplicates(self) -> None:
        spectrum = Spectrum(
            "s",
            wavelength=np.array([3.0, 1.0, 2.0, 2.0, np.nan]),
            flux=np.array([30.0, 10.0, 20.0, 21.0, 9.0]),
            uncertainty=np.ones(5),
            redshift=0.1,
            mask=np.array([False, False, False, False, False]),
        ).prepared()
        np.testing.assert_array_equal(spectrum.wavelength, [1.0, 2.0, 3.0])
        np.testing.assert_array_equal(spectrum.flux, [10.0, 20.0, 30.0])

    def test_linear_fits_axis_honors_crpix(self) -> None:
        header = {"CRVAL3": 5000.0, "CRPIX3": 3.0, "CDELT3": 2.0}
        np.testing.assert_array_equal(
            linear_wavelength(header, 5, fits_axis=3),
            [4996.0, 4998.0, 5000.0, 5002.0, 5004.0],
        )


class ModelTests(unittest.TestCase):
    def test_relativistic_doppler_is_reciprocal(self) -> None:
        self.assertAlmostEqual(
            relativistic_doppler(300.0) * relativistic_doppler(-300.0), 1.0, places=12
        )

    def test_excluded_wavelength_window_is_removed_from_fit(self) -> None:
        config = fit_config(max_components=0)
        config["exclude_windows"] = [[4900.0, 4920.0]]
        prepared = prepare_spectrum(blank_spectrum(), config)
        self.assertFalse(
            np.any(
                (prepared.wavelength >= 4900.0)
                & (prepared.wavelength <= 4920.0)
            )
        )

    def test_ar1_noise_correlation_is_estimated_from_continuum(self) -> None:
        config = fit_config(max_components=0)
        config["noise"] = {"model": "ar1", "rho": "auto"}
        wavelength = np.linspace(4750.0, 5150.0, 4001)
        rng = np.random.default_rng(101)
        residual = np.empty(wavelength.size)
        residual[0] = rng.normal()
        for index in range(1, residual.size):
            residual[index] = 0.7 * residual[index - 1] + np.sqrt(1.0 - 0.7**2) * rng.normal()
        spectrum = Spectrum(
            "correlated",
            wavelength,
            residual,
            np.ones_like(wavelength),
            redshift=0.0,
        )
        prepared = prepare_spectrum(spectrum, config)
        self.assertEqual(prepared.noise_model, "ar1")
        self.assertAlmostEqual(prepared.noise_rho, 0.7, delta=0.08)

    def test_ar1_likelihood_matches_closed_form(self) -> None:
        config = fit_config(max_components=0)
        config["noise"] = {"model": "ar1", "rho": 0.5}
        prepared = prepare_spectrum(blank_spectrum(), config)
        model = ModelDefinition(prepared, config, 0)
        params = np.array([0.0])
        expected = -0.5 * (
            prepared.wavelength.size * np.log(2.0 * np.pi)
            + (prepared.wavelength.size - 1) * np.log(1.0 - 0.5**2)
        )
        self.assertAlmostEqual(model.log_likelihood(params), expected)

    def test_auto_marginal_noise_scale_inflates_but_never_shrinks(self) -> None:
        config = fit_config(max_components=0)
        config["noise"] = {"model": "independent", "marginal_scale": "auto"}
        wavelength = np.linspace(4750.0, 5150.0, 4001)
        rng = np.random.default_rng(102)
        high_scatter = Spectrum(
            "high_scatter",
            wavelength,
            rng.normal(0.0, 2.0, wavelength.size),
            np.ones_like(wavelength),
            redshift=0.0,
        )
        low_scatter = Spectrum(
            "low_scatter",
            wavelength,
            rng.normal(0.0, 0.5, wavelength.size),
            np.ones_like(wavelength),
            redshift=0.0,
        )
        self.assertGreater(prepare_spectrum(high_scatter, config).noise_marginal_scale, 1.5)
        self.assertEqual(prepare_spectrum(low_scatter, config).noise_marginal_scale, 1.0)

    def test_line_ratio_and_integrated_flux(self) -> None:
        prepared = prepare_spectrum(blank_spectrum(), fit_config(max_components=1))
        model = ModelDefinition(prepared, fit_config(max_components=1), 1)
        # continuum, velocity, sigma, H-beta flux, [O III] 5007 flux
        params = np.array([0.0, 0.0, 100.0, 100.0, 300.0])
        profile = model.evaluate(params)
        integrated = np.trapz(profile, prepared.wavelength)
        self.assertAlmostEqual(integrated, 500.0, places=4)

    def test_prior_orders_component_velocities(self) -> None:
        config = fit_config(max_components=2)
        prepared = prepare_spectrum(blank_spectrum(), config)
        model = ModelDefinition(prepared, config, 2)
        unit = np.full(model.ndim, 0.5)
        # Layout: cont; v1,s1,f1,f2; v2,s2,f1,f2
        unit[1] = 0.9
        unit[5] = 0.1
        params = model.prior_transform(unit)
        self.assertLess(params[1], params[5])

    def test_three_component_model_is_parameterized_generically(self) -> None:
        config = fit_config(max_components=3)
        prepared = prepare_spectrum(blank_spectrum(), config)
        model = ModelDefinition(prepared, config, 3)
        self.assertEqual(model.n_components, 3)
        self.assertIn("component.3.velocity_kms", model.parameter_names)
        self.assertIn("component.3.oiii5007.flux", model.parameter_names)
        params = model.prior_transform(np.full(model.ndim, 0.5))
        self.assertTrue(np.all(np.isfinite(model.evaluate(params))))

    def test_free_broad_component_is_independent_of_narrow_kinematics(self) -> None:
        config = fit_config(max_components=1)
        config["broad_components"] = [
            {
                "name": "broad_hbeta",
                "line": "hbeta",
                "velocity_kms": [-1000.0, 1000.0],
                "sigma_kms": [600.0, 3000.0],
            }
        ]
        prepared = prepare_spectrum(blank_spectrum(), config)
        model = ModelDefinition(prepared, config, 1)
        self.assertIn("broad.broad_hbeta.velocity_kms", model.parameter_names)
        self.assertIn("broad.broad_hbeta.sigma_kms", model.parameter_names)
        params = model.prior_transform(np.full(model.ndim, 0.5))
        broad_sigma = params[
            model.parameter_names.index("broad.broad_hbeta.sigma_kms")
        ]
        narrow_sigma = params[
            model.parameter_names.index("component.1.sigma_kms")
        ]
        self.assertGreater(broad_sigma, narrow_sigma)
        self.assertTrue(np.all(np.isfinite(model.evaluate(params))))

    def test_resolving_power_lsf_and_quadrature_width(self) -> None:
        center = 5000.0
        instrumental = lsf_sigma_angstrom(
            center, {"model": "resolving_power", "value": 3000.0}
        )
        self.assertAlmostEqual(instrumental, center / 3000.0 * FWHM_TO_SIGMA)
        total = convolved_sigma(
            center,
            0.0,
            100.0,
            0.0,
            {"model": "resolving_power", "value": 3000.0},
            {},
        )
        intrinsic = center * 100.0 / 299_792.458
        self.assertAlmostEqual(total, np.hypot(intrinsic, instrumental))

    def test_polynomial_resolving_power_lsf(self) -> None:
        wavelength = 143000.0
        resolving_power = 4603.0 - 128.0 * 14.3
        sigma = lsf_sigma_angstrom(
            wavelength,
            {
                "model": "polynomial_resolving_power",
                "coefficients": [4603.0, -128.0],
                "scale_angstrom": 1.0e4,
            },
        )
        self.assertAlmostEqual(
            sigma, wavelength / resolving_power * FWHM_TO_SIGMA
        )

    def test_bundled_nirspec_etc_resolving_power_lsf(self) -> None:
        wavelength = 19699.4030166
        sigma = lsf_sigma_angstrom(
            wavelength,
            {"model": "nirspec_etc_resolving_power", "grating": "G235H"},
        )
        resolving_power = wavelength / (sigma / FWHM_TO_SIGMA)
        self.assertAlmostEqual(resolving_power, 2240.452, places=2)

    def test_bundled_nirspec_lsf_rejects_out_of_range_wavelength(self) -> None:
        with self.assertRaisesRegex(ValueError, "outside the bundled NIRSpec"):
            lsf_sigma_angstrom(
                32000.0,
                {"model": "nirspec_etc_resolving_power", "grating": "G235H"},
            )


class FakeSampler:
    def __init__(self, names, loglike, transform):
        self.names = names
        self.loglike = loglike
        self.transform = transform

    def run(self, **kwargs):
        point = self.transform(np.full(len(self.names), 0.5))
        component_numbers = [
            int(name.split(".")[1])
            for name in self.names
            if name.startswith("component.")
        ]
        n_components = max(component_numbers, default=0)
        logz = {0: 0.0, 1: 10.0, 2: 12.0}[n_components]
        return {
            "logz": logz,
            "logzerr": 0.1,
            "ncall": 10,
            "maximum_likelihood": {"point": point, "logl": self.loglike(point)},
            "posterior": {"median": point, "stdev": np.ones_like(point) * 0.1},
        }


class AdaptiveFakeSampler(FakeSampler):
    """Return a moderate third component in pilot mode and reject it tightly."""

    def run(self, **kwargs):
        point = self.transform(np.full(len(self.names), 0.5))
        component_numbers = [
            int(name.split(".")[1])
            for name in self.names
            if name.startswith("component.")
        ]
        n_components = max(component_numbers, default=0)
        tight = int(kwargs["min_num_live_points"]) >= 100
        logz = (
            {0: 0.0, 1: 10.0, 2: 20.0, 3: 20.0}
            if tight
            else {0: 0.0, 1: 10.0, 2: 20.0, 3: 26.0}
        )[n_components]
        return {
            "logz": logz,
            "logzerr": 0.1,
            "ncall": 10,
            "maximum_likelihood": {"point": point, "logl": self.loglike(point)},
            "posterior": {"median": point, "stdev": np.ones_like(point) * 0.1},
        }


class FitTests(unittest.TestCase):
    def test_evidence_selection_stops_after_rejected_model(self) -> None:
        result = fit_spectrum(blank_spectrum(), fit_config(), sampler_factory=FakeSampler)
        self.assertEqual(result["selected_components"], 1)
        self.assertEqual([item["n_components"] for item in result["models"]], [0, 1, 2])
        line_fluxes = result["components"][0]["lines"]
        self.assertAlmostEqual(
            line_fluxes["oiii4959"]["flux"], line_fluxes["oiii5007"]["flux"] / 3.0
        )

    def test_moderate_max_component_selection_is_flagged(self) -> None:
        config = fit_config(max_components=3)
        config["sampling"].update(
            {"min_num_live_points": 40, "min_ess": 40, "dlogz": 2.0}
        )
        config["selection"]["audit"] = {
            "mode": "flag",
            "uncertainty_sigma": 1.0,
            "minimum_margin": 0.5,
            "max_component_decisive_delta_logz": 20.0,
        }
        result = fit_spectrum(
            blank_spectrum(), config, sampler_factory=AdaptiveFakeSampler
        )
        self.assertEqual(result["selected_components"], 3)
        self.assertEqual(result["selection_status"], "convergence_unverified")
        self.assertFalse(result["selection_audit"]["performed"])
        self.assertIn(
            "maximum component count is selected without decisive evidence",
            result["selection_audit"]["trigger_reasons"],
        )

    def test_adaptive_audit_can_overturn_max_component_selection(self) -> None:
        config = fit_config(max_components=3)
        config["sampling"].update(
            {"min_num_live_points": 40, "min_ess": 40, "dlogz": 2.0}
        )
        config["selection"]["audit"] = {
            "mode": "rerun",
            "uncertainty_sigma": 1.0,
            "minimum_margin": 0.5,
            "max_component_decisive_delta_logz": 20.0,
            "sampling": {
                "min_num_live_points": 100,
                "min_ess": 200,
                "dlogz": 0.5,
                "nsteps": 20,
            },
        }
        result = fit_spectrum(
            blank_spectrum(), config, sampler_factory=AdaptiveFakeSampler
        )
        self.assertEqual(result["selection_audit"]["pilot_selected_components"], 3)
        self.assertEqual(result["selected_components"], 2)
        self.assertTrue(result["selection_audit"]["performed"])
        self.assertEqual(result["selection_status"], "accepted_after_audit")
        self.assertEqual(
            result["selection_audit"]["audit_sampling"]["min_num_live_points"],
            100,
        )
        self.assertAlmostEqual(
            result["selection_diagnostics"][-1]["delta_logz"], 0.0
        )


class InjectionTests(unittest.TestCase):
    def test_analytic_injection_is_deterministic_and_records_truth(self) -> None:
        fit = muse_oiii_fit()
        components = [InjectedComponent(25.0, 80.0, 10.0)]
        first, truth = analytic_injection_spectrum(
            "injected", 0.01, components, fit, seed=123
        )
        second, _ = analytic_injection_spectrum(
            "injected", 0.01, components, fit, seed=123
        )
        np.testing.assert_array_equal(first.flux, second.flux)
        self.assertEqual(truth["true_components"], 1)
        self.assertGreater(truth["components"][0]["lines"]["oiii5007"]["flux"], 0)

    def test_real_residual_baseline_excludes_line_zones(self) -> None:
        fit = muse_oiii_fit()
        wavelength = np.arange(4800.0, 5100.1, 1.0)
        rng = np.random.default_rng(9)
        flux = 5.0 + rng.normal(0.0, 0.2, wavelength.size)
        flux[np.argmin(abs(wavelength - 5006.84))] += 100.0
        source = Spectrum(
            "real",
            wavelength,
            flux,
            np.full(wavelength.shape, 0.2),
            redshift=0.0,
        )
        out_wave, baseline, uncertainty, details = real_residual_baseline(
            source, fit, seed=17
        )
        expected_size = int(np.sum((wavelength >= 4920.0) & (wavelength <= 5050.0)))
        self.assertEqual(out_wave.size, expected_size)
        self.assertEqual(baseline.size, expected_size)
        self.assertEqual(uncertainty.size, expected_size)
        self.assertLess(np.max(baseline), 20.0)
        self.assertGreater(details["available_blocks"], 0)

    def test_real_residual_baseline_rejects_mismatched_noise_donor(self) -> None:
        fit = muse_oiii_fit()
        wavelength = np.arange(4800.0, 5100.1, 1.0)
        rng = np.random.default_rng(11)
        source = Spectrum(
            "bad_donor",
            wavelength,
            rng.normal(0.0, 5.0, wavelength.size),
            np.full(wavelength.shape, 0.1),
            redshift=0.0,
        )
        with self.assertRaisesRegex(ValueError, "residual donor is inconsistent"):
            real_residual_baseline(source, fit, seed=3)

    def test_real_residual_baseline_can_calibrate_underestimated_uncertainty(self) -> None:
        fit = muse_oiii_fit()
        wavelength = np.arange(4800.0, 5100.1, 1.0)
        rng = np.random.default_rng(12)
        source = Spectrum(
            "underestimated_uncertainty",
            wavelength,
            rng.normal(0.0, 1.0, wavelength.size),
            np.full(wavelength.shape, 0.5),
            redshift=0.0,
        )
        _, _, uncertainty, details = real_residual_baseline(
            source,
            fit,
            seed=4,
            max_residual_to_uncertainty=3.0,
            calibrate_uncertainty=True,
        )
        self.assertTrue(details["uncertainty_calibration_applied"])
        self.assertGreater(details["uncertainty_calibration_factor"], 1.0)
        self.assertAlmostEqual(
            details[
                "resampled_residual_to_uncertainty_ratio_after_calibration"
            ],
            1.0,
        )
        self.assertAlmostEqual(
            np.median(uncertainty),
            details["resampled_residual_robust_sigma"],
        )

    def test_score_recovery_accepts_exact_fake_recovery(self) -> None:
        fit = muse_oiii_fit(max_components=1)
        spectrum, truth = analytic_injection_spectrum(
            "score", 0.01, [InjectedComponent(0.0, 80.0, 10.0)], fit, seed=5
        )
        result = fit_spectrum(spectrum, fit, sampler_factory=FakeSampler)
        score = score_recovery(result, truth)
        self.assertEqual(score["true_components"], 1)
        self.assertEqual(score["recovered_components"], 1)

    def test_score_recovery_accepts_blank(self) -> None:
        result = {
            "selected_components": 0,
            "selected_logz": 10.0,
            "components": [],
            "models": [
                {
                    "n_components": 0,
                    "logz": 10.0,
                    "parameter_names": [],
                    "posterior_median": [],
                    "posterior_stdev": [],
                    "ncall": 10,
                }
            ],
        }
        truth = {
            "spectrum_id": "blank",
            "stage": "controlled",
            "seed": 1,
            "true_components": 0,
            "components": [],
        }
        score = score_recovery(result, truth)
        self.assertTrue(score["component_count_correct"])
        self.assertEqual(score["component_matching"], "not applicable")

    def test_score_recovery_matches_three_components_independent_of_list_order(self) -> None:
        fit = muse_oiii_fit(max_components=3)
        _, truth = analytic_injection_spectrum(
            "triple",
            0.01,
            [
                InjectedComponent(-250.0, 70.0, 15.0),
                InjectedComponent(0.0, 100.0, 12.0),
                InjectedComponent(250.0, 80.0, 10.0),
            ],
            fit,
            seed=8,
        )
        components = []
        names = []
        medians = []
        stdevs = []
        for component_number, injected in enumerate(truth["components"], start=1):
            names.extend(
                [
                    f"component.{component_number}.velocity_kms",
                    f"component.{component_number}.sigma_kms",
                ]
            )
            medians.extend([injected["velocity_kms"], injected["sigma_kms"]])
            stdevs.extend([5.0, 5.0])
            components.append(
                {
                    "component": component_number,
                    "velocity_kms": injected["velocity_kms"],
                    "sigma_kms": injected["sigma_kms"],
                    "lines": {
                        "oiii5007": {
                            "flux": injected["lines"]["oiii5007"]["flux"],
                            "flux_stdev": injected["lines"]["oiii5007"]["flux"] * 0.05,
                        }
                    },
                }
            )
        result = {
            "selected_components": 3,
            "selected_logz": 100.0,
            "components": list(reversed(components)),
            "models": [
                {
                    "n_components": 3,
                    "logz": 100.0,
                    "parameter_names": names,
                    "posterior_median": medians,
                    "posterior_stdev": stdevs,
                    "ncall": 10,
                }
            ],
        }
        score = score_recovery(result, truth)
        self.assertTrue(score["component_count_correct"])
        self.assertEqual(score["component_matching"], "minimum velocity-width cost assignment")
        self.assertTrue(
            all(item["velocity_error_kms"] == 0.0 for item in score["parameters"])
        )


class ConfigTests(unittest.TestCase):
    def test_json_config_resolves_relative_paths(self) -> None:
        raw = {
            "version": 2,
            "input": {
                "kind": "cube",
                "path": "cube.fits",
                "flux_hdu": 1,
                "redshift": 0.01,
            },
            "fit": {
                "window": [4800, 5100],
                "lines": [{"name": "hbeta", "wavelength": 4861.33}],
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "config.json"
            path.write_text(json.dumps(raw), encoding="utf-8")
            config = load_config(path)
            self.assertEqual(
                config["input"]["path"], str((Path(directory) / "cube.fits").resolve())
            )
            self.assertEqual(config["fit"]["kinematics"]["max_components"], 3)
            self.assertEqual(config["fit"]["selection"]["audit"]["mode"], "flag")

    def test_invalid_selection_audit_mode_is_rejected(self) -> None:
        config = {
            "version": 2,
            "input": {"kind": "cube", "path": "x", "flux_hdu": 1, "redshift": 0},
            "fit": fit_config(),
            "output": {"directory": "x", "plots": "none"},
        }
        config["fit"]["selection"]["audit"] = {"mode": "sometimes"}
        with self.assertRaisesRegex(ConfigError, "audit.mode"):
            validate_config(config)

    def test_locked_line_must_follow_reference(self) -> None:
        config = {
            "version": 2,
            "input": {"kind": "cube", "path": "x", "flux_hdu": 1, "redshift": 0},
            "fit": fit_config(),
            "output": {"directory": "x", "plots": "none"},
        }
        config["fit"]["lines"][0] = {
            "name": "bad",
            "wavelength": 4861,
            "ratio_to": "later",
            "ratio": 0.5,
        }
        with self.assertRaises(ConfigError):
            validate_config(config)

    def test_muse_adapter_defaults_and_requires_bounded_region(self) -> None:
        raw = {
            "version": 2,
            "input": {
                "adapter": "muse",
                "path": "cube.fits",
                "redshift": 0.01,
                "x_range": [10, 15],
                "y_range": [20, 25],
            },
            "fit": {
                "window": [4800, 5100],
                "lines": [{"name": "hbeta", "wavelength": 4861.33}],
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "muse.json"
            path.write_text(json.dumps(raw), encoding="utf-8")
            config = load_config(path)
            self.assertEqual(config["input"]["flux_hdu"], "DATA")
            self.assertEqual(config["input"]["uncertainty_kind"], "variance")
            self.assertEqual(config["fit"]["wavelength_medium"], "air")
            self.assertEqual(config["fit"]["lsf"]["model"], "instrument")

            del raw["input"]["x_range"]
            path.write_text(json.dumps(raw), encoding="utf-8")
            with self.assertRaisesRegex(ConfigError, "require both input.x_range"):
                load_config(path)

    def test_nirspec_adapter_defaults(self) -> None:
        raw = {
            "version": 2,
            "input": {
                "adapter": "nirspec",
                "path": "cube_s3d.fits",
                "redshift": 0.01,
                "x_range": [1, 3],
                "y_range": [2, 4],
            },
            "fit": {
                "window": [18000, 19000],
                "lines": [{"name": "line", "wavelength": 18750.0}],
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "nirspec.json"
            path.write_text(json.dumps(raw), encoding="utf-8")
            config = load_config(path)
            self.assertEqual(config["input"]["flux_hdu"], "SCI")
            self.assertEqual(config["input"]["mask_hdu"], "DQ")
            self.assertEqual(config["input"]["mask_bits"], 1)
            self.assertEqual(config["fit"]["wavelength_medium"], "vacuum")
            self.assertEqual(config["fit"]["lsf"]["model"], "instrument")

    def test_miri_adapter_accepts_segment_glob(self) -> None:
        raw = {
            "version": 2,
            "input": {
                "adapter": "miri",
                "glob": "target_ch*.fits",
                "redshift": 0.01,
                "x_range": [1, 3],
                "y_range": [2, 4],
            },
            "fit": {
                "window": [141000, 145000],
                "lines": [{"name": "nev", "wavelength": 143200.0}],
            },
        }
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "miri.json"
            path.write_text(json.dumps(raw), encoding="utf-8")
            config = load_config(path)
            self.assertNotIn("path", config["input"])
            self.assertTrue(config["input"]["glob"].endswith("target_ch*.fits"))
            self.assertEqual(config["input"]["mask_bits"], 1)
            self.assertEqual(config["fit"]["wavelength_medium"], "vacuum")


class PipelineTests(unittest.TestCase):
    def test_configured_sampler_seed_is_stable_per_spectrum(self) -> None:
        from beat.pipeline import spectrum_seed

        fit = {"sampling": {"seed": 12345}}
        self.assertEqual(spectrum_seed("target-a", fit), spectrum_seed("target-a", fit))
        self.assertNotEqual(spectrum_seed("target-a", fit), spectrum_seed("target-b", fit))
        self.assertIsNone(spectrum_seed("target-a", {"sampling": {}}))

    def test_process_pool_tolerates_unreadable_macos_semaphore_limit(self) -> None:
        import os
        from unittest.mock import patch

        from beat.pipeline import _process_executor

        original = os.sysconf

        def denied(name):
            if name == "SC_SEM_NSEMS_MAX":
                raise PermissionError(1, "not permitted")
            return original(name)

        with patch("os.sysconf", side_effect=denied):
            executor = _process_executor(1)
        executor.shutdown(wait=True)

    def test_atomic_checkpoint_catalog_and_resume(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            spectrum_path = root / "target.txt"
            wavelength = np.linspace(4750.0, 5150.0, 100)
            np.savetxt(
                spectrum_path,
                np.column_stack([wavelength, np.zeros(100), np.ones(100)]),
            )
            config = {
                "version": 2,
                "input": {
                    "kind": "spectrum_files",
                    "files": [str(spectrum_path)],
                    "format": "ascii",
                    "wavelength_column": 0,
                    "flux_column": 1,
                    "uncertainty_column": 2,
                    "uncertainty_kind": "sigma",
                    "redshift": 0.0,
                },
                "fit": fit_config(max_components=0),
                "output": {
                    "directory": str(root / "output"),
                    "resume": True,
                    "plots": "none",
                    "progress_every": 10,
                },
            }
            calls = []

            def fake_worker(spectrum, fit, plots, plots_dir):
                calls.append(spectrum.spectrum_id)
                return {
                    "status": "ok",
                    "spectrum_id": spectrum.spectrum_id,
                    "redshift": spectrum.redshift,
                    "metadata": spectrum.metadata,
                    "n_input_pixels": 100,
                    "n_fit_pixels": 74,
                    "noise_level": 1.0,
                    "selected_components": 0,
                    "selected_logz": 0.0,
                    "selection_status": "convergence_unverified",
                    "selection_reasons": ["test review reason"],
                    "selection_diagnostics": [
                        {
                            "reference_components": 0,
                            "candidate_components": 1,
                            "delta_logz": 5.1,
                            "combined_logz_error": 0.8,
                            "threshold": 5.0,
                            "distance_from_threshold": 0.1,
                        }
                    ],
                    "selection_audit": {"performed": False},
                    "components": [],
                    "models": [
                        {
                            "n_components": 0,
                            "logz": 0.0,
                            "logz_error": 0.1,
                            "parameter_names": ["continuum_c0"],
                            "posterior_median": [0.0],
                            "posterior_stdev": [0.1],
                        }
                    ],
                }

            with patch("beat.pipeline._worker", side_effect=fake_worker):
                first = run_pipeline(config, workers=1)
                second = run_pipeline(config, workers=1)
            self.assertEqual(first["completed_this_run"], 1)
            self.assertEqual(
                first["selection_status_counts"], {"convergence_unverified": 1}
            )
            self.assertEqual(first["selection_review_rows"], 1)
            self.assertEqual(second["resumed"], 1)
            self.assertEqual(calls, ["target"])
            self.assertTrue((root / "output" / "catalog.csv").exists())
            review = (root / "output" / "selection_review.csv").read_text(
                encoding="utf-8"
            )
            self.assertIn("test review reason", review)
            self.assertEqual(
                len(list((root / "output" / "results").rglob("*.json"))), 1
            )


if __name__ == "__main__":
    unittest.main()
