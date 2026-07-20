"""Injection/recovery utilities for scientific validation of BEAT.

The functions in this module keep the injected truth independent of the
sampler. They support analytic Gaussian noise and block-resampled residuals
from real spectra, while using the same documented Gaussian LSF convention as
the fitting model.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import permutations
from typing import Any, Iterable

import numpy as np

from .fitting import selected_model
from .model import (
    SQRT_2PI,
    convolved_sigma,
    gaussian_integrated,
    observed_center,
    robust_sigma,
)
from .spectrum import Spectrum


@dataclass(frozen=True)
class InjectedComponent:
    """Truth parameters for one kinematic component."""

    velocity_kms: float
    sigma_kms: float
    peak_snr: float
    flux_scale: float = 1.0
    line_flux_scales: dict[str, float] = field(default_factory=dict)


@dataclass(frozen=True)
class InjectedBroadComponent:
    """Truth parameters for one permitted-line broad component."""

    name: str
    line: str
    velocity_kms: float
    sigma_kms: float
    peak_snr: float


def muse_oiii_fit(
    max_components: int = 2,
    min_num_live_points: int = 80,
    min_ess: int = 120,
    dlogz: float = 1.0,
) -> dict[str, Any]:
    """Return the MUSE [O III] doublet validation model."""
    return {
        "frame": "rest",
        "wavelength_medium": "air",
        "window": [4920.0, 5050.0],
        "minimum_valid_pixels": 50,
        "continuum": {
            "degree": 1,
            "windows": [[4920.0, 4940.0], [5030.0, 5050.0]],
            "prior_width_rms": 10.0,
        },
        "kinematics": {
            "max_components": int(max_components),
            "velocity_kms": [-600.0, 600.0],
            "sigma_kms": [30.0, 500.0],
        },
        "lines": [
            {"name": "oiii5007", "wavelength": 5006.84},
            {
                "name": "oiii4959",
                "wavelength": 4958.92,
                "ratio_to": "oiii5007",
                "ratio": 0.33557,
            },
        ],
        "lsf": {"model": "instrument"},
        "flux_prior": {"min_snr": 0.05, "max_signal_factor": 30.0},
        "selection": {"delta_logz": 5.0, "stop_when_not_preferred": True},
        "sampling": {
            "min_num_live_points": int(min_num_live_points),
            "min_ess": int(min_ess),
            "dlogz": float(dlogz),
            "show_status": False,
            "stepsampler": "slice",
            "nsteps": 10,
        },
    }


def muse_hbeta_oiii_fit(
    max_components: int = 2,
    min_num_live_points: int = 80,
    min_ess: int = 120,
    dlogz: float = 1.0,
) -> dict[str, Any]:
    """Return the MUSE H-beta+[O III] validation model."""
    fit = muse_oiii_fit(
        max_components=max_components,
        min_num_live_points=min_num_live_points,
        min_ess=min_ess,
        dlogz=dlogz,
    )
    fit["window"] = [4800.0, 5100.0]
    fit["minimum_valid_pixels"] = 100
    fit["continuum"]["windows"] = [[4800.0, 4830.0], [5070.0, 5100.0]]
    fit["lines"] = [
        {"name": "hbeta", "wavelength": 4861.33},
        {"name": "oiii5007", "wavelength": 5006.84},
        {
            "name": "oiii4959",
            "wavelength": 4958.92,
            "ratio_to": "oiii5007",
            "ratio": 0.33557,
        },
    ]
    return fit


def muse_halpha_nii_fit(
    max_components: int = 2,
    include_broad_halpha: bool = False,
    min_num_live_points: int = 80,
    min_ess: int = 120,
    dlogz: float = 1.0,
) -> dict[str, Any]:
    """Return the MUSE H-alpha+[N II] validation model."""
    fit: dict[str, Any] = {
        "frame": "rest",
        "wavelength_medium": "air",
        "window": [6480.0, 6655.0],
        "minimum_valid_pixels": 70,
        "continuum": {
            "degree": 1,
            "windows": [[6480.0, 6515.0], [6625.0, 6655.0]],
            "prior_width_rms": 10.0,
        },
        "kinematics": {
            "max_components": int(max_components),
            "velocity_kms": [-600.0, 600.0],
            "sigma_kms": [30.0, 500.0],
        },
        "lines": [
            {"name": "halpha", "wavelength": 6562.80},
            {"name": "nii6583", "wavelength": 6583.45},
            {
                "name": "nii6548",
                "wavelength": 6548.05,
                "ratio_to": "nii6583",
                "ratio": 0.335,
            },
        ],
        "lsf": {"model": "instrument"},
        "flux_prior": {"min_snr": 0.05, "max_signal_factor": 40.0},
        "selection": {"delta_logz": 5.0, "stop_when_not_preferred": True},
        "sampling": {
            "min_num_live_points": int(min_num_live_points),
            "min_ess": int(min_ess),
            "dlogz": float(dlogz),
            "show_status": False,
            "stepsampler": "slice",
            "nsteps": 12,
        },
    }
    if include_broad_halpha:
        fit["broad_components"] = [
            {
                "name": "broad_halpha",
                "line": "halpha",
                "velocity_kms": [-1500.0, 1500.0],
                "sigma_kms": [600.0, 3500.0],
            }
        ]
    return fit


def muse_metadata(resolving_power: float = 3027.0) -> dict[str, Any]:
    """Minimal MUSE metadata needed by the instrument LSF model."""
    return {
        "adapter": "muse",
        "instrument": "MUSE",
        "wavelength_medium": "air",
        "instrument_lsf": {
            "model": "resolving_power",
            "value": float(resolving_power),
            "source": "injection/recovery MUSE approximation",
            "approximation": True,
        },
    }


def muse_wavelength_grid(
    redshift: float,
    step_angstrom: float = 1.25,
    window: tuple[float, float] | list[float] = (4920.0, 5050.0),
) -> np.ndarray:
    """Create a MUSE-like observed wavelength grid covering the fit window."""
    lo, hi = np.asarray(window, dtype=float) * (1.0 + float(redshift))
    start = np.floor(lo / step_angstrom) * step_angstrom
    stop = np.ceil(hi / step_angstrom) * step_angstrom
    return np.arange(start, stop + 0.5 * step_angstrom, step_angstrom)


def inject_emission_lines(
    wavelength: np.ndarray,
    baseline_flux: np.ndarray,
    uncertainty: np.ndarray,
    redshift: float,
    components: Iterable[InjectedComponent],
    fit: dict[str, Any],
    metadata: dict[str, Any],
    reference_line: str | None = None,
    independent_line_ratios: dict[str, float] | None = None,
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Add line components and return the injected flux and serialisable truth."""
    wavelength = np.asarray(wavelength, dtype=float)
    flux = np.asarray(baseline_flux, dtype=float).copy()
    uncertainty = np.asarray(uncertainty, dtype=float)
    truth: list[dict[str, Any]] = []
    lines = {line["name"]: line for line in fit["lines"]}
    ordered_lines = list(fit["lines"])
    if reference_line is None:
        reference_line = "oiii5007" if "oiii5007" in lines else ordered_lines[0]["name"]
    if reference_line not in lines or "ratio_to" in lines[reference_line]:
        raise ValueError("reference_line must name a free line in the fitted complex")
    ratios = dict(independent_line_ratios or {})
    if "hbeta" in lines and "hbeta" not in ratios and reference_line == "oiii5007":
        ratios["hbeta"] = 0.25
    ratios.setdefault(reference_line, 1.0)
    components = sorted(components, key=lambda item: item.velocity_kms)
    for component_number, component in enumerate(components, start=1):
        reference = lines[reference_line]
        reference_center = observed_center(
            float(reference["wavelength"]), redshift, component.velocity_kms
        )
        local_noise = float(np.interp(reference_center, wavelength, uncertainty))
        reference_sigma = convolved_sigma(
            float(reference["wavelength"]),
            redshift,
            component.sigma_kms,
            component.velocity_kms,
            fit.get("lsf"),
            metadata,
        )
        reference_flux = (
            component.peak_snr
            * component.flux_scale
            * local_noise
            * reference_sigma
            * SQRT_2PI
        )
        line_fluxes: dict[str, float] = {}
        for line in ordered_lines:
            name = line["name"]
            if "ratio_to" in line:
                line_fluxes[name] = line_fluxes[line["ratio_to"]] * float(line["ratio"])
            else:
                scale = component.line_flux_scales.get(name, ratios.get(name, 1.0))
                line_fluxes[name] = reference_flux * float(scale)
        line_truth: dict[str, Any] = {}
        for name in line_fluxes:
            rest = float(lines[name]["wavelength"])
            center = observed_center(rest, redshift, component.velocity_kms)
            sigma = convolved_sigma(
                rest,
                redshift,
                component.sigma_kms,
                component.velocity_kms,
                fit.get("lsf"),
                metadata,
            )
            flux += gaussian_integrated(wavelength, center, sigma, line_fluxes[name])
            line_truth[name] = {
                "flux": float(line_fluxes[name]),
                "observed_center_angstrom": float(center),
                "convolved_sigma_angstrom": float(sigma),
            }
        truth.append(
            {
                "component": component_number,
                "velocity_kms": float(component.velocity_kms),
                "sigma_kms": float(component.sigma_kms),
                "peak_snr": float(component.peak_snr * component.flux_scale),
                "lines": line_truth,
            }
        )
    return flux, truth


def inject_broad_lines(
    wavelength: np.ndarray,
    flux: np.ndarray,
    uncertainty: np.ndarray,
    redshift: float,
    components: Iterable[InjectedBroadComponent],
    fit: dict[str, Any],
    metadata: dict[str, Any],
) -> tuple[np.ndarray, list[dict[str, Any]]]:
    """Add independent broad permitted-line components."""
    output = np.asarray(flux, dtype=float).copy()
    lines = {line["name"]: line for line in fit["lines"]}
    truth: list[dict[str, Any]] = []
    for component in components:
        line = lines[component.line]
        rest = float(line["wavelength"])
        center = observed_center(rest, redshift, component.velocity_kms)
        sigma = convolved_sigma(
            rest,
            redshift,
            component.sigma_kms,
            component.velocity_kms,
            fit.get("lsf"),
            metadata,
        )
        local_noise = float(np.interp(center, wavelength, uncertainty))
        integrated_flux = component.peak_snr * local_noise * sigma * SQRT_2PI
        output += gaussian_integrated(wavelength, center, sigma, integrated_flux)
        truth.append(
            {
                "name": component.name,
                "line": component.line,
                "velocity_kms": float(component.velocity_kms),
                "sigma_kms": float(component.sigma_kms),
                "peak_snr": float(component.peak_snr),
                "flux": float(integrated_flux),
                "observed_center_angstrom": float(center),
                "convolved_sigma_angstrom": float(sigma),
            }
        )
    return output, truth


def analytic_injection_spectrum(
    spectrum_id: str,
    redshift: float,
    components: Iterable[InjectedComponent],
    fit: dict[str, Any],
    seed: int,
    noise_sigma: float = 1.0,
    continuum_level: float = 5.0,
    continuum_slope: float = 0.002,
    resolving_power: float = 3027.0,
    reference_line: str | None = None,
    independent_line_ratios: dict[str, float] | None = None,
    broad_components: Iterable[InjectedBroadComponent] = (),
) -> tuple[Spectrum, dict[str, Any]]:
    """Generate one deterministic MUSE-like spectrum with Gaussian noise."""
    rng = np.random.default_rng(int(seed))
    wavelength = muse_wavelength_grid(redshift, window=fit["window"])
    midpoint = float(np.mean(wavelength))
    continuum = continuum_level + continuum_slope * (wavelength - midpoint)
    uncertainty = np.full(wavelength.shape, float(noise_sigma))
    noisy_baseline = continuum + rng.normal(0.0, uncertainty)
    metadata = muse_metadata(resolving_power)
    metadata.update({"noise_source": "analytic_gaussian", "injection_seed": int(seed)})
    injected_flux, component_truth = inject_emission_lines(
        wavelength,
        noisy_baseline,
        uncertainty,
        redshift,
        components,
        fit,
        metadata,
        reference_line=reference_line,
        independent_line_ratios=independent_line_ratios,
    )
    injected_flux, broad_truth = inject_broad_lines(
        wavelength,
        injected_flux,
        uncertainty,
        redshift,
        broad_components,
        fit,
        metadata,
    )
    resolved_reference = reference_line or (
        "oiii5007" if any(line["name"] == "oiii5007" for line in fit["lines"])
        else fit["lines"][0]["name"]
    )
    truth = {
        "spectrum_id": spectrum_id,
        "stage": "controlled",
        "seed": int(seed),
        "redshift": float(redshift),
        "true_components": len(component_truth),
        "components": component_truth,
        "broad_components": broad_truth,
        "reference_line": resolved_reference,
    }
    spectrum = Spectrum(
        spectrum_id=spectrum_id,
        wavelength=wavelength,
        flux=injected_flux,
        uncertainty=uncertainty,
        redshift=redshift,
        metadata=metadata,
    )
    return spectrum, truth


def real_residual_baseline(
    source: Spectrum,
    fit: dict[str, Any],
    seed: int,
    block_length: int = 8,
    max_residual_to_uncertainty: float = 3.0,
    calibrate_uncertainty: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Create a baseline from block-resampled, line-free real residuals.

    Known line regions and strong outliers are excluded from the block pool.
    Residual and uncertainty samples are moved together, preserving short-range
    residual covariance and heteroscedasticity without retaining native lines.

    Validation runners may add rest wavelengths under
    ``fit.injection_residual_mask_lines`` and may set
    ``fit.injection_residual_mask_half_width_angstrom``. These settings affect
    only construction of the injected baseline, not the fitted line model.
    """
    prepared = source.prepared()
    redshift = prepared.redshift
    lo, hi = np.asarray(fit["window"], dtype=float) * (1.0 + redshift)
    selected = (prepared.wavelength >= lo) & (prepared.wavelength <= hi)
    wavelength = prepared.wavelength[selected]
    flux = prepared.flux[selected]
    if prepared.uncertainty is None:
        raise ValueError("Real-residual injections require supplied uncertainties")
    uncertainty = prepared.uncertainty[selected]
    if wavelength.size < int(fit.get("minimum_valid_pixels", 100)):
        raise ValueError("Real spectrum contains too few pixels in validation window")

    line_free = np.ones(wavelength.size, dtype=bool)
    mask_half_width = float(
        fit.get("injection_residual_mask_half_width_angstrom", 18.0)
    )
    if not np.isfinite(mask_half_width) or mask_half_width <= 0:
        raise ValueError("Injection residual line-mask half-width must be positive")
    mask_wavelengths = [float(line["wavelength"]) for line in fit["lines"]]
    mask_wavelengths.extend(
        float(value) for value in fit.get("injection_residual_mask_lines", [])
    )
    for rest_wavelength in mask_wavelengths:
        center = rest_wavelength * (1.0 + redshift)
        line_free &= np.abs(wavelength - center) > mask_half_width
    x = (wavelength - np.mean(wavelength)) / (0.5 * np.ptp(wavelength))
    extraction_degree = 2
    coefficients = np.polynomial.polynomial.polyfit(
        x[line_free], flux[line_free], extraction_degree
    )
    for _ in range(5):
        continuum = np.polynomial.polynomial.polyval(x, coefficients)
        residual = flux - continuum
        scale = robust_sigma(residual[line_free])
        clean = line_free & np.isfinite(residual)
        if np.isfinite(scale) and scale > 0:
            clean &= np.abs(residual - np.median(residual[line_free])) < 5.0 * scale
        coefficients = np.polynomial.polynomial.polyfit(
            x[clean], flux[clean], extraction_degree
        )
    continuum = np.polynomial.polynomial.polyval(x, coefficients)
    residual = flux - continuum
    scale = robust_sigma(residual[clean])
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("Could not estimate real residual scale")

    uncertainty_check = clean & np.isfinite(uncertainty) & (uncertainty > 0)
    typical_uncertainty = float(np.median(uncertainty[uncertainty_check]))
    residual_to_uncertainty = float(scale / typical_uncertainty)
    if (
        not np.isfinite(residual_to_uncertainty)
        or residual_to_uncertainty > float(max_residual_to_uncertainty)
    ):
        raise ValueError(
            "Real residual donor is inconsistent with its supplied uncertainty: "
            f"robust residual/median uncertainty={residual_to_uncertainty:.2f}, "
            f"allowed maximum={float(max_residual_to_uncertainty):.2f}. "
            "Select a lower-line-emission donor or widen the native-line masks."
        )

    valid = clean & np.isfinite(uncertainty) & (uncertainty > 0)
    starts = [
        index
        for index in range(0, wavelength.size - block_length + 1)
        if np.all(valid[index : index + block_length])
    ]
    if not starts:
        raise ValueError("No contiguous line-free residual blocks are available")
    rng = np.random.default_rng(int(seed))
    sampled_residual: list[float] = []
    sampled_uncertainty: list[float] = []
    while len(sampled_residual) < wavelength.size:
        start = int(rng.choice(starts))
        sampled_residual.extend(residual[start : start + block_length].tolist())
        sampled_uncertainty.extend(uncertainty[start : start + block_length].tolist())
    sampled_residual_array = np.asarray(sampled_residual[: wavelength.size])
    sampled_uncertainty_array = np.asarray(sampled_uncertainty[: wavelength.size])
    resampled_residual_scale = float(robust_sigma(sampled_residual_array))
    resampled_median_uncertainty = float(np.median(sampled_uncertainty_array))
    resampled_residual_to_uncertainty = (
        resampled_residual_scale / resampled_median_uncertainty
    )
    uncertainty_calibration_factor = (
        max(1.0, resampled_residual_to_uncertainty)
        if calibrate_uncertainty
        else 1.0
    )
    sampled_uncertainty_array = (
        sampled_uncertainty_array * uncertainty_calibration_factor
    )
    injected_degree = int(fit.get("continuum", {}).get("degree", 1))
    injected_coefficients = np.polynomial.polynomial.polyfit(
        x[clean], flux[clean], injected_degree
    )
    injected_continuum = np.polynomial.polynomial.polyval(x, injected_coefficients)
    baseline = injected_continuum + sampled_residual_array
    details = {
        "noise_source": "block_resampled_real_residual",
        "source_spectrum_id": source.spectrum_id,
        "block_length": int(block_length),
        "available_blocks": len(starts),
        "source_residual_robust_sigma": float(scale),
        "source_median_uncertainty": typical_uncertainty,
        "source_residual_to_uncertainty_ratio": residual_to_uncertainty,
        "masked_line_wavelengths_rest_angstrom": mask_wavelengths,
        "line_mask_half_width_angstrom": mask_half_width,
        "maximum_allowed_residual_to_uncertainty_ratio": float(
            max_residual_to_uncertainty
        ),
        "resampled_residual_robust_sigma": resampled_residual_scale,
        "resampled_median_uncertainty_before_calibration": (
            resampled_median_uncertainty
        ),
        "resampled_residual_to_uncertainty_ratio_before_calibration": (
            resampled_residual_to_uncertainty
        ),
        "uncertainty_calibration_applied": bool(calibrate_uncertainty),
        "uncertainty_calibration_factor": float(uncertainty_calibration_factor),
        "resampled_residual_to_uncertainty_ratio_after_calibration": float(
            resampled_residual_to_uncertainty / uncertainty_calibration_factor
        ),
        "residual_extraction_continuum_degree": extraction_degree,
        "injected_continuum_degree": injected_degree,
    }
    return wavelength, baseline, sampled_uncertainty_array, details


def real_noise_injection_spectrum(
    spectrum_id: str,
    source: Spectrum,
    components: Iterable[InjectedComponent],
    fit: dict[str, Any],
    seed: int,
    reference_line: str | None = None,
    independent_line_ratios: dict[str, float] | None = None,
    broad_components: Iterable[InjectedBroadComponent] = (),
    max_residual_to_uncertainty: float = 3.0,
    calibrate_uncertainty: bool = False,
) -> tuple[Spectrum, dict[str, Any]]:
    """Inject known lines into a block-resampled real MUSE residual baseline."""
    wavelength, baseline, uncertainty, details = real_residual_baseline(
        source,
        fit,
        seed,
        max_residual_to_uncertainty=max_residual_to_uncertainty,
        calibrate_uncertainty=calibrate_uncertainty,
    )
    metadata = dict(source.metadata)
    metadata.update(details)
    metadata["injection_seed"] = int(seed)
    injected_flux, component_truth = inject_emission_lines(
        wavelength,
        baseline,
        uncertainty,
        source.redshift,
        components,
        fit,
        metadata,
        reference_line=reference_line,
        independent_line_ratios=independent_line_ratios,
    )
    injected_flux, broad_truth = inject_broad_lines(
        wavelength,
        injected_flux,
        uncertainty,
        source.redshift,
        broad_components,
        fit,
        metadata,
    )
    resolved_reference = reference_line or (
        "oiii5007" if any(line["name"] == "oiii5007" for line in fit["lines"])
        else fit["lines"][0]["name"]
    )
    truth = {
        "spectrum_id": spectrum_id,
        "stage": "real_noise",
        "seed": int(seed),
        "redshift": float(source.redshift),
        "true_components": len(component_truth),
        "components": component_truth,
        "broad_components": broad_truth,
        "reference_line": resolved_reference,
        "noise_details": details,
    }
    return (
        Spectrum(
            spectrum_id=spectrum_id,
            wavelength=wavelength,
            flux=injected_flux,
            uncertainty=uncertainty,
            redshift=source.redshift,
            metadata=metadata,
        ),
        truth,
    )


def score_recovery(result: dict[str, Any], truth: dict[str, Any]) -> dict[str, Any]:
    """Compare a selected BEAT result with injected truth."""
    selected = selected_model(result)
    medians = dict(zip(selected["parameter_names"], selected["posterior_median"]))
    stdevs = dict(zip(selected["parameter_names"], selected["posterior_stdev"]))
    true_count = int(truth["true_components"])
    recovered_count = int(result["selected_components"])
    score: dict[str, Any] = {
        "spectrum_id": truth["spectrum_id"],
        "stage": truth["stage"],
        "seed": truth["seed"],
        "true_components": true_count,
        "recovered_components": recovered_count,
        "component_count_correct": recovered_count == true_count,
        "selected_logz": float(result["selected_logz"]),
        "total_ncall": int(sum(model.get("ncall", 0) for model in result["models"])),
        "peak_snr_primary": (
            None if not truth["components"] else truth["components"][0]["peak_snr"]
        ),
        "true_separation_kms": (
            None
            if true_count < 2
            else truth["components"][-1]["velocity_kms"]
            - truth["components"][0]["velocity_kms"]
        ),
        "parameters": [],
    }
    logz_by_count = {int(model["n_components"]): float(model["logz"]) for model in result["models"]}
    score["logz_by_components"] = logz_by_count
    if recovered_count != true_count:
        return score
    if true_count == 0:
        score["component_matching"] = "not applicable"
        return score
    recovered = list(result.get("components", []))
    reference_line = truth.get("reference_line", "oiii5007")
    injected_components = list(truth["components"])
    if true_count <= 8:
        velocity_scale = max(
            100.0,
            float(np.ptp([item["velocity_kms"] for item in injected_components])),
        )

        def assignment_cost(order: tuple[int, ...]) -> float:
            cost = 0.0
            for injected, recovered_index in zip(injected_components, order):
                fitted = recovered[recovered_index]
                cost += abs(
                    float(fitted["velocity_kms"]) - float(injected["velocity_kms"])
                ) / velocity_scale
                cost += 0.25 * abs(
                    np.log(float(fitted["sigma_kms"]) / float(injected["sigma_kms"]))
                )
            return float(cost)

        best_order = min(permutations(range(true_count)), key=assignment_cost)
        matched = [
            (injected, recovered[recovered_index])
            for injected, recovered_index in zip(injected_components, best_order)
        ]
        matching_method = "minimum velocity-width cost assignment"
    else:
        matched = list(
            zip(
                sorted(injected_components, key=lambda item: item["velocity_kms"]),
                sorted(recovered, key=lambda item: item["velocity_kms"]),
            )
        )
        matching_method = "velocity-order assignment"
    score["component_matching"] = matching_method
    for index, (injected, fitted) in enumerate(matched, start=1):
        recovered_component = int(fitted.get("component", index))
        prefix = f"component.{recovered_component}"
        velocity_stdev = float(stdevs.get(f"{prefix}.velocity_kms", np.nan))
        sigma_stdev = float(stdevs.get(f"{prefix}.sigma_kms", np.nan))
        truth_flux = float(injected["lines"][reference_line]["flux"])
        fitted_flux = float(fitted["lines"][reference_line]["flux"])
        flux_stdev = float(fitted["lines"][reference_line]["flux_stdev"])
        velocity_error = float(fitted["velocity_kms"] - injected["velocity_kms"])
        sigma_error = float(fitted["sigma_kms"] - injected["sigma_kms"])
        flux_fractional_error = float((fitted_flux - truth_flux) / truth_flux)
        score["parameters"].append(
            {
                "component": index,
                "recovered_component": recovered_component,
                "true_velocity_kms": float(injected["velocity_kms"]),
                "recovered_velocity_kms": float(fitted["velocity_kms"]),
                "velocity_error_kms": velocity_error,
                "velocity_posterior_stdev_kms": velocity_stdev,
                "velocity_within_1sigma": bool(
                    np.isfinite(velocity_stdev) and abs(velocity_error) <= velocity_stdev
                ),
                "true_sigma_kms": float(injected["sigma_kms"]),
                "recovered_sigma_kms": float(fitted["sigma_kms"]),
                "sigma_error_kms": sigma_error,
                "sigma_fractional_error": sigma_error / float(injected["sigma_kms"]),
                "sigma_posterior_stdev_kms": sigma_stdev,
                "sigma_within_1sigma": bool(
                    np.isfinite(sigma_stdev) and abs(sigma_error) <= sigma_stdev
                ),
                "true_oiii5007_flux": truth_flux,
                "recovered_oiii5007_flux": fitted_flux,
                "reference_line": reference_line,
                "true_reference_flux": truth_flux,
                "recovered_reference_flux": fitted_flux,
                "flux_fractional_error": flux_fractional_error,
                "flux_posterior_stdev": flux_stdev,
                "flux_within_1sigma": bool(
                    np.isfinite(flux_stdev)
                    and abs(fitted_flux - truth_flux) <= flux_stdev
                ),
            }
        )
    return score
