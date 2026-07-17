"""Dataset-independent emission-line model and prior transforms."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .config import ConfigError
from .lsf import lsf_sigma_angstrom
from .spectrum import Spectrum


C_KMS = 299_792.458
SQRT_2PI = float(np.sqrt(2.0 * np.pi))


def robust_sigma(values: np.ndarray) -> float:
    """Return 1.4826 times the median absolute deviation."""
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return float("nan")
    median = np.median(values)
    return float(1.4826 * np.median(np.abs(values - median)))


def relativistic_doppler(velocity_kms: float | np.ndarray) -> float | np.ndarray:
    """Convert line-of-sight velocity to the relativistic wavelength factor."""
    beta = np.asarray(velocity_kms) / C_KMS
    if np.any(np.abs(beta) >= 1):
        raise ValueError("Velocity must have magnitude below the speed of light")
    factor = np.sqrt((1.0 + beta) / (1.0 - beta))
    return float(factor) if factor.ndim == 0 else factor


def observed_interval(interval: list[float], frame: str, redshift: float) -> tuple[float, float]:
    lo, hi = float(interval[0]), float(interval[1])
    if frame == "rest":
        return lo * (1.0 + redshift), hi * (1.0 + redshift)
    return lo, hi


def observed_center(rest_wavelength: float, redshift: float, velocity_kms: float) -> float:
    return float(rest_wavelength * (1.0 + redshift) * relativistic_doppler(velocity_kms))


def observed_sigma(
    rest_wavelength: float,
    redshift: float,
    sigma_kms: float,
    velocity_kms: float = 0.0,
) -> float:
    """Convert Gaussian velocity sigma to wavelength sigma in Angstrom."""
    center = observed_center(rest_wavelength, redshift, velocity_kms)
    return float(center * sigma_kms / C_KMS)


def convolved_sigma(
    rest_wavelength: float,
    redshift: float,
    sigma_kms: float,
    velocity_kms: float,
    lsf_config: dict[str, Any] | None,
    metadata: dict[str, Any] | None,
) -> float:
    """Combine intrinsic and Gaussian instrumental sigmas in quadrature."""
    center = observed_center(rest_wavelength, redshift, velocity_kms)
    intrinsic = observed_sigma(rest_wavelength, redshift, sigma_kms, velocity_kms)
    instrumental = lsf_sigma_angstrom(center, lsf_config, metadata)
    return float(np.hypot(intrinsic, instrumental))


def gaussian_integrated(
    wavelength: np.ndarray, center: float, sigma: float, integrated_flux: float
) -> np.ndarray:
    """Evaluate a Gaussian normalized to integrated flux."""
    if sigma <= 0:
        raise ValueError("Gaussian sigma must be positive")
    amplitude = integrated_flux / (sigma * SQRT_2PI)
    return amplitude * np.exp(-0.5 * ((wavelength - center) / sigma) ** 2)


def _sigma_clip_polyfit(
    x: np.ndarray,
    y: np.ndarray,
    degree: int,
    iterations: int = 5,
) -> tuple[np.ndarray, float]:
    valid = np.isfinite(x) & np.isfinite(y)
    coefficients = np.zeros(degree + 1)
    for _ in range(iterations):
        if valid.sum() < degree + 2:
            break
        coefficients = np.polynomial.polynomial.polyfit(x[valid], y[valid], degree)
        residual = y - np.polynomial.polynomial.polyval(x, coefficients)
        scatter = robust_sigma(residual[valid])
        if not np.isfinite(scatter) or scatter <= 0:
            scatter = float(np.std(residual[valid]))
        if not np.isfinite(scatter) or scatter <= 0:
            break
        new_valid = valid & (np.abs(residual) <= 4.0 * scatter)
        if np.array_equal(valid, new_valid):
            break
        valid = new_valid
    residual = y[valid] - np.polynomial.polynomial.polyval(x[valid], coefficients)
    scatter = robust_sigma(residual)
    if not np.isfinite(scatter) or scatter <= 0:
        scatter = float(np.std(residual))
    return coefficients, scatter


@dataclass
class PreparedSpectrum:
    spectrum: Spectrum
    wavelength: np.ndarray
    flux: np.ndarray
    uncertainty: np.ndarray
    x_normalized: np.ndarray
    continuum_center: np.ndarray
    continuum_half_width: np.ndarray
    noise_level: float
    noise_model: str
    noise_rho: float
    noise_marginal_scale: float
    positive_peak: float
    default_flux_bounds: tuple[float, float]


def prepare_spectrum(spectrum: Spectrum, fit: dict[str, Any]) -> PreparedSpectrum:
    """Clean and crop a spectrum, estimate noise, and construct continuum priors."""
    spectrum = spectrum.prepared()
    frame = fit.get("frame", "rest")
    lo, hi = observed_interval(fit["window"], frame, spectrum.redshift)
    selected = (spectrum.wavelength >= lo) & (spectrum.wavelength <= hi)
    for excluded in fit.get("exclude_windows", []):
        excluded_lo, excluded_hi = observed_interval(
            excluded, frame, spectrum.redshift
        )
        selected &= ~(
            (spectrum.wavelength >= excluded_lo)
            & (spectrum.wavelength <= excluded_hi)
        )
    minimum = int(fit.get("minimum_valid_pixels", 20))
    if selected.sum() < minimum:
        raise ValueError(
            f"{spectrum.spectrum_id}: fitting window contains {selected.sum()} valid "
            f"pixels; at least {minimum} are required"
        )

    outside = []
    for line in fit["lines"]:
        center = float(line["wavelength"]) * (1.0 + spectrum.redshift)
        if not lo <= center <= hi:
            outside.append(str(line["name"]))
        for excluded in fit.get("exclude_windows", []):
            excluded_lo, excluded_hi = observed_interval(
                excluded, frame, spectrum.redshift
            )
            if excluded_lo <= center <= excluded_hi:
                raise ValueError(
                    f"{spectrum.spectrum_id}: systemic center for {line['name']} "
                    "falls inside fit.exclude_windows"
                )
    if outside:
        raise ValueError(
            f"{spectrum.spectrum_id}: systemic line center(s) outside the fit window: "
            + ", ".join(outside)
        )

    wavelength = spectrum.wavelength[selected]
    flux = spectrum.flux[selected]
    uncertainty = None if spectrum.uncertainty is None else spectrum.uncertainty[selected]
    midpoint = 0.5 * (lo + hi)
    half_range = 0.5 * (hi - lo)
    x_normalized = (wavelength - midpoint) / half_range

    continuum_config = fit.get("continuum", {})
    degree = int(continuum_config.get("degree", 1))
    continuum_mask = np.zeros(wavelength.size, dtype=bool)
    for window in continuum_config.get("windows", []):
        window_lo, window_hi = observed_interval(window, frame, spectrum.redshift)
        continuum_mask |= (wavelength >= window_lo) & (wavelength <= window_hi)
    if continuum_config.get("windows") and continuum_mask.sum() < max(8, degree + 2):
        raise ValueError(
            f"{spectrum.spectrum_id}: configured continuum windows contain too few pixels"
        )
    if not continuum_config.get("windows"):
        continuum_mask[:] = True

    coefficients, continuum_scatter = _sigma_clip_polyfit(
        x_normalized[continuum_mask], flux[continuum_mask], degree
    )
    if uncertainty is None:
        if not np.isfinite(continuum_scatter) or continuum_scatter <= 0:
            raise ValueError(
                f"{spectrum.spectrum_id}: could not estimate a positive uncertainty"
            )
        uncertainty = np.full(flux.shape, continuum_scatter)
    else:
        uncertainty = np.asarray(uncertainty, dtype=float)

    noise_config = fit.get("noise", {})
    configured_marginal_scale = noise_config.get("marginal_scale", 1.0)
    if configured_marginal_scale == "auto":
        continuum_model = np.polynomial.polynomial.polyval(
            x_normalized, coefficients
        )
        normalized_continuum_residual = (
            flux[continuum_mask] - continuum_model[continuum_mask]
        ) / uncertainty[continuum_mask]
        noise_marginal_scale = max(
            1.0, float(robust_sigma(normalized_continuum_residual))
        )
    else:
        noise_marginal_scale = float(configured_marginal_scale)
    uncertainty = uncertainty * noise_marginal_scale

    noise_level = float(np.median(uncertainty[np.isfinite(uncertainty)]))
    if not np.isfinite(noise_level) or noise_level <= 0:
        raise ValueError(f"{spectrum.spectrum_id}: invalid uncertainty values")
    if not np.isfinite(continuum_scatter) or continuum_scatter <= 0:
        continuum_scatter = noise_level

    noise_model = str(noise_config.get("model", "independent"))
    noise_rho = 0.0
    if noise_model == "ar1":
        configured_rho = noise_config.get("rho", "auto")
        if configured_rho == "auto":
            continuum_model = np.polynomial.polynomial.polyval(
                x_normalized, coefficients
            )
            normalized = (flux - continuum_model) / uncertainty
            normalized = normalized - np.median(normalized[continuum_mask])
            spacing = np.diff(wavelength)
            typical_spacing = float(np.median(spacing[spacing > 0]))
            adjacent = (
                continuum_mask[:-1]
                & continuum_mask[1:]
                & np.isfinite(normalized[:-1])
                & np.isfinite(normalized[1:])
                & (spacing <= 1.5 * typical_spacing)
            )
            if adjacent.sum() < 8:
                raise ValueError(
                    f"{spectrum.spectrum_id}: too few adjacent continuum pixels "
                    "to estimate AR(1) noise correlation"
                )
            left = normalized[:-1][adjacent]
            right = normalized[1:][adjacent]
            denominator = float(np.dot(left, left))
            if not np.isfinite(denominator) or denominator <= 0:
                raise ValueError(
                    f"{spectrum.spectrum_id}: could not estimate AR(1) noise correlation"
                )
            noise_rho = float(np.clip(np.dot(left, right) / denominator, -0.95, 0.95))
        else:
            noise_rho = float(configured_rho)

    prior_width = float(continuum_config.get("prior_width_rms", 10.0))
    coefficient_half_width = np.full(degree + 1, prior_width * continuum_scatter)
    # Higher polynomial terms can move by the same total flux scale over x=[-1, 1].
    coefficient_half_width[0] = max(
        coefficient_half_width[0], 0.1 * max(abs(coefficients[0]), noise_level)
    )

    kin = fit["kinematics"]
    sigma_min, sigma_max = map(float, kin["sigma_kms"])
    rest_wavelengths = [float(line["wavelength"]) for line in fit["lines"]]
    lsf_config = fit.get("lsf", {"model": "none"})
    min_sigma_angstrom = min(
        convolved_sigma(
            rest, spectrum.redshift, sigma_min, 0.0, lsf_config, spectrum.metadata
        )
        for rest in rest_wavelengths
    )
    max_sigma_angstrom = max(
        convolved_sigma(
            rest, spectrum.redshift, sigma_max, 0.0, lsf_config, spectrum.metadata
        )
        for rest in rest_wavelengths
    )
    baseline = np.polynomial.polynomial.polyval(x_normalized, coefficients)
    positive_peak = float(np.max(np.maximum(flux - baseline, 0.0)))
    flux_prior = fit.get("flux_prior", {})
    lower = (
        float(flux_prior.get("min_snr", 0.1))
        * noise_level
        * min_sigma_angstrom
        * SQRT_2PI
    )
    upper = (
        float(flux_prior.get("max_signal_factor", 20.0))
        * max(positive_peak, 5.0 * noise_level)
        * max_sigma_angstrom
        * SQRT_2PI
    )
    lower = max(lower, np.finfo(float).tiny)
    upper = max(upper, lower * 100.0)

    return PreparedSpectrum(
        spectrum=spectrum,
        wavelength=wavelength,
        flux=flux,
        uncertainty=uncertainty,
        x_normalized=x_normalized,
        continuum_center=coefficients,
        continuum_half_width=coefficient_half_width,
        noise_level=noise_level,
        noise_model=noise_model,
        noise_rho=noise_rho,
        noise_marginal_scale=noise_marginal_scale,
        positive_peak=positive_peak,
        default_flux_bounds=(lower, upper),
    )


class ModelDefinition:
    """Parameter layout, prior transform, and likelihood for one component count."""

    def __init__(self, prepared: PreparedSpectrum, fit: dict[str, Any], n_components: int):
        self.prepared = prepared
        self.fit = fit
        self.n_components = int(n_components)
        self.lines = fit["lines"]
        self.line_by_name = {line["name"]: line for line in self.lines}
        self.free_lines = [line for line in self.lines if "ratio_to" not in line]
        self.fixed_components = fit.get("fixed_components", [])
        self.broad_components = fit.get("broad_components", [])
        self.lsf_config = fit.get("lsf", {"model": "none"})
        self.continuum_count = int(fit.get("continuum", {}).get("degree", 1)) + 1
        self.parameter_names = self._parameter_names()
        self._validate_fixed_components()
        self._validate_broad_components()

    def _validate_fixed_components(self) -> None:
        names: set[str] = set()
        for component in self.fixed_components:
            name = component.get("name")
            if not name or name in names:
                raise ConfigError("Every fixed component needs a unique name")
            names.add(name)
            if component.get("line") not in self.line_by_name:
                raise ConfigError(f"Fixed component {name} refers to an unknown line")
            if float(component.get("sigma_kms", 0)) <= 0:
                raise ConfigError(f"Fixed component {name} needs positive sigma_kms")

    def _validate_broad_components(self) -> None:
        names: set[str] = set()
        for component in self.broad_components:
            name = component.get("name")
            if not name or name in names:
                raise ConfigError("Every broad component needs a unique name")
            names.add(name)
            if component.get("line") not in self.line_by_name:
                raise ConfigError(f"Broad component {name} refers to an unknown line")

    def _parameter_names(self) -> list[str]:
        names = [f"continuum_c{index}" for index in range(self.continuum_count)]
        names.extend(f"fixed.{item['name']}.flux" for item in self.fixed_components)
        for item in self.broad_components:
            names.extend(
                [
                    f"broad.{item['name']}.velocity_kms",
                    f"broad.{item['name']}.sigma_kms",
                    f"broad.{item['name']}.flux",
                ]
            )
        for component in range(1, self.n_components + 1):
            names.extend(
                [f"component.{component}.velocity_kms", f"component.{component}.sigma_kms"]
            )
            names.extend(
                f"component.{component}.{line['name']}.flux" for line in self.free_lines
            )
        return names

    @property
    def ndim(self) -> int:
        return len(self.parameter_names)

    def _flux_bounds(self, spec: dict[str, Any]) -> tuple[float, float]:
        values = spec.get("flux_bounds")
        if values is None:
            if "sigma_kms" not in spec or "line" not in spec:
                return self.prepared.default_flux_bounds
            line = self.line_by_name[spec["line"]]
            sigma_value = spec["sigma_kms"]
            if isinstance(sigma_value, (list, tuple)):
                sigma_value = max(map(float, sigma_value))
            velocity_value = spec.get("velocity_kms", 0.0)
            if isinstance(velocity_value, (list, tuple)):
                velocity_value = 0.5 * sum(map(float, velocity_value))
            sigma_angstrom = convolved_sigma(
                float(line["wavelength"]),
                self.prepared.spectrum.redshift,
                float(sigma_value),
                float(velocity_value),
                self.lsf_config,
                self.prepared.spectrum.metadata,
            )
            flux_prior = self.fit.get("flux_prior", {})
            lo = (
                float(flux_prior.get("min_snr", 0.1))
                * self.prepared.noise_level
                * sigma_angstrom
                * SQRT_2PI
            )
            hi = (
                float(flux_prior.get("max_signal_factor", 20.0))
                * max(self.prepared.positive_peak, 5.0 * self.prepared.noise_level)
                * sigma_angstrom
                * SQRT_2PI
            )
            return max(lo, np.finfo(float).tiny), max(hi, lo * 100.0)
        lo, hi = float(values[0]), float(values[1])
        if lo <= 0 or hi <= lo:
            raise ConfigError("flux_bounds must be two increasing positive values")
        return lo, hi

    @staticmethod
    def _log_uniform(unit_value: float, bounds: tuple[float, float]) -> float:
        lo, hi = bounds
        return float(lo * np.exp(float(unit_value) * np.log(hi / lo)))

    def prior_transform(self, unit_cube: np.ndarray) -> np.ndarray:
        unit_cube = np.asarray(unit_cube, dtype=float)
        if unit_cube.size != self.ndim:
            raise ValueError(f"Expected {self.ndim} prior coordinates")
        params = np.empty_like(unit_cube)
        index = 0
        for coefficient in range(self.continuum_count):
            center = self.prepared.continuum_center[coefficient]
            half_width = self.prepared.continuum_half_width[coefficient]
            params[index] = center + (2.0 * unit_cube[index] - 1.0) * half_width
            index += 1

        for component in self.fixed_components:
            params[index] = self._log_uniform(
                unit_cube[index], self._flux_bounds(component)
            )
            index += 1

        for component in self.broad_components:
            velocity_lo, velocity_hi = map(
                float, component.get("velocity_kms", [-2000.0, 2000.0])
            )
            params[index] = velocity_lo + unit_cube[index] * (
                velocity_hi - velocity_lo
            )
            params[index + 1] = self._log_uniform(
                unit_cube[index + 1], tuple(map(float, component["sigma_kms"]))
            )
            params[index + 2] = self._log_uniform(
                unit_cube[index + 2], self._flux_bounds(component)
            )
            index += 3

        velocity_lo, velocity_hi = map(float, self.fit["kinematics"]["velocity_kms"])
        sigma_lo, sigma_hi = map(float, self.fit["kinematics"]["sigma_kms"])
        step = 2 + len(self.free_lines)
        velocity_units = [unit_cube[index + comp * step] for comp in range(self.n_components)]
        velocities = np.sort(
            velocity_lo + np.asarray(velocity_units) * (velocity_hi - velocity_lo)
        )
        for component in range(self.n_components):
            params[index] = velocities[component]
            params[index + 1] = self._log_uniform(
                unit_cube[index + 1], (sigma_lo, sigma_hi)
            )
            for line_index, line in enumerate(self.free_lines):
                params[index + 2 + line_index] = self._log_uniform(
                    unit_cube[index + 2 + line_index], self._flux_bounds(line)
                )
            index += step
        return params

    def evaluate(self, params: np.ndarray) -> np.ndarray:
        params = np.asarray(params, dtype=float)
        if params.size != self.ndim:
            raise ValueError(f"Expected {self.ndim} model parameters")
        p = self.prepared
        index = self.continuum_count
        model = np.polynomial.polynomial.polyval(
            p.x_normalized, params[: self.continuum_count]
        )

        for component in self.fixed_components:
            line = self.line_by_name[component["line"]]
            velocity = float(component.get("velocity_kms", 0.0))
            sigma_kms = float(component["sigma_kms"])
            model += gaussian_integrated(
                p.wavelength,
                observed_center(float(line["wavelength"]), p.spectrum.redshift, velocity),
                convolved_sigma(
                    float(line["wavelength"]),
                    p.spectrum.redshift,
                    sigma_kms,
                    velocity,
                    self.lsf_config,
                    p.spectrum.metadata,
                ),
                params[index],
            )
            index += 1

        for component in self.broad_components:
            line = self.line_by_name[component["line"]]
            velocity, sigma_kms, flux = params[index : index + 3]
            rest = float(line["wavelength"])
            model += gaussian_integrated(
                p.wavelength,
                observed_center(rest, p.spectrum.redshift, velocity),
                convolved_sigma(
                    rest,
                    p.spectrum.redshift,
                    sigma_kms,
                    velocity,
                    self.lsf_config,
                    p.spectrum.metadata,
                ),
                flux,
            )
            index += 3

        step = 2 + len(self.free_lines)
        for _ in range(self.n_components):
            velocity, sigma_kms = params[index], params[index + 1]
            free_fluxes = {
                line["name"]: params[index + 2 + line_index]
                for line_index, line in enumerate(self.free_lines)
            }
            all_fluxes: dict[str, float] = {}
            for line in self.lines:
                if "ratio_to" in line:
                    all_fluxes[line["name"]] = (
                        all_fluxes[line["ratio_to"]] * float(line["ratio"])
                    )
                else:
                    all_fluxes[line["name"]] = free_fluxes[line["name"]]
                rest = float(line["wavelength"])
                model += gaussian_integrated(
                    p.wavelength,
                    observed_center(rest, p.spectrum.redshift, velocity),
                    convolved_sigma(
                        rest,
                        p.spectrum.redshift,
                        sigma_kms,
                        velocity,
                        self.lsf_config,
                        p.spectrum.metadata,
                    ),
                    all_fluxes[line["name"]],
                )
            index += step
        return model

    def log_likelihood(self, params: np.ndarray) -> float:
        residual = (self.prepared.flux - self.evaluate(params)) / self.prepared.uncertainty
        log_variance = 2.0 * np.log(self.prepared.uncertainty)
        if self.prepared.noise_model == "independent":
            return float(
                -0.5
                * np.sum(residual**2 + np.log(2.0 * np.pi) + log_variance)
            )

        # A stationary AR(1) process on an irregular wavelength grid is the
        # exponential-covariance Markov model.  Scaling rho by the local pixel
        # separation avoids correlating across gaps in masked spectra.
        spacing = np.diff(self.prepared.wavelength)
        typical_spacing = float(np.median(spacing[spacing > 0]))
        phi = np.sign(self.prepared.noise_rho) * np.abs(
            self.prepared.noise_rho
        ) ** (spacing / typical_spacing)
        innovation_variance = np.maximum(1.0 - phi**2, np.finfo(float).eps)
        quadratic = residual[0] ** 2 + np.sum(
            (residual[1:] - phi * residual[:-1]) ** 2 / innovation_variance
        )
        log_determinant_correlation = float(np.sum(np.log(innovation_variance)))
        normalization = (
            residual.size * np.log(2.0 * np.pi)
            + np.sum(log_variance)
            + log_determinant_correlation
        )
        return float(-0.5 * (quadratic + normalization))

    def component_models(self, params: np.ndarray) -> list[np.ndarray]:
        """Return line-only arrays for diagnostic plots."""
        params = np.asarray(params, dtype=float)
        p = self.prepared
        models: list[np.ndarray] = []
        index = self.continuum_count
        for component in self.fixed_components:
            line = self.line_by_name[component["line"]]
            rest = float(line["wavelength"])
            models.append(
                gaussian_integrated(
                    p.wavelength,
                    observed_center(rest, p.spectrum.redshift, float(component.get("velocity_kms", 0))),
                    convolved_sigma(
                        rest,
                        p.spectrum.redshift,
                        float(component["sigma_kms"]),
                        float(component.get("velocity_kms", 0)),
                        self.lsf_config,
                        p.spectrum.metadata,
                    ),
                    params[index],
                )
            )
            index += 1

        for component in self.broad_components:
            line = self.line_by_name[component["line"]]
            rest = float(line["wavelength"])
            velocity, sigma_kms, flux = params[index : index + 3]
            models.append(
                gaussian_integrated(
                    p.wavelength,
                    observed_center(rest, p.spectrum.redshift, velocity),
                    convolved_sigma(
                        rest,
                        p.spectrum.redshift,
                        sigma_kms,
                        velocity,
                        self.lsf_config,
                        p.spectrum.metadata,
                    ),
                    flux,
                )
            )
            index += 3

        step = 2 + len(self.free_lines)
        for _ in range(self.n_components):
            velocity, sigma_kms = params[index], params[index + 1]
            free_fluxes = {
                line["name"]: params[index + 2 + line_index]
                for line_index, line in enumerate(self.free_lines)
            }
            component_model = np.zeros_like(p.wavelength)
            all_fluxes: dict[str, float] = {}
            for line in self.lines:
                if "ratio_to" in line:
                    all_fluxes[line["name"]] = all_fluxes[line["ratio_to"]] * float(line["ratio"])
                else:
                    all_fluxes[line["name"]] = free_fluxes[line["name"]]
                rest = float(line["wavelength"])
                component_model += gaussian_integrated(
                    p.wavelength,
                    observed_center(rest, p.spectrum.redshift, velocity),
                    convolved_sigma(
                        rest,
                        p.spectrum.redshift,
                        sigma_kms,
                        velocity,
                        self.lsf_config,
                        p.spectrum.metadata,
                    ),
                    all_fluxes[line["name"]],
                )
            models.append(component_model)
            index += step
        return models
