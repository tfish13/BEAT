"""Instrumental line-spread-function models.

BEAT currently treats the LSF as Gaussian and combines its sigma with the
intrinsic Gaussian line sigma in quadrature. More complex non-Gaussian kernels
can be added later without changing the instrument-adapter interface.
"""

from __future__ import annotations

from functools import lru_cache
from importlib import resources
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits


FWHM_TO_SIGMA = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))


class LSFError(ValueError):
    """Raised for a missing or invalid LSF model."""


def validate_lsf_config(config: dict[str, Any], instrument_available: bool) -> None:
    model = str(config.get("model", "none")).lower()
    allowed = {
        "none",
        "instrument",
        "resolving_power",
        "polynomial_resolving_power",
        "nirspec_etc_resolving_power",
        "constant_fwhm_angstrom",
        "polynomial_fwhm_angstrom",
        "table",
    }
    if model not in allowed:
        raise LSFError(f"Unknown LSF model {model!r}; choose from {', '.join(sorted(allowed))}")
    if model == "instrument" and not instrument_available:
        raise LSFError("LSF model 'instrument' requires a named input adapter")
    if model == "resolving_power" and float(config.get("value", 0)) <= 0:
        raise LSFError("resolving_power LSF requires a positive value")
    if model == "polynomial_resolving_power":
        coefficients = config.get("coefficients")
        if not isinstance(coefficients, list) or not coefficients:
            raise LSFError("polynomial_resolving_power requires a coefficient list")
        if float(config.get("scale_angstrom", 1.0)) == 0:
            raise LSFError("polynomial resolving-power scale_angstrom cannot be zero")
    if model == "nirspec_etc_resolving_power":
        grating = str(config.get("grating", "")).strip().upper()
        if grating not in {"G235H", "G395H"}:
            raise LSFError(
                "nirspec_etc_resolving_power currently supports G235H and G395H"
            )
    if model == "constant_fwhm_angstrom" and float(config.get("value", 0)) <= 0:
        raise LSFError("constant_fwhm_angstrom LSF requires a positive value")
    if model == "polynomial_fwhm_angstrom":
        coefficients = config.get("coefficients")
        if not isinstance(coefficients, list) or not coefficients:
            raise LSFError("polynomial_fwhm_angstrom requires a coefficient list")
        if float(config.get("scale_angstrom", 1.0)) == 0:
            raise LSFError("polynomial LSF scale_angstrom cannot be zero")
    if model == "table":
        if not config.get("path"):
            raise LSFError("table LSF requires path")
        if config.get("quantity", "fwhm_angstrom") not in {
            "fwhm_angstrom",
            "sigma_angstrom",
            "resolving_power",
        }:
            raise LSFError(
                "table LSF quantity must be fwhm_angstrom, sigma_angstrom, or resolving_power"
            )


def _resolved_config(
    config: dict[str, Any] | None, metadata: dict[str, Any] | None
) -> dict[str, Any]:
    config = dict(config or {"model": "none"})
    if config.get("model", "none") != "instrument":
        return config
    if not metadata or not metadata.get("instrument_lsf"):
        raise LSFError("Spectrum metadata does not provide an instrument LSF")
    return dict(metadata["instrument_lsf"])


@lru_cache(maxsize=32)
def _load_table(path: str) -> tuple[np.ndarray, np.ndarray]:
    values = np.loadtxt(Path(path), dtype=float)
    if values.ndim != 2 or values.shape[1] < 2:
        raise LSFError("LSF table must contain at least two numeric columns")
    order = np.argsort(values[:, 0])
    wavelength = values[order, 0]
    profile = values[order, 1]
    if wavelength.size < 2 or np.any(np.diff(wavelength) <= 0):
        raise LSFError("LSF table wavelengths must be unique and increasing")
    return wavelength, profile


@lru_cache(maxsize=8)
def _load_nirspec_etc_table(grating: str) -> tuple[np.ndarray, np.ndarray]:
    """Load an STScI ETC NIRSpec dispersion/resolution table bundled with BEAT."""
    grating = grating.strip().lower()
    if grating not in {"g235h", "g395h"}:
        raise LSFError(f"No bundled NIRSpec ETC resolution table for {grating!r}")
    filename = f"jwst_nirspec_{grating}_disp.fits"
    resource = resources.files("beat").joinpath("calibration", filename)
    with resources.as_file(resource) as path:
        try:
            values = fits.getdata(path, ext=1)
            wavelength = np.asarray(values["WAVELENGTH"], dtype=float) * 1.0e4
            resolving_power = np.asarray(values["R"], dtype=float)
        except (OSError, KeyError, ValueError) as exc:
            raise LSFError(f"Could not read bundled NIRSpec table {filename}") from exc
    if (
        wavelength.size < 2
        or wavelength.shape != resolving_power.shape
        or np.any(np.diff(wavelength) <= 0)
        or np.any(~np.isfinite(resolving_power))
        or np.any(resolving_power <= 0)
    ):
        raise LSFError(f"Bundled NIRSpec table {filename} is invalid")
    return wavelength, resolving_power


def lsf_sigma_angstrom(
    wavelength_angstrom: float | np.ndarray,
    config: dict[str, Any] | None,
    metadata: dict[str, Any] | None = None,
) -> float | np.ndarray:
    """Evaluate Gaussian instrumental sigma at one or more wavelengths."""
    wavelength = np.asarray(wavelength_angstrom, dtype=float)
    resolved = _resolved_config(config, metadata)
    model = str(resolved.get("model", "none")).lower()
    if model == "none":
        sigma = np.zeros_like(wavelength)
    elif model == "resolving_power":
        sigma = wavelength / float(resolved["value"]) * FWHM_TO_SIGMA
    elif model == "polynomial_resolving_power":
        reference = float(resolved.get("reference_angstrom", 0.0))
        scale = float(resolved.get("scale_angstrom", 1.0))
        coordinate = (wavelength - reference) / scale
        resolving_power = np.polynomial.polynomial.polyval(
            coordinate, resolved["coefficients"]
        )
        if np.any(resolving_power <= 0):
            raise LSFError("Polynomial LSF produced non-positive resolving power")
        sigma = wavelength / resolving_power * FWHM_TO_SIGMA
    elif model == "nirspec_etc_resolving_power":
        table_wave, resolving_power_table = _load_nirspec_etc_table(
            str(resolved["grating"])
        )
        if np.any(wavelength < table_wave[0]) or np.any(wavelength > table_wave[-1]):
            raise LSFError(
                "Requested wavelength lies outside the bundled NIRSpec resolution table"
            )
        resolving_power = np.interp(wavelength, table_wave, resolving_power_table)
        sigma = wavelength / resolving_power * FWHM_TO_SIGMA
    elif model == "constant_fwhm_angstrom":
        sigma = np.full_like(wavelength, float(resolved["value"]) * FWHM_TO_SIGMA)
    elif model == "polynomial_fwhm_angstrom":
        reference = float(resolved.get("reference_angstrom", 0.0))
        scale = float(resolved.get("scale_angstrom", 1.0))
        coordinate = (wavelength - reference) / scale
        fwhm = np.polynomial.polynomial.polyval(coordinate, resolved["coefficients"])
        sigma = fwhm * FWHM_TO_SIGMA
    elif model == "table":
        table_wave, table_value = _load_table(str(resolved["path"]))
        if np.any(wavelength < table_wave[0]) or np.any(wavelength > table_wave[-1]):
            raise LSFError("Requested wavelength lies outside the LSF table")
        profile = np.interp(wavelength, table_wave, table_value)
        quantity = resolved.get("quantity", "fwhm_angstrom")
        if quantity == "fwhm_angstrom":
            sigma = profile * FWHM_TO_SIGMA
        elif quantity == "sigma_angstrom":
            sigma = profile
        else:
            sigma = wavelength / profile * FWHM_TO_SIGMA
    else:
        raise LSFError(f"Unsupported resolved LSF model: {model}")
    if np.any(~np.isfinite(sigma)) or np.any(sigma < 0):
        raise LSFError("LSF produced a non-finite or negative sigma")
    return float(sigma) if sigma.ndim == 0 else sigma
