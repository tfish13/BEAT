"""The dataset-independent spectrum object used by the fitter."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np


@dataclass
class Spectrum:
    """A single target or cube spaxel with wavelength in Angstrom."""

    spectrum_id: str
    wavelength: np.ndarray
    flux: np.ndarray
    uncertainty: np.ndarray | None
    redshift: float
    mask: np.ndarray | None = None
    metadata: dict[str, Any] = field(default_factory=dict)

    def prepared(self) -> "Spectrum":
        """Return a sorted, finite copy with duplicate wavelengths removed."""
        if not np.isfinite(self.redshift) or self.redshift <= -1:
            raise ValueError(f"{self.spectrum_id}: redshift must be finite and greater than -1")
        wavelength = np.asarray(self.wavelength, dtype=float).reshape(-1)
        flux = np.asarray(self.flux, dtype=float).reshape(-1)
        if wavelength.shape != flux.shape:
            raise ValueError(f"{self.spectrum_id}: wavelength and flux shapes differ")
        if wavelength.size < 2:
            raise ValueError(f"{self.spectrum_id}: spectrum has fewer than two pixels")

        uncertainty = None
        if self.uncertainty is not None:
            uncertainty = np.asarray(self.uncertainty, dtype=float).reshape(-1)
            if uncertainty.shape != wavelength.shape:
                raise ValueError(
                    f"{self.spectrum_id}: uncertainty shape differs from flux"
                )

        valid = np.isfinite(wavelength) & np.isfinite(flux)
        if uncertainty is not None:
            valid &= np.isfinite(uncertainty) & (uncertainty > 0)
        if self.mask is not None:
            mask = np.asarray(self.mask, dtype=bool).reshape(-1)
            if mask.shape != wavelength.shape:
                raise ValueError(f"{self.spectrum_id}: mask shape differs from flux")
            valid &= ~mask

        wavelength = wavelength[valid]
        flux = flux[valid]
        if uncertainty is not None:
            uncertainty = uncertainty[valid]

        order = np.argsort(wavelength)
        wavelength = wavelength[order]
        flux = flux[order]
        if uncertainty is not None:
            uncertainty = uncertainty[order]

        unique = np.concatenate(([True], np.diff(wavelength) > 0))
        wavelength = wavelength[unique]
        flux = flux[unique]
        if uncertainty is not None:
            uncertainty = uncertainty[unique]
        if wavelength.size < 2:
            raise ValueError(
                f"{self.spectrum_id}: fewer than two valid unique wavelength pixels"
            )

        return Spectrum(
            spectrum_id=str(self.spectrum_id),
            wavelength=wavelength,
            flux=flux,
            uncertainty=uncertainty,
            redshift=float(self.redshift),
            mask=None,
            metadata=dict(self.metadata),
        )
