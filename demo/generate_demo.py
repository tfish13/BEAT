#!/usr/bin/env python3
"""Generate a deterministic [O III] doublet spectrum for the BEAT smoke test."""

from pathlib import Path

import numpy as np


def gaussian(wavelength, center, sigma, integrated_flux):
    amplitude = integrated_flux / (sigma * np.sqrt(2.0 * np.pi))
    return amplitude * np.exp(-0.5 * ((wavelength - center) / sigma) ** 2)


redshift = 0.01
rng = np.random.default_rng(20260717)
wavelength = np.arange(4900.0, 5120.0, 0.5)
uncertainty = np.full(wavelength.shape, 1.0)
flux = 5.0 + 0.002 * (wavelength - wavelength.mean())

velocity_kms = 80.0
doppler = np.sqrt((1.0 + velocity_kms / 299792.458) / (1.0 - velocity_kms / 299792.458))
sigma_kms = 110.0
center_5007 = 5006.84 * (1.0 + redshift) * doppler
center_4959 = 4958.92 * (1.0 + redshift) * doppler
sigma_5007 = center_5007 * sigma_kms / 299792.458
sigma_4959 = center_4959 * sigma_kms / 299792.458
flux += gaussian(wavelength, center_5007, sigma_5007, 120.0)
flux += gaussian(wavelength, center_4959, sigma_4959, 120.0 * 0.33557)
flux += rng.normal(0.0, uncertainty)

output = Path(__file__).with_name("demo_spectrum.txt")
np.savetxt(
    output,
    np.column_stack([wavelength, flux, uncertainty]),
    header="wavelength_angstrom flux_density uncertainty_1sigma",
)
print(f"Wrote {output}")
