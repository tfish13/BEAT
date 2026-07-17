# JWST/NIRSpec IFU adapter

The `nirspec` adapter targets calibrated JWST Stage-3 `s3d` IFU cubes. It
requires `TELESCOP=JWST`, `INSTRUME=NIRSPEC`, `EXP_TYPE=NRS_IFU`, and matching
three-dimensional `SCI`, `ERR`, and `DQ` extensions.

## Safe subset workflow

Both `x_range` and `y_range` are mandatory zero-based, half-open bounds unless
`allow_full_cube: true` is deliberately selected. The example configurations
use only 5-by-5 regions. SCI and ERR remain memory-mapped. The unsigned DQ cube
is read in its raw FITS representation so it can also remain memory-mapped;
the adapter tests the low-order `DO_NOT_USE` bit separately for each selected
spaxel and never constructs a full-cube Boolean mask.

## Wavelength, flux, and uncertainty

NIRSpec wavelengths are vacuum wavelengths. BEAT converts the cube's micron
axis to Angstrom. It uses the `WCS-TABLE` wavelength vector when present and
otherwise uses the linear SCI-extension WCS. This supports both current files:
the NGC 4151 pipeline-2.0.1 product contains a wavelength table, while the
IC 5063 pipeline-1.20.2 product uses the linear axis.

SCI and ERR are retained in their calibrated native `MJy/sr` units. Therefore,
BEAT's integrated line fluxes are initially reported in `MJy/sr Angstrom`.
This avoids silently assuming a point-source aperture or converting surface
brightness into per-spaxel flux. The result metadata records `PIXAR_SR` and
`PIXAR_A2`, allowing an explicit surface-brightness or aperture conversion in
a later, documented analysis step.

## Line-spread function and caveats

G235H and G395H use the wavelength-dependent dispersion/resolution FITS tables
published by STScI and used by the JWST Exposure Time Calculator. The tables
are bundled with BEAT so the LSF is reproducible offline, and their SHA-256
checksums are written to spectrum metadata. STScI defines the tabulated
resolving power using a 2.2-pixel resolution element, representative of a
fully illuminated IFU aperture. BEAT converts `R(lambda)` to a Gaussian
equivalent with `FWHM=lambda/R` and combines that width in quadrature with the
intrinsic Gaussian line width.

At the adopted redshifts, the tables give `R=2240.45` at NGC 4151 [Si VI]
1.9634 micron and `R=2811.60` at IC 5063 Br-alpha. These differ materially
from a scalar `R=2700`. Wavelengths outside a bundled table fail explicitly
rather than being extrapolated. Other NIRSpec dispersers retain nominal
fallbacks (about 2700 for high, 1000 for medium, and 100 for the prism) until
their tables are packaged.

This remains an approximation: source extent within an IFU slice and a
non-Gaussian or empirically measured LSF can change the effective profile.
The cube's CRDS disperser reference is recorded separately because it
calibrates the wavelength mapping; it is not treated as an empirical LSF.

High-resolution IFU products can contain detector-gap wavelengths, and the
pipeline's cube resampling can introduce correlated structure in
individual-spaxel spectra. DQ masking handles unusable samples but cannot make
the remaining pixels statistically independent. For correlated products, use
`fit.noise: {model: ar1, rho: auto, marginal_scale: auto}`. BEAT conservatively
inflates underestimated formal errors, estimates rho from adjacent pixels
inside the configured continuum windows, and uses the corresponding Gaussian
AR(1) likelihood. Both resolved values are recorded in each result. This
first-order model should still be checked against residual diagnostics and
small-aperture extractions.

Known nearby lines should not be left for the component model to absorb. The
NGC 4151 example uses a broad continuum window for stable noise estimation and
`fit.exclude_windows` around H2 1-0 S(3), whose wavelength otherwise overlaps
the allowed high-velocity [Si VI] range.

## Injection/recovery validation

`validation/run_nirspec_injection_recovery.py` runs controlled Gaussian-noise
and block-resampled real-residual cases on the actual G235H and G395H cube
sampling. It keeps the automatically selected high-S/N science spaxels
separate from low-emission noise donors. The complete 16-case standard core
recovers all component counts but has only four cases per count class. A
separate 40-case, two-donor matrix with empirical marginal-error calibration
and AR(1) noise recovers 20/20 blanks and 20/20 S/N=10 singles with no evidence
flags. This passes the predeclared blank/single gates but does not power the
double/triple classes. See `validation/NIRSPEC_INJECTION_RECOVERY.md`.
