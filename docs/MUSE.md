# MUSE adapter

The `muse` adapter targets ESO Phase-3 `SCIENCE.CUBE.IFS` products with named
`DATA` and `STAT` extensions. It validates the instrument and product type,
the three-dimensional extension layout, matching variance shape, and the
`AWAV` air-wavelength axis before yielding any spectra.

## Safe subset workflow

MUSE cubes are large. Adapter configurations must supply both `x_range` and
`y_range`, expressed as zero-based, half-open bounds. A range of `[82, 87]`
contains five pixels: 82 through 86. BEAT opens FITS files with memory mapping
and materializes only one selected spaxel spectrum at a time.

Start with a 5-by-5 or 10-by-10 region, run `beat validate`, inspect diagnostic
plots, and only then expand the region. A full cube is possible only with the
explicit `input.allow_full_cube: true` override. That override should normally
be paired with a batch-computing plan rather than a laptop run.

The two `examples/muse_*.local.yaml` files select 25 spaxels near the target
centres in the local NGC 2992 and NGC 3393 products.

## Units and wavelength convention

Adapter defaults are:

- `DATA` flux and `STAT` variance extensions;
- a flux-density multiplier of `1e-20`, matching these Phase-3 products;
- air wavelengths in Angstrom from the `DATA` extension WCS; and
- target names from the primary `OBJECT` header.

Configured line wavelengths must therefore use the air convention. The raw
FITS `BUNIT` strings are preserved in each result's metadata.

## Line-spread function

The model treats the LSF as Gaussian and adds instrumental and intrinsic
sigmas in quadrature. `fit.lsf.model: instrument` uses the scalar `SPEC_RES`
resolving power recorded in the Phase-3 primary header. BEAT marks this as an
approximation and preserves the referenced `LSF_PROFILE` calibration filename
and checksum when available.

MUSE resolution varies with wavelength, so publication analyses should use a
wavelength-dependent calibration when available. Override the adapter default
with `polynomial_fwhm_angstrom` or `table`; see the configuration reference.
The scalar-header mode is suitable for initial fitting and adapter regression,
not a claim that the true MUSE LSF is constant.

The H-alpha alpha-test configurations instead use the representative UDF-10
WFM Gaussian approximation

`FWHM(lambda) = 5.866e-8 lambda^2 - 9.187e-4 lambda + 6.040 Angstrom`.

This empirical relation supports initial WFM validation but is not a
substitute for a configuration-matched `LSF_PROFILE`. Its identity and
approximation status should remain part of run provenance.
