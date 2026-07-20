# Configuration reference

All wavelengths supplied to the model are Angstrom after `wavelength_unit`
conversion. Unknown keys are retained in the archived manifest, but misspelled
core keys may therefore not have an effect; use `beat validate` and review its
reported wavelength range.

## `input`

Common scaling fields:

- `wavelength_unit`: `angstrom`, `nm`, `micron`, or `m`;
- `flux_scale`: multiplicative flux-density scale (default 1);
- `uncertainty_scale`: multiplicative uncertainty scale (defaults to
  `flux_scale`); and
- `uncertainty_kind`: `sigma`, `variance`, or `inverse_variance`.

For `kind: cube`, required fields are `path`, `flux_hdu`, and `redshift`.
Optional fields include `uncertainty_hdu`, `spectral_axis` (default 0),
`wavelength_hdu`, `wavelength_header_hdu`, `wavelength_fits_axis`,
`wavelength_log10`, `mask_hdu`, `x_range`, and `y_range`.

`input.adapter: muse` fills the standard Phase-3 MUSE cube defaults, including
`kind`, named DATA/STAT HDUs, variance handling, air-wavelength WCS, and flux
scale. Both `x_range` and `y_range` are then mandatory. Set
`allow_full_cube: true` only to make a deliberate full-cube run.

`input.adapter: nirspec` supports JWST Stage-3 NIRSpec IFU `s3d` products. It
defaults to SCI flux, ERR sigma uncertainty, DQ `DO_NOT_USE` masking, micron to
Angstrom conversion, vacuum wavelengths, native MJy/sr units, and target names
from `TARGNAME`. Bounded `x_range` and `y_range` selections are mandatory by
default. If a WCS-TABLE wavelength vector exists it takes precedence over the
linear SCI-header wavelength axis.

`input.adapter: miri` accepts `glob` in place of `path` for the 12 MIRI MRS
band cubes. It selects the one segment containing the observed fit window.
Use `segment: 1A` through `4C` (or forms such as `3-medium`) when overlapping
bands make selection ambiguous. SCI/ERR/DQ, vacuum micron wavelengths, native
MJy/sr units, and bounded spatial selections follow the JWST IFU conventions
described above.

For `kind: survey_table`, required fields are `path`, `hdu`, `id_column`,
`wavelength_column`, and `flux_column`, plus either `redshift_column` or
`redshift`. Optional fields include `uncertainty_column`, `mask_column`,
`row_start`, and `row_stop`.

For `kind: spectrum_files`, provide `glob` or `files` and set `format` to
`fits_table`, `fits_image`, or `ascii`. FITS-image wavelength axes use the same
linear-WCS fields as a cube. A redshift is read from `redshift_header` or the
fixed `redshift` value. `id_header` overrides a filename-derived ID.

## `fit`

- `frame`: whether `window` and continuum windows are `rest` or `observed`;
- `window`: two wavelength bounds;
- `exclude_windows`: optional wavelength intervals, in the configured frame,
  removed before fitting (useful for known contaminating lines or artifacts);
- `minimum_valid_pixels`: fail a spectrum with fewer valid pixels (default 20);
- `continuum.degree`: 0, 1, or 2;
- `continuum.windows`: optional intervals used for initial continuum and noise
  estimates; if too few pixels are present, the full fitting window is used;
- `continuum.prior_width_rms`: half-width of coefficient priors in robust-RMS
  units;
- `noise.model`: `independent` (default) or `ar1`; and
- `noise.rho`: for `ar1`, either `auto` (estimate lag-one correlation from
  adjacent pixels in the configured continuum windows) or a fixed value from
  -0.95 through 0.95. The resolved coefficient is stored in each result; and
- `noise.marginal_scale`: a positive uncertainty multiplier or `auto`. Auto
  matches the robust continuum-residual scatter when formal uncertainties are
  underestimated and never reduces them. The applied factor is serialized;
- `kinematics.max_components`: largest narrow-component model;
- `kinematics.velocity_kms`: uniform velocity bounds;
- `kinematics.sigma_kms`: log-uniform dispersion bounds;
- `lines`: arbitrary rest wavelengths and optional fixed ratios;
- `wavelength_medium`: `air` or `vacuum` (MUSE AWAV products require `air`);
- `sampling.seed`: optional non-negative base seed. Batch runs derive and
  record a deterministic per-spectrum sampler seed from this value and the
  spectrum ID, independent of worker scheduling;
- `lsf.model`: `none`, `instrument`, `resolving_power`,
  `polynomial_resolving_power`, `nirspec_etc_resolving_power`,
  `constant_fwhm_angstrom`,
  `polynomial_fwhm_angstrom`, or `table`;
- `fixed_components`: always-present, fixed-kinematic, free-flux lines;
- `broad_components`: always-present permitted-line Gaussians with their own
  free velocity, log-uniform sigma, and free integrated flux;
- `flux_prior.min_snr`: automatic lower integrated-flux scale;
- `flux_prior.max_signal_factor`: automatic upper scale;
- `selection.delta_logz`: required evidence improvement;
- `selection.stop_when_not_preferred`: stop after the first rejected larger
  model (default true); and
- `sampling`: UltraNest `min_num_live_points`, `min_ess`, `dlogz`, `show_status`,
  `stepsampler` (`slice` or `none`), and optional `nsteps`.

Line `ratio` is defined as `F(this line) / F(ratio_to)`. The referenced line
must appear earlier in the list. Add positive `flux_bounds: [min, max]` to a
free line or fixed component to replace automatic bounds.

A broad component is independent of the shared narrow-line kinematics:

```yaml
broad_components:
  - name: broad_halpha
    line: halpha
    velocity_kms: [-1500, 1500]
    sigma_kms: [600, 3500]
```

To test whether broad emission is required, compare otherwise identical runs
with and without this block. Broad components are present in every tested
narrow-component-count model when configured.

LSF resolving power uses `value`. Constant FWHM uses `value` in Angstrom.
Polynomial FWHM uses an increasing-power `coefficients` list and optional
`reference_angstrom` and `scale_angstrom`. A two-column table uses `path` and
`quantity` (`fwhm_angstrom`, `sigma_angstrom`, or `resolving_power`). Tables may
not be extrapolated. Instrumental and intrinsic Gaussian sigmas are combined
in quadrature and both are reported alongside the convolved sigma.

## `output`

- `directory`: output root;
- `resume`: reuse matching atomic results;
- `plots`: `none`, `selected`, or `all`;
- `workers`: process count; and
- `progress_every`: reporting interval.
