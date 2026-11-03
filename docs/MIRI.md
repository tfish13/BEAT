# JWST/MIRI MRS adapter

The `miri` adapter accepts a glob covering the 12 calibrated Stage-3 MRS band
cubes: four channels, each with SHORT (A), MEDIUM (B), and LONG (C) bands. It
validates JWST MIRI `MIR_MRS` IFUCubeModel products and their SCI, ERR, and DQ
extensions.

## Segment selection

BEAT does not concatenate the 12 cubes by spatial pixel index. MIRI channel
fields of view, spaxel sizes, spatial array shapes, wavelength sampling, and
spectral resolution differ. Instead, BEAT transforms the configured fitting
window into the observed frame, inspects all matching cube coverages, and
opens only the segment that contains the complete window.

If two overlapping bands both contain the window, validation stops and asks
for `input.segment`. Accepted forms include `3B`, `3-medium`, and equivalent
channel/band combinations. If no segment contains the complete window, narrow
the fitting window or divide the analysis into separate scientifically
independent runs. The selected segment and selection method are recorded in
every result.

## Memory, masks, and units

Both spatial ranges are mandatory unless `allow_full_cube: true` is explicitly
set. SCI, ERR, and raw DQ storage remain memory-mapped. Only the `DO_NOT_USE`
DQ bit for each selected spaxel is materialized.

Vacuum wavelengths are converted from micron to Angstrom. SCI and ERR remain
in native `MJy/sr`, so integrated line results have units of
`MJy/sr Angstrom`. Results retain the channel-dependent `PIXAR_SR` and
`PIXAR_A2` values. Any conversion to summed aperture flux must explicitly use
the relevant channel's pixel area rather than one value for all 12 cubes.

## LSF

The default instrumental LSF uses STScI's approximate in-flight relation

`R(lambda) = 4603 - 128 lambda`,

with wavelength in micron, and combines its Gaussian sigma with the intrinsic
line sigma in quadrature. It is marked as an approximation and accompanied by
the pipeline version, CRDS context, and cube, distortion, fringe, and photom
reference provenance. A release-quality analysis may override it with a more
appropriate wavelength-dependent calibration table.

### Validated alpha domain

The bounded `miri-bounded-v1` gate passes product routing in all 12 MRS
sub-bands for both supplied targets and controlled recovery tests in 2A, 3B,
and 4C. Those bands span 1.76, 2.10, and 3.67 pixels per
Gaussian-equivalent LSF FWHM. The tested Gaussian approximation recovers all
24 reference component counts and both declared shifted-wing stress cases,
with no evidence flags.

This supports isolated narrow singles at effective peak S/N=15, resolved
singles at S/N=10, and resolved doubles separated by at least three
instrumental FWHM with the weaker component at S/N>=10. Intrinsic-width
recovery is validated only when intrinsic sigma is at least 0.75 times the
instrumental sigma. Weaker or closer components, intrinsically narrow widths
in undersampled regions, other sub-bands, and more asymmetric profiles require
review or a tabulated/empirical LSF. See
[`validation/MIRI_BOUNDED_VALIDATION.md`](../validation/MIRI_BOUNDED_VALIDATION.md).

## Examples

For diagnostic MIRI runs, `output.posteriors: all` saves compressed weighted
posterior draws for every fitted component count under `posteriors/`, with
their paths recorded in each result JSON. The default is `none`; `selected`
saves only the chosen count. These are output controls and do not change the
likelihood, priors, evidence, or component selection. Each result also reports
`component_quality_diagnostics` (peak relative to continuum-window residual
RMS, integrated flux relative to posterior standard deviation, and intrinsic
width-bound proximity). These fields are flags for review, not rejection
rules. Duplicate YAML keys are rejected at configuration load time.

An optional independent-seed check is available without changing the selected
component count:

```yaml
fit:
  sampling: {seed: 5728, min_num_live_points: 200, min_ess: 400, dlogz: 0.5}
  selection:
    stop_when_not_preferred: false
    convergence_audit:
      mode: flag
      models: [2, 3]
      repeats: 1  # one additional independent seed per requested count
      evidence_sigma: 3
      minimum_logz_difference: 2
      parameter_sigma: 3
      velocity_floor_kms: 10
      width_floor_kms: 10
      flux_floor_fraction: 0.1
output:
  posteriors: all
```

The audit compares ln Z using the quadrature sum of sampler-reported errors
and compares velocity/width/line-flux posteriors after uncertainty-aware component
matching. If either discrepancy exceeds its configured threshold, the
`convergence_audit` record flags it and `selection_status` becomes
`convergence_unverified`; `selected_components` is **not** changed. All runs,
seeds, comparisons, and (when requested) weighted posterior paths are retained.
`models: selected` or `all` may replace the explicit list. A fixed sampling
seed is required. This can substantially increase runtime and is intended for
targeted quality checks, not a default survey-wide setting.

The IC 5063 and NGC 4151 examples request a 5-by-5 region around a provisional
[Ne V] 14.32 micron window. Both should automatically select segment 3B. The
rest wavelength and continuum windows remain pilot choices to confirm before
scientific fitting.
