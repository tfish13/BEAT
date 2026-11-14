# ALMA frequency-cube adapter (initial support)

`input.adapter: alma` reads a standard ALMA image cube with a linear `FREQ`
axis in Hz, `Jy/beam` data, `RESTFRQ`/`RESTFREQ`, and either 3-D data or a 4-D
cube with one `STOKES` plane. The adapter selects that plane without copying
the full cube, skips non-finite spaxels, and records the rest frequency,
spectral frame, channel spacing, and beam metadata. Use bounded `x_range` and
`y_range` for initial tests; a full 540×540 cube is costly.

The user-facing line and fitting-window fields can be in GHz:

```yaml
version: 2
input:
  adapter: alma
  path: /path/to/target.cube.I.pbcor.fits
  redshift: 0.00935  # replace with your adopted systemic redshift
  x_range: [300, 302]
  y_range: [270, 272]
output:
  directory: alma_co32_output
fit:
  frequency_window_ghz: [342.0, 343.2]  # observed frequency, increasing
  continuum:
    degree: 1
    frequency_windows_ghz:
      - [342.0, 342.2]
      - [343.0, 343.2]
  lines:
    - {name: CO32, rest_frequency_ghz: 345.79599}
  kinematics:
    max_components: 2
    velocity_kms: [-400, 400]
    sigma_kms: [30, 350]
```

`fit.exclude_frequency_windows_ghz` is also available. Do not mix a GHz
setting with its Ångström counterpart (`window`, `continuum.windows`, or
`exclude_windows`). GHz intervals are increasing even though the internally
converted vacuum-wavelength intervals reverse their endpoints. `fit.frame`
defaults to `observed` for ALMA and must remain so for these frequency
windows. The configured line frequencies are **rest** frequencies. Verify the
redshift and line-free windows for each target before interpreting results.

BEAT retains its existing wavelength-based model internally: it converts each
frequency channel exactly using λ=c/ν and sorts the pixels by wavelength.
The data ordinate stays in `Jy/beam`; the top plot axis shows the corresponding
observed frequency in GHz. The raw `flux` parameter is therefore an integral
in `Jy/beam Å`, not `Jy km/s`. ALMA catalogue rows additionally give
`flux_jy_beam_kms` (and its posterior standard deviation), converted locally
at the fitted line centroid using the **radio velocity convention**. This is
the narrow-line approximation to the integral in `Jy/beam km/s`; it does not
convert a beam measurement to a total-source flux.

The default ALMA LSF is `none` because the effective channel response depends
on the observation and reduction (including possible Hanning smoothing), and
the cube header alone is not a calibrated LSF. BEAT currently evaluates the
Gaussian model at channel centers rather than integrating it over channels.
Do not infer intrinsic widths near the channel resolution until the effective
spectral response has been supplied and validated. The optional AR(1) noise
model is not yet validated for this cube's inter-channel covariance, and BEAT
does not model primary-beam uncertainty. A primary-beam-corrected
cube can have strongly position-dependent noise: inspect the continuum RMS
and use genuinely line-free windows. This initial ALMA path has an ingestion
and one-spaxel smoke test, not a powered ALMA science-validation gate.
