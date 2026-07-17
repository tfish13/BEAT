# BEAT v2 alpha

BEAT (Bayesian Evidence Analysis Tool) fits zero or more shared kinematic
Gaussian components to astronomical emission lines and selects the component
count with Bayesian evidence from UltraNest.

This v2 alpha separates the scientific model from the input dataset. The same
fitter can consume:

- a three-dimensional FITS datacube for one target;
- a FITS binary table containing one vector spectrum per survey target; or
- a directory of FITS or ASCII one-dimensional spectrum files.

The package is a development foundation for validation before a public
astronomy release. It should not yet be treated as a scientifically certified
release; see [the audit and validation roadmap](docs/V1_AUDIT_AND_ROADMAP.md).

## Install and run

Use Python 3.9 or newer in a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e .
beat validate examples/cube.yaml
beat run examples/cube.yaml --workers 4
```

`beat validate` opens the first spectrum, applies the fit window, checks its
uncertainties, and reports the automatically scaled flux prior. Run it before
every new dataset. `--config-only` validates without opening the input file.

Copy the nearest configuration in [`examples/`](examples/) and edit it; users
do not edit Python code.

## Scientific model

For every target or spaxel, BEAT tries models containing 0 through
`max_components` narrow kinematic components. Each component has:

- one velocity offset and one Gaussian velocity dispersion in km/s;
- one integrated flux for each free line; and
- any number of lines tied to earlier lines by a fixed flux ratio.

Centroids and widths are transformed to observed Angstrom independently for
each line using the spectrum's redshift. A polynomial continuum of degree 0,
1, or 2 is fit simultaneously and is present in every compared model. The
likelihood currently assumes independent Gaussian pixel uncertainties.

The next component is selected when its evidence exceeds the best simpler
model by `selection.delta_logz`. Evidence depends on the prior volume, so the
velocity, width, continuum, and flux priors are scientific configuration—not
mere optimizer settings.

Every fit also records an evidence-reliability status. The default `flag` mode
marks threshold-overlap cases as `ambiguous` and moderately supported
maximum-component selections as `convergence_unverified` without increasing
runtime. Optional `rerun` mode repeats only triggered spectra with a tighter
sampling profile and retains both pilot and audit evidence. See
[`docs/ADAPTIVE_SELECTION.md`](docs/ADAPTIVE_SELECTION.md).

The default integrated-flux prior is scaled per spectrum from the supplied or
continuum-estimated noise and the observed signal. Absolute `flux_bounds` may
instead be placed on any free line when a survey needs one fixed prior across
all targets.

```yaml
lines:
  - name: oiii5007
    wavelength: 5006.84
    # flux_bounds: [1.0e-18, 1.0e-12]
  - name: oiii4959
    wavelength: 4958.92
    ratio_to: oiii5007
    ratio: 0.33557  # F(4959) / F(5007)
```

Always state whether wavelengths are air or vacuum and use one convention for
both the data calibration and configured rest wavelengths.

## Inputs

### Datacube

Set `input.kind: cube`, the flux and uncertainty HDUs, and the NumPy spectral
axis. BEAT supports sigma, variance, and inverse-variance arrays, an optional
2-D or 3-D mask, a 1-D wavelength HDU, or a standard linear FITS wavelength
axis using `CRVALn`, `CRPIXn`, and `CDELTn`/`CDn_n`.

Spatial `x_range` and `y_range` values are zero-based, half-open Python ranges.
For example, `[80, 90]` fits pixels 80 through 89. Results record both these
coordinates and one-based FITS coordinates.

For ESO Phase-3 MUSE products, set `input.adapter: muse`. The adapter supplies
the DATA/STAT layout and unit defaults, validates the product, records LSF
provenance, and requires a bounded spatial subset unless `allow_full_cube` is
explicitly enabled. See [the MUSE guide](docs/MUSE.md).

For JWST Stage-3 NIRSpec IFU cubes, set `input.adapter: nirspec`. It handles
SCI/ERR/DQ, vacuum micron wavelength axes (including WCS-TABLE products),
native MJy/sr calibration, and bundled wavelength-dependent STScI G235H/G395H
resolution tables. See the
[NIRSpec guide](docs/NIRSPEC.md).

For JWST MIRI MRS observations, `input.adapter: miri` accepts a glob spanning
the 12 channel/band cubes and selects the one containing the complete observed
fit window. It preserves the distinct spatial grids and refuses ambiguous
band overlaps without an explicit segment choice. See the
[MIRI guide](docs/MIRI.md).

### Survey FITS table

Set `input.kind: survey_table`. Each table row represents one target and the
wavelength, flux, uncertainty, and optional mask columns contain arrays.
Redshift can come from a scalar column per row or one fixed configuration
value. The parent process streams rows and keeps only a bounded number of
spectra queued to workers.

### One spectrum per file

Set `input.kind: spectrum_files`, plus `glob` or `files`. Supported formats are
`fits_table`, `fits_image`, and `ascii`. A collection currently uses one fixed
redshift or a FITS header keyword per file. For a heterogeneous ASCII survey,
convert it to the survey-table layout or split it into redshift-homogeneous
runs.

See [the full configuration reference](docs/CONFIG_REFERENCE.md).

## Output and resuming

BEAT writes:

- `results/<shard>/<spectrum-id>-<hash>.json`: one atomic checkpoint per
  spectrum (hash-sharded so survey directories remain manageable);
- `catalog.csv`: selected parameters, posterior standard deviations, and all
  attempted evidence values;
- `plots/`: selected or all model diagnostics when enabled;
- `failures.csv`: per-spectrum errors that do not abort the survey; and
- `run_manifest.json`: expanded configuration, version, platform, timestamps,
  and a configuration fingerprint.

With `resume: true`, a result is skipped only when both its spectrum ID and
configuration fingerprint match. Worker processes never write a shared
catalogue.

## Fixed-kinematic components

A broad or otherwise pre-characterized line can be present in every evidence
model with fixed velocity and width but free flux:

```yaml
fixed_components:
  - name: broad_hbeta
    line: hbeta
    velocity_kms: 0.0
    sigma_kms: 1800.0
    # flux_bounds: [1.0e-17, 1.0e-12]
```

This is the v2 replacement for v1.1 `prefit_instructions`. It deliberately
does not pretend that fixed kinematics have zero uncertainty; use this feature
only when that conditional model is scientifically intended.

## Free broad permitted-line components

Broad Balmer or Paschen emission can instead have kinematics independent of
the narrow-line family:

```yaml
broad_components:
  - name: broad_halpha
    line: halpha
    velocity_kms: [-1500, 1500]
    sigma_kms: [600, 3500]
```

Velocity is uniform over the configured bounds; sigma and integrated flux use
log-uniform priors. Compare otherwise identical configurations with and
without the broad block when broad-component presence is not known in advance.

## Development checks

```bash
python -m unittest discover -s tests -v
```

Before publishing, expand the instrument-specific regression and
injection/recovery grids, add posterior predictive checks, and compare with
independent fitters on representative data.

The first two-stage MUSE injection/recovery framework and its recorded
26-spectrum standard grid are now available; see
[the injection/recovery guide](docs/INJECTION_RECOVERY.md). The current result
is a beta validation baseline, not final scientific certification.

The bounded multi-instrument alpha run, H-alpha+[N II] injection pilot, and
NGC 1365 broad-H-alpha comparison are summarized in
[`validation/ALPHA_TEST_REPORT.md`](validation/ALPHA_TEST_REPORT.md). That
report records passed adapter/model checks and the scientific gates that still
fail; the project remains `2.0.0a1`.

The subsequent NGC 2992 audit found that the apparent H-alpha false component
used a residual donor contaminated by native line wings. Injection utilities
now reject residual donors whose robust scatter is more than three times their
supplied uncertainty. The corrected initial 12-case grid passes all
provisional gates, but remains too small for release-level rate claims; see
[`validation/NGC2992_FALSE_POSITIVE_AUDIT.md`](validation/NGC2992_FALSE_POSITIVE_AUDIT.md).

The subsequent zero-through-three-component H-alpha pilot confirms that the
generic three-component engine runs without a core rewrite, but scientific
selection is not yet accepted: it recovers 7/8 counts and fails the
provisional triple-recovery and median velocity-error gates. See
[`validation/HALPHA_THREE_COMPONENT_PILOT.md`](validation/HALPHA_THREE_COMPONENT_PILOT.md).

Controlled Gaussian-noise anchors further show that the tested equal triple
is identifiable at 300 km/s adjacent spacing, whereas a moderate unequal
triple at the same spacing is not stable under tight evidence; all tested
patterns recover three at 400 km/s. These conditional boundaries are reported
in [`validation/HALPHA_CONTROLLED_TRIPLE_CALIBRATION.md`](validation/HALPHA_CONTROLLED_TRIPLE_CALIBRATION.md).

The NIRSpec validation uses bundled STScI G235H/G395H `R(lambda)` tables. Its
16-case standard core recovers every zero-through-three-component count. A
separate powered residual-noise gate exposed strong cube-pixel covariance;
with empirical marginal-error calibration and the new opt-in AR(1) likelihood,
it recovers 20/20 blanks and 20/20 S/N=10 singles with no evidence flags.
Double/triple completeness is still underpowered; see
[`validation/NIRSPEC_INJECTION_RECOVERY.md`](validation/NIRSPEC_INJECTION_RECOVERY.md).
