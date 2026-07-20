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
The default likelihood assumes independent Gaussian pixel uncertainties.
For correlated products, an opt-in AR(1) residual-noise model and an
inflation-only empirical marginal-error calibration are available through the
configuration schema.

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

The bounded multi-instrument alpha work is summarized in
[`validation/ALPHA_TEST_REPORT.md`](validation/ALPHA_TEST_REPORT.md). The
powered MUSE and NIRSpec gates, bounded MIRI gate, nine-spaxel broad-H-alpha
study, and synthetic 1D operational regression are complete within their
documented supported domains. Version `2.0.0a2` is the frozen astronomer-pilot
build; it is not yet a scientifically certified final release.

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

The frozen powered MUSE red gate uses real residuals from NGC 2992 and
NGC 3393. In its documented supported domain it recovers 0/60 blank false
positives and 20/20 singles, doubles, and triples, with no evidence flags and
median absolute velocity, width, and flux errors of 4.91 km/s, 5.29%, and
4.41%. A 16-case H-beta+[O III] check exercises the wavelength-dependent LSF
near 5000 A. See
[`validation/MUSE_POWERED_VALIDATION.md`](validation/MUSE_POWERED_VALIDATION.md).

The NIRSpec validation uses bundled STScI G235H/G395H `R(lambda)` tables. Its
16-case standard core recovers every zero-through-three-component count. A
separate powered residual-noise gate exposed strong cube-pixel covariance;
with empirical marginal-error calibration and the new opt-in AR(1) likelihood,
it recovers 20/20 blanks and 20/20 S/N=10 singles with no evidence flags.
The powered double/triple expansion recovers 18/20 doubles and 17/20 triples
after targeted standard-profile audits, passing the predeclared count and
median-parameter-error thresholds. A paired weak-component experiment closes
the scoped G235H reliability gate with a conservative supported alpha domain:
effective component S/N>=10. All eight tested S/N=10 double/triple anchors are
correct and accepted under standard sampling; lower-S/N components are
explicitly exploratory. See
[`validation/NIRSPEC_G235H_WEAK_COMPONENT_BOUNDARY.md`](validation/NIRSPEC_G235H_WEAK_COMPONENT_BOUNDARY.md).

The bounded MIRI MRS gate selects the correct product in 24/24 checks spanning
all 12 sub-bands for both supplied targets. On actual 2A, 3B, and 4C wavelength
grids, it recovers all 24 reference blank/single/double counts and both
controlled non-Gaussian profile-mismatch cases, with no evidence flags. The
Gaussian-equivalent LSF is accepted only within the documented sampling,
S/N, separation, and profile-shape limits; see
[`validation/MIRI_BOUNDED_VALIDATION.md`](validation/MIRI_BOUNDED_VALIDATION.md).

The NGC 1365 broad-H-alpha gate compares narrow-only, one-broad, and two-broad
models across all nine nuclear spaxels under linear and quadratic continua.
One broad component is required in every spaxel, while a second Gaussian is
intermittent and fails the frozen spatial and continuum-stability criteria.
Routine fitting should use one broad component; detailed BLR work may use a
flexible asymmetric profile, without automatically assigning physical meaning
to multiple broad Gaussians. See
[`validation/NGC1365_BROAD_STABILITY.md`](validation/NGC1365_BROAD_STABILITY.md).

The production 1D survey path has also passed a 512-row synthetic SDSS-like
regression: 500 valid spectra reached one deterministic catalog, 12 malformed
rows were isolated without losing successful outputs, a four-worker run was
interrupted and resumed from 32 atomic checkpoints, and a one-worker override
left every successful hash unchanged. A small real H-beta+[O III] set
exercises selected component counts zero, one, and two. This is operational
validation, not real-SDSS scientific completeness; see
[`validation/SURVEY_1D_REGRESSION.md`](validation/SURVEY_1D_REGRESSION.md).
