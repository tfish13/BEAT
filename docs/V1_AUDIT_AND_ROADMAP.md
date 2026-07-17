# v1.1 audit and release roadmap

The supplied `beat_v1.1.py` establishes the Bayesian-evidence workflow, but it
contains dataset assumptions and correctness risks that should be removed
before distribution.

## High-impact v1.1 findings

- The reader assumes flux in HDU 1, uncertainty in HDU 2, a three-dimensional
  `(wavelength, y, x)` array, and only `CRVAL3` plus `CD3_3`. It ignores
  `CRPIX3`, other uncertainty conventions, masks, non-Angstrom axes, and 1-D
  survey layouts.
- Fitting and continuum windows mix wavelength indices with observed-frame
  wavelengths. A configuration can silently select the wrong pixels.
- All lines use the reference line's wavelength bounds; `minwave` and
  `wave_range` on H-beta in the example are not used.
- A shared width in Angstrom does not represent shared velocity dispersion for
  separated transitions. v2 uses km/s and converts each line independently.
- Prefit model construction consumes one parameter per component while its
  prior and parameter count allocate three. Result-writing indices are also
  inconsistent. This can attach parameters to the wrong physical quantities.
- Locked fluxes depend on dictionary order and on earlier entries. v2 validates
  the dependency and defines ratio direction explicitly.
- Independent workers read and append the same CSV file, which permits lost
  updates or corruption. The in-memory catalogue is never safely reconciled.
- `imap_unordered` results are not consumed, so worker exceptions may be hidden
  from the main process.
- Result columns labelled as posterior uncertainties are created but never
  populated; maximum-likelihood points are written instead of the documented
  posterior median.
- The residual early-stop test runs only before component 2 and divides by one
  scalar continuum RMS rather than the supplied pixel uncertainties.
- Global arrays make the likelihood hard to test, reuse for 1-D inputs, or run
  safely outside the original fork-based multiprocessing assumptions.
- Invalid or zero uncertainties, too-short continuum windows, duplicate
  wavelengths, decreasing axes, and insufficient fit pixels are not validated
  before nested sampling.

The v2 alpha addresses these structural items with a `Spectrum` interface,
strict configuration checks, physical velocity parameters, a simultaneous
continuum, per-spectrum immutable model state, explicit failure reporting, and
atomic checkpoints.

## Scientific validation required before a public release

1. Expand the initial two-stage MUSE injection/recovery grid across additional
   sampling, line ratios, continuum shapes, redshifts, spaxels, and noise seeds.
   The recorded 26-case baseline now measures component-count confusion, bias,
   and posterior-stdev coverage; see `docs/INJECTION_RECOVERY.md`.
2. Verify evidence stability against prior widths, live-point counts,
   `dlogz`, and step-sampler settings. Publish recommended settings and a
   convergence policy.
3. Replace approximation LSFs with calibrated wavelength- and, where needed,
   spaxel-dependent products. BEAT now convolves intrinsic widths with Gaussian
   instrument LSFs and reports intrinsic, instrumental, and convolved widths;
   MUSE, NIRSpec, and MIRI adapters currently label approximation modes.
4. Decide how to handle resampling covariance, sky residuals, bad-pixel bit
   masks, and underestimated variances. The current likelihood is diagonal
   Gaussian.
5. Test air-versus-vacuum wavelength conventions and relativistic velocity
   definitions against each supported instrument pipeline.
6. Compare posterior summaries and evidence decisions with at least one
   independent fitter on MUSE/KMOS cubes and a representative 1-D survey.
7. Add spatial strategies for cubes (binning, neighboring-spaxel priors, or
   post-fit regularization) as optional layers; independent spaxels should
   remain available as the reproducible baseline.
8. Add provenance tests, API documentation, citation metadata, a code of
   conduct, contributor guide, versioned example data, and continuous
   integration across supported Python versions.

## Suggested release sequence

- `2.0.0a1`: architecture and internal synthetic tests (this deliverable);
- `2.0.0b1`: instrument adapters, LSF support, injection/recovery report, and
  astronomer pilot feedback;
- `2.0.0rc1`: frozen configuration schema, reproducibility audit, tutorials,
  packaging, and archival DOI; and
- `2.0.0`: only after scientific acceptance criteria and regression datasets
  are documented and passing.

Current beta status: MUSE, NIRSpec, and MIRI adapters, Gaussian LSF support,
the initial injection/recovery report, and tabulated G235H/G395H resolution
curves are implemented. The corrected MUSE donor grids supersede the original
contaminated case. NIRSpec's standard core recovers 16/16 counts, and its
AR(1)-aware residual-noise matrix passes the powered blank/single scope at
40/40. Powered NIRSpec double/triple recovery, MIRI profile validation, and
astronomer pilot feedback remain before tagging `2.0.0b1`.
