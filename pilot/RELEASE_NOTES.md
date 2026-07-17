# BEAT 2.0.0a1 pilot release notes

Release date: 2026-07-17

## Included

- MUSE Phase-3, JWST/NIRSpec Stage-3 IFU, and JWST/MIRI MRS adapters;
- bounded cube processing and resumable multi-spectrum pipelines;
- FITS survey-table, FITS 1D, and ASCII spectrum inputs;
- zero through three shared narrow kinematic components;
- tied line ratios and independent broad permitted-line components;
- UltraNest evidence selection with ambiguity/convergence diagnostics;
- Gaussian, resolving-power, polynomial, tabulated, and instrument LSF modes;
- bundled STScI G235H/G395H wavelength-dependent resolution tables;
- empirical marginal uncertainty calibration and AR(1) correlated-noise mode;
- fit-window exclusions for known contaminating lines or artifacts; and
- machine-readable results, catalogues, manifests, and diagnostic plots.

## Validation status

- NIRSpec standard core: 16/16 component counts recovered.
- NIRSpec powered blank/single scope: 40/40 recovered, including 0/20 false
  blank detections and 20/20 S/N=10 singles.
- The worst NIRSpec evidence margins remain correct under tighter sampling.
- MUSE and broad H-alpha pilots establish functionality but retain unresolved
  three-component completeness and physical-profile questions.

## Known limitations

- The configuration schema is not frozen until `2.0.0rc1`.
- NIRSpec double/triple completeness is not yet powered to 20 cases per class.
- MIRI's Gaussian-equivalent resolution model does not yet validate
  non-Gaussian or undersampled line profiles across all bands.
- The representative MUSE LSF is not a dataset-matched `LSF_PROFILE`.
- AR(1) is a first-order approximation, not a full spectral covariance matrix.
- Broad-line decompositions can absorb continuum or calibration structure and
  require scientific review.
- Full cubes can be computationally expensive; pilot configurations must use
  bounded spatial regions.

## Feedback priority

The most valuable pilot results are real-world adapter compatibility,
configuration friction, unexpected component selections, runtime scaling,
and residual structures that are not captured by the current LSF or noise
models.
