# BEAT 2.0.0a1 — astronomer pilot

This is an alpha prerelease of BEAT (Bayesian Evidence Analysis Tool) for
astronomer evaluation. It supports bounded MUSE, JWST/NIRSpec IFU, and
JWST/MIRI MRS cubes, collections of 1D spectra, and survey-style FITS tables.

## Pilot package

Download `BEAT-2.0.0a1-astronomer-pilot-20260717.zip`. It contains:

- an installable pure-Python wheel;
- a matching source snapshot and 44 automated tests;
- conda and pip dependency files;
- a verified synthetic [O III] smoke test;
- portable MUSE, NIRSpec, MIRI, ASCII, and survey-table templates;
- configuration and instrument documentation;
- validation and release notes; and
- a structured pilot-feedback form.

Archive SHA-256:

`58f84f8ab261e78cf44c1335ea57c321f1687fab50f810ab2d675ee01abad166`

## Validation status

- NIRSpec standard core: 16/16 component counts recovered.
- Powered NIRSpec blank/single scope: 20/20 blanks and 20/20 S/N=10 singles
  recovered, with no evidence flags.
- The packaged-wheel smoke test selects the injected one-component model with
  accepted evidence status.
- All 44 packaged source tests pass.

## Important limitations

- This is alpha software; the configuration schema is not frozen.
- NIRSpec double/triple completeness is not yet powered to 20 cases per class.
- MIRI non-Gaussian/undersampled profile validation is incomplete.
- The representative MUSE LSF is not a dataset-matched calibration file.
- AR(1) is a first-order covariance approximation.
- Broad-line decompositions require scientific review.

Please read the included `README.md`, `RELEASE_NOTES.md`, and validation reports
before scientific use, and return pilot feedback using `PILOT_FEEDBACK.md`.
