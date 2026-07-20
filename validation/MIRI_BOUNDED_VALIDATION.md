# Bounded JWST/MIRI MRS validation

Date: 2026-07-18  
BEAT version: `2.0.0a1`  
Matrix version: `miri-bounded-v1`

## Outcome

Every predeclared bounded MIRI alpha gate passes. Automatic wavelength routing
selected the intended MRS cube in 24/24 checks across all 12 sub-bands for
both IC 5063 and NGC 4151. The 24-case reference injection matrix recovered
every blank, single, and double count with no evidence-reliability flags. Both
controlled non-Gaussian profile-mismatch cases also recovered the correct
single-component model without flags.

| Gate | Result | Threshold | Status |
|---|---:|---:|---|
| Segment selection | 24/24 | 24/24 | pass |
| Blank false positives | 0/6 | 0 | pass |
| Single-component recovery | 12/12 | 12/12 | pass |
| Double-component recovery | 6/6 | at least 5/6 | pass |
| Reference evidence flags | 0/24 | 0 | pass |
| Median velocity error | 0.030 resolution FWHM | <=0.25 | pass |
| Median integrated-flux error | 3.61% | <=15% | pass |
| Median resolved-width error | 11.95% | <=25% | pass |
| Profile-mismatch count | 2/2 | 2/2 | pass |
| Mismatch velocity error | 0.038 resolution FWHM | <=0.25 | pass |
| Mismatch integrated-flux error | 9.86% | <=15% | pass |
| Mismatch intrinsic-width error | 7.14% | <=30% | pass |

No standard-profile audits were triggered because no screening result had an
incorrect count or reliability flag. The summed sampler runtime was 505.2
seconds.

## Bounded matrix

The reference grid uses the actual spectral wavelength arrays and adapter/LSF
metadata from the supplied Stage-3 products. It tests both targets in:

- 2A, with 1.76 pixels per Gaussian-equivalent LSF FWHM (`undersampled`);
- 3B, with 2.10 pixels per FWHM (`borderline`); and
- 4C, with 3.67 pixels per FWHM (`well sampled`).

Each target/segment combination contains one blank, an intrinsic-narrow
S/N=15 pixel-centered single, a resolved S/N=10 half-pixel-centered single,
and a resolved double separated by three resolution FWHM. The double's
components have effective peak S/N=15 and 10. Independent Gaussian noise is
used so the experiment isolates sampling, component selection, and LSF-model
effects rather than the covariance or residual structure of a particular
spaxel.

The two mismatch cases are IC 5063 2A and 4C resolved S/N=15 singles. Their
injected instrumental profile contains an 85% nominal Gaussian core and a
15% wing with 1.8 times the nominal LSF sigma, shifted by +0.5 nominal LSF
sigma. This is a controlled stress profile, not an empirical measurement of
the on-orbit MRS line-spread function.

## Supported alpha statement

Within the tested 2A, 3B, and 4C grids, the Gaussian-equivalent MIRI LSF is
adequate for component selection under these conditions:

- isolated narrow singles at effective peak S/N=15;
- resolved singles at effective peak S/N=10;
- two resolved components separated by at least three instrumental FWHM,
  with the weaker component at effective peak S/N>=10; and
- profiles no more extreme than the declared 15% shifted-wing stress test.

The aggregate velocity and flux recovery gates pass. Intrinsic-width recovery
is supported only for components with intrinsic sigma at least 0.75 times the
Gaussian-equivalent instrumental sigma. Widths of intrinsically narrow lines
in undersampled regions are not validated. Weak components below effective
S/N=10, closer blends, other sub-bands, and more asymmetric profiles remain
review-required rather than accepted by extrapolation.

The other nine MRS sub-bands pass automatic product selection only; their
scientific recovery was intentionally not expanded into a powered campaign.
Actual MRS profile residuals should still be inspected for scientific work.
If a line is narrower, weaker, or more asymmetric than this supported domain,
use a tabulated/empirical LSF when available or treat its decomposition as
diagnostic.

## Reproducibility and freeze

- Predeclared design: `validation/MIRI_BOUNDED_GATE_PLAN.md`
- Runner: `validation/run_miri_bounded_validation.py`
- Results: `validation/miri_bounded_validation/`
- Frozen hashes: `validation/miri_bounded_validation/frozen_gate_manifest.json`

The frozen manifest covers all 26 case checkpoints, the 24-selection audit,
the matrix definition, scores, and summary. Do not add cases unless a
predeclared gate later fails or an unexplained evidence flag appears inside
the supported domain. The complete current suite passed 57/57 tests on
2026-07-18; release closeout must record the then-current final count rather
than treating 57 as permanently fixed.
