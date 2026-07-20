# NIRSpec G235H/G395H injection/recovery validation

Date: 2026-07-17  
BEAT version: `2.0.0a1`

## Outcome

The NIRSpec gate now has three distinct results. First, the complete 16-case
standard-profile core matrix recovers all zero-, one-, two-, and
three-component truth cases, with all selections accepted. Its median absolute
errors are 6.38 km/s in velocity, 8.25% in intrinsic width, and 6.34% in
integrated flux. Four cases per count class are scientifically encouraging but
not statistically powered.

Second, a 20-case boundary matrix identifies conditional limits. Equal pairs
separated by 150 km/s remain single components under tight sampling in both
G235H and G395H. A G235H equal triple with 300 km/s adjacent spacing requires
the standard sampler to recover three components. S/N=5 and sigma=250 km/s
singles, 300 km/s 1:0.25 pairs, and 400 km/s 1:0.5:0.25 triples pass in both
gratings. These are calibration anchors, not universal completeness limits.

Third, a powered blank/single residual-noise matrix exposed and corrected a
likelihood misspecification. Treating resampled cube pixels as independent
gave six false positives in 16 blanks and recovered only 13/16 singles. Merely
calibrating the marginal uncertainty improved this to 15/16 blanks and 14/16
singles but still failed the predeclared gates. The residuals have strong
short-range covariance, especially in G235H. An explicit AR(1) likelihood,
with rho estimated only from configured line-free continuum windows, yields:

| Scoped gate | Result | Threshold | Status |
|---|---:|---:|---|
| Blank false-positive rate | 0/20 (0%) | <=5% | pass |
| S/N=10 single recovery | 20/20 (100%) | >=90% | pass |
| Ambiguous/unverified selections | 0/40 | 0 | pass |
| Median absolute velocity error | 5.28 km/s | <=10 km/s | pass |
| Median absolute width error | 10.04% | <=15% | pass |
| Median absolute flux error | 6.46% | <=10% | pass |

This is a statistically powered pass for the blank and single classes only.
There are no double or triple cases in that 40-case matrix.

Fourth, the powered double/triple expansion contains 20 cases in each class,
balanced across G235H and G395H, primary and alternate residual donors, 300--
500 km/s separations, and equal through 1:0.25 component ratios. The expanded
screening profile recovered 16/20 doubles and 11/20 triples. Every incorrect
or non-accepted selection was then rerun with the standard profile, giving the
audited result:

| Powered double/triple gate | Audited result | Threshold | Status |
|---|---:|---:|---|
| Double recovery | 18/20 (90%) | >=80% | pass |
| Triple recovery | 17/20 (85%) | >=70% | pass |
| Ambiguous/unverified selections | 3/40 | 0 | **fail** |
| Median absolute velocity error | 6.00 km/s | <=10 km/s | pass |
| Median absolute width error | 8.95% | <=15% | pass |
| Median absolute flux error | 6.19% | <=10% | pass |

The component-count and parameter-accuracy gates pass, but the predeclared
overall gate remains failed because three standard-profile evidence
comparisons overlap the selection threshold. The threshold was not relaxed
after seeing the results. Final statuses are 23 `accepted`, 14
`accepted_after_audit`, and 3 `ambiguous`.

All 20 G395H counts are correct after audit. The five remaining count failures
are all G235H weak-component cases: the alternate donor loses a 1:0.25 double
at both 300 and 400 km/s; the primary donor loses 1:0.5:0.25 triples at 400
and 500 km/s; and the alternate donor loses the weak triple at 500 km/s. A
correctly recovered alternate-donor weak triple at 400 km/s is nevertheless
ambiguous (`delta ln Z=5.44`), while the corresponding primary case is an
ambiguous two-component selection (`delta ln Z=4.08`). The 300 km/s
alternate-donor 1:0.25 double is also ambiguous (`delta ln Z=4.59`). These are
conditional completeness/reliability boundaries, not software exceptions.

Screening required 2169 seconds (36.2 minutes); the 17 targeted standard
audits required 7284 seconds (121.4 minutes). This supports the intended
operating pattern: screen large samples, then tightly rerun only flagged or
scientifically critical spectra.

## Correlated-noise audit

The marginal uncertainty calibration never reduces an uncertainty. It matches
the robust scatter of each block-resampled residual realization and records
the applied factor. G235H factors span 1.31-2.41 and fitted AR(1) coefficients
span 0.52-0.81. G395H is less correlated: its fitted rho spans 0.07-0.46, and
many realizations require no marginal inflation.

Before AR(1), tight reruns confirm a persistent G235H false blank and a false
second component in an injected single. After AR(1), the weakest evidence
margins are audited with the standard profile: the blank remains zero and the
single remains one, both with accepted status. The audited single has a
13.1 km/s velocity error and a 27.8% width error. Individual parameter outliers
therefore remain possible even though the powered-set medians pass.

The AR(1) option is explicit (`fit.noise.model: ar1`, `rho: auto`) and its
estimated coefficient is serialized in every result. It is a first-order
covariance model, not a full treatment of wavelength-dependent cube covariance.

## LSF calibration

BEAT bundles the G235H and G395H resolution tables published by STScI for the
JWST ETC. They assume a 2.2-pixel fully illuminated aperture and are converted
to a Gaussian-equivalent width using `FWHM=lambda/R`.

| Dataset/line | Observed wavelength | Tabulated R | Gaussian sigma |
|---|---:|---:|---:|
| NGC 4151 [Si VI] | 19699.40 A | 2240.45 | 3.734 A |
| IC 5063 Br-alpha | 40982.55 A | 2811.60 | 6.190 A |

This nominal curve is not an empirical, spaxel-specific LSF. Illumination,
source extent, resampling, and non-Gaussian profile structure remain caveats.

## Donors and reproducibility

Science and noise donors are distinct. NGC 4151 uses primary `(28, 38)` and
alternate `(27, 38)` donors; IC 5063 uses `(41, 55)` and `(42, 54)`.
Coordinates are zero-based. The NGC 4151 nuclear science spaxel was excluded
as a donor because its masked residual scatter is 7.6 times its formal
uncertainty.

- Runner: `validation/run_nirspec_injection_recovery.py`
- Standard core: `validation/nirspec_injection_recovery_standard/`
- Boundary matrix and tight audit: `validation/nirspec_boundary_expansion/`
  and `validation/nirspec_boundary_audit/`
- Uncalibrated donor control: `validation/nirspec_donor_replicates/`
- Marginally calibrated control and tight audit:
  `validation/nirspec_donor_replicates_calibrated/` and
  `validation/nirspec_donor_replicates_calibrated_audit/`
- Final AR(1) matrix and worst-margin audit:
  `validation/nirspec_donor_replicates_ar1/` and
  `validation/nirspec_donor_replicates_ar1_audit/`
- Powered double/triple matrix and targeted standard audit:
  `validation/nirspec_powered_components_ar1/` and
  `validation/nirspec_powered_components_ar1_audit/`. The reproducible merged
  result is `validation/nirspec_powered_components_ar1/audited_summary.json`.
- STScI source: [NIRSpec dispersers and filters](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-dispersers-and-filters)

The powered NIRSpec count-accuracy gate is met. A paired 32-case G235H boundary
experiment subsequently showed that recovery at weak-component S/N=3.75--6.25
is donor- and geometry-dependent. The supported alpha domain therefore
requires effective component S/N>=10. All eight tested S/N=10 double/triple
anchors recover the correct count with accepted standard-profile evidence and
no unexplained flags. Automatic rerun was exercised on all three earlier
threshold-overlap cases: two flags resolved and the remaining S/N=5 triple is
explicitly outside the supported domain. See
`NIRSPEC_G235H_WEAK_COMPONENT_BOUNDARY.md`. MIRI profile-shape validation
remains separate.
