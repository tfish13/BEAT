# Powered MUSE narrow-line injection/recovery validation

Date: 2026-07-17  
BEAT version: `2.0.0a1`  
Matrix version: `muse-powered-v1`

## Outcome

The powered MUSE red gate passes every predeclared criterion in its documented
supported domain. The frozen gate contains 120 H-alpha+[N II] real-residual
cases balanced across NGC 2992 and NGC 3393 donors.

| Supported red gate | Result | Threshold | Status |
|---|---:|---:|---|
| Blank false positives | 0/60 (0%) | <=5% | pass |
| One-sided 95% upper false-positive bound | 4.87% | <=5% | pass |
| Single-component recovery | 20/20 (100%) | >=90% | pass |
| Double-component recovery | 20/20 (100%) | >=80% | pass |
| Triple-component recovery | 20/20 (100%) | >=70% | pass |
| In-domain evidence flags | 0/120 | 0 | pass |
| Median absolute velocity error | 4.91 km/s | <=10 km/s | pass |
| Median absolute width error | 5.29% | <=15% | pass |
| Median absolute flux error | 4.41% | <=10% | pass |

This resolves the earlier pilot failures. Triple recovery improves from 1/2
at the old mixed 300 km/s boundary to 20/20 in the supported 400--500 km/s
domain. The median absolute velocity error improves from 11.20 to 4.91 km/s.
The evidence and parameter thresholds were not changed after inspecting the
powered results.

The 20 cases per nonblank class meet the predeclared alpha minimum and support
the reported point-estimate gates. They do not make completeness exactly 100%
for the underlying population. The blank class was expanded to 60 because
0/20 false positives would not place a one-sided 95% upper bound below 5%; at
0/60 the exact upper bound is 4.87%.

## Supported red domain

The result applies to the tested MUSE WFM-like wavelength grid, empirical
Gaussian LSF approximation, real-residual noise treatment, priors, and line
configurations below.

- H-alpha+[N II] with independent H-alpha and [N II] 6583 fluxes per
  kinematic component and fixed [N II] 6548/6583 ratio 0.335.
- AR(1) likelihood with `rho: auto`, automatic marginal scaling, and
  block-resampled line-free residuals.
- Singles with effective H-alpha peak S/N=10--15 and intrinsic sigma=50--180
  km/s.
- Doubles in four tested families: 300 km/s/equal flux/sigma=(60,80) km/s;
  400 km/s/ratio=0.67/sigma=(80,120) km/s; 500 km/s/ratio=0.5/
  sigma=(120,160) km/s; and 400 km/s/equal flux/sigma=(100,60) km/s. The
  weakest component has effective peak S/N>=10.
- Triples with adjacent spacing 400--500 km/s, flux patterns from
  1:0.75:0.5 through 1:0.5:0.4, intrinsic sigma=60--160 km/s, and weakest
  effective peak S/N>=10.
- Two donors per class split evenly between NGC 2992 `(96,144)` and NGC 3393
  `(224,96)`; coordinates are zero-based.

The representative LSF is the UDF-10 WFM polynomial. Its Gaussian FWHM is
approximately 2.54 A at H-alpha/[N II]. It remains an approximation rather
than a dataset-matched `LSF_PROFILE`.

## Mapped red boundary

The complete red boundary matrix contains 124 cases. Four cases exercise a
300 km/s, 0.5-ratio double with intrinsic sigma=(160,80) km/s. Both NGC 2992
realizations recover two components, but both NGC 3393 realizations underfit
as one; one of the latter remains evidence-ambiguous after standard sampling.

The entire geometry, including its successful donor realizations, is excluded
from the supported domain. This exclusion was applied as a scientific
identifiability boundary, not by removing only failed noise realizations. The
full matrix therefore reports 22/24 doubles and one explained flag, while the
frozen supported subset reports 20/20 doubles and no flags.

## Blue-wavelength LSF check

Sixteen additional H-beta+[O III] real-residual cases exercise the same LSF
polynomial near 4861--5007 A, where its Gaussian FWHM is 2.96--2.91 A rather
than the 2.54 A red value.

| Blue diagnostic | Result |
|---|---:|
| Blanks | 4/4 correct |
| Singles | 4/4 correct |
| Doubles | 4/4 correct |
| Triples | 4/4 correct |
| Median absolute velocity error | 4.38 km/s |
| Median absolute width error | 6.63% |
| Median absolute flux error | 3.68% |

One correctly counted NGC 3393 400 km/s triple remains
`convergence_unverified` because the maximum tested component count was
selected without decisive evidence. This flag is explained and retained. The
blue matrix is an LSF/wavelength diagnostic, not a powered per-class
acceptance claim.

## Sampling and reproducibility

All cases first use the screening profile (`40` live points, `40` ESS,
`dlogz=5`). A standard independent rerun (`100` live points, `200` ESS,
`dlogz=0.5`, 20 slice steps) is triggered when the screening count differs
from truth or its evidence diagnostics are flagged. Validation truth triggers
are used only to decide which synthetic case receives an audit; they do not
force a component count or alter the evidence threshold.

Fifteen of 124 red cases and three of 16 blue cases required standard audits.
The summed sampler runtime is 21,188 seconds (5.89 CPU-hours); the two donor
streams were run concurrently with one large cube per process.

- Runner: `validation/run_muse_powered_validation.py`
- Checkpoints and aggregate products: `validation/muse_powered_validation/`
- Frozen supported-case hashes:
  `validation/muse_powered_validation/frozen_red_gate_manifest.json`
- H-beta+[O III] fit factory: `beat.injection.muse_hbeta_oiii_fit`

The frozen manifest contains 120 supported checkpoint hashes and four excluded
boundary checkpoint hashes. Per the freeze decision, no cases should be added
unless a predeclared gate later fails or an unexplained flag is found inside
the supported domain.
