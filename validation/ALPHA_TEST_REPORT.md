# BEAT 2.0 alpha-test report

Date: 2026-07-18  
Version tested: `2.0.0a1`

## Outcome

The bounded alpha workflow now runs on the supplied MUSE, NIRSpec, and MIRI
cubes without loading full cubes into memory. Instrument-adapter smoke tests,
automatic region selection, representative Gaussian LSF handling, tied
H-alpha+[N II] fitting, and independent broad permitted-line components are
operational.

Production fitting now preserves the delta ln Z component-selection rule while
adding evidence diagnostics and trust statuses. Full-cube/survey `flag` mode
does not add sampler work; optional `rerun` mode repeats only threshold-overlap
or moderately supported maximum-component selections with a tighter profile
and retains both results.

Replaying the policy against stored calibration results flags the known loose-
sampler false third component and the marginal delta ln Z=5.23 triple, accepts
their tight two-component results, accepts the decisive delta ln Z=280 triple,
and flags a delta ln Z=15 maximum-component result for review. This validates
the intended convergence triage but does not remove empirical completeness
limits for blended or weak injected components.

This is not yet a complete `2.0.0b1` release because astronomer-pilot and
formal closeout work remain. The synthetic 1D survey workflow gate now passes;
real-SDSS scientific completeness remains explicitly unvalidated. The broad-line alpha recommendation is now
closed: one broad H-alpha component is the routine default, while multiple
broad Gaussians are not automatically physical. The bounded MIRI profile and
undersampling gate also passes within its documented domain. The MUSE
narrow-line gate is now powered and frozen: its supported red domain has 0/60
blank false positives, 20/20 correct singles, doubles, and triples, no evidence
flags, and passing median parameter errors. The original contaminated-donor
and undersized H-alpha pilots are retained below as the audit trail they
supersede. NIRSpec's covariance-aware powered gates also pass within their
documented domains; G235H weak components require effective S/N>=10, with
lower-S/N threshold overlaps treated as explained out-of-domain boundaries.

## Inputs and adopted redshifts

| Target | Instrument | Adopted redshift | Basis |
|---|---|---:|---|
| NGC 2992 | MUSE | 0.00771 | Marinucci et al. 2022, citing Keel 1996 |
| NGC 3393 | MUSE | 0.0125 | Contini 2012 |
| NGC 1365 | MUSE | 0.00546 | Venturi et al. 2018 |
| NGC 4151 | NIRSpec/MIRI | 0.003326 | literature/AGN Variability Archive value |
| IC 5063 | NIRSpec/MIRI | 0.01135 | Travascio et al. 2021 |

The NGC 1365 systemic convention is close to, but not identical with, the
1671 km/s value used by Lena et al. (2016). Publication analyses should state
the adopted convention explicitly and may refit the velocity zero point.

## Automatic bounded-region selection

Regions were ranked from continuum-subtracted line-slab signal with propagated
uncertainty in 3-by-3 spaxel boxes. Bounds below are zero-based and half-open.
They are test selections, not scientific aperture definitions.

| Dataset/diagnostic | High-S/N x range | High-S/N y range |
|---|---:|---:|
| NGC 2992 MUSE H-alpha+[N II] | [162, 165] | [150, 153] |
| NGC 3393 MUSE H-alpha+[N II] | [147, 150] | [167, 170] |
| NGC 1365 MUSE broad H-alpha+[N II] | [164, 167] | [154, 157] |
| NGC 4151 NIRSpec Pa-alpha | [15, 18] | [26, 29] |
| NGC 4151 NIRSpec [Si VI] 1.963 micron | [27, 30] | [15, 18] |
| IC 5063 NIRSpec Br-alpha | [43, 46] | [43, 46] |
| NGC 4151 MIRI [Ne V] 14.32 micron, 3B | [31, 34] | [16, 19] |
| IC 5063 MIRI [Ne V] 14.32 micron, 3B | [26, 29] | [25, 28] |

The selection manifest also records candidate regions for NIRSpec H2 and
[Mg IV], and for MIRI [Ne VI] 7.65, [Ne II] 12.81, [Ne III] 15.56, and
[O IV] 25.89 micron. The MIRI candidates span segments 2A through 4C and are
therefore suitable for later wavelength-dependent resolution tests.

## LSF conventions

The H-alpha pilot uses the representative UDF-10 WFM Gaussian relation

`FWHM(lambda) = 5.866e-8 lambda^2 - 9.187e-4 lambda + 6.040 Angstrom`.

This is an empirical WFM approximation, not a dataset-matched MUSE
`LSF_PROFILE`. NIRSpec G235H and G395H now use the STScI ETC tabulated
`R(lambda)` curves under a 2.2-pixel, fully illuminated-aperture convention.
MIRI uses the documented approximate in-flight relation
`R(lambda_micron) = 4603 - 128 lambda_micron`. Every result records the
approximation and instrument/pipeline provenance.

## H-alpha+[N II] injection/recovery pilot

The narrow model fits independent H-alpha and [N II] 6583 fluxes while fixing
`F([N II] 6548)/F([N II] 6583)=0.335`. The original pilot and corrected
12-case grid compared zero, one, and two narrow components with delta ln Z=5.

| Check | Result | Provisional gate | Status |
|---|---:|---:|---|
| Controlled blank false positive | no | <=5% rate | pilot pass; rate not yet measurable |
| Controlled S/N=10 single count | 1 -> 1 | >=90% | pilot pass |
| Controlled 300 km/s pair | 2 -> 2 | >=80% | pilot pass |
| Real NGC 3393 residual 300 km/s pair | 2 -> 2 | >=80% | pilot pass |
| Original NGC 2992 residual single | 1 -> 2 | n/a | invalid donor; superseded |
| Original median absolute velocity error | 11.52 km/s | n/a | superseded donor set |
| Original median absolute width fractional error | 0.166 | n/a | superseded donor set |
| Original median absolute flux fractional error | 0.129 | n/a | superseded donor set |

Only five cases were present in the original pilot. The NGC 2992 audit found
that the donor's line-free residual scatter was 14.20 times its uncertainty
because native H-alpha/[N II] wings contaminated the residual pool. Corrected
low-emission donors have residual/uncertainty ratios of 0.94-1.25. Ten of ten
independent/alternate corrected single-component spectra recover one
component. One loose-sampler repeat crossed the evidence threshold, but three
tighter repeats consistently recovered one component. See
`NGC2992_FALSE_POSITIVE_AUDIT.md`.

The corrected initial acceptance grid contains four blank, four S/N 10 single,
and four 300 km/s double cases split across NGC 2992 and NGC 3393 donors. It
recovers 12/12 component counts, with median velocity error 7.56 km/s, width
fractional error 0.108, and flux fractional error 0.067. These pass the
provisional thresholds, but four realizations per scenario are not sufficient
to measure a 5% false-positive rate.

| Corrected-grid check | Result | Provisional gate | Initial status |
|---|---:|---:|---|
| Blank false positives | 0/4 | <=5% | pass, undersized |
| S/N 10 single count | 4/4 | >=90% | pass, undersized |
| 300 km/s double count | 4/4 | >=80% | pass, undersized |
| Median absolute velocity error | 7.56 km/s | <=10 km/s | pass |
| Median absolute width fractional error | 0.108 | <=0.15 | pass |
| Median absolute flux fractional error | 0.067 | <=0.10 | pass |

### Three-component expansion

The subsequent `three-component-v1` pilot allows zero through three narrow
components for every case and adds a true triple at velocities -300, 0, and
+300 km/s with peak S/N 20, 14, and 10. It uses one realization for each of
four scenarios in both corrected donors.

| Expanded-grid check | Result | Provisional gate | Status |
|---|---:|---:|---|
| Blank false positives | 0/2 | <=5% | point estimate passes; undersized |
| S/N 10 single count | 2/2 | >=90% | point estimate passes; undersized |
| 300 km/s double count | 2/2 | >=80% | pass after one tight rerun; undersized |
| Three-component count | 1/2 | >=70% | fail |
| False third components after audit | 0/6 simpler cases | <=5% | point estimate passes; undersized |
| Median absolute velocity error | 11.20 km/s | <=10 km/s | fail |
| Median absolute width fractional error | 0.110 | <=0.15 | pass |
| Median absolute flux fractional error | 0.095 | <=0.10 | pass |

The NGC 2992 triple remains a two-component recovery under tight sampling;
its central and red injections merge. An NGC 3393 double initially selected a
false third component, but tight sampling reduced delta ln Z(3-2) from 9.55
to 2.77 and restored the correct pair. Tight audit runtimes were approximately
12-14 minutes per failed pilot case. See `HALPHA_THREE_COMPONENT_PILOT.md`.

A controlled Gaussian-noise follow-up separates model capability from real
residual structure. Equal triples recover as two at 200 km/s adjacent spacing
and as three decisively at 300 and 400 km/s. The 300 km/s 1:0.7:0.5 case flips
from a marginal pilot three-component result to two under tight sampling. At
400 km/s, all tested flux patterns through 1:0.5:0.25 recover three, although
the faintest component's flux and width errors are approximately 21% and 15%.
These six anchors are not a powered completeness grid. See
`HALPHA_CONTROLLED_TRIPLE_CALIBRATION.md`.

### Powered MUSE gate

The frozen `muse-powered-v1` real-residual matrix supersedes the undersized
H-alpha count and parameter gates. Its supported red domain contains 60
blanks and 20 cases in each nonblank class, balanced across NGC 2992 and
NGC 3393. It recovers 60/60 blanks, 20/20 singles, 20/20 doubles, and 20/20
triples with no evidence flags. The blank point estimate is 0%, with a 4.87%
one-sided 95% upper bound. Median absolute errors are 4.91 km/s in velocity,
5.29% in intrinsic width, and 4.41% in integrated flux. Every predeclared
powered red gate passes.

The complete boundary matrix retains four 300 km/s, 0.5-ratio,
sigma=(160,80) km/s doubles. Both NGC 3393 realizations underfit, one
ambiguously; the full geometry is therefore excluded from the supported
domain. A separate 16-case H-beta+[O III] wavelength/LSF check recovers all
counts and has median velocity error 4.38 km/s, but one correct NGC 3393
triple remains convergence-unverified. See `MUSE_POWERED_VALIDATION.md`.

## Broad H-alpha validation

BEAT now supports broad permitted-line components with velocity, sigma, and
flux parameters independent of the shared narrow-line family.

For a controlled narrow plus broad injection (`sigma_broad=1200 km/s`):

- broad width was recovered to 3.3%;
- broad integrated flux was recovered to 0.1%;
- broad velocity error was 34.5 km/s; and
- the broad+narrow model was preferred over a narrow-only model by
  delta ln Z=36.45.

For the automatically selected NGC 1365 nuclear spaxel:

- a one-broad-Gaussian model is overwhelmingly preferred to a narrow-only
  model in the same narrow pilot window (delta ln Z approximately 8794);
- the recovered one-Gaussian broad sigma is approximately 1160 km/s, close to
  the published combined broad-line dispersion of 1181 km/s;
- on an identical wider 6420-6710 Angstrom window, two broad Gaussians are
  preferred to one by delta ln Z=341.15; and
- the two-broad model reduces normalized residual RMS from 2.62 to 1.89.

The initial two-Gaussian decomposition is not physically accepted. The frozen
follow-up contains 54 fits across all nine nuclear spaxels, comparing
narrow-only, one-broad, and two-broad models under both linear and quadratic
continua, wider broad priors, and the validated MUSE AR(1) likelihood.

One broad component is preferred in 9/9 spaxels under both continua. A second
broad Gaussian is preferred in only 2/9 linear and 3/9 quadratic fits. Its
velocity sign agrees across continua in only 6/9 spaxels; its median
cross-continuum velocity change is 772 km/s; its flux fraction changes by
0.154 in the median; and seven of 18 wing widths approach the lower width-prior
boundary. The central-spaxel delta ln Z(two minus one) changes from +341 in
the old independent-pixel pilot to -1.4 with AR(1) and a linear continuum,
and +5.2 with a quadratic continuum.

The study therefore recommends one broad component for routine fitting and a
flexible asymmetric profile for detailed BLR analysis. Multiple broad
Gaussians must not be assigned a direct physical interpretation automatically.
See `NGC1365_BROAD_STABILITY.md`.

## JWST bounded smoke tests

All four selected real spectra completed with the Stage-3 cube adapters:

| Dataset/line | Segment | Selected narrow components | LSF approximation |
|---|---|---:|---|
| NGC 4151 NIRSpec [Si VI] | G235H/F170LP | 2 | tabulated R=2240 at the line |
| IC 5063 NIRSpec Br-alpha | G395H/F290LP | 2 | tabulated R=2812 at the line |
| NGC 4151 MIRI [Ne V] | 3B | 2 | R(lambda), R=2764 at the line |
| IC 5063 MIRI [Ne V] | 3B | 2 | R(lambda), R=2749 at the line |

NIRSpec [Si VI] is visibly asymmetric/double-peaked. MIRI [Ne V] shows
structured residuals around its undersampled profile. These observations make
useful diagnostic targets, but the selected counts are not truth labels.

### Bounded MIRI profile and sampling gate

The predeclared MIRI matrix checks automatic segment routing across all 12 MRS
sub-bands for both targets, then performs controlled injections on actual 2A,
3B, and 4C wavelength grids. These represent undersampled (1.76 pixels per
LSF FWHM), borderline (2.10), and well-sampled (3.67) configurations.

| Gate | Result | Threshold | Status |
|---|---:|---:|---|
| Segment selection | 24/24 | 24/24 | pass |
| Blank false positives | 0/6 | 0 | pass |
| Singles | 12/12 | 12/12 | pass |
| Doubles | 6/6 | at least 5/6 | pass |
| Reference evidence flags | 0/24 | 0 | pass |
| Median velocity error | 0.030 resolution FWHM | <=0.25 | pass |
| Median flux error | 3.61% | <=15% | pass |
| Median resolved-width error | 11.95% | <=25% | pass |
| Shifted-wing mismatch cases | 2/2 | 2/2 | pass |

The Gaussian-equivalent LSF is therefore adequate for component selection in
the tested domain: isolated S/N=15 narrow singles, S/N=10 resolved singles,
and resolved doubles separated by three instrumental FWHM with a weakest
component at S/N=10. The 15% shifted-wing stress profile also passes in 2A and
4C. This does not validate intrinsically narrow width recovery in undersampled
regions, weak components below S/N=10, closer blends, all profiles, or
scientific recovery in the other nine sub-bands. See
`MIRI_BOUNDED_VALIDATION.md`.

## NIRSpec injection/recovery validation

The G235H/G395H validation uses actual cube wavelength grids, bundled STScI
resolution tables, separate science/noise spaxels, and up to three narrow
components. The complete standard-profile core recovers 16/16 counts (four in
each class), with median errors of 6.38 km/s, 8.25% in width, and 6.34% in
flux. Boundary tests show that equal 150 km/s pairs remain unresolved and that
a G235H equal triple at 300 km/s adjacent spacing needs standard sampling.

Replicating donors and noise seeds exposed a more important issue: an
independent-pixel likelihood produced six false positives in 16 blanks. A
marginal uncertainty rescale alone reduced but did not eliminate the problem.
The residuals have strong lag-one covariance, particularly in G235H. BEAT now
offers an explicit AR(1) likelihood with rho estimated from configured
line-free continuum windows. With this model and recorded marginal-error
calibration, the powered scoped matrix gives:

| Gate | Result | Threshold | Status |
|---|---:|---:|---|
| Blank false positives | 0/20 | <=5% | pass |
| S/N=10 singles | 20/20 | >=90% | pass |
| Evidence flags | 0/40 | 0 | pass |
| Median velocity error | 5.28 km/s | <=10 km/s | pass |
| Median width error | 10.04% | <=15% | pass |
| Median flux error | 6.46% | <=10% | pass |

The two worst blank/single evidence margins remain correct and accepted under standard
sampling. One audited single has a 13.1 km/s velocity error and 27.8% width
error, so individual outliers remain even though aggregate medians pass.

The powered double/triple expansion balances 20 cases per class across both
gratings, two residual donors per target, 300--500 km/s separations, and equal
through 1:0.25 ratios. Replacing every incorrect or non-accepted screening
result with its standard-profile audit gives:

| Gate | Audited result | Threshold | Status |
|---|---:|---:|---|
| Doubles | 18/20 | >=80% | pass |
| Triples | 17/20 | >=70% | pass |
| Evidence flags | 3/40 | 0 | **fail** |
| Median velocity error | 6.00 km/s | <=10 km/s | pass |
| Median width error | 8.95% | <=15% | pass |
| Median flux error | 6.19% | <=10% | pass |

All G395H counts are correct after audit. A subsequent paired G235H boundary
matrix maps the weak-component transition and supports a conservative alpha
domain of effective component S/N>=10. All eight tested S/N=10 double/triple
anchors recover the correct count under standard sampling with no evidence
flags. Automatic rerun was exercised on the three earlier threshold overlaps;
two resolved, while the remaining S/N=5 triple is an explained out-of-domain
ambiguity. The scoped G235H reliability gate is therefore closed without
relaxing the evidence threshold. See `NIRSPEC_G235H_WEAK_COMPONENT_BOUNDARY.md`.

## Alpha decision and next work

Passed:

- bounded reading and fitting for all three adapters;
- automatic compact-region selection;
- representative wavelength-dependent MUSE/MIRI LSF evaluation;
- bounded MIRI routing, sampling, and controlled profile-mismatch validation;
- the frozen powered MUSE H-alpha+[N II] false-positive, completeness, and
  parameter-error gate, plus the blue H-beta+[O III] LSF check;
- tabulated wavelength-dependent NIRSpec G235H/G395H LSF evaluation, the
  complete standard core, the powered covariance-aware blank/single gate, and
  powered double/triple count-accuracy gates;
- H-alpha+[N II] tied-ratio inference;
- independent broad permitted-line modeling;
- detection of the known NGC 1365 BLR requirement;
- the frozen nine-spaxel one-broad routine-model recommendation; and
- the 512-row synthetic SDSS-like streaming, multiprocessing,
  interruption/resume, partial-failure, and deterministic-provenance gate.

The next release work should therefore be:

1. Run the astronomer pilot and resolve blocker or scientifically dangerous
   findings.
2. Perform the formal alpha closeout and clean-environment package checks.

## Reproducibility products

- `alpha_targets/selection_manifest.json`: all ranked regions and line/segment
  coverage;
- `halpha_alpha_results/`: injected truth, complete inference results,
  same-window BLR comparisons, summaries, and plots;
- `ngc2992_false_positive_audit/`: corrected donors, sampler/continuum audits,
  and the invalid-original-case diagnosis;
- `halpha_acceptance_grid/`: the corrected 12-case starting grid and automatic
  evidence-ambiguity rerun policy;
- `halpha_acceptance_grid_3component/`: the eight-case zero-through-three
  component pilot, including retained pilot and tight audit results;
- `halpha_three_component_calibration/`: controlled separation/flux-ratio
  anchors and retained pilot/tight transition evidence;
- `muse_powered_validation/`: frozen powered H-alpha+[N II] gate, complete red
  boundary matrix, blue H-beta+[O III] LSF checks, and checkpoint hashes;
- `adaptive_selection_replay.json`: production trust-status replay against
  stored false-third, threshold-overlap, and decisive-three-component cases;
- `jwst_alpha_results/`: NIRSpec/MIRI real-spectrum results and plots;
- `miri_bounded_validation/`: all-band routing audit, 2A/3B/4C controlled
  recovery and profile-mismatch results, and frozen checkpoint hashes;
- `ngc1365_broad_stability/`: complete nine-spaxel/two-continuum screening
  matrix, modeling recommendation, and frozen checkpoint hashes;
- `survey_1d_regression/`: generated SDSS-like tables, interrupted/resumed
  workflow outputs, partial-failure catalogs, real-inference smoke results,
  and frozen hashes;
- `nirspec_injection_recovery/`: NIRSpec injected truth, spectra, fits,
  summaries, and plots;
- `nirspec_injection_recovery_audit/`: tighter NGC 4151 blank/triple reruns;
- `nirspec_injection_recovery_standard/`: complete 16-case standard core;
- `nirspec_boundary_expansion/` and `nirspec_boundary_audit/`: resolution and
  weak-component boundary anchors;
- `nirspec_donor_replicates*/`: independent-noise control, marginal-error
  control, final 40-case AR(1) matrix, and worst-margin standard audit;
- `prepare_alpha_targets.py`, `run_halpha_alpha.py`, `run_jwst_alpha.py`,
  `run_nirspec_injection_recovery.py`, and the NGC 1365 refinement runners; and
- `examples/*.local.yaml`: bounded, runnable configurations.

## Literature and calibration references

- [NGC 1365 ionized-gas kinematics and broad-line decomposition](https://academic.oup.com/mnras/article/459/4/4485/2624045)
- [NGC 1365 MAGNUM/MUSE context and adopted redshift](https://www.frontiersin.org/journals/astronomy-and-space-sciences/articles/10.3389/fspas.2017.00046/full)
- [NGC 2992 redshift context](https://doi.org/10.1093/mnras/stac1381)
- [NGC 3393 redshift context](https://academic.oup.com/mnras/article/425/2/1205/1190869)
- [IC 5063 redshift context](https://openaccess.inaf.it/entities/publication/b52d7a81-dffa-42ac-8aa6-790dc5435e10)
- [MUSE empirical LSF polynomial](https://eso.org/public/archives/releases/sciencepapers/eso1738/eso1738a.pdf)
- [NIRSpec dispersion/resolution tables](https://jwst-docs.stsci.edu/jwst-near-infrared-spectrograph/nirspec-instrumentation/nirspec-dispersers-and-filters)
- [MIRI MRS wavelength-dependent resolving power](https://jwst-docs.stsci.edu/jwst-mid-infrared-instrument/miri-observing-modes/miri-medium-resolution-spectroscopy)
