# BEAT 2.0 alpha-test report

Date: 2026-07-17  
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

This is not yet a `2.0.0b1` scientific acceptance result. The original
H-alpha pilot used high-line-emission science spaxels as residual-noise donors;
the NGC 2992 donor was contaminated by native line wings and is invalid for
false-positive scoring. After separating science and noise-donor apertures,
the corrected initial 12-case H-alpha grid passes every provisional gate.
That historical grid only allowed two components. A new eight-case pilot that
allows up to three recovers 7/8 component counts and fails the provisional
triple-recovery and median velocity-error gates. The samples remain too small
for beta-level rate claims. NIRSpec's 16-case standard core now recovers all
zero-through-three-component counts. A separate covariance-aware residual-noise
matrix passes the statistically powered blank/single scope at 40/40, but the
NIRSpec double/triple classes remain underpowered.

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

The two-Gaussian decomposition is not yet physically accepted. Its broad core
has sigma about 1069 km/s, while the second red wing has sigma about 3471 km/s
and velocity about +1134 km/s. The wing approaches the upper width prior and
may absorb BLR asymmetry, continuum mismatch, Fe emission, or calibration
structure. It requires prior expansion and continuum/profile stress tests
across the full nuclear 3-by-3 region.

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
good later injection targets, but the selected counts are not truth labels.
The MIRI residuals motivate empirical or non-Gaussian LSF tests in addition to
the current Gaussian-equivalent resolving power.

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

The two worst evidence margins remain correct and accepted under standard
sampling. One audited single has a 13.1 km/s velocity error and 27.8% width
error, so individual outliers remain even though aggregate medians pass. This
is a powered blank/single result, not powered double/triple completeness. See
`NIRSPEC_INJECTION_RECOVERY.md` for the full control sequence and provenance.

## Alpha decision and next work

Passed:

- bounded reading and fitting for all three adapters;
- automatic compact-region selection;
- representative wavelength-dependent MUSE/MIRI LSF evaluation;
- tabulated wavelength-dependent NIRSpec G235H/G395H LSF evaluation, the
  complete standard core, and the powered covariance-aware blank/single gate;
- H-alpha+[N II] tied-ratio inference;
- independent broad permitted-line modeling; and
- detection of the known NGC 1365 BLR requirement.

Not yet passed:

- three-component H-alpha recovery and median velocity-error tolerances;
- statistically powered H-alpha recovery tolerances (both grids are
  undersized);
- statistically meaningful MUSE false-positive/completeness rates;
- powered NIRSpec double/triple recovery grids;
- MIRI profile-shape/undersampling validation; and
- a stable physical two-component BLR interpretation.

The next release work should therefore be:

1. Add controlled and real-residual three-component calibration cases across
   flux ratio, separation, and donor before scaling to at least 20 independent
   realizations per accepted gate.
2. Expand NIRSpec doubles/triples to 20 cases per class across G235H/G395H
   seeds, widths, flux ratios, separations, and residual donors.
3. Run MIRI injections in representative 2A, 3B/3C, and 4C bands, including
   undersampling and non-Gaussian/profile-mismatch cases.
4. Stress-test one versus two broad H-alpha components over all nine selected
   NGC 1365 nuclear spaxels with wider priors and alternative continua.
5. Freeze acceptance datasets only after those diagnostics are reviewed by
   astronomer pilot users.

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
- `adaptive_selection_replay.json`: production trust-status replay against
  stored false-third, threshold-overlap, and decisive-three-component cases;
- `jwst_alpha_results/`: NIRSpec/MIRI real-spectrum results and plots;
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
