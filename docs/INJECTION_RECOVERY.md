# Injection/recovery validation

BEAT includes a reproducible two-stage MUSE [O III] doublet validation
experiment. It is intended to measure scientific identifiability, bias, false
components, and runtime—not merely whether the program executes.

## Experiment design

The controlled stage generates spectra with 1.25 Angstrom sampling, a linear
continuum, known Gaussian uncertainties, deterministic random seeds, and a
Gaussian MUSE LSF approximation with resolving power 3027. It covers blank,
single-component, broad, velocity-offset, unequal-flux double, and equal-flux
double cases over peak S/N 3–20 and velocity separations 75–300 km/s.

The real-noise stage begins with selected NGC 2992 and NGC 3393 MUSE spaxels.
It masks the native [O III] line regions, removes a sigma-clipped quadratic
trend, and resamples contiguous residual and uncertainty blocks together.
Artificial lines are then placed on a continuum from the same family fitted by
BEAT. This preserves non-Gaussian residual amplitudes, heteroscedasticity, and
short-range structure without retaining the native emission lines.

BEAT is not told the injected component count. Every spectrum compares zero,
one, and two components using `delta_logz: 5`. UltraNest uses the slice sampler
because its region sampler was pathologically slow for the tight high-S/N
two-component alternative during the recorded pilot audit.

## Running the experiment

From the project root, using an environment that contains BEAT's dependencies:

```bash
PYTHONPATH=src python validation/run_injection_recovery.py \
  --profile standard \
  --output validation/injection_recovery_results
```

The runner checkpoints each case. Rerunning the same command resumes completed
cases and rebuilds the aggregate tables, plots, and report. Use `--no-resume`
only when intentionally repeating every nested-sampling run. `--profile pilot`
runs five calibration cases. Repeatable `--only CASE_NAME` arguments support
targeted convergence audits.

## Standard-grid result

The recorded 26-case grid found:

- 0/4 false positives in blank spectra;
- 9/10 correct single-component counts;
- 4/4 correct counts for two components separated by at least 300 km/s with
  primary peak S/N at least 10;
- failure to separate every tested 150 km/s pair at `delta_logz: 5`; and
- one real-noise NGC 3393 single-component injection incorrectly selected as
  two components, confirmed at tighter sampler settings.

For correctly counted cases, the median absolute velocity error was 4.5 km/s,
the median absolute width error was 12.2%, and the median absolute [O III] 5007
flux error was 5.9%. These values describe this small diagnostic grid and must
not be advertised as population-wide completeness or accuracy.

The full recorded report, individual spectra, injected truth, sampler results,
analysis tables, and plots live under
`validation/injection_recovery_results/`.

## Interpretation

The 150 km/s result is an identifiability limit for this particular line width,
S/N, evidence threshold, prior volume, sampling, and LSF—not a universal hard
resolution limit. Similarly, the NGC 3393 false component shows that structured
residuals can mimic the tied doublet pattern closely enough to overcome the
current evidence threshold.

Before 2.0.0, expand the experiment across more noise seeds and spaxels,
wavelength-dependent calibrated LSFs, asymmetric profiles, different width and
flux-ratio grids, correlated-noise models, and expanded NIRSpec/MIRI cases.
Equal-tailed posterior intervals should replace posterior-standard-deviation
coverage for the final calibration audit.

## H-alpha and broad-line alpha pilot

`validation/run_halpha_alpha.py` repeats the controlled/real-residual design
for H-alpha+[N II], using independent H-alpha and [N II] 6583 fluxes and a
fixed [N II] 6548/6583 ratio. It also compares a known broad-H-alpha injection
with and without an independent broad component and fits the selected NGC 1365
nuclear spaxel.

The recorded pilot passes the controlled blank, S/N=10 single, and 300 km/s
pair checks. The apparent real-residual NGC 2992 false component was later
shown to use an invalid high-line-emission residual donor: native line wings
made the residual scatter 14.2 times the supplied uncertainty. Corrected
donors and the initial 12-case acceptance grid recover all component counts
and pass the provisional parameter-error thresholds, though the grid remains
too small for rate claims. The same-window NGC 1365 refinement strongly prefers two broad
Gaussians over one, but the second component approaches its width prior and is
therefore a profile-asymmetry diagnostic rather than an accepted physical
decomposition. See `validation/ALPHA_TEST_REPORT.md` for exact values and
next steps.

The historical `three-component-v1` H-alpha pilot allows zero through three narrow
components in every fit. It recovers 7/8 counts: both blanks, both singles,
both doubles, and one of two triples. The missed NGC 2992 triple remains an
underfit under tighter sampling, while an apparent false third component in an
NGC 3393 double disappears under tighter sampling. The triple-recovery and
median velocity-error gates failed in that undersized pilot. See
`validation/HALPHA_THREE_COMPONENT_PILOT.md` and the machine-readable results
under `validation/halpha_acceptance_grid_3component/`.

The controlled follow-up maps six three-component anchors in Gaussian noise.
For the tested S/N and widths, equal components are unresolved at 200 km/s
adjacent spacing and decisive at 300 km/s. A 300 km/s 1:0.7:0.5 triple is not
stable under tight evidence, while 400 km/s recovers equal, moderate, and
1:0.5:0.25 patterns in the tested seeds. These are conditional calibration
points, not a universal spectral-resolution rule. See
`validation/HALPHA_CONTROLLED_TRIPLE_CALIBRATION.md`.

The frozen powered follow-up supersedes those provisional H-alpha gates. Its
supported red domain contains 60 blanks and 20 singles, doubles, and triples;
all 120 counts are correct, no evidence flags remain, and the median absolute
velocity, width, and flux errors are 4.91 km/s, 5.29%, and 4.41%. The exact
one-sided 95% upper false-positive bound is 4.87%. The complete boundary matrix
also retains a donor-dependent 300 km/s broad-primary/0.5-ratio double outside
the supported domain. A 16-case H-beta+[O III] check exercises the LSF near
5000 A. See `validation/MUSE_POWERED_VALIDATION.md` and the frozen hashes under
`validation/muse_powered_validation/`.

## NIRSpec G235H/G395H pilot

The NIRSpec runner uses the actual wavelength grids from the supplied NGC 4151
G235H and IC 5063 G395H cubes, together with the bundled STScI `R(lambda)`
tables. Controlled spectra use Gaussian noise scaled to the local cube
uncertainty. Real-noise spectra use contiguous blocks from separately selected
low-emission spaxels; target-line regions are excluded and donors must have
robust residual scatter no more than three times their formal uncertainty.

```bash
PYTHONPATH=src python validation/run_nirspec_injection_recovery.py --profile pilot
```

The complete standard core recovers 16/16 component counts, with four cases in
each zero-through-three class. Boundary tests confirm that equal 150 km/s pairs
are not separable and that some G235H triples require standard sampling.

A donor-replicate control then showed that independent-pixel likelihoods fail
on correlated cube residuals. After empirical marginal-error calibration and
an AR(1) likelihood estimated from line-free continuum pixels, the final
40-case matrix recovers 20/20 blanks and 20/20 S/N=10 singles, with no evidence
flags. Median errors are 5.28 km/s, 10.04% in intrinsic width, and 6.46% in
integrated flux. The subsequent powered double/triple expansion passes its
predeclared count gates after audit. G235H weak-component reliability is
supported at effective component S/N>=10; lower-S/N cases remain exploratory.
See `validation/NIRSPEC_INJECTION_RECOVERY.md`.
