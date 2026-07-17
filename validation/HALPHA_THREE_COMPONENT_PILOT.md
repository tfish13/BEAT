# H-alpha three-component acceptance pilot

Date: 2026-07-16  
Version tested: `2.0.0a1`  
Grid version: `three-component-v1`

## Purpose

This pilot tests the existing generic `max_components: 3` implementation. It
is separate from the earlier 12-case grid, which only allowed zero through two
components and therefore could not measure false third components or recovery
of a true third component.

Every case in this pilot allows zero, one, two, or three narrow kinematic
components. The model retains independent H-alpha and [N II] 6583 fluxes,
fixes [N II] 6548/6583 to 0.335, uses the representative MUSE WFM Gaussian LSF,
and selects an added component at delta ln Z greater than 5.

## Design

One residual realization was run for each scenario and each uncertainty-
consistent donor in NGC 2992 and NGC 3393:

| Scenario | Injected components |
|---|---|
| Blank | none |
| Single | v=0 km/s, sigma=80 km/s, peak S/N=10 |
| Double | v=-150,+150 km/s, sigma=80 km/s, relative peaks 1:0.5 |
| Triple | v=-300,0,+300 km/s; sigma=70,100,80 km/s; peak S/N=20,14,10 |

Recovered and injected components are paired by a minimum joint velocity-width
cost rather than by list position. A known injection misclassification or an
evidence comparison overlapping the selection threshold automatically triggers
a tighter convergence run.

## Results

| Scenario | Correct/total | Overfit | Underfit |
|---|---:|---:|---:|
| Blank | 2/2 | 0 | 0 |
| Single | 2/2 | 0 | 0 |
| Double | 2/2 | 0 | 0 |
| Triple | 1/2 | 0 | 1 |

Overall component-count accuracy is 7/8. Among correctly counted nonblank
cases, the median absolute velocity error is 11.20 km/s, median absolute width
fractional error is 0.110, and median absolute flux fractional error is 0.095.
The width and flux gates pass. The provisional triple-recovery gate and
10 km/s median velocity-error gate fail. Two samples per scenario are not
enough for rate claims even where the point estimate passes.

## Tight convergence audits

The NGC 2992 triple injection selected two components in both profiles. The
pilot delta ln Z(3-2) was +1.68; under tighter sampling it became -3.47. The
selected two-component solution recovers the blue component but combines the
central and red components into a broad sigma approximately 269 km/s feature.
The full pilot plus tight run took 822 seconds.

The NGC 3393 double injection initially selected a false third component with
pilot delta ln Z(3-2)=+9.55. Tighter sampling reduced this to +2.77, below the
selection threshold, and restored the correct two-component result. The full
pilot plus tight run took 740 seconds. Thus the apparent overfit was a loose-
evidence convergence failure, not stable evidence for residual structure.

## Interpretation and next action

No core rewrite was required: construction, priors, likelihood evaluation,
serialization, and plotting all operated with three components. The current
scientific selection calibration is not yet accepted for three components.

Do not launch the large realization grid with the current design yet. First:

1. add controlled Gaussian-noise triple injections to separate model behavior
   from real-residual identifiability;
2. sample triple flux ratios and adjacent separations, including an equal-
   strength high-S/N reference case;
3. use multiple valid residual donors per target;
4. retain automatic tight reruns for misclassifications, while recording pilot
   and tight results independently; and
5. keep tight three-component evidence as a validation or selected-spectrum
   mode, because it is too expensive for brute-force full-cube fitting.

Machine-readable results are under
`validation/halpha_acceptance_grid_3component/`.

## Controlled-noise follow-up

The subsequent analytic-noise calibration confirms that the engine recovers
three equal components decisively at 300 km/s adjacent spacing, while an equal
triple at 200 km/s is recovered as two. At 300 km/s, changing the flux pattern
to 1:0.7:0.5 produces a marginal pilot pass that becomes a two-component result
under tight sampling. At 400 km/s, equal, moderate, and 1:0.5:0.25 patterns all
recover three in the tested seed. See
`HALPHA_CONTROLLED_TRIPLE_CALIBRATION.md`.
