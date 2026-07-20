# NGC 1365 nine-spaxel broad-H-alpha stability study

Date: 2026-07-18  
BEAT version: `2.0.0a1`  
Matrix version: `ngc1365-broad-stability-v1`

## Outcome

One broad H-alpha component is required in all nine nuclear spaxels under
both linear and quadratic continua. A second broad Gaussian is not a stable
descriptive or physical component: it is preferred in only 2/9 linear and
3/9 quadratic fits and fails the predeclared spatial-coherence,
continuum-robustness, velocity-sign, flux-fraction, width-stability, and
prior-boundary criteria.

The alpha recommendation is therefore:

> Use one broad component for routine BEAT fitting. For detailed BLR work,
> use a flexible asymmetric broad-line profile. Do not automatically assign
> multiple broad Gaussians a direct physical interpretation.

## Frozen design

The matrix contains 54 fits: nine zero-based spaxels (`x=164..166`,
`y=154..156`), each fit with narrow-only, one-broad, and two-broad models
under degree-1 and degree-2 continua. Every fit uses the same 6420--6710
Angstrom rest-frame window, sidebands, H-alpha+[N II] narrow model,
representative MUSE Gaussian LSF, at most two narrow components, and the
validated AR(1) likelihood with automatic marginal-error scaling.

Broad priors were expanded relative to the initial pilot. The one-broad model
allows sigma 400--6000 km/s. The two-broad model uses a 400--1800 km/s core
and a separately labeled 1800--6000 km/s wing, preventing label exchange.

## Results

| Criterion | Linear | Quadratic | Required | Status |
|---|---:|---:|---:|---|
| One broad preferred over narrow-only | 9/9 | 9/9 | at least 7/9 each | pass |
| Two broad preferred over one | 2/9 | 3/9 | at least 7/9 each | fail |
| Redward wing velocity | 7/9 | 6/9 | at least 7/9 each | fail |
| Wing velocity MAD | 372 km/s | 989 km/s | <=300 km/s each | fail |
| Wing width fractional MAD | 10.2% | 32.0% | <=25% each | fail |

Across matched spaxels, the continuum choices agree on the wing-velocity sign
in only 6/9 cases, below the required 8/9. The median absolute velocity change
is 772 km/s, versus the 250 km/s limit. The median wing-flux-fraction change
is 0.154, versus the 0.10 limit. Width changes between continua pass in the
median (19.2% versus 25%), but seven of 18 wing solutions lie within 5% of the
lower wing-width prior boundary. Passing one width metric does not rescue the
failed combined stability gate.

The original central-spaxel, independent-pixel pilot favored two broad
Gaussians by delta ln Z=341. With the validated AR(1) likelihood, wider priors,
and the linear continuum, the same central spaxel instead gives delta ln Z
approximately -1.4. Under a quadratic continuum it gives approximately +5.2.
This reversal demonstrates that the original second Gaussian was strongly
sensitive to noise and continuum assumptions.

## Audit stopping decision

The screening matrix consumed 7,157 sampler-seconds (1.99 CPU-hours) and
already failed the second-component stability gate in several independent
ways. Standard audits took approximately 17--20 minutes per one-broad fit.
Because marginal evidence refinements could not plausibly repair the failed
spatial and methodological criteria, the screening matrix was declared
terminal and the remaining audits were cancelled.

Three one-broad standard checkpoints completed during that stopping decision.
None has a completed paired two-broad standard checkpoint. They are retained
for provenance but excluded from every reported evidence comparison and
stability metric. This avoids mixing evidence estimates from different
sampler profiles.

## Reproducibility and freeze

- Predeclared plan and stopping amendments:
  `validation/NGC1365_BROAD_STABILITY_GATE_PLAN.md`
- Runner: `validation/run_ngc1365_broad_stability.py`
- Complete screening summary:
  `validation/ngc1365_broad_stability/summary.json`
- Frozen hashes:
  `validation/ngc1365_broad_stability/frozen_gate_manifest.json`

The frozen product contains all 54 screening checkpoints, the matrix
definition, summary, and the three explicitly excluded standard checkpoints.
No additional broad-Gaussian audits are required for the alpha modeling
recommendation.
