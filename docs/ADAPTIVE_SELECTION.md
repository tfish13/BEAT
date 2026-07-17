# Adaptive evidence selection

BEAT still selects component count using the configured Bayesian-evidence
threshold. For candidate model `N` and the currently preferred simpler model
`K`, the candidate is selected when

`ln Z(N) - ln Z(K) > selection.delta_logz`.

The default threshold is 5. Adaptive selection does not change this scientific
rule. It records whether the numerical evidence estimate is reliable enough to
use without further convergence work.

## Result statuses

Every successful result contains `selection_status`, `selection_reasons`,
`selection_diagnostics`, and `selection_audit`.

- `accepted`: no configured reliability trigger fired.
- `ambiguous`: the selection threshold lies within the configured evidence-
  uncertainty margin.
- `convergence_unverified`: a reliability trigger fired but no resolving audit
  was completed, or an evidence uncertainty was non-finite.
- `accepted_after_audit`: a configured tight rerun completed and its final
  evidence comparison no longer overlaps the threshold.

The catalogue repeats the status, reasons, and whether an audit was performed.
The run manifest counts each status across the fitted dataset. Every batch run
also writes `selection_review.csv`, containing only `ambiguous` and
`convergence_unverified` spectra plus the evidence comparison closest to the
threshold. This is the candidate list for bounded follow-up.

## Configuration

```yaml
fit:
  kinematics:
    max_components: 3
  selection:
    delta_logz: 5.0
    stop_when_not_preferred: true
    audit:
      mode: flag  # none, flag, or rerun
      uncertainty_sigma: 1.0
      minimum_margin: 0.5
      max_component_decisive_delta_logz: 20.0
      sampling:
        min_num_live_points: 100
        min_ess: 200
        dlogz: 0.5
        nsteps: 20
```

`uncertainty_sigma` and `minimum_margin` define the larger of two margins around
`delta_logz`. If the evidence difference falls inside that margin, the result
is ambiguous.

Selecting the configured maximum multi-component count (two or more) can
deserve extra scrutiny even when it does not overlap the threshold. When its
supporting delta ln Z is below `max_component_decisive_delta_logz`, `flag` mode marks it
`convergence_unverified`, while `rerun` mode repeats the full model sequence
with the audit sampling profile. Set this option to `null` to disable that
trigger.

The audit profile cannot make minimum live points, effective samples, slice
steps, or the evidence tolerance less rigorous than the initial profile.
Completed reruns retain the pilot component count, evidence, models,
diagnostics, and sampling settings in `selection_audit`.

## Recommended operating modes

Use `flag` for full cubes and large surveys. It adds negligible fitting cost
and lets users select a manageable candidate set for review or later reruns.

Use `rerun` for bounded apertures, individual spectra, or a curated candidate
set. A tight three-component audit took approximately 12-21 minutes in the
recorded MUSE calibration cases, so automatic reruns across every cube spaxel
are not recommended.

`none` retains the evidence diagnostics but does not claim that triggered
selections were convergence-audited.

## Validation replay and limitation

The recorded replay under `validation/adaptive_selection_replay.json` applies
this policy to stored MUSE calibration evidence. It flags the known loose-
sampler false third component, marks the delta ln Z=5.23 case ambiguous,
accepts their tight two-component results, accepts a delta ln Z approximately
280 triple without audit, and flags a delta ln Z approximately 15 maximum-
component selection for review.

Adaptive evidence assessment cannot reveal a component that the data do not
identify. For example, a true injected triple whose three-component evidence
is below the two-component evidence remains a scientifically valid underfit.
Injection/recovery grids are still required to map those completeness limits.
