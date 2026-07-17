# NGC 2992 H-alpha false-positive audit

Date: 2026-07-16  
BEAT version: `2.0.0a1`

## Finding

The original NGC 2992 H-alpha single-component “false positive” was not a
valid noise-recovery case. Its residual donor was the automatically selected
high-H-alpha science spaxel `(x, y)=(162, 150)`. Native H-alpha and [N II]
wings extended into the nominal line-free pool, producing a robust residual
scatter 14.20 times the median uncertainty. The two-component model was then
overwhelmingly preferred (`delta ln Z(2-1)=622.95`) because it was fitting
unmasked astrophysical structure.

The roles of science aperture and noise donor are now separated. Corrected
donors are selected for low line emission and residual/uncertainty consistency:

| Donor label | x | y | Residual/uncertainty ratio |
|---|---:|---:|---:|
| S/N approximately 5 | 96 | 144 | 0.94 |
| S/N approximately 10 | 156 | 35 | 1.25 |
| S/N approximately 20 | 216 | 90 | 1.17 |

BEAT's injection utility now rejects a donor when its robust line-free
residual exceeds three times its median supplied uncertainty. This check fails
before nested sampling and is covered by an automated test.

## Corrected audit

- Eight independent block-resampled residual realizations from the primary
  corrected donor recovered 8/8 single components.
- Two alternate donors recovered 2/2 singles.
- Continuum degrees 0 and 2 both recovered one component when applied to the
  same degree-1 generated spectrum.
- Three loose pilot-sampler runs of one borderline spectrum produced
  `delta ln Z(2-1)=4.03, 0.65, 5.23`; the last barely selected two.
- Three tighter repeats of that identical spectrum produced
  `delta ln Z(2-1)=1.24, 1.44, 1.93`, all selecting one component.

The remaining loose-sampler flip is therefore an evidence-convergence issue,
not persistent evidence for an additional astrophysical component. The
acceptance-grid runner now reruns a case at tighter settings whenever the
reported combined log-evidence uncertainty overlaps the `delta ln Z=5`
selection boundary.

## Initial expanded grid

The corrected grid uses one uncertainty-consistent donor in NGC 2992 and one
in NGC 3393. Two realizations per target were run for each scenario:

| Scenario | Correct/total | Result |
|---|---:|---|
| Blank | 4/4 | no false component |
| S/N 10 single | 4/4 | pass |
| 300 km/s double | 4/4 | pass |

For correctly counted nonblank cases:

- median absolute velocity error: 7.56 km/s;
- median absolute width fractional error: 0.108; and
- median absolute flux fractional error: 0.067.

These values pass the provisional thresholds of 10 km/s, 0.15, and 0.10,
respectively. No one of the 12 cases triggered an ambiguity rerun.

## Interpretation

The prior alpha report's NGC 2992 false-positive finding is superseded and
must not be included in a scientific false-positive rate. The corrected
initial grid passes every provisional gate, but four realizations per scenario
cannot establish a <=5% false-positive rate or stable completeness. Expand to
at least 20 independent realizations per scenario, preferably 50 for the blank
rate, before beta acceptance.

This conclusion applies to the historical grid capped at two components. The
later zero-through-three component pilot preserves the donor correction but
adds a harder selection question. Its NGC 2992 triple injection is recovered
as two components even under tighter sampling; see
`HALPHA_THREE_COMPONENT_PILOT.md`.

Reproducibility products are under:

- `validation/ngc2992_false_positive_audit/`;
- `validation/halpha_acceptance_grid/`; and
- `validation/halpha_acceptance_grid_3component/`.
