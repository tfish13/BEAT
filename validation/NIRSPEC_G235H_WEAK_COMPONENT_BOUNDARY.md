# NIRSpec G235H weak-component reliability boundary

Date: 2026-07-17  
BEAT version: `2.0.0a1`

## Decision

For the tested NGC 4151 G235H/F170LP [Si VI] configuration, BEAT's supported
alpha domain requires an effective peak S/N of at least 10 for every component
whose separate recovery is to be interpreted scientifically. Components below
S/N=10 are not claimed to be complete or reliably separable; their fits are
exploratory and require explicit review.

This is deliberately more conservative than an S/N=4--5 cutoff. Recovery near
that boundary depends on residual donor and component geometry, and remains
non-monotonic even when the same residual realization is used at every injected
S/N. Effective component S/N alone is therefore not a universal selection
function.

Within the stated domain, all eight tested S/N=10 anchors recover the injected
component count under the standard profile and have accepted evidence status:

| Truth family | Geometry | Residual donors | Correct/total | Evidence flags |
|---|---|---:|---:|---:|
| Double | 300 km/s separation | 2 | 2/2 | 0 |
| Double | 400 km/s separation | 2 | 2/2 | 0 |
| Triple | 400 km/s adjacent spacing | 2 | 2/2 | 0 |
| Triple | 500 km/s adjacent spacing | 2 | 2/2 | 0 |
| **Total** | | | **8/8** | **0** |

The alternate-donor 400 km/s triple selected two components during screening,
but recovered three with accepted status under the standard profile. This is
an expected screen-to-standard transition, not an unresolved evidence flag.

## Boundary experiment

The paired 32-case matrix crosses:

- weak-component effective peak S/N = 3.75, 5.0, 6.25, and 10.0;
- primary and alternate real-residual donors;
- doubles separated by 300 and 400 km/s; and
- triples with adjacent spacing of 400 and 500 km/s.

For each donor/family/geometry tuple, all four S/N levels use the same residual
realization. The strongest component has peak S/N=15 for doubles and 20 for
triples; the middle triple component has peak S/N=10. The line width is 80
km/s, the wavelength grid is taken from the real cube, the tabulated G235H LSF
is applied, and the likelihood uses empirically calibrated marginal errors plus
the configured AR(1) covariance model.

Across the full screening matrix, 19/32 component counts were correct. The
S/N=3.75--6.25 results include underfits, overfits, ambiguous selections, and
convergence-unverified selections. No single cutoff in that interval clears
both donors and all geometries. These failures map the unsupported boundary;
they are not software exceptions.

An earlier partial matrix used a different residual realization at every S/N
and is superseded. Only `nirspec_g235h_weak_boundary_paired/` is used for the
boundary decision.

## Automatic-rerun exercise

Automatic `rerun` mode was exercised on the three ambiguous cases inherited
from the powered double/triple gate. It retained the pilot results, triggered
on the documented evidence diagnostics, and independently reran all models:

| Case | Weak S/N | Final count | Truth | Final status |
|---|---:|---:|---:|---|
| Alternate donor, double, 300 km/s | 3.75 | 1 | 2 | `accepted_after_audit` |
| Primary donor, triple, 400 km/s | 5.0 | 2 | 3 | `ambiguous` |
| Alternate donor, triple, 400 km/s | 5.0 | 3 | 3 | `accepted_after_audit` |

Two numerical flags resolved. The remaining flag is explained by the mapped
weak-component boundary and lies outside the supported S/N>=10 domain. The
rerun does not force a higher component count, and the evidence threshold was
not changed after inspecting the results.

Automatic rerun can respond to evidence or convergence diagnostics, but it
cannot recognize an accepted underfit when truth is unknown. Production users
must therefore apply the supported-domain rule to scientific component claims.

## Scope and provenance

The S/N threshold is the injection runner's effective line-peak S/N relative
to the calibrated per-pixel uncertainty. It is conditional on the tested line,
G235H wavelength/LSF, priors, 80 km/s intrinsic widths, separations, real
residual donors, and AR(1) treatment. It is not yet a guarantee for arbitrary
NIRSpec targets, gratings, extended-source morphologies, or covariance models.

- Runner: `validation/run_nirspec_injection_recovery.py`
- Paired boundary matrix: `validation/nirspec_g235h_weak_boundary_paired/`
- Ambiguous-case automatic reruns:
  `validation/nirspec_g235h_ambiguous_autorerun/`
- Standard S/N=10 transition audit:
  `validation/nirspec_g235h_weak_boundary_snr10_audit/`

## Gate status

The scoped G235H weak-component reliability gate is closed for alpha testing:
the supported domain is documented, all tested cases within it have correct
component counts and accepted evidence status, and all flags outside it have a
documented boundary interpretation. Expanding that domain below S/N=10 would
require additional donors and a predeclared conditional selection function.
