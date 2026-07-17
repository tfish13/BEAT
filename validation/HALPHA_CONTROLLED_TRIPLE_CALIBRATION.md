# Controlled H-alpha three-component calibration

Date: 2026-07-16  
Version tested: `2.0.0a1`  
Calibration version: `controlled-triple-v1`

## Purpose

This experiment isolates three-component identifiability from real MUSE
residual structure. It uses analytic Gaussian noise, deterministic seeds, the
representative MUSE WFM Gaussian LSF, and the same H-alpha+[N II] model and
delta ln Z=5 selection rule as the real-residual pilot.

Each injection has component velocity dispersions 70, 100, and 80 km/s and a
primary H-alpha peak S/N of 20. The tested `separation` is the adjacent spacing:
the component velocities are `-separation`, 0, and `+separation` km/s.

## Planned matrix and anchor strategy

The resumable runner defines nine cells:

- adjacent separations 200, 300, and 400 km/s; and
- flux patterns 1:1:1 (`equal`), 1:0.7:0.5 (`moderate`), and 1:0.5:0.25
  (`weak`).

Six monotonic anchors were run. Cells that are both less separated and weaker
than a failed anchor were not spent automatically; those implications require
confirmation with additional noise realizations before being treated as
measured completeness.

## Results

| Adjacent separation | Flux pattern | Active profile | Recovered | delta ln Z(3-2) | Runtime |
|---:|---:|---|---:|---:|---:|
| 200 km/s | 1:1:1 | pilot | 2 | -2.37 | 253 s |
| 300 km/s | 1:1:1 | pilot | 3 | +87.46 | 244 s |
| 300 km/s | 1:0.7:0.5 | tight | 2 | -0.29 | 1259 s |
| 400 km/s | 1:1:1 | pilot | 3 | +279.98 | 257 s |
| 400 km/s | 1:0.7:0.5 | pilot | 3 | +89.67 | 269 s |
| 400 km/s | 1:0.5:0.25 | pilot | 3 | +14.78 | 218 s |

The 300 km/s moderate case initially selected three at pilot precision with
delta ln Z(3-2)=+5.23 and combined evidence error 1.98. The tight profile
changed the comparison to -0.29 with combined error 1.29 and selected two.
Both profiles are retained in the case checkpoint. A nominal result just over
the evidence threshold must therefore not be treated as a stable detection.

At 400 km/s the weak pattern recovers all three counts, including the peak
S/N=5 third component. Its recovered third-component velocity is within
10.5 km/s, but its flux is high by 20.6%, width is high by 15.4%, and posterior
width uncertainty is large. Count recovery does not imply precision recovery
at that boundary.

## Current scientific interpretation

For this single deterministic noise realization, line widths, LSF, priors,
and S/N:

- 200 km/s adjacent spacing is unresolved even for equal components;
- 300 km/s is decisive for equal components but not for 1:0.7:0.5;
- 400 km/s recovers all three tested flux patterns; and
- a transition result close to delta ln Z=5 requires a tighter audit.

These are calibration anchors, not universal resolution limits. They must not
be generalized across line width, S/N, LSF, residual covariance, or target.

## Next calibration step

Before a powered grid:

1. repeat the 300 km/s equal and 400 km/s weak boundary cells across noise
   seeds;
2. add 350 km/s moderate and weak cells to refine the transition;
3. repeat accepted controlled cells with valid real-residual donors; and
4. retain staged screening, because tight boundary fits take approximately
   21 minutes on the test laptop.

Machine-readable checkpoints, manifest, and summary are under
`validation/halpha_three_component_calibration/`.
