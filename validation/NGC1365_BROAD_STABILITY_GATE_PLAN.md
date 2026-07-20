# Predeclared NGC 1365 nine-spaxel broad-H-alpha stability gate

Date frozen: 2026-07-18  
Matrix version: `ngc1365-broad-stability-v1`

## Scientific question

The earlier central-spaxel pilot strongly preferred two broad Gaussians over
one, but its second Gaussian approached the upper width prior. This bounded
gate asks whether that red wing is spatially coherent and robust to continuum
description. It is a modeling diagnostic, not a claim that individual broad
Gaussians are distinct physical BLR components.

## Frozen matrix

Use the automatically selected 3-by-3 nuclear region, with zero-based
coordinates `x=164..166` and `y=154..156`. Every spaxel is fit on the same
rest-frame 6420--6710 Angstrom window with H-alpha and the fixed-ratio [N II]
doublet, the representative MUSE Gaussian LSF, at most two shared narrow
kinematic components, and the validated AR(1) residual-noise treatment with
automatic inflation-only marginal scaling.

For every spaxel, compare six otherwise identical configurations:

- narrow-only, one broad H-alpha Gaussian, and two broad H-alpha Gaussians;
- each under a degree-1 and a degree-2 continuum; and
- the same sidebands, 6420--6470 and 6660--6710 Angstrom, for both continua.

The one-broad prior is velocity `[-2500, 2500]` km/s and intrinsic sigma
`[400, 6000]` km/s. The two-broad model uses a core with velocity
`[-2000, 2000]` km/s and sigma `[400, 1800]` km/s, plus a wing with velocity
`[-3500, 3500]` km/s and sigma `[1800, 6000]` km/s. The non-overlapping width
priors prevent label exchange. Broad-component flux priors remain automatic
and positive.

All 54 fits use the screening profile: 40 live points, 60 ESS, `dlogz=2`, and
12 slice steps. The component-count evidence threshold remains delta ln Z=5.
Checkpoints are deterministic and resumable.

## Predeclared interpretation criteria

Broad H-alpha is operationally required if one broad Gaussian exceeds the
narrow-only evidence by delta ln Z>=5 in at least 7/9 spaxels under each
continuum.

A second broad Gaussian may be called a stable descriptive requirement only
if all of the following hold:

- it exceeds the one-broad evidence by delta ln Z>=5 in at least 7/9 spaxels
  under each continuum;
- the wing velocity is redward (`>0` km/s) in at least 7/9 spaxels under each
  continuum;
- the continuum choices agree on the wing-velocity sign in at least 8/9
  spaxels;
- the median absolute cross-continuum wing-velocity change is <=250 km/s;
- the median cross-continuum wing-width change is <=25%;
- the median absolute cross-continuum wing-flux-fraction change is <=0.10;
- within each continuum, wing-velocity MAD is <=300 km/s and wing-width
  fractional MAD is <=25%; and
- no accepted wing posterior median lies within 5% of either width-prior
  boundary.

Here wing flux fraction is `F_wing/(F_core+F_wing)`. A failed stability
criterion is a successful diagnostic outcome: routine fitting should then use
one broad component, while detailed BLR work should use a flexible asymmetric
profile and should not physically interpret multiple broad Gaussians.

The screen matrix will not be enlarged merely to improve the outcome.

### Final audit stopping rule

After the complete screening matrix, the physical-stability criteria failed
in multiple independent ways. Standard fits required approximately 17
minutes each and could not plausibly overturn the failed spatial, continuum,
sign, flux-fraction, width, and prior-boundary criteria. On 2026-07-18 the
screening matrix was therefore declared terminal for the alpha modeling
recommendation and all remaining audits were cancelled.

Two one-broad standard checkpoints completed during the stopping decision,
but neither has a completed paired two-broad standard fit. They are retained
as provenance and explicitly excluded from all evidence comparisons. All
reported model differences and stability metrics use the internally
consistent 54-fit screening matrix. Disfavored or marginal comparisons are
not audited merely to measure an already negative physical-stability
conclusion more precisely.
