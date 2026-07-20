# Predeclared bounded MIRI alpha gate

Date frozen: 2026-07-18  
Matrix version: `miri-bounded-v1`

## Scope

This gate is intentionally bounded. It tests segment routing across all 12 MRS
sub-bands for both supplied targets, then tests component selection and
parameter recovery in representative 2A, 3B, and 4C cubes. It is not a powered
completeness survey across all MRS configurations.

Measured Gaussian-equivalent LSF sampling is 1.76 pixels/FWHM in 2A
(`undersampled`), 2.10 in 3B (`borderline`), and 3.67 in 4C (`well sampled`).
The values are identical in the two supplied reductions because their spectral
grids match, although their pipeline/CRDS versions differ.

## Frozen matrix

For each of NGC 4151 and IC 5063 in 2A, 3B, and 4C:

- one blank;
- one narrow S/N=15 single centered on a spectral pixel;
- one resolved S/N=10 single centered halfway between pixels; and
- one double separated by three Gaussian-equivalent resolution FWHM, with
  effective component S/N=15 and 10.

Two additional IC 5063 singles, in 2A and 4C, replace the instrumental Gaussian
with an explicitly non-Gaussian stress profile: 85% nominal core plus a 15%
broader, shifted wing. This is a controlled mismatch, not an empirical claim
that it exactly represents the on-orbit MRS LSF.

## Predeclared thresholds

- automatic wavelength selection chooses the expected segment in 24/24
  target/sub-band checks;
- blank false positives: 0/6;
- reference single recovery: 12/12;
- reference double recovery: at least 5/6;
- unexplained reference evidence flags: 0;
- median absolute velocity error <=0.25 instrumental FWHM;
- median absolute flux error <=15%;
- median absolute intrinsic-width error <=25% for components whose intrinsic
  sigma is at least 0.75 times the Gaussian-equivalent instrumental sigma;
- mismatch cases select one component in 2/2, with median velocity error
  <=0.25 instrumental FWHM, flux error <=15%, and intrinsic-width error <=30%.

Every reference failure or evidence flag receives one independent standard
audit without changing `delta ln Z=5`. No cases will be added merely to improve
the result. If the mismatch gate fails, the affected scope will not be called
Gaussian-adequate; it will require an empirical/tabulated profile or explicit
review language.
