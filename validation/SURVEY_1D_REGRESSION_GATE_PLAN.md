# Predeclared synthetic 1D survey regression gate

Date frozen: 2026-07-20  
Matrix version: `survey-1d-synthetic-v1`

## Scope

This is an engineering acceptance test for the production `survey_table`
workflow, not a scientific completeness study of the SDSS QSO population. It
uses deterministic SDSS-like logarithmic wavelength grids around rest-frame
H-beta+[O III] and requires no downloaded public sample.

## Large workflow matrix

The FITS binary table contains 512 rows of 321-pixel spectra:

- 500 valid rows with unique deterministic target IDs, redshifts, wavelength,
  flux, inverse variance, and injected class metadata;
- four all-NaN flux rows;
- four zero-inverse-variance rows; and
- four non-finite-redshift rows.

The workflow configuration uses a real, inexpensive UltraNest continuum-only
fit (`max_components: 0`) so hundreds of production worker tasks can exercise
streaming, bounded process submission, checkpoints, catalogs, and failures
without launching a scientific sampling campaign.

Run with four workers, terminate the complete process group after at least 32
successful checkpoints, then resume. After completion, run once with a
one-worker CLI override to confirm that operational worker count does not
alter the scientific fingerprint or successful checkpoints.

## Small real-inference matrix

A separate four-row table contains two blanks, one one-component spectrum,
and one well-separated two-component spectrum. It runs the real H-beta+[O
III] model with tied [O III] 4959/5007 flux and 40 live points/40 ESS. A
planned 12-row repetition was stopped after the first blank/double/blank
checkpoints showed that the nonzero fits were the expensive step. The single
was then run with a max-one configuration after its max-two run exceeded the
15-minute cutoff; the already completed double retains the max-two path. This
is a fitter-to-result/catalog path check, not a powered recovery claim.

## Acceptance gates and stopping rule

- interruption occurs after 32--499 successful checkpoints;
- every pre-interruption successful checkpoint hash is unchanged after resume;
- final workflow catalog has 500 unique deterministic IDs;
- failures table has exactly the 12 declared malformed IDs while all valid
  outputs survive;
- resumed run reports all preexisting valid checkpoints and creates one
  coherent catalog;
- config fingerprint is stable across interruption, resume, and worker-count
  override;
- final manifest records the effective CLI worker count;
- every successful result records the deterministic per-spectrum sampler seed;
- measured peak process-tree RSS, when available, is <=1.5 GiB;
- small real-inference checkpoints contain zero failures and at least one
  selected result in each zero-, one-, and two-component class; and
- full current unit/integration suite passes.

Stop after these gates. Do not expand into SDSS population completeness or
parameter-error calibration. The supported statement is operational only:
the 1D survey workflow is validated on a synthetic SDSS-like regression
dataset; real-SDSS scientific completeness remains unestablished.
