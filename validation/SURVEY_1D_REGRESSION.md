# End-to-end synthetic 1D survey regression

Date: 2026-07-20  
BEAT version: `2.0.0a1`  
Matrix version: `survey-1d-synthetic-v1`

## Outcome

Every declared operational gate passes. BEAT streamed a 512-row SDSS-like
FITS table through the production CLI, performed real UltraNest fits with four
workers, survived intentional interruption and resume, preserved deterministic
checkpoint provenance, isolated 12 malformed spectra, generated coherent
catalog and failure products, and exercised real H-beta+[O III] zero-, one-,
and two-component inference.

The supported alpha statement is deliberately limited:

> The 1D survey workflow is operationally validated on a synthetic SDSS-like
> regression dataset. Scientific completeness for a real SDSS QSO population
> has not been established.

## Large workflow result

| Gate | Result | Status |
|---|---:|---|
| Input rows | 512 | complete |
| Valid catalog rows | 500 | pass |
| Unique deterministic IDs | 500/500 | pass |
| Declared malformed failures | 12/12 | pass |
| Mid-run interruption | 32 successful checkpoints | pass |
| Interrupted checkpoint hashes preserved | 32/32 | pass |
| Four-worker resume | 468 new, 32 resumed, 12 failed | pass |
| One-worker resume | 0 new, 500 resumed, 12 failed | pass |
| Successful hashes unchanged by worker override | 500/500 | pass |
| Scientific config fingerprint stable | identical | pass |
| Effective CLI workers recorded | 4, then 1 | pass |
| Deterministic sampler seeds recorded | 500/500 | pass |

The intentional interruption occurred after 5.66 seconds. The four-worker
resume completed in 56.0 seconds. The subsequent one-worker resume took 1.01
seconds and did not modify any successful checkpoint. CLI exit code 1 on the
completed runs is expected because BEAT reports partial failures while
preserving successful results.

The table reader keeps the FITS table memory-mapped, copies only one vector
cell at a time, and the process-pool scheduler holds at most twice the worker
count in pending tasks. Catalog creation uses two streaming passes over
checkpoint files. Process-tree RSS could not be sampled because `psutil` is
not installed in the test environment; the structural streaming path was
exercised, but an empirical peak-RSS number is not claimed.

## Malformed inputs and failure isolation

Four rows each contained all-NaN flux, zero inverse variance, or non-finite
redshift. All 12 appear once in `failures.csv`, while all 500 valid spectra
remain in `catalog.csv`.

The test exposed and fixed an error-reporting defect for spectra with no valid
pixels. Those rows previously leaked an internal boolean-index `IndexError`.
They now produce explicit `ValueError` messages stating that fewer than two
valid wavelength/flux pixels remain. Non-finite redshifts likewise produce a
specific redshift-validation message.

## Multiprocessing and deterministic provenance

The restricted macOS execution environment denied Python's
`SC_SEM_NSEMS_MAX` preflight query before worker creation. BEAT now treats a
permission-denied semaphore limit as indeterminate (`-1`), matching Python's
existing behavior for unavailable limits, while allowing real worker creation
to proceed. Actual four-process execution then passed.

When `fit.sampling.seed` is configured, BEAT derives a stable 32-bit seed from
the base seed and spectrum ID. This seed is recorded in every success or
failure checkpoint. Worker count remains operational rather than scientific,
so changing `--workers` does not change the configuration fingerprint; the
effective override is now recorded in `run_manifest.json`.

## Real H-beta+[O III] inference

The real inference smoke set uses logarithmic observed-frame wavelength grids,
individual redshifts, H-beta and [O III] 5007 fluxes, a fixed [O III]
4959/5007 ratio, and real UltraNest sampling.

| Spectrum | Injected class | Selected components |
|---|---|---:|
| `SDSS-SCIENCE-000` | blank | 0 |
| `SDSS-SCIENCE-001` | single | 1 |
| `SDSS-SCIENCE-002` | double | 2 |
| `SDSS-SCIENCE-003` | blank | 0 |

A planned 12-case repetition was stopped when the max-two single fit crossed
the 15-minute cutoff. The single was rerun with a max-one configuration and
finished in 5.15 seconds; the completed double retains the real max-two path.
This result confirms all three fitter/result branches but is explicitly not a
completeness or parameter-recovery claim.

## Reproducibility and freeze

- Predeclared scope and stopping rule:
  `validation/SURVEY_1D_REGRESSION_GATE_PLAN.md`
- Generator and orchestrator:
  `validation/run_survey_1d_regression.py`
- Aggregate result: `validation/survey_1d_regression/summary.json`
- Workflow catalog: `validation/survey_1d_regression/workflow_output/catalog.csv`
- Failure catalog: `validation/survey_1d_regression/workflow_output/failures.csv`
- Frozen hashes: `validation/survey_1d_regression/frozen_gate_manifest.json`

No public SDSS spectra were used. A representative public QSO sample remains
appropriate for beta scientific validation or the astronomer pilot, but is
not required to close this operational alpha gate. The complete current suite
passed 66/66 tests on 2026-07-20; formal release closeout must record the
then-current final count rather than treating 66 as permanently fixed.
