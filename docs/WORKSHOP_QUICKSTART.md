# BEAT at the ESAC workshop: first 20 minutes

Start with the small synthetic [O III] demonstration. It tests installation
and output creation without a FITS download; it does **not** validate a
scientific result on real data. BEAT v2 is an alpha/pilot package, so inspect
the fit diagnostics and validation limits before interpreting new targets.

## Install and check

Follow the [README installation commands](../README.md#install-and-try-the-small-demo)
for macOS/Linux or Windows PowerShell. Run all commands below from the cloned
`BEAT` repository root with its virtual environment active:

```text
beat --version
python demo/generate_demo.py
beat validate demo/demo.yaml
beat run demo/demo.yaml
```

The version command should print `BEAT 2.0.0a2`. The generator writes
`demo/demo_spectrum.txt`; `beat validate` opens it and reports the fit pixel
count, wavelength range, and noise/flux-prior summary. The run uses one worker
and tutorial-only, reduced UltraNest sampling in `demo/demo.yaml`. It should
write:

- `demo/beat_demo_output/plots/`: a PNG with data, model, components, and
  residuals;
- `demo/beat_demo_output/catalog.csv`: one row with the selected component
  count, fitted parameters, and `model.N.logz`/`model.N.logz_error` for each
  attempted count;
- `demo/beat_demo_output/results/<shard>/`: the full per-spectrum JSON result;
- `demo/beat_demo_output/run_manifest.json`: expanded configuration and run
  provenance; and
- `demo/beat_demo_output/failures.csv`: failures, if any.

An existing output can be resumed. To force a new fit, use
`beat run demo/demo.yaml --no-resume`; this is unnecessary for the first run.
`beat validate CONFIG.yaml --config-only` checks syntax/configuration without
opening a dataset, whereas plain `beat validate CONFIG.yaml` checks the first
spectrum. Always use plain validation before a new real-data run.

## Try the supplied NGC 5728 MIRI exercise

The real-data exercise uses two **separately supplied, working** YAML files
and JWST/MIRI MRS cubes: H₂ 0–0 S(1) near 17.035 µm in 3C, and [Ne III] near
15.56 µm in 3C for the supplied NGC 5728 fitting window. Although the
[Ne III] line lies near the 3B/3C overlap, that full observed window extends
beyond 3B. Check `input.segment` against each supplied YAML and cube rather
than assuming the same segment for every [Ne III] analysis. The repository
does not contain those large FITS cubes or
invent a replacement scientific configuration. BEAT's existing
[`miri` adapter](MIRI.md) accepts these cubes, selects the segment containing
the observed fit window (or uses `input.segment` if specified), and records
the selected product in the result.

In your copy of each supplied YAML, update its existing `input.path` (one
cube) **or** `input.glob` (pattern matching MIRI cubes) to point to your
local FITS files. Quote paths containing spaces or glob characters. Relative
input paths are resolved from the YAML file's directory; an absolute path is
easiest when the YAML and FITS files are in different locations. Choose a
fresh `output.directory` for each experiment so precomputed results remain
untouched.

Keep the exercise short by editing only the spatial selection. `x_range` and
`y_range` are zero-based, half-open: `x_range: [44, 45]` and
`y_range: [17, 18]` select just spaxel (44, 17); `[44, 46]` selects x=44 and
45. Pick a spaxel inside the cube's valid footprint. Avoid
`allow_full_cube: true` during the workshop: a whole-cube nested-sampling run
can take hours, well beyond the session.

From the repository root, run each supplied configuration as follows (replace
the placeholder with its actual filename):

```text
beat validate /path/to/your-ngc5728.yaml
beat run /path/to/your-ngc5728.yaml --workers 1
```

`--workers 1` overrides `output.workers` in the YAML for this invocation;
it does not modify the file. If you received **precomputed** output folders,
open their `plots/` PNGs and `catalog.csv` directly—there is no need to rerun
the cube. Use `run_manifest.json` to confirm the configuration, code version,
and input product that produced them. In the catalog, compare
`selected_components`, `selection_status`, attempted `model.N.logz` values,
component velocities/widths, and line fluxes. A flagged selection calls for
review; a neighboring spaxel need not have the same component count. See
[the full configuration reference](CONFIG_REFERENCE.md) for field details.

## Read a YAML at a glance

- `version`: configuration schema (`2`).
- `input`: FITS path or glob, adapter, redshift, and spatial ranges.
- `fit`: wavelength frame/window, lines, continuum, noise model, instrumental
  LSF, kinematic priors, evidence selection, and sampler settings. These are
  scientific choices; do not change them merely to make a workshop run fast.
- `output`: destination, resume/plot options, and worker count.

## If something goes wrong

- `beat: command not found`: activate `.venv`, then check `python -m pip show
  beat-astro` and rerun `python -m pip install -e .` from the repository root.
- `pip install -e .` says editable installs require `setup.py` or `setup.cfg`
  despite finding `pyproject.toml`: the virtual environment has an older `pip`.
  Run `python -m pip install --upgrade pip setuptools wheel` with `.venv`
  active, then retry `python -m pip install -e .`. Do not add a `setup.py`.
- PowerShell blocks `Activate.ps1`: in that PowerShell window, run
  `Set-ExecutionPolicy -Scope Process RemoteSigned`, then activate again.
- `pip install -e .` fails while compiling UltraNest: a native compiler may be
  missing or, on macOS, the Xcode/Command Line Tools license may need to be
  accepted. If you already use Conda, a prebuilt UltraNest package is an
  alternative: from the repository root, run
  `conda create -n beat-workshop -c conda-forge python=3.11 ultranest`,
  `conda activate beat-workshop`, then `python -m pip install -e .`.
  See [UltraNest installation](https://johannesbuchner.github.io/UltraNest/installation.html).
- `demo_spectrum.txt` missing: run `python demo/generate_demo.py` from the
  repository root before validation.
- FITS path not found or MIRI segment ambiguous: check the path/glob in the
  YAML, run `beat validate` again, and set the supplied configuration's
  `input.segment` only when the adapter requests it.
- Too few valid pixels: choose an in-footprint spaxel, or inspect DQ and the
  configured fitting window; do not bypass the error blindly.
- A run is unexpectedly slow: first verify the spatial ranges and worker
  count. The demo is the quick installation test; real MIRI fits are not.
- No plot: check `output.plots` is `selected` or `all`, and inspect
  `failures.csv` and the command's final summary.

For installation help during the workshop, bring the exact command, error
message, Python version (`python --version`), and YAML path.
