#!/usr/bin/env python3
"""Generate and run the bounded synthetic 1D survey regression."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import signal
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits

from beat.model import C_KMS, SQRT_2PI, gaussian_integrated, relativistic_doppler


PROJECT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT / "validation" / "survey_1d_regression"
PYTHON = Path("/Users/tfischer/opt/anaconda3/envs/astroenv/bin/python")
MATRIX_VERSION = "survey-1d-synthetic-v1"
N_VALID = 500
N_MALFORMED = 12
N_PIXELS = 321
SCIENCE_CLASSES = ("blank", "single", "double")
INTERRUPT_MINIMUM = 32


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identifiers() -> tuple[list[str], dict[str, str]]:
    valid = [f"SDSS-SYNTH-{index:06d}" for index in range(N_VALID)]
    malformed = {}
    for kind in ("nan-flux", "zero-ivar", "nan-redshift"):
        for index in range(4):
            name = f"SDSS-MALFORMED-{kind.upper()}-{index:02d}"
            malformed[name] = kind
    return valid, malformed


def rest_grid(redshift: float) -> np.ndarray:
    observed_lo = 4750.0 * (1.0 + redshift)
    observed_hi = 5150.0 * (1.0 + redshift)
    return np.power(10.0, np.linspace(np.log10(observed_lo), np.log10(observed_hi), N_PIXELS))


def add_component(
    wavelength: np.ndarray,
    flux: np.ndarray,
    redshift: float,
    velocity: float,
    sigma_kms: float,
    peak_snr: float,
    noise: float,
) -> None:
    for rest, scale in ((5006.84, 1.0), (4958.92, 0.33557), (4861.33, 0.55)):
        center = rest * (1.0 + redshift) * relativistic_doppler(velocity)
        sigma = center * sigma_kms / C_KMS
        integrated = peak_snr * noise * sigma * SQRT_2PI * scale
        flux += gaussian_integrated(wavelength, center, sigma, integrated)


def synthetic_spectrum(index: int, target_id: str, science: bool = False) -> tuple:
    redshift = 0.05 + 0.45 * ((index % 97) / 96.0)
    wavelength = rest_grid(redshift)
    rng = np.random.default_rng(71_000 + index + (10_000 if science else 0))
    noise = 1.0
    midpoint = wavelength.mean()
    flux = 10.0 + 2.0e-4 * (wavelength - midpoint) + rng.normal(0.0, noise, wavelength.size)
    truth = SCIENCE_CLASSES[index % 3] if science else ("blank", "single", "double")[index % 3]
    if truth == "single":
        add_component(wavelength, flux, redshift, 0.0, 80.0, 18.0 if science else 8.0, noise)
    elif truth == "double":
        add_component(wavelength, flux, redshift, -250.0, 70.0, 20.0 if science else 9.0, noise)
        add_component(wavelength, flux, redshift, 250.0, 90.0, 15.0 if science else 7.0, noise)
    return target_id, redshift, wavelength, flux, np.full(wavelength.size, noise**-2), truth


def write_table(path: Path, rows: list[tuple]) -> None:
    columns = [
        fits.Column(name="TARGETID", format="40A", array=[row[0] for row in rows]),
        fits.Column(name="Z", format="D", array=[row[1] for row in rows]),
        fits.Column(name="WAVE", format=f"{N_PIXELS}D", array=[row[2] for row in rows]),
        fits.Column(name="FLUX", format=f"{N_PIXELS}D", array=[row[3] for row in rows]),
        fits.Column(name="IVAR", format=f"{N_PIXELS}D", array=[row[4] for row in rows]),
        fits.Column(name="TRUTH", format="12A", array=[row[5] for row in rows]),
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    fits.HDUList([
        fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="SPECTRA")
    ]).writeto(path, overwrite=True)


def base_input(path: Path) -> dict[str, Any]:
    return {
        "kind": "survey_table", "path": str(path), "hdu": "SPECTRA",
        "id_column": "TARGETID", "redshift_column": "Z",
        "wavelength_column": "WAVE", "flux_column": "FLUX",
        "uncertainty_column": "IVAR", "uncertainty_kind": "inverse_variance",
        "wavelength_unit": "angstrom",
    }


def fit_config(max_components: int, seed: int) -> dict[str, Any]:
    return {
        "frame": "rest", "wavelength_medium": "vacuum", "window": [4750.0, 5150.0],
        "minimum_valid_pixels": 250,
        "continuum": {"degree": 0, "windows": [[4750.0, 4810.0], [5080.0, 5150.0]]},
        "noise": {"model": "independent"},
        "kinematics": {"max_components": max_components, "velocity_kms": [-600.0, 600.0], "sigma_kms": [30.0, 300.0]},
        "lines": [
            {"name": "hbeta", "wavelength": 4861.33},
            {"name": "oiii5007", "wavelength": 5006.84},
            {"name": "oiii4959", "wavelength": 4958.92, "ratio_to": "oiii5007", "ratio": 0.33557},
        ],
        "lsf": {"model": "none"},
        "selection": {"delta_logz": 5.0, "stop_when_not_preferred": True, "audit": {"mode": "none"}},
        "sampling": {"min_num_live_points": 40, "min_ess": 40, "dlogz": 5.0, "show_status": False, "stepsampler": "none", "seed": seed},
    }


def config(input_path: Path, output_path: Path, max_components: int, seed: int) -> dict[str, Any]:
    return {
        "version": 2, "input": base_input(input_path),
        "fit": fit_config(max_components, seed),
        "output": {"directory": str(output_path), "resume": True, "plots": "none", "workers": 4, "progress_every": 25},
    }


def generate() -> None:
    valid, malformed = identifiers()
    rows = [synthetic_spectrum(index, target_id) for index, target_id in enumerate(valid)]
    for offset, (target_id, kind) in enumerate(malformed.items()):
        row = list(synthetic_spectrum(N_VALID + offset, target_id))
        if kind == "nan-flux":
            row[3] = np.full(N_PIXELS, np.nan)
        elif kind == "zero-ivar":
            row[4] = np.zeros(N_PIXELS)
        else:
            row[1] = np.nan
        row[5] = kind
        rows.append(tuple(row))
    workflow_table = OUTPUT / "data" / "sdss_like_workflow.fits"
    write_table(workflow_table, rows)
    science_rows = [
        synthetic_spectrum(index, f"SDSS-SCIENCE-{index:03d}", science=True)
        for index in range(4)
    ]
    science_table = OUTPUT / "data" / "sdss_like_science.fits"
    write_table(science_table, science_rows)
    write_json(OUTPUT / "workflow_config.json", config(workflow_table, OUTPUT / "workflow_output", 0, 24_681))
    science_config = config(science_table, OUTPUT / "science_output", 2, 13_579)
    science_config["output"]["workers"] = 2
    write_json(OUTPUT / "science_config.json", science_config)
    single_config = config(science_table, OUTPUT / "science_single_output", 1, 13_579)
    single_config["input"]["row_start"] = 1
    single_config["input"]["row_stop"] = 2
    single_config["output"]["workers"] = 1
    write_json(OUTPUT / "science_single_config.json", single_config)
    pilot_config = config(science_table, OUTPUT / "science_pilot_output", 2, 13_579)
    pilot_config["input"]["row_stop"] = 1
    pilot_config["output"]["workers"] = 1
    write_json(OUTPUT / "science_pilot_config.json", pilot_config)
    write_json(OUTPUT / "dataset_manifest.json", {
        "matrix_version": MATRIX_VERSION, "workflow_rows": len(rows),
        "valid_ids": valid, "malformed_ids": malformed,
        "science_truth": {row[0]: row[5] for row in science_rows},
        "workflow_table_sha256": sha256(workflow_table),
        "science_table_sha256": sha256(science_table),
    })
    print(f"generated {len(rows)} workflow rows and {len(science_rows)} science rows")


def result_files(output: Path) -> list[Path]:
    return sorted((output / "results").rglob("*.json")) if (output / "results").exists() else []


def command(config_path: Path, workers: int) -> list[str]:
    return [str(PYTHON), "-m", "beat.cli", "run", str(config_path), "--workers", str(workers)]


def environment() -> dict[str, str]:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(PROJECT / "src")
    env["MPLCONFIGDIR"] = "/tmp/beat-mpl"
    return env


def process_tree_rss(pid: int) -> int | None:
    try:
        import psutil
        process = psutil.Process(pid)
        return process.memory_info().rss + sum(
            child.memory_info().rss for child in process.children(recursive=True)
            if child.is_running()
        )
    except Exception:
        return None


def interrupt() -> None:
    config_path = OUTPUT / "workflow_config.json"
    output = OUTPUT / "workflow_output"
    log_path = OUTPUT / "workflow_interrupt.log"
    started = time.time()
    peak_rss = 0
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command(config_path, 4), cwd=PROJECT, env=environment(),
            stdout=log, stderr=subprocess.STDOUT, start_new_session=True,
        )
        while process.poll() is None:
            count = len(result_files(output))
            rss = process_tree_rss(process.pid)
            if rss is not None:
                peak_rss = max(peak_rss, rss)
            if count >= INTERRUPT_MINIMUM:
                os.killpg(process.pid, signal.SIGTERM)
                break
            if time.time() - started > 900:
                os.killpg(process.pid, signal.SIGTERM)
                raise RuntimeError("workflow interruption phase exceeded 15 minutes")
            time.sleep(0.2)
        process.wait(timeout=30)
    files = result_files(output)
    snapshot = {str(path.relative_to(output)): sha256(path) for path in files}
    write_json(OUTPUT / "interruption.json", {
        "checkpoint_count": len(files), "checkpoint_hashes": snapshot,
        "returncode": process.returncode, "peak_process_tree_rss_bytes": peak_rss or None,
        "runtime_seconds": time.time() - started,
    })
    print(f"interrupted after {len(files)} successful checkpoints")


def run_command(config_path: Path, workers: int, log_name: str, timeout: int = 1200) -> dict[str, Any]:
    started = time.time()
    completed = subprocess.run(
        command(config_path, workers), cwd=PROJECT, env=environment(),
        capture_output=True, text=True, timeout=timeout,
    )
    (OUTPUT / log_name).write_text(completed.stdout + completed.stderr, encoding="utf-8")
    return {"returncode": completed.returncode, "runtime_seconds": time.time() - started}


def resume() -> None:
    interrupted = json.loads((OUTPUT / "interruption.json").read_text(encoding="utf-8"))
    first = run_command(OUTPUT / "workflow_config.json", 4, "workflow_resume.log")
    output = OUTPUT / "workflow_output"
    preserved = {
        relative: sha256(output / relative) == expected
        for relative, expected in interrupted["checkpoint_hashes"].items()
    }
    manifest_four = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    hashes_before_one = {str(path.relative_to(output)): sha256(path) for path in result_files(output)}
    second = run_command(OUTPUT / "workflow_config.json", 1, "workflow_one_worker_resume.log")
    hashes_after_one = {str(path.relative_to(output)): sha256(path) for path in result_files(output)}
    manifest_one = json.loads((output / "run_manifest.json").read_text(encoding="utf-8"))
    write_json(OUTPUT / "resume_checks.json", {
        "four_worker_run": first, "one_worker_run": second,
        "preinterruption_checkpoints_preserved": all(preserved.values()),
        "preservation_details": preserved,
        "successful_hashes_unchanged_by_one_worker_resume": hashes_before_one == hashes_after_one,
        "four_worker_config_hash": manifest_four["config_hash"],
        "one_worker_config_hash": manifest_one["config_hash"],
        "four_worker_effective_workers": manifest_four["execution"]["workers"],
        "one_worker_effective_workers": manifest_one["execution"]["workers"],
        "four_worker_summary": manifest_four.get("summary"),
        "one_worker_summary": manifest_one.get("summary"),
    })
    print(f"resume return codes: workers4={first['returncode']} workers1={second['returncode']}")


def science(pilot: bool = False) -> None:
    name = "science_pilot_config.json" if pilot else "science_config.json"
    log = "science_pilot.log" if pilot else "science.log"
    result = run_command(OUTPUT / name, 1 if pilot else 2, log, timeout=1800)
    write_json(OUTPUT / ("science_pilot_run.json" if pilot else "science_run.json"), result)
    print(json.dumps(result, indent=2))


def science_single() -> None:
    result = run_command(OUTPUT / "science_single_config.json", 1, "science_single.log", timeout=1200)
    write_json(OUTPUT / "science_single_run.json", result)
    print(json.dumps(result, indent=2))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def summarize() -> dict[str, Any]:
    dataset = json.loads((OUTPUT / "dataset_manifest.json").read_text(encoding="utf-8"))
    interrupted = json.loads((OUTPUT / "interruption.json").read_text(encoding="utf-8"))
    resumed = json.loads((OUTPUT / "resume_checks.json").read_text(encoding="utf-8"))
    workflow = OUTPUT / "workflow_output"
    workflow_manifest = json.loads((workflow / "run_manifest.json").read_text(encoding="utf-8"))
    catalog = read_csv(workflow / "catalog.csv")
    failures = read_csv(workflow / "failures.csv")
    result_payloads = [json.loads(path.read_text(encoding="utf-8")) for path in result_files(workflow)]
    science_output = OUTPUT / "science_output"
    single_output = OUTPUT / "science_single_output"
    science_payloads = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in result_files(science_output) + result_files(single_output)
    ]
    science_failures = list((science_output / "failure_details").rglob("*.json")) if (science_output / "failure_details").exists() else []
    science_failures += list((single_output / "failure_details").rglob("*.json")) if (single_output / "failure_details").exists() else []
    science_counts = sorted({int(payload["selected_components"]) for payload in science_payloads})
    malformed_expected = set(dataset["malformed_ids"])
    gates = {
        "interrupted_midstream": INTERRUPT_MINIMUM <= interrupted["checkpoint_count"] < N_VALID,
        "preinterruption_hashes_preserved": resumed["preinterruption_checkpoints_preserved"],
        "catalog_count": len(catalog) == N_VALID,
        "catalog_ids_unique": len({row["spectrum_id"] for row in catalog}) == N_VALID,
        "malformed_failures_exact": {row["spectrum_id"] for row in failures} == malformed_expected,
        "valid_checkpoints_preserved_on_worker_override": resumed["successful_hashes_unchanged_by_one_worker_resume"],
        "fingerprint_stable": resumed["four_worker_config_hash"] == resumed["one_worker_config_hash"] == workflow_manifest["config_hash"],
        "worker_override_recorded": resumed["four_worker_effective_workers"] == 4 and resumed["one_worker_effective_workers"] == 1,
        "deterministic_seeds_recorded": all(payload.get("sampler_seed") is not None for payload in result_payloads),
        "memory_bounded": interrupted["peak_process_tree_rss_bytes"] is None or interrupted["peak_process_tree_rss_bytes"] <= int(1.5 * 1024**3),
        "science_results_complete": len({payload["spectrum_id"] for payload in science_payloads}) == 4 and not science_failures,
        "science_exercises_all_component_counts": science_counts == [0, 1, 2],
    }
    summary = {
        "matrix_version": MATRIX_VERSION, "gates": gates,
        "all_gates_pass": all(gates.values()),
        "workflow": {
            "catalog_rows": len(catalog), "failure_rows": len(failures),
            "interrupted_checkpoint_count": interrupted["checkpoint_count"],
            "peak_process_tree_rss_bytes": interrupted["peak_process_tree_rss_bytes"],
            "config_hash": workflow_manifest["config_hash"],
        },
        "science": {"result_rows": len(science_payloads), "failure_rows": len(science_failures), "selected_component_counts": science_counts},
        "supported_statement": "The 1D survey workflow is operationally validated on a synthetic SDSS-like regression dataset; scientific completeness for a real SDSS QSO population is not established.",
    }
    write_json(OUTPUT / "summary.json", summary)
    if summary["all_gates_pass"]:
        files = [
            OUTPUT / "dataset_manifest.json", OUTPUT / "workflow_config.json",
            OUTPUT / "science_config.json", OUTPUT / "science_single_config.json", OUTPUT / "interruption.json",
            OUTPUT / "resume_checks.json", OUTPUT / "summary.json",
            workflow / "run_manifest.json", workflow / "catalog.csv", workflow / "failures.csv",
            single_output / "run_manifest.json", single_output / "catalog.csv",
        ] + result_files(workflow) + result_files(science_output) + result_files(single_output)
        write_json(OUTPUT / "frozen_gate_manifest.json", {
            "matrix_version": MATRIX_VERSION, "status": "passed_and_frozen",
            "files": {str(path.relative_to(PROJECT)): sha256(path) for path in files},
        })
    print(json.dumps(summary, indent=2))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("phase", choices=("generate", "interrupt", "resume", "science-pilot", "science", "science-single", "summarize"))
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    if args.phase == "generate": generate()
    elif args.phase == "interrupt": interrupt()
    elif args.phase == "resume": resume()
    elif args.phase == "science-pilot": science(True)
    elif args.phase == "science": science(False)
    elif args.phase == "science-single": science_single()
    else: summarize()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
