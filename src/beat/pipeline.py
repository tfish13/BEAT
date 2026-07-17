"""Resumable, race-free batch execution for BEAT datasets."""

from __future__ import annotations

import csv
import hashlib
import json
import os
import platform
import re
import traceback
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator

import numpy as np

from . import __version__
from .config import public_config
from .data import iter_spectra
from .fitting import fit_spectrum, make_diagnostic_plot, selected_model
from .spectrum import Spectrum


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(f"Cannot serialize {type(value).__name__}")


def config_fingerprint(config: dict[str, Any]) -> str:
    # Operational choices (worker count, plots, resume mode, output path) do not
    # change the scientific result and therefore must not invalidate checkpoints.
    public = public_config(config)
    scientific = {
        "version": public.get("version"),
        "input": public.get("input"),
        "fit": public.get("fit"),
    }
    payload = json.dumps(
        scientific, sort_keys=True, separators=(",", ":"), default=_json_default
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def safe_stem(identifier: str) -> str:
    readable = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(identifier)).strip("._")
    readable = readable[:80] or "spectrum"
    digest = hashlib.sha1(str(identifier).encode("utf-8")).hexdigest()[:10]
    return f"{readable}-{digest}"


def _atomic_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + f".tmp-{os.getpid()}")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )
    os.replace(temporary, path)


def _result_path(results_dir: Path, identifier: str) -> Path:
    digest = hashlib.sha1(str(identifier).encode("utf-8")).hexdigest()
    return results_dir / digest[:2] / f"{safe_stem(identifier)}.json"


def _matching_completed(path: Path, fingerprint: str) -> bool:
    if not path.exists():
        return False
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return False
    return payload.get("config_hash") == fingerprint and payload.get("status") == "ok"


def _worker(
    spectrum: Spectrum,
    fit: dict[str, Any],
    plots: str,
    plots_dir: str,
) -> dict[str, Any]:
    try:
        result = fit_spectrum(spectrum, fit)
        if plots == "selected":
            plot_path = Path(plots_dir) / f"{safe_stem(spectrum.spectrum_id)}.png"
            make_diagnostic_plot(spectrum, fit, result, plot_path)
            result["plot"] = str(plot_path)
        elif plots == "all":
            paths = []
            for model in result["models"]:
                n_components = int(model["n_components"])
                plot_path = (
                    Path(plots_dir)
                    / f"{safe_stem(spectrum.spectrum_id)}-n{n_components}.png"
                )
                make_diagnostic_plot(
                    spectrum, fit, result, plot_path, n_components=n_components
                )
                paths.append(str(plot_path))
            result["plots"] = paths
        return result
    except Exception as exc:  # return failures to the parent process
        return {
            "status": "failed",
            "spectrum_id": spectrum.spectrum_id,
            "redshift": spectrum.redshift,
            "metadata": spectrum.metadata,
            "error_type": type(exc).__name__,
            "error": str(exc),
            "traceback": traceback.format_exc(),
        }


def flatten_result(result: dict[str, Any]) -> dict[str, Any]:
    """Flatten the selected model into one catalogue row."""
    row: dict[str, Any] = {
        "spectrum_id": result["spectrum_id"],
        "redshift": result["redshift"],
        "selected_components": result["selected_components"],
        "selected_logz": result["selected_logz"],
        "selection_status": result.get("selection_status", "not_recorded"),
        "selection_reasons": "; ".join(result.get("selection_reasons", [])),
        "selection_audit_performed": result.get("selection_audit", {}).get(
            "performed", False
        ),
        "selection_audit_reasons": "; ".join(
            result.get("selection_audit", {}).get("trigger_reasons", [])
        ),
        "noise_level": result["noise_level"],
        "n_input_pixels": result["n_input_pixels"],
        "n_fit_pixels": result["n_fit_pixels"],
    }
    for key, value in result.get("metadata", {}).items():
        if np.isscalar(value) or value is None:
            row[f"meta.{key}"] = value
    for model in result["models"]:
        n_components = model["n_components"]
        row[f"model.{n_components}.logz"] = model["logz"]
        row[f"model.{n_components}.logz_error"] = model["logz_error"]

    chosen = selected_model(result)
    for name, median, stdev in zip(
        chosen["parameter_names"],
        chosen["posterior_median"],
        chosen["posterior_stdev"],
    ):
        row[name] = median
        row[f"{name}.stdev"] = stdev
    for component in result.get("components", []):
        prefix = f"component.{component['component']}"
        for line_name, values in component["lines"].items():
            line_prefix = f"{prefix}.{line_name}"
            row[f"{line_prefix}.flux"] = values["flux"]
            row[f"{line_prefix}.flux.stdev"] = values["flux_stdev"]
            row[f"{line_prefix}.observed_center_angstrom"] = values[
                "observed_center_angstrom"
            ]
            row[f"{line_prefix}.observed_sigma_angstrom"] = values[
                "observed_sigma_angstrom"
            ]
            for field in (
                "intrinsic_sigma_angstrom",
                "lsf_sigma_angstrom",
                "convolved_sigma_angstrom",
            ):
                if field in values:
                    row[f"{line_prefix}.{field}"] = values[field]
    return row


def _matching_json_payloads(directory: Path, fingerprint: str) -> Iterator[dict[str, Any]]:
    for path in sorted(directory.rglob("*.json")):
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if payload.get("config_hash") == fingerprint:
            yield payload


def _write_catalog(output_dir: Path, fingerprint: str) -> int:
    # Two streaming passes avoid keeping a large survey catalogue in memory.
    columns: list[str] = []
    known_columns: set[str] = set()
    count = 0
    for payload in _matching_json_payloads(output_dir / "results", fingerprint):
        if payload.get("status") != "ok":
            continue
        row = flatten_result(payload)
        count += 1
        for column in row:
            if column not in known_columns:
                known_columns.add(column)
                columns.append(column)

    path = output_dir / "catalog.csv"
    temporary = path.with_suffix(".csv.tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        if columns:
            writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            for payload in _matching_json_payloads(output_dir / "results", fingerprint):
                if payload.get("status") == "ok":
                    writer.writerow(flatten_result(payload))
    os.replace(temporary, path)
    return count


def _write_failures(output_dir: Path, fingerprint: str) -> int:
    path = output_dir / "failures.csv"
    temporary = path.with_suffix(".csv.tmp")
    columns = ["spectrum_id", "redshift", "error_type", "error"]
    count = 0
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        for payload in _matching_json_payloads(output_dir / "failure_details", fingerprint):
            writer.writerow(payload)
            count += 1
    os.replace(temporary, path)
    return count


def _selection_summary(output_dir: Path, fingerprint: str) -> tuple[dict[str, int], int]:
    counts: dict[str, int] = {}
    audits_performed = 0
    for payload in _matching_json_payloads(output_dir / "results", fingerprint):
        if payload.get("status") != "ok":
            continue
        status = str(payload.get("selection_status", "not_recorded"))
        counts[status] = counts.get(status, 0) + 1
        audits_performed += int(
            bool(payload.get("selection_audit", {}).get("performed", False))
        )
    return counts, audits_performed


def _write_selection_review(output_dir: Path, fingerprint: str) -> int:
    """Write a compact candidate list for selections needing convergence review."""
    path = output_dir / "selection_review.csv"
    temporary = path.with_suffix(".csv.tmp")
    columns = [
        "spectrum_id",
        "redshift",
        "selected_components",
        "selection_status",
        "selection_reasons",
        "audit_performed",
        "reference_components",
        "candidate_components",
        "delta_logz",
        "combined_logz_error",
        "threshold",
    ]
    count = 0
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for payload in _matching_json_payloads(output_dir / "results", fingerprint):
            if payload.get("status") != "ok" or payload.get("selection_status") not in {
                "ambiguous",
                "convergence_unverified",
            }:
                continue
            comparisons = payload.get("selection_diagnostics", [])
            closest = (
                min(
                    comparisons,
                    key=lambda item: float(item.get("distance_from_threshold", np.inf)),
                )
                if comparisons
                else {}
            )
            writer.writerow(
                {
                    "spectrum_id": payload["spectrum_id"],
                    "redshift": payload["redshift"],
                    "selected_components": payload["selected_components"],
                    "selection_status": payload["selection_status"],
                    "selection_reasons": "; ".join(
                        payload.get("selection_reasons", [])
                    ),
                    "audit_performed": payload.get("selection_audit", {}).get(
                        "performed", False
                    ),
                    "reference_components": closest.get("reference_components"),
                    "candidate_components": closest.get("candidate_components"),
                    "delta_logz": closest.get("delta_logz"),
                    "combined_logz_error": closest.get("combined_logz_error"),
                    "threshold": closest.get("threshold"),
                }
            )
            count += 1
    os.replace(temporary, path)
    return count


def _persist_result(
    result: dict[str, Any],
    output_dir: Path,
    fingerprint: str,
) -> None:
    result["config_hash"] = fingerprint
    result_path = _result_path(output_dir / "results", result["spectrum_id"])
    failure_path = _result_path(output_dir / "failure_details", result["spectrum_id"])
    if result.get("status") == "ok":
        result_path.parent.mkdir(parents=True, exist_ok=True)
        _atomic_json(result_path, result)
        failure_path.unlink(missing_ok=True)
    else:
        failure_path.parent.mkdir(parents=True, exist_ok=True)
        _atomic_json(failure_path, result)
        result_path.unlink(missing_ok=True)


def _iter_unique_spectra(config: dict[str, Any]) -> Iterator[Spectrum]:
    identifiers: set[str] = set()
    for spectrum in iter_spectra(config):
        identifier = str(spectrum.spectrum_id)
        if identifier in identifiers:
            raise ValueError(f"Duplicate spectrum_id in input: {identifier}")
        identifiers.add(identifier)
        yield spectrum


def run_pipeline(config: dict[str, Any], workers: int | None = None) -> dict[str, Any]:
    """Fit the configured dataset, writing atomic per-spectrum checkpoints."""
    output_config = config["output"]
    output_dir = Path(output_config["directory"])
    results_dir = output_dir / "results"
    plots_dir = output_dir / "plots"
    failure_details_dir = output_dir / "failure_details"
    results_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)
    failure_details_dir.mkdir(parents=True, exist_ok=True)

    fingerprint = config_fingerprint(config)
    started = datetime.now(timezone.utc).isoformat()
    manifest = {
        "beat_version": __version__,
        "config_hash": fingerprint,
        "started_utc": started,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "config": public_config(config),
    }
    _atomic_json(output_dir / "run_manifest.json", manifest)

    if workers is None:
        workers = int(output_config.get("workers", max(1, (os.cpu_count() or 2) - 1)))
    if workers < 1:
        raise ValueError("workers must be at least one")
    resume = bool(output_config.get("resume", True))
    plots = output_config.get("plots", "selected")
    progress_every = max(1, int(output_config.get("progress_every", 10)))
    completed = 0
    failed = 0
    resumed = 0
    discovered = 0

    spectra = _iter_unique_spectra(config)

    def should_skip(spectrum: Spectrum) -> bool:
        return resume and _matching_completed(
            _result_path(results_dir, spectrum.spectrum_id), fingerprint
        )

    if workers == 1:
        for spectrum in spectra:
            discovered += 1
            if should_skip(spectrum):
                resumed += 1
                continue
            result = _worker(spectrum, config["fit"], plots, str(plots_dir))
            _persist_result(result, output_dir, fingerprint)
            if result.get("status") == "ok":
                completed += 1
            else:
                failed += 1
            processed = completed + failed
            if processed % progress_every == 0:
                print(f"Processed {processed} spectrum/spectra ({failed} failed)")
    else:
        pending: dict[Any, str] = {}
        exhausted = False

        def submit_one(executor: ProcessPoolExecutor) -> bool:
            nonlocal discovered, resumed, exhausted
            if exhausted:
                return False
            for spectrum in spectra:
                discovered += 1
                if should_skip(spectrum):
                    resumed += 1
                    continue
                future = executor.submit(
                    _worker, spectrum, config["fit"], plots, str(plots_dir)
                )
                pending[future] = spectrum.spectrum_id
                return True
            exhausted = True
            return False

        with ProcessPoolExecutor(max_workers=workers) as executor:
            for _ in range(2 * workers):
                if not submit_one(executor):
                    break
            while pending:
                done, _ = wait(pending, return_when=FIRST_COMPLETED)
                for future in done:
                    identifier = pending.pop(future)
                    try:
                        result = future.result()
                    except Exception as exc:  # executor-level failure
                        result = {
                            "status": "failed",
                            "spectrum_id": identifier,
                            "redshift": np.nan,
                            "error_type": type(exc).__name__,
                            "error": str(exc),
                            "traceback": traceback.format_exc(),
                        }
                    _persist_result(result, output_dir, fingerprint)
                    if result.get("status") == "ok":
                        completed += 1
                    else:
                        failed += 1
                    processed = completed + failed
                    if processed % progress_every == 0:
                        print(
                            f"Processed {processed} spectrum/spectra "
                            f"({failed} failed)"
                        )
                    submit_one(executor)

    catalog_rows = _write_catalog(output_dir, fingerprint)
    failure_rows = _write_failures(output_dir, fingerprint)
    selection_status_counts, selection_audits_performed = _selection_summary(
        output_dir, fingerprint
    )
    selection_review_rows = _write_selection_review(output_dir, fingerprint)
    manifest["finished_utc"] = datetime.now(timezone.utc).isoformat()
    manifest["summary"] = {
        "discovered": discovered,
        "completed_this_run": completed,
        "resumed": resumed,
        "failed_this_run": failed,
        "failure_rows": failure_rows,
        "catalog_rows": catalog_rows,
        "selection_status_counts": selection_status_counts,
        "selection_audits_performed": selection_audits_performed,
        "selection_review_rows": selection_review_rows,
    }
    _atomic_json(output_dir / "run_manifest.json", manifest)
    return manifest["summary"]
