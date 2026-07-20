#!/usr/bin/env python3
"""Bounded MIRI MRS segment, sampling, and LSF-profile alpha validation."""

from __future__ import annotations

import argparse
import copy
import glob
import hashlib
import json
import time
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits

from beat.config import load_config
from beat.data import iter_spectra, linear_wavelength
from beat.fitting import fit_spectrum
from beat.injection import InjectedComponent, inject_emission_lines, score_recovery
from beat.lsf import lsf_sigma_angstrom
from beat.model import C_KMS, SQRT_2PI, gaussian_integrated
from beat.spectrum import Spectrum


PROJECT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT / "validation" / "miri_bounded_validation"
MATRIX_VERSION = "miri-bounded-v1"
SEGMENTS = ("1A", "1B", "1C", "2A", "2B", "2C", "3A", "3B", "3C", "4A", "4B", "4C")
REPRESENTATIVE_SEGMENTS = {
    "2A": "undersampled",
    "3B": "borderline",
    "4C": "well_sampled",
}
TARGETS = {
    "ic5063": {
        "config": "miri_ic5063.local.yaml",
        "glob": "/Users/tfischer/research/data/ic5063/jwst/jw02004-o003_t001_miri_ch*.fits",
    },
    "ngc4151": {
        "config": "miri_ngc4151.local.yaml",
        "glob": "/Users/tfischer/research/data/ngc4151/jwst/miri/jw02773-o002_t002_miri_ch*.fits",
    },
}
THRESHOLDS = {
    "segment_selection_correct_min": 24,
    "blank_false_positives_max": 0,
    "single_correct_min": 12,
    "double_correct_min": 5,
    "reference_evidence_flags_max": 0,
    "median_absolute_velocity_resolution_fraction_max": 0.25,
    "median_absolute_flux_fractional_error_max": 0.15,
    "median_absolute_resolved_width_fractional_error_max": 0.25,
    "mismatch_single_correct_min": 2,
    "mismatch_velocity_resolution_fraction_max": 0.25,
    "mismatch_flux_fractional_error_max": 0.15,
    "mismatch_width_fractional_error_max": 0.30,
}


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def segment_from_headers(primary: Any) -> str:
    return str(primary["CHANNEL"]) + {
        "SHORT": "A", "MEDIUM": "B", "LONG": "C"
    }[str(primary["BAND"]).strip().upper()]


def segment_headers(target: str) -> dict[str, dict[str, Any]]:
    output = {}
    for path_text in sorted(glob.glob(TARGETS[target]["glob"])):
        path = Path(path_text)
        primary = fits.getheader(path, 0)
        science = fits.getheader(path, "SCI")
        segment = segment_from_headers(primary)
        wavelength = linear_wavelength(
            science, int(science["NAXIS3"]), fits_axis=3,
            unit=science.get("CUNIT3", "micron"),
        )
        center = float(np.median(wavelength))
        delta = float(np.median(np.diff(wavelength)))
        resolving_power = 4603.0 - 128.0 * (center / 1.0e4)
        fwhm = center / resolving_power
        output[segment] = {
            "path": str(path),
            "segment": segment,
            "coverage_angstrom": [float(np.min(wavelength)), float(np.max(wavelength))],
            "center_angstrom": center,
            "spectral_step_angstrom": delta,
            "resolving_power": resolving_power,
            "lsf_fwhm_angstrom": fwhm,
            "pixels_per_lsf_fwhm": fwhm / abs(delta),
            "spatial_shape_yx": [int(science["NAXIS2"]), int(science["NAXIS1"])],
            "spectral_pixels": int(science["NAXIS3"]),
            "pipeline_version": primary.get("CAL_VER"),
            "crds_context": primary.get("CRDS_CTX"),
        }
    return output


def validation_fit(rest_center: float, rest_half_width: float, profile: str) -> dict[str, Any]:
    if profile == "screening":
        sampling = {"min_num_live_points": 40, "min_ess": 40, "dlogz": 5.0, "nsteps": 10}
    elif profile == "standard":
        sampling = {"min_num_live_points": 100, "min_ess": 200, "dlogz": 0.5, "nsteps": 20}
    else:  # pragma: no cover
        raise ValueError(profile)
    return {
        "frame": "rest",
        "wavelength_medium": "vacuum",
        "window": [rest_center - rest_half_width, rest_center + rest_half_width],
        "minimum_valid_pixels": 50,
        "continuum": {
            "degree": 1,
            "windows": [
                [rest_center - rest_half_width, rest_center - 0.55 * rest_half_width],
                [rest_center + 0.55 * rest_half_width, rest_center + rest_half_width],
            ],
            "prior_width_rms": 10.0,
        },
        "noise": {"model": "independent"},
        "kinematics": {
            "max_components": 2,
            "velocity_kms": [-1200.0, 1200.0],
            "sigma_kms": [10.0, 800.0],
        },
        "lines": [{"name": "validation_line", "wavelength": rest_center}],
        "lsf": {"model": "instrument"},
        "flux_prior": {"min_snr": 0.05, "max_signal_factor": 40.0},
        "selection": {
            "delta_logz": 5.0,
            "stop_when_not_preferred": True,
            "audit": {
                "mode": "flag",
                "sampling": {
                    "min_num_live_points": 100, "min_ess": 200,
                    "dlogz": 0.5, "stepsampler": "slice", "nsteps": 20,
                },
            },
        },
        "sampling": {
            **sampling, "show_status": False,
            "stepsampler": "slice",
        },
    }


def selection_config(target: str, details: dict[str, Any]) -> dict[str, Any]:
    config = load_config(PROJECT / "examples" / TARGETS[target]["config"])
    source = config["input"]
    source.pop("path", None)
    source.pop("segment", None)
    ny, nx = details["spatial_shape_yx"]
    source["x_range"] = [nx // 2, nx // 2 + 1]
    source["y_range"] = [ny // 2, ny // 2 + 1]
    observed_center = details["center_angstrom"]
    rest_center = observed_center / (1.0 + float(source["redshift"]))
    config["fit"] = validation_fit(rest_center, 50.0, "screening")
    return config


def audit_segment_selection() -> dict[str, Any]:
    checks = []
    inventories = {}
    for target in TARGETS:
        headers = segment_headers(target)
        inventories[target] = headers
        for segment in SEGMENTS:
            config = selection_config(target, headers[segment])
            spectrum = next(iter_spectra(config))
            selected = spectrum.metadata.get("segment")
            checks.append({
                "target": target,
                "expected_segment": segment,
                "selected_segment": selected,
                "correct": selected == segment,
                "selection_method": spectrum.metadata.get("segment_selection"),
                "input_file": spectrum.metadata.get("input_file"),
                **headers[segment],
            })
            print(f"segment {target} {segment}: selected={selected}", flush=True)
    result = {
        "matrix_version": MATRIX_VERSION,
        "n_checks": len(checks),
        "correct": int(sum(item["correct"] for item in checks)),
        "checks": checks,
        "inventories": inventories,
    }
    write_json(OUTPUT / "segment_selection.json", result)
    return result


def definitions() -> list[dict[str, Any]]:
    cases = []
    for target_index, target in enumerate(TARGETS):
        for segment_index, (segment, sampling_class) in enumerate(REPRESENTATIVE_SEGMENTS.items()):
            for kind_index, kind in enumerate(("blank", "single_narrow_pixel", "single_resolved_half", "double_resolved")):
                cases.append({
                    "name": f"{target}_{segment}_{kind}",
                    "target": target,
                    "segment": segment,
                    "sampling_class": sampling_class,
                    "kind": kind,
                    "profile_mismatch": False,
                    "noise_seed": 90_000 + target_index * 10_000 + segment_index * 1000 + kind_index * 100,
                    "sampler_seed": 4_000_000 + target_index * 10_000 + segment_index * 1000 + kind_index * 100,
                })
    for segment_index, segment in enumerate(("2A", "4C")):
        cases.append({
            "name": f"ic5063_{segment}_single_winged_profile_mismatch",
            "target": "ic5063",
            "segment": segment,
            "sampling_class": REPRESENTATIVE_SEGMENTS[segment],
            "kind": "single_winged_profile_mismatch",
            "profile_mismatch": True,
            "noise_seed": 120_000 + segment_index * 1000,
            "sampler_seed": 4_120_000 + segment_index * 1000,
        })
    return cases


def load_segment_source(target: str, segment: str, details: dict[str, Any]):
    config = selection_config(target, details)
    config["input"]["segment"] = segment
    return next(iter_spectra(config))


def build_case(case: dict[str, Any], source: Spectrum) -> tuple[Spectrum, dict[str, Any], dict[str, Any]]:
    wavelength = np.asarray(source.wavelength, dtype=float)
    delta = float(np.median(np.diff(wavelength)))
    center_index = wavelength.size // 2
    half_phase = case["kind"] == "single_resolved_half"
    observed_center = float(
        0.5 * (wavelength[center_index] + wavelength[center_index + 1])
        if half_phase else wavelength[center_index]
    )
    half_width_observed = 35.0 * abs(delta)
    selected = np.abs(wavelength - observed_center) <= half_width_observed
    wavelength = wavelength[selected]
    redshift = float(source.redshift)
    rest_center = observed_center / (1.0 + redshift)
    rest_half_width = half_width_observed / (1.0 + redshift)
    fit = validation_fit(rest_center, rest_half_width, "screening")
    metadata = dict(source.metadata)
    lsf_sigma = float(lsf_sigma_angstrom(observed_center, fit["lsf"], metadata))
    lsf_sigma_kms = lsf_sigma / observed_center * C_KMS
    resolution_fwhm_kms = 2.354820045 * lsf_sigma_kms
    pixels_per_lsf_fwhm = 2.354820045 * lsf_sigma / abs(delta)
    narrow_sigma = max(10.0, 0.25 * lsf_sigma_kms)
    resolved_sigma = max(40.0, lsf_sigma_kms)
    components: list[InjectedComponent]
    if case["kind"] == "blank":
        components = []
    elif case["kind"] == "single_narrow_pixel":
        components = [InjectedComponent(0.0, narrow_sigma, 15.0)]
    elif case["kind"] in {"single_resolved_half", "single_winged_profile_mismatch"}:
        components = [InjectedComponent(0.0, resolved_sigma, 15.0 if case["profile_mismatch"] else 10.0)]
    else:
        separation = 3.0 * resolution_fwhm_kms
        components = [
            InjectedComponent(-separation / 2.0, resolved_sigma, 15.0),
            InjectedComponent(separation / 2.0, resolved_sigma, 15.0, flux_scale=2.0 / 3.0),
        ]
    rng = np.random.default_rng(int(case["noise_seed"]))
    uncertainty = np.ones(wavelength.size, dtype=float)
    midpoint = float(np.mean(wavelength))
    baseline = 10.0 + 0.001 * (wavelength - midpoint) + rng.normal(0.0, 1.0, wavelength.size)
    if not case["profile_mismatch"]:
        flux, component_truth = inject_emission_lines(
            wavelength, baseline, uncertainty, redshift, components, fit, metadata,
            reference_line="validation_line",
        )
        profile = {"kind": "gaussian_equivalent_reference"}
    else:
        component = components[0]
        total_flux = component.peak_snr * lsf_sigma * SQRT_2PI
        core_sigma = np.hypot(lsf_sigma, observed_center * component.sigma_kms / C_KMS)
        wing_sigma = np.hypot(1.8 * lsf_sigma, observed_center * component.sigma_kms / C_KMS)
        wing_shift = 0.5 * lsf_sigma
        flux = baseline.copy()
        flux += gaussian_integrated(wavelength, observed_center, core_sigma, 0.85 * total_flux)
        flux += gaussian_integrated(wavelength, observed_center + wing_shift, wing_sigma, 0.15 * total_flux)
        component_truth = [{
            "component": 1,
            "velocity_kms": 0.0,
            "sigma_kms": float(component.sigma_kms),
            "peak_snr": float(component.peak_snr),
            "lines": {"validation_line": {
                "flux": float(total_flux),
                "observed_center_angstrom": observed_center,
                "convolved_sigma_angstrom": core_sigma,
            }},
        }]
        profile = {
            "kind": "controlled_winged_mismatch",
            "core_flux_fraction": 0.85,
            "wing_flux_fraction": 0.15,
            "wing_lsf_sigma_factor": 1.8,
            "wing_shift_lsf_sigma": 0.5,
        }
    truth = {
        "spectrum_id": case["name"],
        "stage": "controlled_actual_miri_grid",
        "seed": int(case["noise_seed"]),
        "redshift": redshift,
        "true_components": len(component_truth),
        "components": component_truth,
        "broad_components": [],
        "reference_line": "validation_line",
        "injected_profile": profile,
    }
    spectrum = Spectrum(
        spectrum_id=case["name"], wavelength=wavelength, flux=flux,
        uncertainty=uncertainty, redshift=redshift, metadata=metadata,
    )
    diagnostics = {
        "observed_center_angstrom": observed_center,
        "spectral_step_angstrom": delta,
        "lsf_sigma_angstrom": lsf_sigma,
        "lsf_sigma_kms": lsf_sigma_kms,
        "resolution_fwhm_kms": resolution_fwhm_kms,
        "pixels_per_lsf_fwhm": pixels_per_lsf_fwhm,
        "intrinsic_sigma_kms": None if not components else components[0].sigma_kms,
        "intrinsic_to_lsf_sigma": None if not components else components[0].sigma_kms / lsf_sigma_kms,
    }
    return spectrum, truth, {"fit": fit, "diagnostics": diagnostics}


def run_fit(spectrum: Spectrum, fit: dict[str, Any], seed: int) -> tuple[dict[str, Any], float, str, str]:
    np.random.seed(int(seed))
    stdout, stderr = StringIO(), StringIO()
    start = time.perf_counter()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        result = fit_spectrum(spectrum, fit)
    return result, time.perf_counter() - start, stdout.getvalue(), stderr.getvalue()


def is_reliable(result: dict[str, Any]) -> bool:
    return result.get("selection_status", "accepted") in {"accepted", "accepted_after_audit"}


def execute_case(case: dict[str, Any], source: Spectrum) -> dict[str, Any]:
    spectrum, truth, build = build_case(case, source)
    result, screen_runtime, stdout, stderr = run_fit(
        spectrum, copy.deepcopy(build["fit"]), case["sampler_seed"]
    )
    screen_score = score_recovery(result, truth)
    reasons = []
    if not screen_score["component_count_correct"]:
        reasons.append("screening component count differs from injected truth")
    if not is_reliable(result):
        reasons.extend(result.get("selection_reasons", ["screening reliability flag"]))
    standard_runtime = 0.0
    screening_result = result
    if reasons:
        standard_fit = validation_fit(
            build["fit"]["lines"][0]["wavelength"],
            0.5 * np.ptp(build["fit"]["window"]), "standard",
        )
        result, standard_runtime, std_out, std_err = run_fit(
            spectrum, standard_fit, case["sampler_seed"] + 500_000
        )
        stdout += "\n--- standard audit ---\n" + std_out
        stderr += "\n--- standard audit ---\n" + std_err
    score = score_recovery(result, truth)
    diagnostics = build["diagnostics"]
    for parameter in score.get("parameters", []):
        parameter["absolute_velocity_resolution_fraction"] = (
            abs(parameter["velocity_error_kms"]) / diagnostics["resolution_fwhm_kms"]
        )
    score.update({
        "target": case["target"], "segment": case["segment"],
        "sampling_class": case["sampling_class"], "kind": case["kind"],
        "profile_mismatch": case["profile_mismatch"],
        "selection_status": result.get("selection_status", "accepted"),
        "selection_reasons": result.get("selection_reasons", []),
        "standard_audit_performed": bool(reasons),
        "screening_runtime_seconds": screen_runtime,
        "standard_runtime_seconds": standard_runtime,
        "runtime_seconds": screen_runtime + standard_runtime,
        **diagnostics,
    })
    return {
        "matrix_version": MATRIX_VERSION,
        "case": case,
        "truth": truth,
        "fit": build["fit"],
        "screening_result": screening_result,
        "audit_reasons": reasons,
        "result": result,
        "score": score,
        "spectrum": {
            "wavelength_angstrom": spectrum.wavelength.tolist(),
            "flux": spectrum.flux.tolist(),
            "uncertainty": spectrum.uncertainty.tolist(),
        },
        "sampler_stdout_tail": stdout[-4000:],
        "sampler_stderr_tail": stderr[-4000:],
    }


def median(parameters: list[dict[str, Any]], key: str) -> float | None:
    values = [abs(float(item[key])) for item in parameters if np.isfinite(item.get(key, np.nan))]
    return float(np.median(values)) if values else None


def summarize(scores: list[dict[str, Any]], selection: dict[str, Any] | None) -> dict[str, Any]:
    reference = [score for score in scores if not score["profile_mismatch"]]
    mismatch = [score for score in scores if score["profile_mismatch"]]
    ref_parameters = [item for score in reference for item in score.get("parameters", [])]
    resolved_parameters = [
        item for score in reference if (score.get("intrinsic_to_lsf_sigma") or 0.0) >= 0.75
        for item in score.get("parameters", [])
    ]
    mismatch_parameters = [item for score in mismatch for item in score.get("parameters", [])]
    by_kind = {}
    for kind in ("blank", "single_narrow_pixel", "single_resolved_half", "double_resolved"):
        subset = [score for score in reference if score["kind"] == kind]
        by_kind[kind] = {
            "n": len(subset),
            "correct": int(sum(score["component_count_correct"] for score in subset)),
        }
    by_sampling = {}
    for sampling_class in REPRESENTATIVE_SEGMENTS.values():
        subset = [score for score in reference if score["sampling_class"] == sampling_class]
        by_sampling[sampling_class] = {
            "n": len(subset),
            "correct": int(sum(score["component_count_correct"] for score in subset)),
            "pixels_per_lsf_fwhm": sorted({round(score["pixels_per_lsf_fwhm"], 3) for score in subset}),
        }
    reference_flags = int(sum(not is_reliable_score(score) for score in reference))
    mismatch_velocity = median(mismatch_parameters, "absolute_velocity_resolution_fraction")
    mismatch_flux = median(mismatch_parameters, "flux_fractional_error")
    mismatch_width = median(mismatch_parameters, "sigma_fractional_error")
    metrics = {
        "segment_selection": {
            "n": 0 if selection is None else selection["n_checks"],
            "correct": 0 if selection is None else selection["correct"],
        },
        "reference": {
            "n": len(reference), "by_kind": by_kind, "by_sampling": by_sampling,
            "evidence_flags": reference_flags,
            "median_absolute_velocity_resolution_fraction": median(ref_parameters, "absolute_velocity_resolution_fraction"),
            "median_absolute_flux_fractional_error": median(ref_parameters, "flux_fractional_error"),
            "median_absolute_resolved_width_fractional_error": median(resolved_parameters, "sigma_fractional_error"),
        },
        "profile_mismatch": {
            "n": len(mismatch),
            "correct": int(sum(score["component_count_correct"] for score in mismatch)),
            "evidence_flags": int(sum(not is_reliable_score(score) for score in mismatch)),
            "median_absolute_velocity_resolution_fraction": mismatch_velocity,
            "median_absolute_flux_fractional_error": mismatch_flux,
            "median_absolute_width_fractional_error": mismatch_width,
        },
    }
    blank = by_kind["blank"]
    singles_correct = by_kind["single_narrow_pixel"]["correct"] + by_kind["single_resolved_half"]["correct"]
    doubles = by_kind["double_resolved"]
    gates = {
        "segment_selection": metrics["segment_selection"]["correct"] >= THRESHOLDS["segment_selection_correct_min"],
        "blank_false_positives": (blank["n"] - blank["correct"]) <= THRESHOLDS["blank_false_positives_max"],
        "single_recovery": singles_correct >= THRESHOLDS["single_correct_min"],
        "double_recovery": doubles["correct"] >= THRESHOLDS["double_correct_min"],
        "reference_evidence_reliability": reference_flags <= THRESHOLDS["reference_evidence_flags_max"],
        "reference_velocity": metric_pass(metrics["reference"]["median_absolute_velocity_resolution_fraction"], THRESHOLDS["median_absolute_velocity_resolution_fraction_max"]),
        "reference_flux": metric_pass(metrics["reference"]["median_absolute_flux_fractional_error"], THRESHOLDS["median_absolute_flux_fractional_error_max"]),
        "reference_resolved_width": metric_pass(metrics["reference"]["median_absolute_resolved_width_fractional_error"], THRESHOLDS["median_absolute_resolved_width_fractional_error_max"]),
        "mismatch_count": metrics["profile_mismatch"]["correct"] >= THRESHOLDS["mismatch_single_correct_min"],
        "mismatch_velocity": metric_pass(mismatch_velocity, THRESHOLDS["mismatch_velocity_resolution_fraction_max"]),
        "mismatch_flux": metric_pass(mismatch_flux, THRESHOLDS["mismatch_flux_fractional_error_max"]),
        "mismatch_width": metric_pass(mismatch_width, THRESHOLDS["mismatch_width_fractional_error_max"]),
    }
    return {
        "matrix_version": MATRIX_VERSION,
        "n_cases": len(scores),
        "metrics": metrics,
        "acceptance_thresholds": THRESHOLDS,
        "gate_results": gates,
        "all_bounded_gates_pass": all(gates.values()),
        "failed_gates": [name for name, passed in gates.items() if not passed],
        "runtime_seconds": float(sum(score["runtime_seconds"] for score in scores)),
    }


def metric_pass(value: float | None, threshold: float) -> bool:
    return value is not None and value <= threshold


def is_reliable_score(score: dict[str, Any]) -> bool:
    return score.get("selection_status", "accepted") in {"accepted", "accepted_after_audit"}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inventory-only", action="store_true")
    parser.add_argument("--injections-only", action="store_true")
    parser.add_argument("--only", action="append", default=[])
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    if args.inventory_only and args.injections_only:
        parser.error("choose at most one phase-only option")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    selection_path = OUTPUT / "segment_selection.json"
    selection = None
    if not args.injections_only:
        selection = audit_segment_selection()
    elif selection_path.exists():
        selection = json.loads(selection_path.read_text(encoding="utf-8"))
    if args.inventory_only:
        return 0
    cases = definitions()
    names = {case["name"] for case in cases}
    unknown = sorted(set(args.only) - names)
    if unknown:
        parser.error("unknown --only case(s): " + ", ".join(unknown))
    active = [case for case in cases if not args.only or case["name"] in args.only]
    headers = {target: segment_headers(target) for target in TARGETS}
    sources = {}
    for case in active:
        key = (case["target"], case["segment"])
        if key not in sources:
            sources[key] = load_segment_source(case["target"], case["segment"], headers[case["target"]][case["segment"]])
    for index, case in enumerate(active, start=1):
        checkpoint = OUTPUT / "cases" / case["name"] / "result.json"
        if checkpoint.exists() and not args.no_resume:
            print(f"[{index}/{len(active)}] resumed {case['name']}", flush=True)
            continue
        payload = execute_case(case, sources[(case["target"], case["segment"])])
        write_json(checkpoint, payload)
        score = payload["score"]
        print(
            f"[{index}/{len(active)}] {case['name']}: {score['true_components']} -> "
            f"{score['recovered_components']} status={score['selection_status']} "
            f"audit={score['standard_audit_performed']} ({score['runtime_seconds']:.1f}s)",
            flush=True,
        )
    scores = []
    for case in cases:
        checkpoint = OUTPUT / "cases" / case["name"] / "result.json"
        if checkpoint.exists():
            scores.append(json.loads(checkpoint.read_text(encoding="utf-8"))["score"])
    summary = summarize(scores, selection)
    write_json(OUTPUT / "summary.json", summary)
    write_json(OUTPUT / "scores.json", scores)
    write_json(OUTPUT / "manifest.json", {
        "matrix_version": MATRIX_VERSION,
        "representative_segments": REPRESENTATIVE_SEGMENTS,
        "targets": TARGETS,
        "case_definitions": cases,
        "thresholds": THRESHOLDS,
        "profile_mismatch_statement": "85% nominal Gaussian core plus 15% wing with 1.8x LSF sigma and +0.5 LSF-sigma shift; controlled stress profile, not an empirical MRS calibration",
    })
    if summary["all_bounded_gates_pass"] and len(scores) == len(cases):
        frozen_files = [
            OUTPUT / "manifest.json",
            OUTPUT / "scores.json",
            OUTPUT / "segment_selection.json",
            OUTPUT / "summary.json",
            *(OUTPUT / "cases" / case["name"] / "result.json" for case in cases),
        ]
        write_json(OUTPUT / "frozen_gate_manifest.json", {
            "matrix_version": MATRIX_VERSION,
            "status": "passed_and_frozen",
            "case_count": len(cases),
            "segment_selection_checks": selection["n_checks"] if selection else 0,
            "files": {
                str(path.relative_to(PROJECT)): sha256(path)
                for path in frozen_files
            },
            "expansion_rule": (
                "Do not add cases unless a predeclared gate later fails or an "
                "unexplained evidence flag appears inside the supported domain."
            ),
        })
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
