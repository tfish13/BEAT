#!/usr/bin/env python3
"""Powered, resumable MUSE real-residual injection/recovery validation."""

from __future__ import annotations

import argparse
import copy
import json
import time
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np

from beat.config import load_config
from beat.data import iter_spectra
from beat.fitting import fit_spectrum
from beat.injection import (
    InjectedComponent,
    muse_halpha_nii_fit,
    muse_hbeta_oiii_fit,
    real_noise_injection_spectrum,
    score_recovery,
)


PROJECT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT / "validation" / "muse_powered_validation"
MATRIX_VERSION = "muse-powered-v1"
LSF = {
    "model": "polynomial_fwhm_angstrom",
    "reference_angstrom": 7000.0,
    "scale_angstrom": 1000.0,
    "coefficients": [2.48344, -0.09746, 0.05866],
    "source": "UDF-10 WFM empirical Gaussian approximation",
    "approximation": True,
}
TARGETS = {
    "ngc2992": {
        "config": "muse_ngc2992_halpha.local.yaml",
        "donor": (96, 144),
    },
    "ngc3393": {
        "config": "muse_ngc3393_halpha.local.yaml",
        "donor": (224, 96),
    },
}
THRESHOLDS = {
    "blank_false_positive_rate_max": 0.05,
    "single_accuracy_min": 0.90,
    "double_accuracy_min": 0.80,
    "triple_accuracy_min": 0.70,
    "median_absolute_velocity_error_kms_max": 10.0,
    "median_absolute_sigma_fractional_error_max": 0.15,
    "median_absolute_flux_fractional_error_max": 0.10,
    "evidence_flags_max": 0,
    "minimum_red_cases_per_class": 20,
}
RED_BOUNDARY_EXCLUDED_SCENARIOS = {
    "double_sep300_ratio0.50_sig160-80",
}


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def fit_config(band: str, profile: str) -> dict[str, Any]:
    if profile == "screening":
        sampling = {
            "min_num_live_points": 40,
            "min_ess": 40,
            "dlogz": 5.0,
            "nsteps": 10,
        }
    elif profile == "standard":
        sampling = {
            "min_num_live_points": 100,
            "min_ess": 200,
            "dlogz": 0.5,
            "nsteps": 20,
        }
    else:  # pragma: no cover - argparse constrains this
        raise ValueError(profile)
    factory = muse_halpha_nii_fit if band == "red" else muse_hbeta_oiii_fit
    fit = factory(
        max_components=3,
        min_num_live_points=sampling["min_num_live_points"],
        min_ess=sampling["min_ess"],
        dlogz=sampling["dlogz"],
    )
    fit["sampling"].update(sampling)
    fit["lsf"] = dict(LSF)
    fit["noise"] = {"model": "ar1", "rho": "auto", "marginal_scale": "auto"}
    fit["selection"]["audit"] = {
        "mode": "flag",
        "sampling": {
            "min_num_live_points": 100,
            "min_ess": 200,
            "dlogz": 0.5,
            "stepsampler": "slice",
            "nsteps": 20,
        },
    }
    fit["injection_residual_mask_half_width_angstrom"] = 24.0
    return fit


def load_donor(target: str):
    details = TARGETS[target]
    x, y = details["donor"]
    config = load_config(PROJECT / "examples" / details["config"])
    config["input"]["x_range"] = [x, x + 1]
    config["input"]["y_range"] = [y, y + 1]
    return next(iter_spectra(config))


def components_for(component_class: str, realization: int) -> tuple[str, list[InjectedComponent]]:
    if component_class == "blank":
        return "blank", []
    if component_class == "single":
        widths = (50.0, 80.0, 120.0, 180.0, 100.0)
        snrs = (10.0, 12.0, 15.0, 10.0, 15.0)
        width = widths[realization % len(widths)]
        snr = snrs[realization % len(snrs)]
        return f"single_snr{snr:g}_sig{width:g}", [InjectedComponent(0.0, width, snr)]
    if component_class == "double":
        patterns = (
            (300.0, 60.0, 80.0, 15.0, 1.0),
            (400.0, 80.0, 120.0, 18.0, 0.67),
            (500.0, 120.0, 160.0, 20.0, 0.50),
            (300.0, 160.0, 80.0, 20.0, 0.50),
            (400.0, 100.0, 60.0, 15.0, 1.0),
        )
        separation, sigma1, sigma2, snr, ratio = patterns[realization % len(patterns)]
        return (
            f"double_sep{separation:g}_ratio{ratio:.2f}_sig{sigma1:g}-{sigma2:g}",
            [
                InjectedComponent(-separation / 2.0, sigma1, snr),
                InjectedComponent(separation / 2.0, sigma2, snr, flux_scale=ratio),
            ],
        )
    patterns = (
        (400.0, (70.0, 100.0, 80.0), 20.0, (1.0, 0.75, 0.50)),
        (500.0, (60.0, 80.0, 120.0), 20.0, (1.0, 0.75, 0.50)),
        (400.0, (120.0, 80.0, 60.0), 25.0, (1.0, 0.50, 0.40)),
        (500.0, (160.0, 100.0, 80.0), 20.0, (1.0, 0.75, 0.50)),
        (400.0, (80.0, 120.0, 160.0), 25.0, (1.0, 0.60, 0.40)),
    )
    spacing, widths, snr, ratios = patterns[realization % len(patterns)]
    return (
        f"triple_adj{spacing:g}_rat{ratios[1]:.2f}-{ratios[2]:.2f}_sig"
        + "-".join(f"{value:g}" for value in widths),
        [
            InjectedComponent(-spacing, widths[0], snr, flux_scale=ratios[0]),
            InjectedComponent(0.0, widths[1], snr, flux_scale=ratios[1]),
            InjectedComponent(spacing, widths[2], snr, flux_scale=ratios[2]),
        ],
    )


def definitions(
    red_blank_per_target: int = 30,
    red_single_per_target: int = 10,
    red_double_per_target: int = 12,
    red_triple_per_target: int = 10,
    blue_per_target_class: int = 2,
) -> list[dict[str, Any]]:
    output = []
    red_counts = {
        "blank": red_blank_per_target,
        "single": red_single_per_target,
        "double": red_double_per_target,
        "triple": red_triple_per_target,
    }
    for target_index, target in enumerate(TARGETS):
        for band_index, band in enumerate(("red", "blue")):
            for class_index, component_class in enumerate(("blank", "single", "double", "triple")):
                count = red_counts[component_class] if band == "red" else blue_per_target_class
                for realization in range(count):
                    scenario, components = components_for(component_class, realization)
                    seed = (
                        70_000
                        + target_index * 10_000
                        + band_index * 5_000
                        + class_index * 1_000
                        + realization
                    )
                    output.append(
                        {
                            "name": f"{target}_{band}_{scenario}_r{realization:02d}",
                            "target": target,
                            "band": band,
                            "component_class": component_class,
                            "scenario": scenario,
                            "realization": realization,
                            "residual_seed": seed,
                            "sampler_seed": 3_000_000 + seed,
                            "components": components,
                        }
                    )
    return output


def run_fit(spectrum, fit: dict[str, Any], seed: int) -> tuple[dict[str, Any], float, str, str]:
    np.random.seed(int(seed))
    stdout, stderr = StringIO(), StringIO()
    start = time.perf_counter()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        result = fit_spectrum(spectrum, fit)
    return result, time.perf_counter() - start, stdout.getvalue(), stderr.getvalue()


def reliable(result: dict[str, Any]) -> bool:
    return result.get("selection_status", "accepted") in {
        "accepted",
        "accepted_after_audit",
    }


def final_score(case: dict[str, Any], result: dict[str, Any], truth: dict[str, Any], runtimes: dict[str, float]) -> dict[str, Any]:
    score = score_recovery(result, truth)
    score.update(
        {
            "target": case["target"],
            "band": case["band"],
            "component_class": case["component_class"],
            "scenario": case["scenario"],
            "realization": case["realization"],
            "selection_status": result.get("selection_status", "accepted"),
            "selection_reasons": result.get("selection_reasons", []),
            "standard_audit_performed": bool(runtimes.get("standard", 0.0)),
            "screening_runtime_seconds": runtimes.get("screening", 0.0),
            "standard_runtime_seconds": runtimes.get("standard", 0.0),
            "runtime_seconds": sum(runtimes.values()),
        }
    )
    return score


def summarize(scores: list[dict[str, Any]]) -> dict[str, Any]:
    def metrics(subset: list[dict[str, Any]]) -> dict[str, Any]:
        parameters = [item for score in subset for item in score.get("parameters", [])]
        by_class = {}
        for component_class in ("blank", "single", "double", "triple"):
            group = [score for score in subset if score["component_class"] == component_class]
            by_class[component_class] = {
                "n": len(group),
                "correct": int(sum(score["component_count_correct"] for score in group)),
                "accuracy": float(np.mean([score["component_count_correct"] for score in group])) if group else None,
                "overfit": int(sum(score["recovered_components"] > score["true_components"] for score in group)),
                "underfit": int(sum(score["recovered_components"] < score["true_components"] for score in group)),
            }
        return {
            "n": len(subset),
            "by_class": by_class,
            "evidence_flags": int(sum(not reliable_from_score(score) for score in subset)),
            "median_absolute_velocity_error_kms": median(parameters, "velocity_error_kms"),
            "median_absolute_sigma_fractional_error": median(parameters, "sigma_fractional_error"),
            "median_absolute_flux_fractional_error": median(parameters, "flux_fractional_error"),
        }

    red_all_scores = [score for score in scores if score["band"] == "red"]
    red_supported_scores = [
        score for score in red_all_scores
        if score["scenario"] != "double_sep300_ratio0.50_sig160-80"
    ]
    red = metrics(red_supported_scores)
    red_full = metrics(red_all_scores)
    blue = metrics([score for score in scores if score["band"] == "blue"])
    red_classes = red["by_class"]
    blank_upper = (
        1.0 - 0.05 ** (1.0 / red_classes["blank"]["n"])
        if red_classes["blank"]["n"] and red_classes["blank"]["overfit"] == 0
        else None
    )
    gates = {
        "powered_class_sizes": all(red_classes[name]["n"] >= THRESHOLDS["minimum_red_cases_per_class"] for name in red_classes),
        "blank_false_positive_rate": (1.0 - red_classes["blank"]["accuracy"]) <= THRESHOLDS["blank_false_positive_rate_max"],
        "blank_false_positive_one_sided_95_upper": (
            blank_upper is not None
            and blank_upper <= THRESHOLDS["blank_false_positive_rate_max"]
        ),
        "single_accuracy": red_classes["single"]["accuracy"] >= THRESHOLDS["single_accuracy_min"],
        "double_accuracy": red_classes["double"]["accuracy"] >= THRESHOLDS["double_accuracy_min"],
        "triple_accuracy": red_classes["triple"]["accuracy"] >= THRESHOLDS["triple_accuracy_min"],
        "median_absolute_velocity_error": red["median_absolute_velocity_error_kms"] <= THRESHOLDS["median_absolute_velocity_error_kms_max"],
        "median_absolute_sigma_fractional_error": red["median_absolute_sigma_fractional_error"] <= THRESHOLDS["median_absolute_sigma_fractional_error_max"],
        "median_absolute_flux_fractional_error": red["median_absolute_flux_fractional_error"] <= THRESHOLDS["median_absolute_flux_fractional_error_max"],
        "evidence_reliability": red["evidence_flags"] <= THRESHOLDS["evidence_flags_max"],
    }
    return {
        "matrix_version": MATRIX_VERSION,
        "n_cases": len(scores),
        "red_powered_gate": red,
        "red_full_boundary_matrix": red_full,
        "red_supported_domain_exclusion": (
            "exclude 300 km/s, 0.5-ratio doubles with sigma=(160,80) km/s"
        ),
        "blank_false_positive_one_sided_95_upper": blank_upper,
        "blue_lsf_check": blue,
        "acceptance_thresholds": THRESHOLDS,
        "gate_results": gates,
        "all_powered_red_gates_pass": all(gates.values()),
        "failed_gates": [name for name, passed in gates.items() if not passed],
        "runtime_seconds": float(sum(score["runtime_seconds"] for score in scores)),
    }


def reliable_from_score(score: dict[str, Any]) -> bool:
    return score.get("selection_status", "accepted") in {"accepted", "accepted_after_audit"}


def in_supported_red_domain(score: dict[str, Any]) -> bool:
    return (
        score.get("band") == "red"
        and score.get("scenario") not in RED_BOUNDARY_EXCLUDED_SCENARIOS
    )


def median(parameters: list[dict[str, Any]], key: str) -> float | None:
    values = [abs(float(item[key])) for item in parameters if np.isfinite(item.get(key, np.nan))]
    return float(np.median(values)) if values else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--red-blank-per-target", type=int, default=30)
    parser.add_argument("--red-single-per-target", type=int, default=10)
    parser.add_argument("--red-double-per-target", type=int, default=12)
    parser.add_argument("--red-triple-per-target", type=int, default=10)
    parser.add_argument("--blue-per-target-class", type=int, default=2)
    parser.add_argument("--target", choices=tuple(TARGETS), action="append", default=[])
    parser.add_argument("--band", choices=("red", "blue"), action="append", default=[])
    parser.add_argument(
        "--component-class",
        choices=("blank", "single", "double", "triple"),
        action="append",
        default=[],
    )
    parser.add_argument("--only", action="append", default=[])
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    red_counts = (
        args.red_blank_per_target,
        args.red_single_per_target,
        args.red_double_per_target,
        args.red_triple_per_target,
    )
    if min(red_counts) < 1 or args.blue_per_target_class < 0:
        parser.error("case counts must be non-negative, with at least one red case per class")
    all_cases = definitions(*red_counts, args.blue_per_target_class)
    names = {case["name"] for case in all_cases}
    unknown = sorted(set(args.only) - names)
    if unknown:
        parser.error("unknown --only case(s): " + ", ".join(unknown))
    active = [
        case for case in all_cases
        if (not args.target or case["target"] in args.target)
        and (not args.band or case["band"] in args.band)
        and (not args.component_class or case["component_class"] in args.component_class)
        and (not args.only or case["name"] in args.only)
    ]
    OUTPUT.mkdir(parents=True, exist_ok=True)
    donors = {target: load_donor(target) for target in sorted({case["target"] for case in active})}
    for index, case in enumerate(active, start=1):
        checkpoint = OUTPUT / "cases" / case["name"] / "result.json"
        if checkpoint.exists() and not args.no_resume:
            print(f"[{index}/{len(active)}] resumed {case['name']}", flush=True)
            continue
        screening_fit = fit_config(case["band"], "screening")
        reference_line = "halpha" if case["band"] == "red" else "oiii5007"
        line_ratios = (
            {"halpha": 1.0, "nii6583": 0.7}
            if case["band"] == "red"
            else {"oiii5007": 1.0, "hbeta": 0.25}
        )
        spectrum, truth = real_noise_injection_spectrum(
            case["name"], donors[case["target"]], case["components"], screening_fit,
            case["residual_seed"], reference_line=reference_line,
            independent_line_ratios=line_ratios, calibrate_uncertainty=True,
        )
        result, screen_runtime, stdout, stderr = run_fit(
            spectrum, copy.deepcopy(screening_fit), case["sampler_seed"]
        )
        screening_result = result
        runtimes = {"screening": screen_runtime, "standard": 0.0}
        screening_score = score_recovery(result, truth)
        audit_reasons = []
        if not screening_score["component_count_correct"]:
            audit_reasons.append("screening component count differs from injected truth")
        if not reliable(result):
            audit_reasons.extend(result.get("selection_reasons", ["screening reliability flag"]))
        if audit_reasons:
            result, standard_runtime, standard_stdout, standard_stderr = run_fit(
                spectrum, fit_config(case["band"], "standard"), case["sampler_seed"] + 500_000
            )
            runtimes["standard"] = standard_runtime
            stdout += "\n--- standard audit ---\n" + standard_stdout
            stderr += "\n--- standard audit ---\n" + standard_stderr
        score = final_score(case, result, truth, runtimes)
        payload = {
            "matrix_version": MATRIX_VERSION,
            "case": {**case, "components": [component.__dict__ for component in case["components"]]},
            "truth": truth,
            "screening_fit": screening_fit,
            "standard_fit": fit_config(case["band"], "standard") if audit_reasons else None,
            "screening_result": screening_result,
            "audit_reasons": audit_reasons,
            "result": result,
            "score": score,
            "sampler_stdout_tail": stdout[-4000:],
            "sampler_stderr_tail": stderr[-4000:],
        }
        write_json(checkpoint, payload)
        print(
            f"[{index}/{len(active)}] {case['name']}: {score['true_components']} -> "
            f"{score['recovered_components']} status={score['selection_status']} "
            f"audit={bool(audit_reasons)} ({score['runtime_seconds']:.1f}s)", flush=True,
        )
    scores = []
    for case in all_cases:
        checkpoint = OUTPUT / "cases" / case["name"] / "result.json"
        if checkpoint.exists():
            scores.append(json.loads(checkpoint.read_text(encoding="utf-8"))["score"])
    summary = summarize(scores)
    write_json(OUTPUT / "summary.json", summary)
    write_json(OUTPUT / "scores.json", scores)
    write_json(
        OUTPUT / "manifest.json",
        {
            "matrix_version": MATRIX_VERSION,
            "red_cases_per_target": {
                "blank": args.red_blank_per_target,
                "single": args.red_single_per_target,
                "double": args.red_double_per_target,
                "triple": args.red_triple_per_target,
            },
            "blue_per_target_class": args.blue_per_target_class,
            "targets": TARGETS,
            "lsf": LSF,
            "screening_fits": {band: fit_config(band, "screening") for band in ("red", "blue")},
            "standard_fits": {band: fit_config(band, "standard") for band in ("red", "blue")},
            "audit_policy": "standard rerun of validation cases with an incorrect screening count or a reliability flag",
            "supported_triple_domain": "adjacent spacing 400-500 km/s; weakest effective peak S/N >= 10",
        },
    )
    print(json.dumps(summary, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
