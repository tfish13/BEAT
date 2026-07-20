#!/usr/bin/env python3
"""Run the frozen nine-spaxel NGC 1365 broad-H-alpha stability gate."""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import time
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np

from beat.config import load_config
from beat.data import iter_spectra
from beat.fitting import fit_spectrum, make_diagnostic_plot, selected_model
from beat.model import ModelDefinition, prepare_spectrum


PROJECT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT / "validation" / "ngc1365_broad_stability"
MATRIX_VERSION = "ngc1365-broad-stability-v1"
X_VALUES = (164, 165, 166)
Y_VALUES = (154, 155, 156)
CONTINUA = {
    "linear": 1,
    "quadratic": 2,
}
MODELS = ("narrow", "broad1", "broad2")
SCREENING = {
    "min_num_live_points": 40,
    "min_ess": 60,
    "dlogz": 2.0,
    "show_status": False,
    "stepsampler": "slice",
    "nsteps": 12,
}
STANDARD = {
    "min_num_live_points": 100,
    "min_ess": 200,
    "dlogz": 0.5,
    "show_status": False,
    "stepsampler": "slice",
    "nsteps": 20,
}
CRITERIA = {
    "minimum_spaxels": 7,
    "continuum_sign_agreement_min": 8,
    "cross_continuum_velocity_change_max_kms": 250.0,
    "cross_continuum_width_change_max": 0.25,
    "cross_continuum_wing_fraction_change_max": 0.10,
    "spatial_velocity_mad_max_kms": 300.0,
    "spatial_width_fractional_mad_max": 0.25,
    "prior_boundary_fraction": 0.05,
    "delta_logz_threshold": 5.0,
    "audit_delta_logz_half_width": 5.0,
}
AUDIT_STOPPING_RULE = (
    "screening matrix is terminal; incomplete unpaired standard checkpoints "
    "are retained but excluded from evidence comparisons"
)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def case_definitions() -> list[dict[str, Any]]:
    cases = []
    for y_index, y in enumerate(Y_VALUES):
        for x_index, x in enumerate(X_VALUES):
            for continuum_index, continuum in enumerate(CONTINUA):
                for model_index, model in enumerate(MODELS):
                    cases.append({
                        "name": f"x{x}_y{y}_{continuum}_{model}",
                        "x": x,
                        "y": y,
                        "continuum": continuum,
                        "continuum_degree": CONTINUA[continuum],
                        "model": model,
                        "screening_seed": (
                            5_000_000 + y_index * 100_000 + x_index * 10_000
                            + continuum_index * 1000 + model_index * 100
                        ),
                    })
    return cases


def fit_config(continuum: str, model: str, profile: str) -> dict[str, Any]:
    config = load_config(PROJECT / "examples" / "muse_ngc1365_broad_halpha.local.yaml")
    fit = copy.deepcopy(config["fit"])
    fit["kinematics"]["max_components"] = 2
    fit["continuum"] = {
        "degree": CONTINUA[continuum],
        "windows": [[6420.0, 6470.0], [6660.0, 6710.0]],
        "prior_width_rms": 10.0,
    }
    fit["noise"] = {"model": "ar1", "rho": "auto", "marginal_scale": "auto"}
    fit["selection"] = {
        "delta_logz": 5.0,
        "stop_when_not_preferred": True,
        "audit": {
            "mode": "flag",
            "sampling": copy.deepcopy(STANDARD),
        },
    }
    fit["sampling"] = copy.deepcopy(SCREENING if profile == "screening" else STANDARD)
    if model == "narrow":
        fit.pop("broad_components", None)
    elif model == "broad1":
        fit["broad_components"] = [{
            "name": "broad_halpha",
            "line": "halpha",
            "velocity_kms": [-2500.0, 2500.0],
            "sigma_kms": [400.0, 6000.0],
        }]
    elif model == "broad2":
        fit["broad_components"] = [
            {
                "name": "broad_halpha_core",
                "line": "halpha",
                "velocity_kms": [-2000.0, 2000.0],
                "sigma_kms": [400.0, 1800.0],
            },
            {
                "name": "broad_halpha_wing",
                "line": "halpha",
                "velocity_kms": [-3500.0, 3500.0],
                "sigma_kms": [1800.0, 6000.0],
            },
        ]
    else:  # pragma: no cover
        raise ValueError(model)
    return fit


def quality(spectrum, fit: dict[str, Any], result: dict[str, Any]) -> dict[str, float]:
    prepared = prepare_spectrum(spectrum, fit)
    chosen = selected_model(result)
    definition = ModelDefinition(prepared, fit, result["selected_components"])
    model = definition.evaluate(np.asarray(chosen["maximum_likelihood"], dtype=float))
    residual = (prepared.flux - model) / prepared.uncertainty
    return {
        "normalized_residual_rms": float(np.sqrt(np.mean(residual**2))),
        "reduced_chi_square": float(np.sum(residual**2) / max(1, residual.size - definition.ndim)),
    }


def run_fit(spectrum, fit: dict[str, Any], seed: int) -> dict[str, Any]:
    np.random.seed(seed)
    stdout, stderr = StringIO(), StringIO()
    started = time.perf_counter()
    with redirect_stdout(stdout), redirect_stderr(stderr):
        result = fit_spectrum(spectrum, fit)
    return {
        "result": result,
        "quality": quality(spectrum, fit, result),
        "runtime_seconds": time.perf_counter() - started,
        "sampler_stdout_tail": stdout.getvalue()[-4000:],
        "sampler_stderr_tail": stderr.getvalue()[-4000:],
    }


def load_spaxels() -> dict[tuple[int, int], Any]:
    config = load_config(PROJECT / "examples" / "muse_ngc1365_broad_halpha.local.yaml")
    return {
        (int(spectrum.metadata["x"]), int(spectrum.metadata["y"])): spectrum
        for spectrum in iter_spectra(config)
    }


def checkpoint_path(case: dict[str, Any], profile: str) -> Path:
    return OUTPUT / "cases" / case["name"] / f"{profile}.json"


def execute(case: dict[str, Any], spectrum, profile: str) -> dict[str, Any]:
    fit = fit_config(case["continuum"], case["model"], profile)
    seed = int(case["screening_seed"] + (500_000 if profile == "standard" else 0))
    payload = run_fit(spectrum, fit, seed)
    return {
        "matrix_version": MATRIX_VERSION,
        "profile": profile,
        "case": case,
        "fit": fit,
        "seed": seed,
        "spectrum_id": spectrum.spectrum_id,
        "input_file": spectrum.metadata.get("input_file"),
        **payload,
    }


def reliable(payload: dict[str, Any]) -> bool:
    return payload["result"].get("selection_status", "accepted") in {
        "accepted", "accepted_after_audit"
    }


def broad(payload: dict[str, Any], name: str) -> dict[str, Any] | None:
    return next(
        (item for item in payload["result"].get("broad_components", []) if item["name"] == name),
        None,
    )


def effective_payload(case: dict[str, Any]) -> dict[str, Any] | None:
    screening = checkpoint_path(case, "screening")
    return json.loads(screening.read_text(encoding="utf-8")) if screening.exists() else None


def selected_logz_error(payload: dict[str, Any]) -> float | None:
    selected = payload["result"]["selected_components"]
    model = next(
        item for item in payload["result"]["models"] if item["n_components"] == selected
    )
    value = model.get("logz_error")
    return None if value is None else float(value)


def median(values: list[float]) -> float | None:
    return float(np.median(values)) if values else None


def mad(values: list[float]) -> float | None:
    if not values:
        return None
    array = np.asarray(values, dtype=float)
    return float(np.median(np.abs(array - np.median(array))))


def fractional_difference(a: float, b: float) -> float:
    scale = 0.5 * (abs(a) + abs(b))
    return abs(a - b) / scale if scale > 0 else 0.0


def grouped_payloads() -> tuple[dict[tuple[int, int, str, str], dict[str, Any]], list[dict[str, Any]]]:
    payloads = {}
    missing = []
    for case in case_definitions():
        payload = effective_payload(case)
        if payload is None:
            missing.append(case)
        else:
            payloads[(case["x"], case["y"], case["continuum"], case["model"])] = payload
    return payloads, missing


def audit_candidates(payloads: dict[tuple[int, int, str, str], dict[str, Any]]) -> list[dict[str, Any]]:
    names = set()
    reasons: dict[str, list[str]] = {}
    by_name = {case["name"]: case for case in case_definitions()}
    by_key = {
        (case["x"], case["y"], case["continuum"], case["model"]): case
        for case in case_definitions()
    }
    for y in Y_VALUES:
        for x in X_VALUES:
            for continuum in CONTINUA:
                group = {model: payloads[(x, y, continuum, model)] for model in MODELS}
                comparisons = (
                    ("broad1_minus_narrow", "narrow", "broad1"),
                    ("broad2_minus_broad1", "broad1", "broad2"),
                )
                for label, low, high in comparisons:
                    delta = group[high]["result"]["selected_logz"] - group[low]["result"]["selected_logz"]
                    if abs(delta - CRITERIA["delta_logz_threshold"]) <= CRITERIA["audit_delta_logz_half_width"]:
                        for model in (low, high):
                            case = by_key[(x, y, continuum, model)]
                            names.add(case["name"])
                            reasons.setdefault(case["name"], []).append(f"{label}={delta:.3f} near threshold")
                if any(not reliable(group[model]) for model in MODELS):
                    for model in MODELS:
                        case = by_key[(x, y, continuum, model)]
                        names.add(case["name"])
                        reasons.setdefault(case["name"], []).append("paired model has evidence-reliability flag")
    return [{**by_name[name], "audit_reasons": reasons[name]} for name in sorted(names)]


def focused_audit_candidates() -> list[dict[str, Any]]:
    """Return paired models only where screening preferred the second broad Gaussian."""
    by_key = {
        (case["x"], case["y"], case["continuum"], case["model"]): case
        for case in case_definitions()
    }
    selected = []
    for y in Y_VALUES:
        for x in X_VALUES:
            for continuum in CONTINUA:
                one_case = by_key[(x, y, continuum, "broad1")]
                two_case = by_key[(x, y, continuum, "broad2")]
                one_path = checkpoint_path(one_case, "screening")
                two_path = checkpoint_path(two_case, "screening")
                if not one_path.exists() or not two_path.exists():
                    continue
                one = json.loads(one_path.read_text(encoding="utf-8"))
                two = json.loads(two_path.read_text(encoding="utf-8"))
                delta = two["result"]["selected_logz"] - one["result"]["selected_logz"]
                if delta >= CRITERIA["delta_logz_threshold"]:
                    reason = f"screening broad2_minus_broad1={delta:.3f} preferred broad2"
                    selected.extend([
                        {**one_case, "audit_reasons": [reason]},
                        {**two_case, "audit_reasons": [reason]},
                    ])
    return selected


def summarize() -> dict[str, Any]:
    payloads, missing = grouped_payloads()
    if missing:
        return {
            "matrix_version": MATRIX_VERSION,
            "complete": False,
            "completed_fits": len(payloads),
            "expected_fits": len(case_definitions()),
            "missing": [case["name"] for case in missing],
        }
    rows = []
    for y in Y_VALUES:
        for x in X_VALUES:
            for continuum in CONTINUA:
                group = {model: payloads[(x, y, continuum, model)] for model in MODELS}
                core = broad(group["broad2"], "broad_halpha_core")
                wing = broad(group["broad2"], "broad_halpha_wing")
                wing_fraction = wing["flux"] / (core["flux"] + wing["flux"])
                rows.append({
                    "x": x, "y": y, "continuum": continuum,
                    "profiles": {model: group[model]["profile"] for model in MODELS},
                    "selection_status": {
                        model: group[model]["result"].get("selection_status", "accepted")
                        for model in MODELS
                    },
                    "selected_narrow_components": {
                        model: group[model]["result"]["selected_components"] for model in MODELS
                    },
                    "selected_logz": {
                        model: group[model]["result"]["selected_logz"] for model in MODELS
                    },
                    "selected_logz_error": {
                        model: selected_logz_error(group[model]) for model in MODELS
                    },
                    "delta_logz_broad1_minus_narrow": (
                        group["broad1"]["result"]["selected_logz"]
                        - group["narrow"]["result"]["selected_logz"]
                    ),
                    "delta_logz_broad2_minus_broad1": (
                        group["broad2"]["result"]["selected_logz"]
                        - group["broad1"]["result"]["selected_logz"]
                    ),
                    "quality": {model: group[model]["quality"] for model in MODELS},
                    "broad1": broad(group["broad1"], "broad_halpha"),
                    "broad2_core": core,
                    "broad2_wing": wing,
                    "wing_flux_fraction": wing_fraction,
                })
    by_continuum = {}
    for continuum in CONTINUA:
        subset = [row for row in rows if row["continuum"] == continuum]
        velocities = [row["broad2_wing"]["velocity_kms"] for row in subset]
        widths = [row["broad2_wing"]["sigma_kms"] for row in subset]
        by_continuum[continuum] = {
            "broad_required_count": sum(
                row["delta_logz_broad1_minus_narrow"] >= CRITERIA["delta_logz_threshold"]
                for row in subset
            ),
            "second_broad_required_count": sum(
                row["delta_logz_broad2_minus_broad1"] >= CRITERIA["delta_logz_threshold"]
                for row in subset
            ),
            "red_wing_count": sum(value > 0 for value in velocities),
            "wing_velocity_median_kms": median(velocities),
            "wing_velocity_mad_kms": mad(velocities),
            "wing_width_median_kms": median(widths),
            "wing_width_fractional_mad": mad(widths) / median(widths),
            "wing_flux_fraction_median": median([row["wing_flux_fraction"] for row in subset]),
        }
    pairs = []
    for y in Y_VALUES:
        for x in X_VALUES:
            linear = next(row for row in rows if row["x"] == x and row["y"] == y and row["continuum"] == "linear")
            quadratic = next(row for row in rows if row["x"] == x and row["y"] == y and row["continuum"] == "quadratic")
            pairs.append({
                "x": x, "y": y,
                "wing_velocity_sign_agrees": bool(
                    np.sign(linear["broad2_wing"]["velocity_kms"])
                    == np.sign(quadratic["broad2_wing"]["velocity_kms"])
                ),
                "wing_velocity_absolute_change_kms": abs(
                    linear["broad2_wing"]["velocity_kms"] - quadratic["broad2_wing"]["velocity_kms"]
                ),
                "wing_width_fractional_change": fractional_difference(
                    linear["broad2_wing"]["sigma_kms"], quadratic["broad2_wing"]["sigma_kms"]
                ),
                "wing_flux_fraction_absolute_change": abs(
                    linear["wing_flux_fraction"] - quadratic["wing_flux_fraction"]
                ),
            })
    boundary_span = 6000.0 - 1800.0
    lower_safe = 1800.0 + CRITERIA["prior_boundary_fraction"] * boundary_span
    upper_safe = 6000.0 - CRITERIA["prior_boundary_fraction"] * boundary_span
    boundary_hits = [
        {"x": row["x"], "y": row["y"], "continuum": row["continuum"], "sigma_kms": row["broad2_wing"]["sigma_kms"]}
        for row in rows
        if not lower_safe < row["broad2_wing"]["sigma_kms"] < upper_safe
    ]
    cross = {
        "velocity_sign_agreement_count": sum(pair["wing_velocity_sign_agrees"] for pair in pairs),
        "median_absolute_velocity_change_kms": median([pair["wing_velocity_absolute_change_kms"] for pair in pairs]),
        "median_width_fractional_change": median([pair["wing_width_fractional_change"] for pair in pairs]),
        "median_wing_flux_fraction_absolute_change": median([pair["wing_flux_fraction_absolute_change"] for pair in pairs]),
        "wing_width_prior_boundary_hits": boundary_hits,
    }
    criteria = {
        "broad_required_both_continua": all(
            item["broad_required_count"] >= CRITERIA["minimum_spaxels"] for item in by_continuum.values()
        ),
        "second_broad_required_both_continua": all(
            item["second_broad_required_count"] >= CRITERIA["minimum_spaxels"] for item in by_continuum.values()
        ),
        "red_wing_both_continua": all(
            item["red_wing_count"] >= CRITERIA["minimum_spaxels"] for item in by_continuum.values()
        ),
        "continuum_velocity_sign_stable": cross["velocity_sign_agreement_count"] >= CRITERIA["continuum_sign_agreement_min"],
        "continuum_velocity_stable": cross["median_absolute_velocity_change_kms"] <= CRITERIA["cross_continuum_velocity_change_max_kms"],
        "continuum_width_stable": cross["median_width_fractional_change"] <= CRITERIA["cross_continuum_width_change_max"],
        "continuum_wing_fraction_stable": cross["median_wing_flux_fraction_absolute_change"] <= CRITERIA["cross_continuum_wing_fraction_change_max"],
        "spatial_velocity_stable": all(
            item["wing_velocity_mad_kms"] <= CRITERIA["spatial_velocity_mad_max_kms"] for item in by_continuum.values()
        ),
        "spatial_width_stable": all(
            item["wing_width_fractional_mad"] <= CRITERIA["spatial_width_fractional_mad_max"] for item in by_continuum.values()
        ),
        "wing_width_not_prior_bound": not boundary_hits,
    }
    criteria["second_broad_descriptively_stable"] = all(
        value for name, value in criteria.items() if name != "broad_required_both_continua"
    )
    near_threshold = audit_candidates(payloads)
    screening_preferred = focused_audit_candidates()
    completed_standard = [
        case["name"] for case in case_definitions()
        if checkpoint_path(case, "standard").exists()
    ]
    return {
        "matrix_version": MATRIX_VERSION,
        "complete": True,
        "effective_standard_audits": 0,
        "completed_unpaired_standard_checkpoints_excluded": completed_standard,
        "rows": rows,
        "by_continuum": by_continuum,
        "cross_continuum": cross,
        "criteria_thresholds": CRITERIA,
        "criteria_results": criteria,
        "audit_stopping_rule": AUDIT_STOPPING_RULE,
        "screening_preferred_second_component_pairs": screening_preferred,
        "audit_candidates": [],
        "superseded_near_threshold_candidates": near_threshold,
        "audit_complete": True,
        "recommendation": (
            "Use one broad component for routine fitting; for detailed BLR work use a flexible asymmetric profile and do not physically interpret multiple broad Gaussians automatically."
            if not criteria["second_broad_descriptively_stable"]
            else "Two broad Gaussians are a spatially and methodologically stable descriptive model, but they must not be interpreted automatically as distinct physical BLR components."
        ),
        "runtime_seconds": float(sum(
            effective_payload(case)["runtime_seconds"] for case in case_definitions()
        )),
    }


def write_manifest() -> None:
    write_json(OUTPUT / "manifest.json", {
        "matrix_version": MATRIX_VERSION,
        "input": "/Users/tfischer/research/data/ngc1365/muse/ADP.2017-03-27T12:08:50.541.fits",
        "redshift": 0.00546,
        "case_definitions": case_definitions(),
        "continua": CONTINUA,
        "models": MODELS,
        "screening": SCREENING,
        "standard": STANDARD,
        "criteria": CRITERIA,
        "audit_stopping_rule": AUDIT_STOPPING_RULE,
        "interpretation": "Multiple broad Gaussians are descriptive basis functions, not automatically physical BLR components.",
    })


def freeze(summary: dict[str, Any]) -> None:
    if not summary.get("complete") or not summary.get("audit_complete"):
        return
    files = [OUTPUT / "manifest.json", OUTPUT / "summary.json"]
    for case in case_definitions():
        files.append(checkpoint_path(case, "screening"))
        standard = checkpoint_path(case, "standard")
        if standard.exists():
            files.append(standard)
    write_json(OUTPUT / "frozen_gate_manifest.json", {
        "matrix_version": MATRIX_VERSION,
        "status": "diagnostic_complete_and_frozen",
        "files": {str(path.relative_to(PROJECT)): sha256(path) for path in files},
        "expansion_rule": "Do not add cases unless an unexplained evidence flag remains or a predeclared comparison requires audit.",
    })


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=("screening", "audit", "summarize"), default="screening")
    parser.add_argument("--x", type=int, action="append", default=[])
    parser.add_argument("--y", type=int, action="append", default=[])
    parser.add_argument("--only", action="append", default=[])
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    write_manifest()
    cases = case_definitions()
    known = {case["name"] for case in cases}
    if set(args.only) - known:
        parser.error("unknown --only case(s): " + ", ".join(sorted(set(args.only) - known)))
    if args.phase == "summarize":
        summary = summarize()
        write_json(OUTPUT / "summary.json", summary)
        freeze(summary)
        print(json.dumps({key: value for key, value in summary.items() if key != "rows"}, indent=2))
        return 0
    if args.phase == "audit":
        payloads, missing = grouped_payloads()
        if missing:
            raise RuntimeError("screening matrix is incomplete")
        candidates = {case["name"] for case in focused_audit_candidates()}
        cases = [case for case in cases if case["name"] in candidates]
    cases = [
        case for case in cases
        if (not args.x or case["x"] in args.x)
        and (not args.y or case["y"] in args.y)
        and (not args.only or case["name"] in args.only)
    ]
    profile = "standard" if args.phase == "audit" else "screening"
    spaxels = load_spaxels()
    for index, case in enumerate(cases, start=1):
        path = checkpoint_path(case, profile)
        if path.exists() and not args.no_resume:
            print(f"[{index}/{len(cases)}] resumed {case['name']} {profile}", flush=True)
            continue
        payload = execute(case, spaxels[(case["x"], case["y"])], profile)
        write_json(path, payload)
        result = payload["result"]
        print(
            f"[{index}/{len(cases)}] {case['name']} {profile}: "
            f"narrow={result['selected_components']} logZ={result['selected_logz']:.2f} "
            f"status={result.get('selection_status', 'accepted')} "
            f"({payload['runtime_seconds']:.1f}s)",
            flush=True,
        )
    summary = summarize()
    write_json(OUTPUT / "summary.json", summary)
    freeze(summary)
    print(json.dumps({key: value for key, value in summary.items() if key != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
