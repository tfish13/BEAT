#!/usr/bin/env python3
"""Run bounded NIRSpec G235H/G395H injection/recovery validation."""

from __future__ import annotations

import argparse
import csv
import json
import time
from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
from typing import Any

import numpy as np

from beat import __version__
from beat.config import load_config
from beat.data import iter_spectra
from beat.fitting import fit_spectrum, make_diagnostic_plot
from beat.injection import (
    InjectedComponent,
    inject_emission_lines,
    real_noise_injection_spectrum,
    score_recovery,
)
from beat.lsf import FWHM_TO_SIGMA, lsf_sigma_angstrom
from beat.spectrum import Spectrum


PROJECT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT = PROJECT / "validation" / "nirspec_injection_recovery"

ACCEPTANCE_THRESHOLDS = {
    "blank_false_positive_rate_max": 0.05,
    "single_component_accuracy_min": 0.90,
    "two_component_accuracy_min": 0.80,
    "three_component_accuracy_min": 0.70,
    "median_absolute_velocity_error_kms_max": 10.0,
    "median_absolute_sigma_fractional_error_max": 0.15,
    "median_absolute_flux_fractional_error_max": 0.10,
    "ambiguous_or_unverified_selections_max": 0,
    "minimum_cases_per_count_class": 20,
}

DATASETS = {
    "ngc4151_sivi": {
        "config": "nirspec_ngc4151.local.yaml",
        "science_xy": (28, 16),
        "donor_xy": (28, 38),
        "alternate_donor_xy": (27, 38),
        "line": {"name": "sivi_1963", "wavelength": 19634.1},
        "window": [19200.0, 20050.0],
        "continuum_windows": [[19200.0, 19400.0], [19800.0, 20050.0]],
        "extra_residual_mask_lines": [19575.6],
        "donor_selection": (
            "lowest residual/uncertainty among 81 valid spectra in five bounded "
            "off-nuclear candidate boxes"
        ),
    },
    "ic5063_bralpha": {
        "config": "nirspec_ic5063.local.yaml",
        "science_xy": (44, 44),
        "donor_xy": (41, 55),
        "alternate_donor_xy": (42, 54),
        "line": {"name": "bralpha", "wavelength": 40522.7},
        "window": [39700.0, 41350.0],
        "continuum_windows": [[39700.0, 40000.0], [41000.0, 41350.0]],
        "extra_residual_mask_lines": [],
        "donor_selection": (
            "lowest residual/uncertainty among 256 valid spectra in four bounded "
            "off-nuclear candidate boxes"
        ),
    },
}


def _json_default(value: Any) -> Any:
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, Path):
        return str(value)
    raise TypeError(type(value).__name__)


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, default=_json_default) + "\n",
        encoding="utf-8",
    )


def validation_fit(
    specification: dict[str, Any], profile: str, correlated_noise: bool = False
) -> dict[str, Any]:
    if profile in {"pilot", "expanded"}:
        sampling = {"min_num_live_points": 40, "min_ess": 40, "dlogz": 5.0}
    else:
        sampling = {"min_num_live_points": 100, "min_ess": 200, "dlogz": 0.5}
    nsteps = 10 if profile in {"pilot", "expanded"} else 20
    return {
        "frame": "rest",
        "wavelength_medium": "vacuum",
        "window": specification["window"],
        "minimum_valid_pixels": 100,
        "continuum": {
            "degree": 1,
            "windows": specification["continuum_windows"],
            "prior_width_rms": 10.0,
        },
        "noise": (
            {"model": "ar1", "rho": "auto", "marginal_scale": "auto"}
            if correlated_noise
            else {"model": "independent"}
        ),
        "kinematics": {
            "max_components": 3,
            "velocity_kms": [-1000.0, 1000.0],
            "sigma_kms": [30.0, 1000.0],
        },
        "lines": [specification["line"]],
        "lsf": {"model": "instrument"},
        "flux_prior": {"min_snr": 0.05, "max_signal_factor": 40.0},
        "selection": {
            "delta_logz": 5.0,
            "stop_when_not_preferred": True,
            "audit": {"mode": "flag"},
        },
        "sampling": {
            **sampling,
            "show_status": False,
            "stepsampler": "slice",
            "nsteps": nsteps,
        },
        "injection_residual_mask_lines": specification[
            "extra_residual_mask_lines"
        ],
        "injection_residual_mask_half_width_angstrom": 30.0,
    }


def load_coordinate(config: dict[str, Any], xy: tuple[int, int]) -> Spectrum:
    x, y = xy
    config = json.loads(json.dumps(config, default=_json_default))
    config["input"]["x_range"] = [x, x + 1]
    config["input"]["y_range"] = [y, y + 1]
    return next(iter_spectra(config))


def controlled_spectrum(
    name: str,
    source: Spectrum,
    fit: dict[str, Any],
    components: list[InjectedComponent],
    seed: int,
) -> tuple[Spectrum, dict[str, Any]]:
    prepared = source.prepared()
    lo, hi = np.asarray(fit["window"], dtype=float) * (1.0 + source.redshift)
    selected = (prepared.wavelength >= lo) & (prepared.wavelength <= hi)
    wavelength = prepared.wavelength[selected]
    source_uncertainty = prepared.uncertainty[selected]
    if source_uncertainty is None or wavelength.size < fit["minimum_valid_pixels"]:
        raise ValueError(f"{name}: insufficient source sampling for injection")
    noise_sigma = float(np.median(source_uncertainty))
    uncertainty = np.full(wavelength.shape, noise_sigma)
    rng = np.random.default_rng(seed)
    midpoint = float(np.mean(wavelength))
    baseline = 5.0 * noise_sigma + 1.0e-4 * noise_sigma * (wavelength - midpoint)
    baseline += rng.normal(0.0, noise_sigma, wavelength.size)
    metadata = dict(source.metadata)
    metadata.update({"noise_source": "analytic_gaussian", "injection_seed": seed})
    flux, component_truth = inject_emission_lines(
        wavelength,
        baseline,
        uncertainty,
        source.redshift,
        components,
        fit,
        metadata,
        reference_line=fit["lines"][0]["name"],
    )
    truth = {
        "spectrum_id": name,
        "stage": "controlled",
        "seed": seed,
        "redshift": source.redshift,
        "true_components": len(component_truth),
        "components": component_truth,
        "broad_components": [],
        "reference_line": fit["lines"][0]["name"],
        "sampling_source_spectrum_id": source.spectrum_id,
        "constant_noise_sigma": noise_sigma,
    }
    return Spectrum(name, wavelength, flux, uncertainty, source.redshift, metadata=metadata), truth


def case_definitions(
    label: str, profile: str, matrix: str = "core"
) -> list[dict[str, Any]]:
    base = 11000 if label.startswith("ngc4151") else 21000
    if matrix == "donor_replicates":
        cases = []
        for donor_index, donor in enumerate(("primary", "alternate")):
            for replicate in range(1, 6):
                cases.append(
                    {
                        "name": f"{donor}_real_blank_r{replicate:02d}",
                        "stage": "real_noise",
                        "seed": base + 7000 + 100 * donor_index + replicate,
                        "components": [],
                        "donor": donor,
                    }
                )
                cases.append(
                    {
                        "name": f"{donor}_real_single_snr10_r{replicate:02d}",
                        "stage": "real_noise",
                        "seed": base + 7050 + 100 * donor_index + replicate,
                        "components": [InjectedComponent(0.0, 80.0, 10.0)],
                        "donor": donor,
                    }
                )
        for case in cases:
            case["name"] = f"{label}_{case['name']}"
        return cases
    if matrix == "boundary":
        cases = [
            {"name": "controlled_single_snr5_sigma80", "stage": "controlled", "seed": base + 5005,
             "components": [InjectedComponent(0.0, 80.0, 5.0)]},
            {"name": "controlled_single_snr10_sigma250", "stage": "controlled", "seed": base + 5250,
             "components": [InjectedComponent(0.0, 250.0, 10.0)]},
            {"name": "controlled_double_sep150_equal", "stage": "controlled", "seed": base + 5150,
             "components": [InjectedComponent(-75.0, 80.0, 15.0), InjectedComponent(75.0, 80.0, 15.0)]},
            {"name": "controlled_double_sep300_ratio025", "stage": "controlled", "seed": base + 5300,
             "components": [InjectedComponent(-150.0, 80.0, 15.0), InjectedComponent(150.0, 80.0, 15.0, flux_scale=0.25)]},
            {"name": "controlled_triple_adj300_equal", "stage": "controlled", "seed": base + 5600,
             "components": [InjectedComponent(-300.0, 80.0, 15.0), InjectedComponent(0.0, 80.0, 15.0), InjectedComponent(300.0, 80.0, 15.0)]},
            {"name": "controlled_triple_adj400_weak", "stage": "controlled", "seed": base + 5800,
             "components": [InjectedComponent(-400.0, 80.0, 20.0), InjectedComponent(0.0, 80.0, 20.0, flux_scale=0.5), InjectedComponent(400.0, 80.0, 20.0, flux_scale=0.25)]},
            {"name": "alternate_real_blank", "stage": "real_noise", "seed": base + 6001,
             "components": [], "donor": "alternate"},
            {"name": "alternate_real_single_snr10", "stage": "real_noise", "seed": base + 6010,
             "components": [InjectedComponent(0.0, 80.0, 10.0)], "donor": "alternate"},
            {"name": "alternate_real_double_sep300_ratio025", "stage": "real_noise", "seed": base + 6300,
             "components": [InjectedComponent(-150.0, 80.0, 15.0), InjectedComponent(150.0, 80.0, 15.0, flux_scale=0.25)], "donor": "alternate"},
            {"name": "alternate_real_triple_adj400_weak", "stage": "real_noise", "seed": base + 6800,
             "components": [InjectedComponent(-400.0, 80.0, 20.0), InjectedComponent(0.0, 80.0, 20.0, flux_scale=0.5), InjectedComponent(400.0, 80.0, 20.0, flux_scale=0.25)], "donor": "alternate"},
        ]
        for case in cases:
            case["name"] = f"{label}_{case['name']}"
        return cases
    cases = [
        {"name": "controlled_blank", "stage": "controlled", "seed": base + 1, "components": []},
        {"name": "controlled_single_snr10", "stage": "controlled", "seed": base + 10,
         "components": [InjectedComponent(0.0, 80.0, 10.0)]},
        {"name": "controlled_double_sep300", "stage": "controlled", "seed": base + 300,
         "components": [InjectedComponent(-150.0, 80.0, 15.0), InjectedComponent(150.0, 80.0, 15.0, flux_scale=0.5)]},
        {"name": "controlled_triple_sep400", "stage": "controlled", "seed": base + 400,
         "components": [InjectedComponent(-400.0, 80.0, 20.0), InjectedComponent(0.0, 80.0, 14.0), InjectedComponent(400.0, 80.0, 10.0)]},
        {"name": "real_blank", "stage": "real_noise", "seed": base + 1001, "components": []},
        {"name": "real_single_snr10", "stage": "real_noise", "seed": base + 1010,
         "components": [InjectedComponent(0.0, 80.0, 10.0)]},
        {"name": "real_double_sep300", "stage": "real_noise", "seed": base + 1300,
         "components": [InjectedComponent(-150.0, 80.0, 15.0), InjectedComponent(150.0, 80.0, 15.0, flux_scale=0.5)]},
        {"name": "real_triple_sep400", "stage": "real_noise", "seed": base + 1400,
         "components": [InjectedComponent(-400.0, 80.0, 20.0), InjectedComponent(0.0, 80.0, 14.0), InjectedComponent(400.0, 80.0, 10.0)]},
    ]
    if profile == "pilot":
        keep = {"controlled_blank", "controlled_single_snr10", "controlled_double_sep300", "real_single_snr10"}
        cases = [case for case in cases if case["name"] in keep]
    for case in cases:
        case["name"] = f"{label}_{case['name']}"
    return cases


def summarize(scores: list[dict[str, Any]]) -> dict[str, Any]:
    parameters = [item for score in scores for item in score.get("parameters", [])]
    grouped = {
        count: [score for score in scores if score["true_components"] == count]
        for count in range(4)
    }
    by_truth = {
        str(count): {
            "n": len(grouped[count]),
            "correct": sum(score["component_count_correct"] for score in grouped[count]),
        }
        for count in range(4)
    }
    median_velocity = (
        float(np.median([abs(item["velocity_error_kms"]) for item in parameters]))
        if parameters else None
    )
    median_sigma = (
        float(np.median([abs(item["sigma_fractional_error"]) for item in parameters]))
        if parameters else None
    )
    median_flux = (
        float(np.median([abs(item["flux_fractional_error"]) for item in parameters]))
        if parameters else None
    )
    rates = {
        "blank_false_positive_rate": (
            float(np.mean([score["recovered_components"] > 0 for score in grouped[0]]))
            if grouped[0] else None
        ),
        "single_component_accuracy": (
            float(np.mean([score["component_count_correct"] for score in grouped[1]]))
            if grouped[1] else None
        ),
        "two_component_accuracy": (
            float(np.mean([score["component_count_correct"] for score in grouped[2]]))
            if grouped[2] else None
        ),
        "three_component_accuracy": (
            float(np.mean([score["component_count_correct"] for score in grouped[3]]))
            if grouped[3] else None
        ),
    }
    unreliable = sum(
        score.get("selection_status") != "accepted" for score in scores
    )
    gates = {
        "blank_false_positive_rate": None
        if rates["blank_false_positive_rate"] is None
        else rates["blank_false_positive_rate"]
        <= ACCEPTANCE_THRESHOLDS["blank_false_positive_rate_max"],
        "single_component_accuracy": None
        if rates["single_component_accuracy"] is None
        else rates["single_component_accuracy"]
        >= ACCEPTANCE_THRESHOLDS["single_component_accuracy_min"],
        "two_component_accuracy": None
        if rates["two_component_accuracy"] is None
        else rates["two_component_accuracy"]
        >= ACCEPTANCE_THRESHOLDS["two_component_accuracy_min"],
        "three_component_accuracy": None
        if rates["three_component_accuracy"] is None
        else rates["three_component_accuracy"]
        >= ACCEPTANCE_THRESHOLDS["three_component_accuracy_min"],
        "median_absolute_velocity_error": median_velocity is not None
        and median_velocity
        <= ACCEPTANCE_THRESHOLDS["median_absolute_velocity_error_kms_max"],
        "median_absolute_sigma_fractional_error": median_sigma is not None
        and median_sigma
        <= ACCEPTANCE_THRESHOLDS["median_absolute_sigma_fractional_error_max"],
        "median_absolute_flux_fractional_error": median_flux is not None
        and median_flux
        <= ACCEPTANCE_THRESHOLDS["median_absolute_flux_fractional_error_max"],
        "selection_reliability": unreliable
        <= ACCEPTANCE_THRESHOLDS["ambiguous_or_unverified_selections_max"],
    }
    populated_counts = [count for count in range(4) if grouped[count]]
    powered = bool(populated_counts) and all(
        len(grouped[count])
        >= ACCEPTANCE_THRESHOLDS["minimum_cases_per_count_class"]
        for count in populated_counts
    )
    full_component_count_coverage = len(populated_counts) == 4
    all_gates_pass = all(value for value in gates.values() if value is not None)
    if all_gates_pass and powered and full_component_count_coverage:
        gate_status = "pass"
    elif all_gates_pass and powered:
        gate_status = "pass_scope"
    elif all_gates_pass:
        gate_status = "pass_underpowered"
    else:
        gate_status = "fail"
    return {
        "n_cases": len(scores),
        "component_count_correct": sum(score["component_count_correct"] for score in scores),
        "by_true_component_count": by_truth,
        "median_absolute_velocity_error_kms": median_velocity,
        "median_absolute_sigma_fractional_error": median_sigma,
        "median_absolute_flux_fractional_error": median_flux,
        "rates": rates,
        "acceptance_thresholds": ACCEPTANCE_THRESHOLDS,
        "provisional_gates": gates,
        "provisional_gate_status": gate_status,
        "statistically_powered": powered,
        "full_component_count_coverage": full_component_count_coverage,
        "total_runtime_seconds": float(sum(score["runtime_seconds"] for score in scores)),
        "selection_statuses": {
            status: sum(score.get("selection_status") == status for score in scores)
            for status in sorted({score.get("selection_status") for score in scores})
        },
        "failures": [score["spectrum_id"] for score in scores if not score["component_count_correct"]],
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--profile", choices=["pilot", "expanded", "standard"], default="pilot"
    )
    parser.add_argument(
        "--matrix",
        choices=["core", "boundary", "donor_replicates"],
        default="core",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument(
        "--empirical-uncertainty-calibration",
        action="store_true",
        help=(
            "Inflate underestimated formal uncertainties to the robust scatter "
            "of each block-resampled real-residual realization"
        ),
    )
    parser.add_argument(
        "--correlated-noise",
        action="store_true",
        help="Use an AR(1) likelihood with rho estimated from continuum windows",
    )
    parser.add_argument(
        "--only", action="append", default=[], metavar="CASE_NAME",
        help="Run only an exact case name; repeat for multiple cases",
    )
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    scores: list[dict[str, Any]] = []
    datasets_manifest: dict[str, Any] = {}
    requested = set(args.only)
    cases_by_dataset = {
        label: [
            case for case in case_definitions(label, args.profile, args.matrix)
            if not requested or case["name"] in requested
        ]
        for label in DATASETS
    }
    found = {case["name"] for cases in cases_by_dataset.values() for case in cases}
    if requested - found:
        parser.error("unknown --only case(s): " + ", ".join(sorted(requested - found)))
    all_cases = sum(len(cases) for cases in cases_by_dataset.values())
    completed = 0
    for label, specification in DATASETS.items():
        config = load_config(PROJECT / "examples" / specification["config"])
        science = load_coordinate(config, specification["science_xy"])
        donor = load_coordinate(config, specification["donor_xy"])
        alternate_donor = load_coordinate(
            config, specification["alternate_donor_xy"]
        )
        donors = {"primary": donor, "alternate": alternate_donor}
        fit = validation_fit(specification, args.profile, args.correlated_noise)
        line = specification["line"]
        observed = float(line["wavelength"]) * (1.0 + science.redshift)
        sigma = float(lsf_sigma_angstrom(observed, fit["lsf"], science.metadata))
        datasets_manifest[label] = {
            "config": specification["config"],
            "input_file": science.metadata.get("input_file"),
            "science_spectrum_id": science.spectrum_id,
            "donor_spectrum_id": donor.spectrum_id,
            "alternate_donor_spectrum_id": alternate_donor.spectrum_id,
            "donor_selection": specification["donor_selection"],
            "redshift": science.redshift,
            "grating": science.metadata.get("grating"),
            "filter": science.metadata.get("filter"),
            "pipeline_version": science.metadata.get("pipeline_version"),
            "crds_context": science.metadata.get("crds_context"),
            "instrument_lsf": science.metadata.get("instrument_lsf"),
            "line": line,
            "line_lsf_sigma_angstrom": sigma,
            "line_resolving_power": observed / (sigma / FWHM_TO_SIGMA),
            "fit": fit,
        }
        for case in cases_by_dataset[label]:
            completed += 1
            case_dir = output / "cases" / case["name"]
            checkpoint = case_dir / "recovery.json"
            if not args.no_resume and checkpoint.exists():
                payload = json.loads(checkpoint.read_text(encoding="utf-8"))
                scores.append(payload["score"])
                print(f"[{completed}/{all_cases}] resumed {case['name']}", flush=True)
                continue
            if case["stage"] == "controlled":
                spectrum, truth = controlled_spectrum(
                    case["name"], science, fit, case["components"], case["seed"]
                )
            else:
                case_donor = donors[case.get("donor", "primary")]
                spectrum, truth = real_noise_injection_spectrum(
                    case["name"], case_donor, case["components"], fit, case["seed"],
                    reference_line=line["name"], max_residual_to_uncertainty=3.0,
                    calibrate_uncertainty=args.empirical_uncertainty_calibration,
                )
            case_dir.mkdir(parents=True, exist_ok=True)
            np.savetxt(
                case_dir / "spectrum.csv",
                np.column_stack([spectrum.wavelength, spectrum.flux, spectrum.uncertainty]),
                delimiter=",", header="wavelength_angstrom,flux_density,uncertainty_1sigma", comments="",
            )
            write_json(case_dir / "truth.json", truth)
            np.random.seed(case["seed"] + 900_000)
            stdout, stderr = StringIO(), StringIO()
            start = time.perf_counter()
            with redirect_stdout(stdout), redirect_stderr(stderr):
                result = fit_spectrum(spectrum, fit)
            runtime = time.perf_counter() - start
            score = score_recovery(result, truth)
            score.update({
                "dataset": label,
                "runtime_seconds": runtime,
                "selection_status": result.get("selection_status"),
                "selection_reasons": result.get("selection_reasons", []),
            })
            payload = {
                "truth": truth, "result": result, "score": score,
                "sampler_stdout_tail": stdout.getvalue()[-4000:],
                "sampler_stderr_tail": stderr.getvalue()[-4000:],
            }
            write_json(checkpoint, payload)
            make_diagnostic_plot(spectrum, fit, result, output / "plots" / f"{case['name']}.png")
            scores.append(score)
            print(
                f"[{completed}/{all_cases}] {case['name']}: true={score['true_components']} "
                f"recovered={score['recovered_components']} status={score['selection_status']} "
                f"runtime={runtime:.1f}s", flush=True,
            )
    summary = {
        "beat_version": __version__,
        "profile": args.profile,
        "matrix": args.matrix,
        "empirical_uncertainty_calibration": (
            args.empirical_uncertainty_calibration
        ),
        "correlated_noise": args.correlated_noise,
        "scope": (
            "NIRSpec G235H/G395H controlled and block-resampled-real-residual "
            f"{args.matrix} matrix with {args.profile} profile"
        ),
        "datasets": datasets_manifest,
        "results": summarize(scores),
        "limitations": [
            "The current grid is not statistically powered for completeness or false-positive rates unless the summary explicitly reports statistically_powered=true.",
            "The bundled STScI curves assume a 2.2-pixel fully illuminated aperture and a Gaussian-equivalent FWHM.",
            "Cube resampling correlations are represented only through short residual blocks.",
        ],
    }
    write_json(output / "summary.json", summary)
    write_json(
        output / "experiment_manifest.json",
        {
            "profile": args.profile,
            "matrix": args.matrix,
            "empirical_uncertainty_calibration": (
                args.empirical_uncertainty_calibration
            ),
            "correlated_noise": args.correlated_noise,
            "datasets": datasets_manifest,
        },
    )
    fields = ["spectrum_id", "dataset", "stage", "seed", "true_components", "recovered_components", "component_count_correct", "selection_status", "runtime_seconds", "total_ncall"]
    with (output / "case_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for score in scores:
            writer.writerow({key: score.get(key) for key in fields})
    print(json.dumps(summary["results"], indent=2), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
