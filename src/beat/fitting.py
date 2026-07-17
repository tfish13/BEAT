"""UltraNest orchestration and result serialization for one spectrum."""

from __future__ import annotations

import copy
from pathlib import Path
from typing import Any, Callable

import numpy as np

from .lsf import lsf_sigma_angstrom
from .model import ModelDefinition, observed_center, observed_sigma, prepare_spectrum
from .spectrum import Spectrum


SamplerFactory = Callable[[list[str], Callable[..., float], Callable[..., np.ndarray]], Any]


def _default_sampler_factory(
    names: list[str], loglike: Callable[..., float], transform: Callable[..., np.ndarray]
) -> Any:
    try:
        import ultranest
    except ImportError as exc:  # pragma: no cover - installation error
        raise RuntimeError(
            "UltraNest is required for fitting; install BEAT with `pip install -e .`"
        ) from exc
    return ultranest.ReactiveNestedSampler(
        names,
        loglike,
        transform,
        log_dir=None,
        resume="overwrite",
    )


def _configure_stepsampler(sampler: Any, sampling: dict[str, Any], ndim: int) -> None:
    if sampling.get("stepsampler", "slice") == "none":
        return
    try:
        import ultranest.stepsampler
    except ImportError:  # fake samplers in tests
        return
    sampler.stepsampler = ultranest.stepsampler.SliceSampler(
        nsteps=int(sampling.get("nsteps", max(ndim, 10))),
        generate_direction=ultranest.stepsampler.generate_mixture_random_direction,
    )


def _result_value(result: dict[str, Any], group: str, key: str, fallback: Any) -> Any:
    value = result.get(group, {}).get(key, fallback)
    return np.asarray(value, dtype=float)


def _run_model(
    model: ModelDefinition,
    sampling: dict[str, Any],
    sampler_factory: SamplerFactory,
) -> dict[str, Any]:
    sampler = sampler_factory(
        model.parameter_names, model.log_likelihood, model.prior_transform
    )
    _configure_stepsampler(sampler, sampling, model.ndim)
    run_kwargs = {
        "min_num_live_points": int(sampling.get("min_num_live_points", 200)),
        "min_ess": int(sampling.get("min_ess", 400)),
        "dlogz": float(sampling.get("dlogz", 0.5)),
        "show_status": bool(sampling.get("show_status", False)),
        "viz_callback": False,
    }
    result = sampler.run(**run_kwargs)
    maximum = _result_value(
        result, "maximum_likelihood", "point", np.zeros(model.ndim)
    )
    median = _result_value(result, "posterior", "median", maximum)
    stdev = _result_value(result, "posterior", "stdev", np.full(model.ndim, np.nan))
    return {
        "n_components": model.n_components,
        "logz": float(result["logz"]),
        "logz_error": float(result.get("logzerr", np.nan)),
        "parameter_names": list(model.parameter_names),
        "posterior_median": median.tolist(),
        "posterior_stdev": stdev.tolist(),
        "maximum_likelihood": maximum.tolist(),
        "maximum_log_likelihood": float(
            result.get("maximum_likelihood", {}).get("logl", np.nan)
        ),
        "ncall": int(result.get("ncall", 0)),
    }


def _fit_spectrum_once(
    spectrum: Spectrum,
    fit: dict[str, Any],
    sampler_factory: SamplerFactory | None = None,
) -> dict[str, Any]:
    """Fit zero through N narrow components once and select by evidence."""
    prepared = prepare_spectrum(spectrum, fit)
    factory = sampler_factory or _default_sampler_factory
    sampling = fit.get("sampling", {})
    selection = fit.get("selection", {})
    threshold = float(selection.get("delta_logz", 5.0))
    max_components = int(fit["kinematics"]["max_components"])
    stop_early = bool(selection.get("stop_when_not_preferred", True))
    models: list[dict[str, Any]] = []
    selected_index = 0

    for n_components in range(max_components + 1):
        model = ModelDefinition(prepared, fit, n_components)
        model_result = _run_model(model, sampling, factory)
        models.append(model_result)
        if n_components == 0:
            continue
        improvement = model_result["logz"] - models[selected_index]["logz"]
        if improvement > threshold:
            selected_index = n_components
        elif stop_early:
            break

    selected = models[selected_index]
    medians = dict(zip(selected["parameter_names"], selected["posterior_median"]))
    stdevs = dict(zip(selected["parameter_names"], selected["posterior_stdev"]))
    derived_broad_components = []
    line_by_name = {line["name"]: line for line in fit["lines"]}
    for component in fit.get("broad_components", []):
        prefix = f"broad.{component['name']}"
        velocity = float(medians[f"{prefix}.velocity_kms"])
        sigma_kms = float(medians[f"{prefix}.sigma_kms"])
        flux = float(medians[f"{prefix}.flux"])
        line = line_by_name[component["line"]]
        rest = float(line["wavelength"])
        center = observed_center(rest, prepared.spectrum.redshift, velocity)
        intrinsic_sigma = observed_sigma(
            rest, prepared.spectrum.redshift, sigma_kms, velocity
        )
        instrumental_sigma = float(
            lsf_sigma_angstrom(
                center,
                fit.get("lsf", {"model": "none"}),
                prepared.spectrum.metadata,
            )
        )
        derived_broad_components.append(
            {
                "name": component["name"],
                "line": component["line"],
                "velocity_kms": velocity,
                "velocity_stdev_kms": float(stdevs[f"{prefix}.velocity_kms"]),
                "sigma_kms": sigma_kms,
                "sigma_stdev_kms": float(stdevs[f"{prefix}.sigma_kms"]),
                "flux": flux,
                "flux_stdev": float(stdevs[f"{prefix}.flux"]),
                "observed_center_angstrom": center,
                "intrinsic_sigma_angstrom": intrinsic_sigma,
                "lsf_sigma_angstrom": instrumental_sigma,
                "convolved_sigma_angstrom": float(
                    np.hypot(intrinsic_sigma, instrumental_sigma)
                ),
            }
        )
    derived_components = []
    for component_number in range(1, int(selected["n_components"]) + 1):
        prefix = f"component.{component_number}"
        velocity = float(medians[f"{prefix}.velocity_kms"])
        sigma_kms = float(medians[f"{prefix}.sigma_kms"])
        line_summaries: dict[str, dict[str, float]] = {}
        for line in fit["lines"]:
            if "ratio_to" in line:
                reference = line_summaries[line["ratio_to"]]
                ratio = float(line["ratio"])
                flux = reference["flux"] * ratio
                flux_stdev = reference["flux_stdev"] * ratio
            else:
                flux_name = f"{prefix}.{line['name']}.flux"
                flux = float(medians[flux_name])
                flux_stdev = float(stdevs[flux_name])
            rest = float(line["wavelength"])
            center = observed_center(rest, prepared.spectrum.redshift, velocity)
            intrinsic_sigma = observed_sigma(
                rest, prepared.spectrum.redshift, sigma_kms, velocity
            )
            instrumental_sigma = float(
                lsf_sigma_angstrom(
                    center,
                    fit.get("lsf", {"model": "none"}),
                    prepared.spectrum.metadata,
                )
            )
            total_sigma = float(np.hypot(intrinsic_sigma, instrumental_sigma))
            line_summaries[line["name"]] = {
                "flux": flux,
                "flux_stdev": flux_stdev,
                "observed_center_angstrom": center,
                "intrinsic_sigma_angstrom": intrinsic_sigma,
                "lsf_sigma_angstrom": instrumental_sigma,
                "convolved_sigma_angstrom": total_sigma,
                # Compatibility alias for 2.0.0a1 result readers.
                "observed_sigma_angstrom": total_sigma,
            }
        derived_components.append(
            {
                "component": component_number,
                "velocity_kms": velocity,
                "sigma_kms": sigma_kms,
                "lines": line_summaries,
            }
        )
    return {
        "status": "ok",
        "spectrum_id": prepared.spectrum.spectrum_id,
        "redshift": prepared.spectrum.redshift,
        "metadata": prepared.spectrum.metadata,
        "n_input_pixels": int(prepared.spectrum.wavelength.size),
        "n_fit_pixels": int(prepared.wavelength.size),
        "noise_level": prepared.noise_level,
        "noise_model": prepared.noise_model,
        "noise_rho": prepared.noise_rho,
        "noise_marginal_scale": prepared.noise_marginal_scale,
        "selected_components": int(selected["n_components"]),
        "maximum_components_considered": max_components,
        "selected_logz": float(selected["logz"]),
        "components": derived_components,
        "broad_components": derived_broad_components,
        "models": models,
    }


def _selection_diagnostics(
    result: dict[str, Any], selection: dict[str, Any]
) -> list[dict[str, Any]]:
    """Reconstruct each evidence decision using the same sequential rule."""
    threshold = float(selection.get("delta_logz", 5.0))
    sigma_factor = float(selection.get("audit", {}).get("uncertainty_sigma", 1.0))
    minimum_margin = float(selection.get("audit", {}).get("minimum_margin", 0.5))
    selected_index = 0
    diagnostics = []
    models = result["models"]
    for model in models[1:]:
        reference = models[selected_index]
        improvement = float(model["logz"] - reference["logz"])
        candidate_error = float(model.get("logz_error", np.nan))
        reference_error = float(reference.get("logz_error", np.nan))
        combined_error = float(np.hypot(candidate_error, reference_error))
        finite_error = bool(np.isfinite(combined_error))
        uncertainty_margin = (
            max(minimum_margin, sigma_factor * combined_error)
            if finite_error
            else None
        )
        distance = abs(improvement - threshold)
        overlaps = bool(
            uncertainty_margin is not None and distance <= uncertainty_margin
        )
        preferred = bool(improvement > threshold)
        diagnostics.append(
            {
                "candidate_components": int(model["n_components"]),
                "reference_components": int(reference["n_components"]),
                "delta_logz": improvement,
                "threshold": threshold,
                "combined_logz_error": combined_error if finite_error else None,
                "uncertainty_margin": uncertainty_margin,
                "distance_from_threshold": distance,
                "uncertainty_overlaps_threshold": overlaps,
                "preferred": preferred,
            }
        )
        if preferred:
            selected_index = int(model["n_components"])
    return diagnostics


def _selection_assessment(
    result: dict[str, Any], selection: dict[str, Any], audited: bool = False
) -> dict[str, Any]:
    diagnostics = _selection_diagnostics(result, selection)
    reasons = []
    ambiguous = [
        item for item in diagnostics if item["uncertainty_overlaps_threshold"]
    ]
    if ambiguous:
        reasons.append("evidence uncertainty overlaps the selection threshold")
    nonfinite = [
        item for item in diagnostics if item["combined_logz_error"] is None
    ]
    if nonfinite:
        reasons.append("one or more evidence uncertainties are not finite")

    audit = selection.get("audit", {})
    decisive = audit.get("max_component_decisive_delta_logz", 20.0)
    max_components = int(
        result.get(
            "maximum_components_considered",
            max(int(model["n_components"]) for model in result["models"]),
        )
    )
    selected_components = int(result["selected_components"])
    max_support = next(
        (
            item["delta_logz"]
            for item in diagnostics
            if item["candidate_components"] == selected_components
            and item["preferred"]
        ),
        None,
    )
    moderate_max_selection = bool(
        decisive is not None
        and selected_components >= 2
        and selected_components == max_components
        and max_support is not None
        and max_support < float(decisive)
    )
    if moderate_max_selection:
        reasons.append(
            "maximum component count is selected without decisive evidence"
        )

    if ambiguous:
        status = "ambiguous"
    elif nonfinite:
        status = "convergence_unverified"
    elif audited:
        status = "accepted_after_audit"
    elif moderate_max_selection:
        status = "convergence_unverified"
    else:
        status = "accepted"
    return {
        "status": status,
        "audit_recommended": bool(ambiguous or nonfinite or moderate_max_selection),
        "reasons": reasons,
        "comparisons": diagnostics,
    }


def _audit_sampling(
    sampling: dict[str, Any], audit_sampling: dict[str, Any]
) -> dict[str, Any]:
    """Merge an audit profile without making minimum rigor settings looser."""
    merged = copy.deepcopy(sampling)
    merged.update(copy.deepcopy(audit_sampling))
    for key, fallback in (
        ("min_num_live_points", 200),
        ("min_ess", 400),
        ("nsteps", 10),
    ):
        merged[key] = max(
            int(sampling.get(key, fallback)),
            int(audit_sampling.get(key, fallback)),
        )
    merged["dlogz"] = min(
        float(sampling.get("dlogz", 0.5)),
        float(audit_sampling.get("dlogz", 0.5)),
    )
    return merged


def fit_spectrum(
    spectrum: Spectrum,
    fit: dict[str, Any],
    sampler_factory: SamplerFactory | None = None,
) -> dict[str, Any]:
    """Fit a spectrum, assess evidence reliability, and optionally audit it."""
    selection = fit.get("selection", {})
    audit = selection.get("audit", {})
    mode = audit.get("mode", "flag")
    pilot = _fit_spectrum_once(spectrum, fit, sampler_factory=sampler_factory)
    pilot_assessment = _selection_assessment(pilot, selection, audited=False)
    result = pilot
    audit_record: dict[str, Any] = {
        "mode": mode,
        "performed": False,
        "trigger_reasons": pilot_assessment["reasons"],
    }

    if mode == "rerun" and pilot_assessment["audit_recommended"]:
        audited_fit = copy.deepcopy(fit)
        audited_fit["sampling"] = _audit_sampling(
            fit.get("sampling", {}), audit.get("sampling", {})
        )
        result = _fit_spectrum_once(
            spectrum, audited_fit, sampler_factory=sampler_factory
        )
        final_assessment = _selection_assessment(result, selection, audited=True)
        audit_record.update(
            {
                "performed": True,
                "pilot_selected_components": pilot["selected_components"],
                "pilot_selected_logz": pilot["selected_logz"],
                "pilot_models": pilot["models"],
                "pilot_selection_diagnostics": pilot_assessment["comparisons"],
                "pilot_sampling": copy.deepcopy(fit.get("sampling", {})),
                "audit_sampling": copy.deepcopy(audited_fit["sampling"]),
            }
        )
    else:
        final_assessment = pilot_assessment
        if mode == "none" and final_assessment["status"] != "accepted":
            final_assessment["status"] = "convergence_unverified"

    result["selection_status"] = final_assessment["status"]
    result["selection_reasons"] = final_assessment["reasons"]
    result["selection_diagnostics"] = final_assessment["comparisons"]
    result["selection_audit"] = audit_record
    return result


def selected_model(result: dict[str, Any]) -> dict[str, Any]:
    target = int(result["selected_components"])
    return next(model for model in result["models"] if model["n_components"] == target)


def make_diagnostic_plot(
    spectrum: Spectrum,
    fit: dict[str, Any],
    result: dict[str, Any],
    output_path: str | Path,
    n_components: int | None = None,
) -> None:
    """Render data, uncertainty, total model, residual, and components."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError as exc:  # pragma: no cover - installation error
        raise RuntimeError("Matplotlib is required to make plots") from exc

    prepared = prepare_spectrum(spectrum, fit)
    if n_components is None:
        n_components = int(result["selected_components"])
    model_result = next(
        item for item in result["models"] if item["n_components"] == n_components
    )
    definition = ModelDefinition(prepared, fit, n_components)
    parameters = np.asarray(model_result["maximum_likelihood"], dtype=float)
    total = definition.evaluate(parameters)
    components = definition.component_models(parameters)
    component_labels = [
        f"fixed {item['name']}" for item in fit.get("fixed_components", [])
    ]
    component_labels.extend(
        f"broad {item['name']}" for item in fit.get("broad_components", [])
    )
    component_labels.extend(
        f"narrow component {index}" for index in range(1, n_components + 1)
    )
    continuum = np.polynomial.polynomial.polyval(
        prepared.x_normalized, parameters[: definition.continuum_count]
    )
    residual = prepared.flux - total

    fig, (ax, residual_ax) = plt.subplots(
        2, 1, figsize=(9, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
    )
    ax.step(prepared.wavelength, prepared.flux, where="mid", color="0.35", lw=1, label="data")
    ax.plot(prepared.wavelength, total, color="black", lw=1.5, label="total model")
    ax.plot(prepared.wavelength, continuum, color="tab:green", lw=1, label="continuum")
    for index, component in enumerate(components):
        ax.plot(
            prepared.wavelength,
            continuum + component,
            lw=0.9,
            alpha=0.8,
            label=component_labels[index],
        )
    ax.fill_between(
        prepared.wavelength,
        prepared.flux - prepared.uncertainty,
        prepared.flux + prepared.uncertainty,
        color="0.75",
        alpha=0.3,
        step="mid",
        label=r"$1\sigma$ uncertainty",
    )
    ax.set_ylabel("Flux density")
    ax.set_title(
        f"{spectrum.spectrum_id}: {n_components} component(s), "
        f"ln Z={model_result['logz']:.2f}"
    )
    ax.legend(fontsize=8, ncol=2)

    residual_ax.axhline(0, color="black", lw=0.8)
    residual_ax.step(prepared.wavelength, residual, where="mid", color="tab:blue", lw=1)
    residual_ax.fill_between(
        prepared.wavelength,
        -prepared.uncertainty,
        prepared.uncertainty,
        color="0.75",
        alpha=0.4,
        step="mid",
    )
    residual_ax.set_xlabel(r"Observed wavelength ($\AA$)")
    residual_ax.set_ylabel("Residual")
    fig.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=160, bbox_inches="tight")
    plt.close(fig)
