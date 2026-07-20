"""Configuration loading and validation for BEAT.

The public configuration remains a plain dictionary so a run can be archived
exactly as the user supplied it.  Validation is deliberately strict: an
expensive nested-sampling run should fail before the first spectrum is fit.
"""

from __future__ import annotations

import copy
import json
import math
import re
from pathlib import Path
from typing import Any

from .adapters import AdapterError, apply_adapter_defaults
from .lsf import LSFError, validate_lsf_config


class ConfigError(ValueError):
    """Raised when a BEAT configuration is internally inconsistent."""


DEFAULTS: dict[str, Any] = {
    "version": 2,
    "output": {
        "directory": "beat_output",
        "resume": True,
        "plots": "selected",  # none, selected, all
    },
    "fit": {
        "frame": "rest",
        "continuum": {"degree": 1, "prior_width_rms": 10.0},
        "noise": {"model": "independent"},
        "kinematics": {
            "max_components": 3,
            "velocity_kms": [-800.0, 800.0],
            "sigma_kms": [30.0, 800.0],
        },
        "selection": {
            "delta_logz": 5.0,
            "stop_when_not_preferred": True,
            "audit": {
                "mode": "flag",
                "uncertainty_sigma": 1.0,
                "minimum_margin": 0.5,
                "max_component_decisive_delta_logz": 20.0,
                "sampling": {
                    "min_num_live_points": 100,
                    "min_ess": 200,
                    "dlogz": 0.5,
                    "nsteps": 20,
                },
            },
        },
        "sampling": {
            "min_num_live_points": 200,
            "min_ess": 400,
            "dlogz": 0.5,
            "show_status": False,
        },
        "flux_prior": {
            "min_snr": 0.1,
            "max_signal_factor": 20.0,
        },
        "minimum_valid_pixels": 20,
    },
}


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _deep_merge(result[key], value)
        else:
            result[key] = copy.deepcopy(value)
    return result


def load_config(path: str | Path) -> dict[str, Any]:
    """Load YAML or JSON, resolve paths relative to the config, and validate."""
    config_path = Path(path).expanduser().resolve()
    if not config_path.exists():
        raise ConfigError(f"Configuration does not exist: {config_path}")

    text = config_path.read_text(encoding="utf-8")
    if config_path.suffix.lower() == ".json":
        raw = json.loads(text)
    else:
        try:
            import yaml
        except ImportError as exc:  # pragma: no cover - installation error
            raise ConfigError("PyYAML is required to read YAML configurations") from exc
        raw = yaml.safe_load(text)

    if not isinstance(raw, dict):
        raise ConfigError("The top level of the configuration must be a mapping")

    config = _deep_merge(DEFAULTS, raw)
    config["_config_path"] = str(config_path)
    config["_config_dir"] = str(config_path.parent)
    try:
        apply_adapter_defaults(config)
    except AdapterError as exc:
        raise ConfigError(str(exc)) from exc
    _resolve_paths(config, config_path.parent)
    validate_config(config)
    return config


def _resolve_paths(config: dict[str, Any], base: Path) -> None:
    input_config = config.get("input", {})
    for key in ("path", "file", "redshift_catalog"):
        value = input_config.get(key)
        if value and not Path(value).expanduser().is_absolute():
            input_config[key] = str((base / value).resolve())

    glob_value = input_config.get("glob")
    if glob_value and not Path(glob_value).expanduser().is_absolute():
        input_config["glob"] = str((base / glob_value).resolve())

    output = config.get("output", {})
    value = output.get("directory")
    if value and not Path(value).expanduser().is_absolute():
        output["directory"] = str((base / value).resolve())

    lsf = config.get("fit", {}).get("lsf", {})
    lsf_path = lsf.get("path")
    if lsf_path and not Path(lsf_path).expanduser().is_absolute():
        lsf["path"] = str((base / lsf_path).resolve())


def _require(mapping: dict[str, Any], key: str, where: str) -> Any:
    if key not in mapping:
        raise ConfigError(f"Missing required setting: {where}.{key}")
    return mapping[key]


def _pair(value: Any, name: str) -> tuple[float, float]:
    if not isinstance(value, (list, tuple)) or len(value) != 2:
        raise ConfigError(f"{name} must contain exactly two numbers")
    try:
        lo, hi = float(value[0]), float(value[1])
    except (TypeError, ValueError) as exc:
        raise ConfigError(f"{name} must contain numbers") from exc
    if not lo < hi:
        raise ConfigError(f"{name} must be increasing")
    return lo, hi


def validate_config(config: dict[str, Any]) -> None:
    """Validate an expanded configuration, raising :class:`ConfigError`."""
    if config.get("version") != 2:
        raise ConfigError("Only configuration version 2 is supported")

    source = _require(config, "input", "root")
    kind = _require(source, "kind", "input")
    if kind not in {"cube", "survey_table", "spectrum_files"}:
        raise ConfigError(
            "input.kind must be cube, survey_table, or spectrum_files"
        )

    if kind in {"cube", "survey_table"}:
        if kind == "cube" and source.get("adapter") == "miri" and source.get("glob"):
            pass
        else:
            _require(source, "path", "input")
    else:
        if not source.get("glob") and not source.get("files"):
            raise ConfigError("spectrum_files input requires input.glob or input.files")

    if kind == "cube":
        _require(source, "flux_hdu", "input")
        if "redshift" not in source:
            raise ConfigError("cube input requires a target redshift")
    elif kind == "survey_table":
        for key in ("hdu", "id_column", "wavelength_column", "flux_column"):
            _require(source, key, "input")
        if "redshift_column" not in source and "redshift" not in source:
            raise ConfigError(
                "survey_table input requires redshift_column or a fixed redshift"
            )
    if "redshift" in source:
        redshift = float(source["redshift"])
        if not math.isfinite(redshift) or redshift <= -1:
            raise ConfigError("input.redshift must be finite and greater than -1")
    adapter = source.get("adapter")
    if adapter in {"miri", "muse", "nirspec"} and not source.get("allow_full_cube", False):
        if source.get("x_range") is None or source.get("y_range") is None:
            raise ConfigError(
                f"{adapter.upper()} runs require both input.x_range and input.y_range; set "
                "allow_full_cube: true only for an intentional full-cube run"
            )
    uncertainty_kind = source.get("uncertainty_kind", "sigma")
    if uncertainty_kind not in {"sigma", "variance", "inverse_variance"}:
        raise ConfigError(
            "input.uncertainty_kind must be sigma, variance, or inverse_variance"
        )

    fit = _require(config, "fit", "root")
    frame = fit.get("frame", "rest")
    if frame not in {"rest", "observed"}:
        raise ConfigError("fit.frame must be rest or observed")
    wavelength_medium = fit.get("wavelength_medium")
    if wavelength_medium is not None and wavelength_medium not in {"air", "vacuum"}:
        raise ConfigError("fit.wavelength_medium must be air or vacuum")
    if adapter == "muse" and wavelength_medium != "air":
        raise ConfigError("MUSE AWAV cubes require fit.wavelength_medium: air")
    if adapter == "nirspec" and wavelength_medium != "vacuum":
        raise ConfigError("JWST/NIRSpec cubes require fit.wavelength_medium: vacuum")
    if adapter == "miri" and wavelength_medium != "vacuum":
        raise ConfigError("JWST/MIRI cubes require fit.wavelength_medium: vacuum")
    _pair(_require(fit, "window", "fit"), "fit.window")
    for index, window in enumerate(fit.get("exclude_windows", [])):
        _pair(window, f"fit.exclude_windows[{index}]")

    continuum = fit.get("continuum", {})
    degree = continuum.get("degree", 1)
    if not isinstance(degree, int) or degree not in {0, 1, 2}:
        raise ConfigError("fit.continuum.degree must be 0, 1, or 2")
    windows = continuum.get("windows", [])
    if windows:
        for index, window in enumerate(windows):
            _pair(window, f"fit.continuum.windows[{index}]")

    noise = fit.get("noise", {})
    if not isinstance(noise, dict):
        raise ConfigError("fit.noise must be a mapping")
    noise_model = noise.get("model", "independent")
    if noise_model not in {"independent", "ar1"}:
        raise ConfigError("fit.noise.model must be independent or ar1")
    if noise_model == "ar1":
        rho = noise.get("rho", "auto")
        if rho != "auto":
            try:
                rho = float(rho)
            except (TypeError, ValueError) as exc:
                raise ConfigError("fit.noise.rho must be auto or a number") from exc
            if not math.isfinite(rho) or not -0.95 <= rho <= 0.95:
                raise ConfigError("fit.noise.rho must lie between -0.95 and 0.95")
    marginal_scale = noise.get("marginal_scale", 1.0)
    if marginal_scale != "auto":
        try:
            marginal_scale = float(marginal_scale)
        except (TypeError, ValueError) as exc:
            raise ConfigError(
                "fit.noise.marginal_scale must be auto or a positive number"
            ) from exc
        if not math.isfinite(marginal_scale) or marginal_scale <= 0:
            raise ConfigError(
                "fit.noise.marginal_scale must be auto or a positive number"
            )

    lines = _require(fit, "lines", "fit")
    if not isinstance(lines, list) or not lines:
        raise ConfigError("fit.lines must be a non-empty list")

    names: list[str] = []
    free_names: set[str] = set()
    for index, line in enumerate(lines):
        if not isinstance(line, dict):
            raise ConfigError(f"fit.lines[{index}] must be a mapping")
        raw_name = _require(line, "name", f"fit.lines[{index}]")
        if not isinstance(raw_name, str) or not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", raw_name):
            raise ConfigError(
                f"fit.lines[{index}].name must be a simple identifier"
            )
        name = raw_name
        if name in names:
            raise ConfigError(f"Duplicate emission-line name: {name}")
        names.append(name)
        wavelength = float(_require(line, "wavelength", f"fit.lines[{index}]"))
        if not math.isfinite(wavelength) or wavelength <= 0:
            raise ConfigError(f"Line {name} has a non-positive wavelength")
        if "ratio_to" in line:
            if line["ratio_to"] not in free_names:
                raise ConfigError(
                    f"Line {name} ratio_to must refer to an earlier free line"
                )
            ratio = float(_require(line, "ratio", f"line {name}"))
            if not math.isfinite(ratio) or ratio <= 0:
                raise ConfigError(f"Line {name} ratio must be positive")
            if "flux_bounds" in line:
                raise ConfigError(f"Locked line {name} cannot define flux_bounds")
        else:
            free_names.add(name)
        if "flux_bounds" in line:
            flux_lo, _ = _pair(line["flux_bounds"], f"line {name}.flux_bounds")
            if flux_lo <= 0:
                raise ConfigError(f"Line {name} flux_bounds must be positive")

    fixed_names: set[str] = set()
    for index, component in enumerate(fit.get("fixed_components", [])):
        where = f"fit.fixed_components[{index}]"
        if not isinstance(component, dict):
            raise ConfigError(f"{where} must be a mapping")
        raw_name = _require(component, "name", where)
        if not isinstance(raw_name, str) or not re.fullmatch(
            r"[A-Za-z][A-Za-z0-9_]*", raw_name
        ):
            raise ConfigError(f"{where}.name must be a simple identifier")
        name = raw_name
        if name in fixed_names:
            raise ConfigError(f"Duplicate fixed-component name: {name}")
        fixed_names.add(name)
        if _require(component, "line", where) not in names:
            raise ConfigError(f"Fixed component {name} refers to an unknown line")
        fixed_sigma = float(_require(component, "sigma_kms", where))
        if not math.isfinite(fixed_sigma) or fixed_sigma <= 0:
            raise ConfigError(f"Fixed component {name} sigma_kms must be positive")
        if "flux_bounds" in component:
            flux_lo, _ = _pair(
                component["flux_bounds"], f"fixed component {name}.flux_bounds"
            )
            if flux_lo <= 0:
                raise ConfigError(
                    f"Fixed component {name} flux_bounds must be positive"
                )

    broad_names: set[str] = set()
    for index, component in enumerate(fit.get("broad_components", [])):
        where = f"fit.broad_components[{index}]"
        if not isinstance(component, dict):
            raise ConfigError(f"{where} must be a mapping")
        raw_name = _require(component, "name", where)
        if not isinstance(raw_name, str) or not re.fullmatch(
            r"[A-Za-z][A-Za-z0-9_]*", raw_name
        ):
            raise ConfigError(f"{where}.name must be a simple identifier")
        name = raw_name
        if name in broad_names:
            raise ConfigError(f"Duplicate broad-component name: {name}")
        broad_names.add(name)
        if _require(component, "line", where) not in names:
            raise ConfigError(f"Broad component {name} refers to an unknown line")
        velocity_lo, velocity_hi = _pair(
            component.get("velocity_kms", [-2000.0, 2000.0]),
            f"broad component {name}.velocity_kms",
        )
        if velocity_lo <= -299_792.458 or velocity_hi >= 299_792.458:
            raise ConfigError(
                f"Broad component {name} velocity bounds exceed the speed of light"
            )
        sigma_lo, sigma_hi = _pair(
            _require(component, "sigma_kms", where),
            f"broad component {name}.sigma_kms",
        )
        if sigma_lo <= 0 or sigma_hi >= 299_792.458:
            raise ConfigError(
                f"Broad component {name} sigma bounds must be positive and subluminal"
            )
        if "flux_bounds" in component:
            flux_lo, _ = _pair(
                component["flux_bounds"], f"broad component {name}.flux_bounds"
            )
            if flux_lo <= 0:
                raise ConfigError(
                    f"Broad component {name} flux_bounds must be positive"
                )

    kin = fit.get("kinematics", {})
    max_components = kin.get("max_components", 3)
    if not isinstance(max_components, int) or max_components < 0:
        raise ConfigError("fit.kinematics.max_components must be a non-negative integer")
    velocity_lo, velocity_hi = _pair(
        kin.get("velocity_kms"), "fit.kinematics.velocity_kms"
    )
    if velocity_lo <= -299_792.458 or velocity_hi >= 299_792.458:
        raise ConfigError("velocity_kms bounds must lie within the speed of light")
    sigma_lo, sigma_hi = _pair(kin.get("sigma_kms"), "fit.kinematics.sigma_kms")
    if sigma_lo <= 0:
        raise ConfigError("The lower sigma_kms bound must be positive")
    if sigma_hi >= 299_792.458:
        raise ConfigError("sigma_kms must be below the speed of light")

    selection = fit.get("selection", {})
    if not isinstance(selection, dict):
        raise ConfigError("fit.selection must be a mapping")
    delta_logz = float(selection.get("delta_logz", 5.0))
    if delta_logz < 0:
        raise ConfigError("fit.selection.delta_logz cannot be negative")
    if not isinstance(selection.get("stop_when_not_preferred", True), bool):
        raise ConfigError("fit.selection.stop_when_not_preferred must be boolean")
    audit = selection.get("audit", {})
    if not isinstance(audit, dict):
        raise ConfigError("fit.selection.audit must be a mapping")
    audit_mode = audit.get("mode", "flag")
    if audit_mode not in {"none", "flag", "rerun"}:
        raise ConfigError("fit.selection.audit.mode must be none, flag, or rerun")
    if float(audit.get("uncertainty_sigma", 1.0)) < 0:
        raise ConfigError("fit.selection.audit.uncertainty_sigma cannot be negative")
    if float(audit.get("minimum_margin", 0.5)) < 0:
        raise ConfigError("fit.selection.audit.minimum_margin cannot be negative")
    decisive = audit.get("max_component_decisive_delta_logz", 20.0)
    if decisive is not None and float(decisive) < delta_logz:
        raise ConfigError(
            "fit.selection.audit.max_component_decisive_delta_logz must be "
            "at least fit.selection.delta_logz or null"
        )
    sampling = fit.get("sampling", {})
    if int(sampling.get("min_num_live_points", 200)) < 40:
        raise ConfigError(
            "fit.sampling.min_num_live_points must be at least 40 for UltraNest"
        )
    if int(sampling.get("min_ess", 400)) < 20:
        raise ConfigError("fit.sampling.min_ess must be at least 20")
    if float(sampling.get("dlogz", 0.5)) <= 0:
        raise ConfigError("fit.sampling.dlogz must be positive")
    if sampling.get("stepsampler", "slice") not in {"slice", "none"}:
        raise ConfigError("fit.sampling.stepsampler must be slice or none")
    if sampling.get("seed") is not None:
        seed = sampling["seed"]
        if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
            raise ConfigError("fit.sampling.seed must be a non-negative integer")
    audit_sampling = audit.get("sampling", {})
    if not isinstance(audit_sampling, dict):
        raise ConfigError("fit.selection.audit.sampling must be a mapping")
    merged_audit_sampling = {**sampling, **audit_sampling}
    if int(merged_audit_sampling.get("min_num_live_points", 200)) < 40:
        raise ConfigError(
            "fit.selection.audit.sampling.min_num_live_points must be at least 40"
        )
    if int(merged_audit_sampling.get("min_ess", 400)) < 20:
        raise ConfigError("fit.selection.audit.sampling.min_ess must be at least 20")
    if float(merged_audit_sampling.get("dlogz", 0.5)) <= 0:
        raise ConfigError("fit.selection.audit.sampling.dlogz must be positive")
    if merged_audit_sampling.get("stepsampler", "slice") not in {"slice", "none"}:
        raise ConfigError(
            "fit.selection.audit.sampling.stepsampler must be slice or none"
        )
    try:
        validate_lsf_config(fit.get("lsf", {"model": "none"}), adapter is not None)
    except LSFError as exc:
        raise ConfigError(str(exc)) from exc

    output = config.get("output", {})
    if output.get("plots", "selected") not in {"none", "selected", "all"}:
        raise ConfigError("output.plots must be none, selected, or all")


def public_config(config: dict[str, Any]) -> dict[str, Any]:
    """Return the serialisable configuration without private loader fields."""
    return {key: value for key, value in config.items() if not key.startswith("_")}
