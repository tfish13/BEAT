"""Command-line interface for BEAT."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from . import __version__
from .config import ConfigError, load_config
from .data import iter_spectra
from .model import prepare_spectrum
from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="beat", description="Bayesian multi-component emission-line fitting"
    )
    parser.add_argument("--version", action="version", version=f"BEAT {__version__}")
    subparsers = parser.add_subparsers(dest="command", required=True)

    validate = subparsers.add_parser(
        "validate", help="validate configuration and inspect the first spectrum"
    )
    validate.add_argument("config", type=Path)
    validate.add_argument(
        "--config-only", action="store_true", help="do not open the input dataset"
    )

    run = subparsers.add_parser("run", help="fit every configured spectrum")
    run.add_argument("config", type=Path)
    run.add_argument("--workers", type=int, default=None)
    run.add_argument(
        "--no-resume", action="store_true", help="refit completed spectrum IDs"
    )
    return parser


def _validate(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    print(f"Configuration valid: {args.config.resolve()}")
    if args.config_only:
        return 0
    try:
        first = next(iter_spectra(config))
    except StopIteration as exc:
        raise ConfigError("The input contains no spectra") from exc
    prepared = prepare_spectrum(first, config["fit"])
    summary = {
        "first_spectrum_id": first.spectrum_id,
        "redshift": first.redshift,
        "input_pixels": int(first.wavelength.size),
        "valid_fit_pixels": int(prepared.wavelength.size),
        "fit_wavelength_angstrom": [
            float(prepared.wavelength.min()),
            float(prepared.wavelength.max()),
        ],
        "median_uncertainty": prepared.noise_level,
        "noise_model": prepared.noise_model,
        "noise_rho": prepared.noise_rho,
        "noise_marginal_scale": prepared.noise_marginal_scale,
        "default_integrated_flux_prior": list(prepared.default_flux_bounds),
        "selection_audit_mode": config["fit"].get("selection", {})
        .get("audit", {})
        .get("mode", "flag"),
    }
    print(json.dumps(summary, indent=2))
    return 0


def _run(args: argparse.Namespace) -> int:
    config = load_config(args.config)
    if args.no_resume:
        config["output"]["resume"] = False
    summary = run_pipeline(config, workers=args.workers)
    print(json.dumps(summary, indent=2))
    return 1 if summary["failed_this_run"] else 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "validate":
            return _validate(args)
        if args.command == "run":
            return _run(args)
    except (ConfigError, ValueError, OSError, RuntimeError) as exc:
        parser.exit(2, f"BEAT error: {exc}\n")
    return 2


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
