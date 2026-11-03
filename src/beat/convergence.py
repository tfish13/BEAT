"""Independent-seed diagnostics; no model-selection decisions live here."""

from __future__ import annotations

import hashlib
import itertools
from typing import Any

import numpy as np


def independent_seed(base_seed: int, identifier: str, count: int, repeat: int) -> int:
    payload = f"{base_seed}:{identifier}:{count}:repeat:{repeat}".encode("utf-8")
    return int.from_bytes(hashlib.sha256(payload).digest()[:4], "big")


def _components(model: dict[str, Any]) -> list[dict[str, Any]]:
    median = dict(zip(model["parameter_names"], model["posterior_median"]))
    stdev = dict(zip(model["parameter_names"], model["posterior_stdev"]))
    components = []
    for i in range(1, int(model["n_components"]) + 1):
        prefix = f"component.{i}."
        components.append({
            "velocity": float(median[prefix + "velocity_kms"]),
            "sigma": float(median[prefix + "sigma_kms"]),
            "velocity_error": float(stdev[prefix + "velocity_kms"]),
            "sigma_error": float(stdev[prefix + "sigma_kms"]),
            "fluxes": {
                name[len(prefix):]: (float(value), float(stdev[name]))
                for name, value in median.items()
                if name.startswith(prefix) and name.endswith(".flux")
            },
        })
    return components


def compare_independent_models(
    first: dict[str, Any], second: dict[str, Any], settings: dict[str, Any]
) -> dict[str, Any]:
    """Compare evidence and matched kinematics with uncertainty-aware floors."""
    if first["n_components"] != second["n_components"]:
        raise ValueError("Independent runs must fit the same component count")
    difference = abs(float(first["logz"]) - float(second["logz"]))
    combined_error = float(np.hypot(first["logz_error"], second["logz_error"]))
    evidence_limit = max(
        float(settings.get("minimum_logz_difference", 2.0)),
        float(settings.get("evidence_sigma", 3.0)) * combined_error,
    )
    evidence_flag = bool(not np.isfinite(combined_error) or difference > evidence_limit)
    left, right = _components(first), _components(second)
    distances: list[dict[str, Any]] = []
    if left:
        # The default maximum of three components makes exhaustive matching
        # cheap, and avoids adding a scipy dependency solely for diagnostics.
        best = None
        for permutation in itertools.permutations(range(len(right))):
            proposed = []
            for i, j in enumerate(permutation):
                a, b = left[i], right[j]
                velocity_error = max(
                    float(settings.get("velocity_floor_kms", 10.0)),
                    float(np.hypot(a["velocity_error"], b["velocity_error"])),
                )
                width_error = max(
                    float(settings.get("width_floor_kms", 10.0)),
                    float(np.hypot(a["sigma_error"], b["sigma_error"])),
                )
                flux_distances = []
                for name in a["fluxes"].keys() & b["fluxes"].keys():
                    av, ae = a["fluxes"][name]
                    bv, be = b["fluxes"][name]
                    flux_error = max(
                        float(np.hypot(ae, be)),
                        float(settings.get("flux_floor_fraction", 0.1))
                        * max(abs(av), abs(bv)),
                    )
                    if flux_error > 0:
                        flux_distances.append(abs(av - bv) / flux_error)
                distance = float(np.sqrt(
                    ((a["velocity"] - b["velocity"]) / velocity_error) ** 2
                    + ((a["sigma"] - b["sigma"]) / width_error) ** 2
                    + max(flux_distances, default=0.0) ** 2
                ))
                proposed.append({
                    "first_component": i + 1,
                    "second_component": j + 1,
                    "standardized_distance": distance,
                    "velocity_difference_kms": b["velocity"] - a["velocity"],
                    "width_difference_kms": b["sigma"] - a["sigma"],
                    "maximum_flux_standardized_distance": max(flux_distances, default=0.0),
                })
            score = sum(item["standardized_distance"] for item in proposed)
            if best is None or score < best[0]:
                best = (score, proposed)
        distances = best[1] if best else []
    largest_distance = max(
        (item["standardized_distance"] for item in distances), default=0.0
    )
    parameter_flag = bool(
        not np.isfinite(largest_distance)
        or largest_distance > float(settings.get("parameter_sigma", 3.0))
    )
    return {
        "n_components": int(first["n_components"]),
        "first_logz": float(first["logz"]),
        "second_logz": float(second["logz"]),
        "absolute_logz_difference": difference,
        "combined_logz_error": combined_error,
        "evidence_disagreement_limit": evidence_limit,
        "evidence_disagreement": evidence_flag,
        "component_matches": distances,
        "maximum_component_distance": largest_distance,
        "decomposition_disagreement": parameter_flag,
        "disagreement": evidence_flag or parameter_flag,
    }
