#!/usr/bin/env python3
"""Freeze the completed powered MUSE red gate and its checkpoint hashes."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

from run_muse_powered_validation import (
    MATRIX_VERSION,
    RED_BOUNDARY_EXCLUDED_SCENARIOS,
    in_supported_red_domain,
)


PROJECT = Path(__file__).resolve().parents[1]
OUTPUT = PROJECT / "validation" / "muse_powered_validation"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> int:
    summary_path = OUTPUT / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if not summary.get("all_powered_red_gates_pass"):
        raise RuntimeError("Refusing to freeze: a powered red gate still fails")
    supported, excluded = [], []
    for path in sorted((OUTPUT / "cases").glob("*/result.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        score = payload["score"]
        if score.get("band") != "red":
            continue
        entry = {
            "spectrum_id": score["spectrum_id"],
            "relative_path": str(path.relative_to(PROJECT)),
            "sha256": digest(path),
            "component_class": score["component_class"],
            "selection_status": score["selection_status"],
            "component_count_correct": score["component_count_correct"],
        }
        (supported if in_supported_red_domain(score) else excluded).append(entry)
    if len(supported) != 120:
        raise RuntimeError(f"Expected 120 supported red cases, found {len(supported)}")
    manifest = {
        "matrix_version": MATRIX_VERSION,
        "frozen_at_utc": datetime.now(timezone.utc).isoformat(),
        "decision": "powered MUSE red gate frozen; do not add cases unless a predeclared gate fails or an unexplained in-domain flag is found",
        "supported_red_case_count": len(supported),
        "excluded_boundary_case_count": len(excluded),
        "excluded_scenarios": sorted(RED_BOUNDARY_EXCLUDED_SCENARIOS),
        "summary_sha256": digest(summary_path),
        "gate_results": summary["gate_results"],
        "supported_red_cases": supported,
        "excluded_boundary_cases": excluded,
    }
    path = OUTPUT / "frozen_red_gate_manifest.json"
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "path": str(path),
        "supported_red_case_count": len(supported),
        "excluded_boundary_case_count": len(excluded),
        "all_powered_red_gates_pass": summary["all_powered_red_gates_pass"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
