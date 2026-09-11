"""Reproduce the Q4 review experiments; no external review values are inputs."""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
import scipy

from drying_model.problem1 import Problem1Parameters, load_chamber_history_csv
from drying_model.problem3 import extend_chamber_history_to_plateau
from drying_model.problem4 import RadiusHistory, density_shrinkage_compatibility, load_radius_history_csv, solve_problem4


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("reports/data/problem4_review_experiments.json"))
    parser.add_argument("--sensitivity-grid", type=int, default=320)
    args = parser.parse_args()
    chamber_path = Path("data/raw/attachment1.csv")
    radius_path = Path("data/raw/attachment2.csv")
    chamber = extend_chamber_history_to_plateau(load_chamber_history_csv(chamber_path), 144 * 3600.0)
    radius = load_radius_history_csv(radius_path)
    parameters = Problem1Parameters()
    results = {
        "schema_version": 1,
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "input_sha256": {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in (chamber_path, radius_path)},
        "threshold": 0.15,
        "default_solver": {"relative_tolerance": 2e-8, "max_step_s": 120.0},
        "sensitivity_radial_intervals": args.sensitivity_grid,
    }

    def save() -> None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(results, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")

    def endpoint(grid, radius_input=radius, params=parameters, rtol=2e-8, step=120.0):
        solution = solve_problem4(chamber, radius_input, grid, [0.0, 144 * 3600.0], parameters=params,
                                  relative_tolerance=rtol, max_step_s=step, maximum_moisture_threshold=0.15)
        return solution

    fixed = RadiusHistory(radius.time_s, np.full_like(radius.radius_m, parameters.radius_m))
    fixed_endpoints = {}
    results["fixed_radius_endpoint_s"] = fixed_endpoints
    for grid in (160, 320, 640, 1280):
        fixed_endpoints[str(grid)] = float(endpoint(grid, fixed).time_s[-1])
        print(f"fixed radius N={grid}: {fixed_endpoints[str(grid)] / 3600:.9f} h", flush=True)
        save()
    fixed_endpoints["richardson_640_1280"] = fixed_endpoints["1280"] + (fixed_endpoints["1280"] - fixed_endpoints["640"]) / 3
    fixed_endpoints["observed_order_320_640_1280"] = float(np.log2((fixed_endpoints["320"] - fixed_endpoints["640"]) / (fixed_endpoints["640"] - fixed_endpoints["1280"])))

    time_checks = []
    results["time_integration"] = time_checks
    for rtol, step in ((2e-8, 120.0), (2e-9, 60.0), (5e-10, 30.0)):
        solution = endpoint(args.sensitivity_grid, rtol=rtol, step=step)
        time_checks.append({"relative_tolerance": rtol, "max_step_s": step, "endpoint_s": float(solution.time_s[-1])})
        print(f"time convergence: {time_checks[-1]}", flush=True)
        save()
    baseline = time_checks[0]["endpoint_s"]
    results["time_integration_endpoint_spread_s"] = max(row["endpoint_s"] for row in time_checks) - min(row["endpoint_s"] for row in time_checks)

    pchip = RadiusHistory(radius.time_s, radius.radius_m, interpolation="pchip")
    pchip_time = float(endpoint(args.sensitivity_grid, pchip).time_s[-1])
    results["radius_interpolation"] = {"linear_endpoint_s": baseline, "pchip_endpoint_s": pchip_time,
                                       "pchip_minus_linear_s": pchip_time - baseline, "observation_times_unchanged": True}
    save()
    print(f"PCHIP minus linear: {pchip_time - baseline:.6f} s", flush=True)

    for name, field in (("mass_transfer", "mass_transfer_coefficient_m_s"), ("heat_transfer", "heat_transfer_coefficient_w_m2_k")):
        rows = []
        results[name + "_sensitivity"] = rows
        for multiplier in (0.8, 1.0, 1.2):
            value = getattr(parameters, field) * multiplier
            time = baseline if multiplier == 1.0 else float(endpoint(args.sensitivity_grid, params=replace(parameters, **{field: value})).time_s[-1])
            rows.append({"multiplier": multiplier, "coefficient": value, "endpoint_s": time, "change_from_baseline_s": time - baseline})
            print(f"{name} x{multiplier}: {time / 3600:.9f} h", flush=True)
            save()

    density_checks = {}
    results["density_shrinkage_compatibility"] = density_checks
    for grid in (640, 1280):
        event_solution = endpoint(grid)
        times = np.append(np.arange(0.0, 49 * 3600.0, 6 * 3600.0), event_solution.time_s[-1])
        full = solve_problem4(chamber, radius, grid, times, relative_tolerance=2e-8, max_step_s=120.0)
        density_checks[str(grid)] = density_shrinkage_compatibility(full, radius)
        density_checks[str(grid)]["center_moisture_at_48h"] = float(full.moisture_concentration[-2, 0])
        print(f"density compatibility N={grid}: final mass ratio {density_checks[str(grid)]['final_mass_ratio']:.9f}", flush=True)
        save()
    save()


if __name__ == "__main__":
    main()
