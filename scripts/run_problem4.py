"""Solve Problem 4, perform grid checks, and create workbook-ready payloads."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from drying_model.problem1 import load_chamber_history_csv
from drying_model.problem3 import extend_chamber_history_to_plateau
from drying_model.problem4 import (
    RadiusHistory,
    load_radius_history_csv,
    moisture_balance_diagnostics,
    sample_moving_solution,
    sample_surface,
    solve_problem4,
    validate_radially_nonincreasing_moisture,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chamber", type=Path, default=Path("data/raw/attachment1.csv"))
    parser.add_argument("--radius", type=Path, default=Path("data/raw/attachment2.csv"))
    parser.add_argument("--output", type=Path, default=Path("tmp/problem4_result.json"))
    parser.add_argument("--figure-dir", type=Path, default=Path("reports/figures"))
    parser.add_argument("--coarse-radial-intervals", type=int, default=640)
    parser.add_argument("--fine-radial-intervals", type=int, default=1280)
    parser.add_argument("--horizon-hours", type=float, default=72.0)
    return parser


def richardson(coarse: np.ndarray, fine: np.ndarray, order: int = 2) -> np.ndarray:
    """Second-order Richardson extrapolation for equal output grids."""
    if coarse.shape != fine.shape:
        raise ValueError("Richardson inputs must have equal shapes")
    return fine + (fine - coarse) / (2.0**order - 1.0)


def interpolate_rows(values: np.ndarray, times: np.ndarray, query: float) -> np.ndarray:
    """Linearly interpolate each spatial column at one time."""
    if query < times[0] or query > times[-1]:
        raise ValueError("Interpolation query lies outside the solution")
    index = int(np.searchsorted(times, query))
    if index == 0:
        return values[0].copy()
    if index == times.size:
        return values[-1].copy()
    if np.isclose(times[index], query, atol=1.0e-10, rtol=0.0):
        return values[index].copy()
    fraction = (query - times[index - 1]) / (times[index] - times[index - 1])
    return values[index - 1] + fraction * (values[index] - values[index - 1])


def save_figures(
    times_s: np.ndarray,
    radius_cm: np.ndarray,
    moisture: np.ndarray,
    surface: np.ndarray,
    radius_history,
    drying_time_s: float,
    figure_dir: Path,
) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    time_h = times_s / 3600.0
    critical_index = int(np.argmin(np.abs(times_s - drying_time_s)))
    selected_hours = [6.0, 12.0, 24.0, 36.0, 48.0]
    selected = [int(np.argmin(np.abs(time_h - value))) for value in selected_hours if value <= time_h[-1]]
    selected.extend([critical_index, len(times_s) - 1])
    selected = list(dict.fromkeys(selected))
    fig, axis = plt.subplots(figsize=(7.4, 4.7), constrained_layout=True)
    for index in selected:
        label = f"{time_h[index]:g} h"
        if np.isclose(times_s[index], drying_time_s, atol=1.0e-6):
            label = "critical threshold"
        elif index == len(times_s) - 1:
            label = "first strict regular time"
        axis.plot(radius_cm, moisture[index], marker="", label=label)
    axis.axhline(0.15, color="#B91C1C", linestyle="--", linewidth=1.0)
    axis.set(xlabel="Fixed physical radius (cm)", ylabel="Moisture concentration (kg/kg)")
    axis.grid(alpha=0.25)
    axis.legend(ncol=2, fontsize=8)
    fig.savefig(figure_dir / "problem4_radial_profiles.png", dpi=220)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.4, 4.7), constrained_layout=True)
    center = moisture[:, 0]
    axis.plot(time_h, center, label="center")
    axis.plot(time_h, surface, label="actual surface")
    axis.axhline(0.15, color="#B91C1C", linestyle="--", linewidth=1.0)
    axis.set(xlabel="Time (h)", ylabel="Moisture concentration (kg/kg)")
    axis.grid(alpha=0.25)
    axis.legend()
    fig.savefig(figure_dir / "problem4_time_histories.png", dpi=220)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.4, 4.7), constrained_layout=True)
    axis.plot(radius_history.time_s / 3600.0, radius_history.radius_m * 100.0, color="#1D4ED8")
    axis.set(xlabel="Time (h)", ylabel="Measured radius (cm)")
    axis.grid(alpha=0.25)
    fig.savefig(figure_dir / "problem4_radius_history.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = build_parser().parse_args()
    chamber = load_chamber_history_csv(args.chamber)
    radius_history = load_radius_history_csv(args.radius)
    horizon_s = min(float(args.horizon_hours * 3600.0), float(radius_history.time_s[-1]))
    extended_chamber = extend_chamber_history_to_plateau(chamber, horizon_s)
    threshold = 0.15
    solver_options = {"relative_tolerance": 2.0e-8, "max_step_s": 120.0}
    endpoint_times: dict[str, float] = {}
    for grid in (args.coarse_radial_intervals, args.fine_radial_intervals):
        probe_times = np.arange(3600.0, horizon_s + 1.0e-9, 3600.0)
        endpoint = solve_problem4(
            extended_chamber,
            radius_history,
            grid,
            probe_times,
            maximum_moisture_threshold=threshold,
            **solver_options,
        )
        endpoint_times[str(grid)] = float(endpoint.time_s[-1])
        print(f"endpoint grid={grid}: {endpoint.time_s[-1]:.6f} s", flush=True)
    drying_time_s = float(
        endpoint_times[str(args.fine_radial_intervals)]
        + (endpoint_times[str(args.fine_radial_intervals)] - endpoint_times[str(args.coarse_radial_intervals)]) / 3.0
    )
    timing_sensitivity: dict[str, float] = {}
    sensitivity_grid = max(160, args.coarse_radial_intervals // 2)
    for name, shift_s in (("half_interval_early", -900.0), ("baseline", 0.0), ("half_interval_late", 900.0)):
        shifted_times = radius_history.time_s.copy()
        shifted_times[1:-1] += shift_s
        shifted_radius = RadiusHistory(shifted_times, radius_history.radius_m.copy())
        sensitivity_solution = solve_problem4(
            extended_chamber,
            shifted_radius,
            sensitivity_grid,
            np.arange(3600.0, horizon_s + 1.0e-9, 3600.0),
            maximum_moisture_threshold=threshold,
            **solver_options,
        )
        timing_sensitivity[name] = float(sensitivity_solution.time_s[-1])
    if not (0.0 < drying_time_s < horizon_s):
        raise RuntimeError("Extrapolated drying endpoint is outside the radius-history horizon")
    strict_time_s = float(np.floor(drying_time_s / 60.0 + 1.0) * 60.0)
    if strict_time_s <= drying_time_s + 1.0e-8:
        strict_time_s += 60.0

    regular_times = np.arange(60.0, strict_time_s + 1.0e-9, 60.0)
    bracket = np.array([drying_time_s - 60.0, drying_time_s + 60.0])
    integration_times = np.unique(np.concatenate([regular_times, bracket]))
    fixed_radius_cm = np.round(np.arange(0.0, 2.0 + 0.0001, 0.1), 10)
    sampled_by_grid: list[np.ndarray] = []
    surface_by_grid: list[np.ndarray] = []
    balance_by_grid: dict[str, dict[str, float]] = {}
    monotonicity_by_grid: dict[str, float] = {}
    diagnostic_times = np.unique(np.concatenate([np.arange(0.0, 601.0), integration_times]))
    for grid in (args.coarse_radial_intervals, args.fine_radial_intervals):
        full = solve_problem4(
            extended_chamber,
            radius_history,
            grid,
            diagnostic_times,
            **solver_options,
        )
        balance_by_grid[str(grid)] = moisture_balance_diagnostics(
            full, extended_chamber, radius_history
        )
        monotonicity_by_grid[str(grid)] = validate_radially_nonincreasing_moisture(full)
        sampled_full = sample_moving_solution(full, radius_history, fixed_radius_cm)
        surface_full = sample_surface(full)
        retained = np.searchsorted(diagnostic_times, integration_times)
        sampled_by_grid.append(sampled_full[retained])
        surface_by_grid.append(surface_full[retained])
        del full
        print(f"full output grid={grid} complete", flush=True)

    extrapolated = richardson(sampled_by_grid[0], sampled_by_grid[1])
    extrapolated_surface = richardson(surface_by_grid[0], surface_by_grid[1])
    critical_fixed = interpolate_rows(extrapolated, integration_times, drying_time_s)
    critical_surface = float(interpolate_rows(extrapolated_surface[:, None], integration_times, drying_time_s)[0])
    strict_index = int(np.where(np.isclose(integration_times, strict_time_s, atol=1.0e-9))[0][0])
    strict_fixed = extrapolated[strict_index]
    strict_surface = float(extrapolated_surface[strict_index])
    regular_before = regular_times[regular_times < drying_time_s - 1.0e-8]
    output_times = np.concatenate([regular_before, [drying_time_s, strict_time_s]])
    output_fixed = np.vstack([
        extrapolated[np.searchsorted(integration_times, regular_before)],
        critical_fixed,
        strict_fixed,
    ])
    output_surface = np.concatenate([
        extrapolated_surface[np.searchsorted(integration_times, regular_before)],
        [critical_surface, strict_surface],
    ])

    endpoint_moisture = float(np.nanmax(critical_fixed))
    strict_maximum = float(np.nanmax(strict_fixed))
    if not np.isclose(endpoint_moisture, threshold, atol=3.0e-4):
        raise RuntimeError(f"Critical moisture does not equal threshold: {endpoint_moisture}")
    if not strict_maximum < threshold:
        raise RuntimeError(f"Strict regular row is not below threshold: {strict_maximum}")
    rounded_fixed = np.round(output_fixed, 4)
    workbook_matrix = [
        [None if not np.isfinite(value) else float(value) for value in row]
        for row in rounded_fixed
    ]
    payload = {
        "time_s": [int(round(t)) if np.isclose(t, round(t), rtol=0.0, atol=1.0e-8) else round(float(t), 6) for t in output_times],
        "radius_cm": fixed_radius_cm.tolist(),
        "moisture_concentration": workbook_matrix,
        "surface_moisture_concentration": np.round(output_surface, 4).tolist(),
        "drying_time_s": round(drying_time_s, 6),
        "drying_time_h": round(drying_time_s / 3600.0, 10),
        "strict_completion_time_s": int(round(strict_time_s)),
        "strict_completion_time_h": round(strict_time_s / 3600.0, 10),
        "threshold": threshold,
        "critical_maximum_moisture_unrounded": endpoint_moisture,
        "strict_completion_maximum_moisture_unrounded": strict_maximum,
        "radius_at_drying_time_cm": float(np.interp(drying_time_s, radius_history.time_s, radius_history.radius_m) * 100.0),
        "grid_endpoint_s": endpoint_times | {"richardson": drying_time_s},
        "outside_radius_cells_are_blank": True,
        "moisture_balance": balance_by_grid,
        "radial_monotonicity": monotonicity_by_grid,
        "radius_timing_sensitivity_s": timing_sensitivity,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":"), allow_nan=False),
        encoding="utf-8",
    )
    save_figures(output_times, fixed_radius_cm, output_fixed, output_surface, radius_history, drying_time_s, args.figure_dir)
    print(f"drying endpoint: {drying_time_s / 3600.0:.6f} h; strict row: {strict_time_s:.0f} s", flush=True)


if __name__ == "__main__":
    main()
