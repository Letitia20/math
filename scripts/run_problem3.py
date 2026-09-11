"""Solve Problem 3, locate the drying endpoint, and export workbook-ready values."""

from __future__ import annotations

import argparse
import gc
import json
from dataclasses import replace
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from drying_model.problem1 import (
    load_chamber_history_csv,
    richardson_extrapolate_solutions,
    sample_solution,
)
from drying_model.problem2 import solve_problem2
from drying_model.problem3 import (
    Problem3Result,
    extend_chamber_history_to_plateau,
    problem3_payload,
    richardson_extrapolate_drying_time,
    solve_problem3,
    truncate_solution_at_threshold,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/raw/attachment1.csv"))
    parser.add_argument("--output", type=Path, default=Path("tmp/problem3_result.json"))
    parser.add_argument("--figure-dir", type=Path, default=Path("reports/figures"))
    parser.add_argument("--coarse-radial-intervals", type=int, default=2560)
    parser.add_argument("--fine-radial-intervals", type=int, default=5120)
    parser.add_argument("--sensitivity-radial-intervals", type=int, default=640)
    parser.add_argument("--horizon-hours", type=float, default=96.0)
    return parser


def save_figures(result: Problem3Result, figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    solution = result.solution
    time_h = solution.time_s / 3600.0
    radius_cm = solution.radius_m * 100.0
    moisture = solution.moisture_concentration

    selected_hours = [6.0, 12.0, 24.0, 36.0, 48.0]
    selected_hours = [hour for hour in selected_hours if hour < time_h[-1]]
    selected_indices = [int(np.argmin(np.abs(time_h - hour))) for hour in selected_hours]
    selected_indices.append(len(time_h) - 1)
    fig, axis = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    for index in selected_indices:
        label = (
            "drying endpoint"
            if index == len(time_h) - 1
            else f"{time_h[index]:g} h"
        )
        axis.plot(radius_cm, moisture[index], label=label)
    axis.axhline(result.threshold, color="#B91C1C", linestyle="--", linewidth=1.0)
    axis.set(
        xlabel="Radius r (cm)",
        ylabel="Moisture concentration (kg/kg)",
    )
    axis.grid(alpha=0.25)
    axis.legend(ncol=2, fontsize=8)
    fig.savefig(figure_dir / "problem3_radial_profiles.png", dpi=220)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.2, 4.6), constrained_layout=True)
    for column, label in ((0, "center"), (10, "mid-radius"), (20, "surface")):
        axis.plot(time_h, moisture[:, column], label=label)
    axis.axhline(result.threshold, color="#B91C1C", linestyle="--", linewidth=1.0)
    axis.set(
        xlabel="Time (h)",
        ylabel="Moisture concentration (kg/kg)",
    )
    axis.grid(alpha=0.25)
    axis.legend()
    fig.savefig(figure_dir / "problem3_time_histories.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = build_parser().parse_args()
    history = load_chamber_history_csv(args.input)
    horizon_s = args.horizon_hours * 3600.0
    solver_options = {"relative_tolerance": 2.0e-9, "max_step_s": 60.0}

    endpoint_results = []
    for radial_intervals in (
        args.coarse_radial_intervals,
        args.fine_radial_intervals,
    ):
        endpoint = solve_problem3(
            history,
            radial_intervals,
            horizon_s,
            output_interval_s=3600.0,
            **solver_options,
        )
        endpoint_results.append(endpoint)
        print(
            f"endpoint grid={radial_intervals}: "
            f"{endpoint.drying_time_s:.6f} s",
            flush=True,
        )
    extrapolated_time = richardson_extrapolate_drying_time(
        endpoint_results[0].drying_time_s,
        endpoint_results[1].drying_time_s,
        order=2,
    )

    common_end_time = extrapolated_time + 60.0
    extended_history = extend_chamber_history_to_plateau(history, common_end_time)
    regular_times = np.arange(60.0, common_end_time + 1.0e-9, 60.0)
    common_times = np.unique(
        np.append(
            regular_times,
            [extrapolated_time - 60.0, extrapolated_time + 60.0],
        )
    )
    requested_radii_cm = np.arange(0.0, 2.0 + 0.05, 0.1)
    sampled_solutions = []
    for radial_intervals in (
        args.coarse_radial_intervals,
        args.fine_radial_intervals,
    ):
        full_solution = solve_problem2(
            extended_history,
            radial_intervals,
            common_times,
            **solver_options,
        )
        sampled_solutions.append(sample_solution(full_solution, requested_radii_cm))
        del full_solution
        gc.collect()
        print(f"sampled field grid={radial_intervals}", flush=True)
    extrapolated_solution = richardson_extrapolate_solutions(
        sampled_solutions[0],
        sampled_solutions[1],
        order=2,
    )
    terminal_solution = truncate_solution_at_threshold(
        extrapolated_solution,
        threshold=0.15,
        output_interval_s=60.0,
    )
    terminal_moisture = terminal_solution.moisture_concentration[-1]
    maximum_index = int(np.argmax(terminal_moisture))
    result = Problem3Result(
        solution=terminal_solution,
        drying_time_s=float(terminal_solution.time_s[-1]),
        maximum_moisture=float(terminal_moisture[maximum_index]),
        maximum_radius_m=float(terminal_solution.radius_m[maximum_index]),
        plateau_temperature_c=float(extended_history.temperature_c[-1]),
        plateau_moisture_concentration=float(
            extended_history.moisture_concentration[-1]
        ),
        threshold=0.15,
        output_interval_s=60.0,
    )

    tail_half_hour = history.time_s >= history.time_s[-1] - 1800.0
    tail_hour = history.time_s >= history.time_s[-1] - 3600.0
    temperature_mean = float(np.mean(history.temperature_c[tail_hour]))
    moisture_mean = float(np.mean(history.moisture_concentration[tail_hour]))
    temperature_std = float(np.std(history.temperature_c[tail_hour], ddof=1))
    moisture_std = float(np.std(history.moisture_concentration[tail_hour], ddof=1))
    scenarios = {
        "last_hour_mean": (temperature_mean, moisture_mean),
        "last_half_hour_mean": (
            float(np.mean(history.temperature_c[tail_half_hour])),
            float(np.mean(history.moisture_concentration[tail_half_hour])),
        ),
        "last_observation": (
            float(history.temperature_c[-1]),
            float(history.moisture_concentration[-1]),
        ),
        "cool_wet_one_standard_deviation": (
            temperature_mean - temperature_std,
            moisture_mean + moisture_std,
        ),
        "warm_dry_one_standard_deviation": (
            temperature_mean + temperature_std,
            max(0.0, moisture_mean - moisture_std),
        ),
    }
    sensitivity: dict[str, dict[str, float]] = {}
    for name, (plateau_temperature, plateau_moisture) in scenarios.items():
        scenario = solve_problem3(
            history,
            args.sensitivity_radial_intervals,
            horizon_s,
            output_interval_s=3600.0,
            plateau_temperature_c=plateau_temperature,
            plateau_moisture_concentration=plateau_moisture,
            **solver_options,
        )
        sensitivity[name] = {
            "plateau_temperature_c": plateau_temperature,
            "plateau_moisture_concentration": plateau_moisture,
            "drying_time_s": scenario.drying_time_s,
            "drying_time_h": scenario.drying_time_s / 3600.0,
        }
        print(
            f"sensitivity {name}: {scenario.drying_time_s / 3600.0:.6f} h",
            flush=True,
        )

    payload = problem3_payload(result)
    payload["grid_endpoint_s"] = {
        str(args.coarse_radial_intervals): endpoint_results[0].drying_time_s,
        str(args.fine_radial_intervals): endpoint_results[1].drying_time_s,
        "richardson": result.drying_time_s,
    }
    payload["sensitivity"] = sensitivity
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
        encoding="utf-8",
    )
    save_figures(result, args.figure_dir)


if __name__ == "__main__":
    main()
