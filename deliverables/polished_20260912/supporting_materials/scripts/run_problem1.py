"""Solve Problem 1 and export workbook-ready values plus diagnostic figures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from drying_model.problem1 import (
    Problem1Parameters,
    load_chamber_history_csv,
    richardson_extrapolate_solutions,
    result_payload,
    sample_solution,
    solve_problem1,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/raw/attachment1.csv"))
    parser.add_argument("--output", type=Path, default=Path("tmp/problem1_result.json"))
    parser.add_argument("--figure-dir", type=Path, default=Path("reports/figures"))
    parser.add_argument("--coarse-radial-intervals", type=int, default=5120)
    parser.add_argument("--fine-radial-intervals", type=int, default=10240)
    return parser


def save_figures(payload: dict[str, list], figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    time = np.asarray(payload["time_s"], dtype=float)
    radius = np.asarray(payload["radius_cm"], dtype=float)
    temperature = np.asarray(payload["temperature_c"], dtype=float)
    moisture = np.asarray(payload["moisture_concentration"], dtype=float)

    profile_times = [100, 300, 600, 900, 1200, 1500, 1800]
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(profile_times)))
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)
    for selected_time, color in zip(profile_times, colors, strict=True):
        row = selected_time - 1
        axes[0].plot(radius, temperature[row], color=color, label=f"{selected_time} s")
        axes[1].plot(radius, moisture[row], color=color, label=f"{selected_time} s")
    axes[0].set(xlabel="Radius r (cm)", ylabel="Temperature (deg C)")
    axes[1].set(xlabel="Radius r (cm)", ylabel="Moisture concentration (kg/kg)")
    for axis in axes:
        axis.grid(alpha=0.25)
    axes[1].legend(ncol=2, fontsize=8)
    fig.savefig(figure_dir / "problem1_radial_profiles.png", dpi=220)
    plt.close(fig)

    locations = [(0, "center"), (10, "mid-radius"), (20, "surface")]
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)
    for column, label in locations:
        axes[0].plot(time, temperature[:, column], label=label)
        axes[1].plot(time, moisture[:, column], label=label)
    axes[0].set(xlabel="Time (s)", ylabel="Temperature (deg C)")
    axes[1].set(xlabel="Time (s)", ylabel="Moisture concentration (kg/kg)")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend()
    fig.savefig(figure_dir / "problem1_time_histories.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = build_parser().parse_args()
    if args.fine_radial_intervals != 2 * args.coarse_radial_intervals:
        raise ValueError("Fine radial grid must have twice as many intervals as the coarse grid")
    history = load_chamber_history_csv(args.input)
    output_times = np.arange(1.0, 1801.0)
    requested_radii_cm = np.arange(0.0, 2.0 + 0.05, 0.1)
    sampled_solutions = []
    for radial_intervals in (
        args.coarse_radial_intervals,
        args.fine_radial_intervals,
    ):
        solution = solve_problem1(
            history,
            radial_intervals,
            output_times,
            parameters=Problem1Parameters(),
            relative_tolerance=5.0e-10,
            max_step_s=1.0,
        )
        sampled_solutions.append(sample_solution(solution, requested_radii_cm))
    extrapolated = richardson_extrapolate_solutions(
        sampled_solutions[0],
        sampled_solutions[1],
        order=2,
    )
    payload = result_payload(extrapolated)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
        encoding="utf-8",
    )
    save_figures(payload, args.figure_dir)


if __name__ == "__main__":
    main()
