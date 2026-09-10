"""Solve Problem 2 and export workbook-ready values plus diagnostic figures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from drying_model.problem1 import load_chamber_history_csv, result_payload, sample_solution
from drying_model.problem2 import diffusivity_q2, solve_problem2


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=Path("data/raw/attachment1.csv"))
    parser.add_argument("--output", type=Path, default=Path("tmp/problem2_result.json"))
    parser.add_argument("--figure-dir", type=Path, default=Path("reports/figures"))
    parser.add_argument("--radial-intervals", type=int, default=640)
    return parser


def save_figures(payload: dict[str, list], figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    time = np.asarray(payload["time_s"], dtype=float)
    radius = np.asarray(payload["radius_cm"], dtype=float)
    temperature = np.asarray(payload["temperature_c"], dtype=float)
    moisture = np.asarray(payload["moisture_concentration"], dtype=float)

    profile_times = [1800, 3600, 5400, 7200, 9000, 10800]
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(profile_times)))
    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.2), constrained_layout=True)
    for selected_time, color in zip(profile_times, colors, strict=True):
        row = selected_time - 1
        label = f"{selected_time / 3600:g} h"
        axes[0].plot(radius, temperature[row], color=color, label=label)
        axes[1].plot(radius, moisture[row], color=color, label=label)
    axes[0].set(xlabel="Radius r (cm)", ylabel="Temperature (deg C)")
    axes[1].set(xlabel="Radius r (cm)", ylabel="Moisture concentration (kg/kg)")
    for axis in axes:
        axis.grid(alpha=0.25)
    axes[1].legend(ncol=2, fontsize=8)
    fig.savefig(figure_dir / "problem2_radial_profiles.png", dpi=220)
    plt.close(fig)

    locations = [(0, "center"), (10, "mid-radius"), (20, "surface")]
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.2), constrained_layout=True)
    for column, label in locations:
        axes[0].plot(time / 3600.0, temperature[:, column], label=label)
        axes[1].plot(time / 3600.0, moisture[:, column], label=label)
        axes[2].plot(
            time / 3600.0,
            diffusivity_q2(moisture[:, column], temperature[:, column]),
            label=label,
        )
    axes[0].set(xlabel="Time (h)", ylabel="Temperature (deg C)")
    axes[1].set(xlabel="Time (h)", ylabel="Moisture concentration (kg/kg)")
    axes[2].set(xlabel="Time (h)", ylabel="Diffusivity (m2/s)")
    for axis in axes:
        axis.grid(alpha=0.25)
        axis.legend()
    fig.savefig(figure_dir / "problem2_time_histories.png", dpi=220)
    plt.close(fig)


def main() -> None:
    args = build_parser().parse_args()
    history = load_chamber_history_csv(args.input)
    solution = solve_problem2(
        history,
        args.radial_intervals,
        np.arange(1.0, 10801.0),
        relative_tolerance=2.0e-10,
        max_step_s=2.0,
    )
    sampled = sample_solution(solution, np.arange(0.0, 2.0 + 0.05, 0.1))
    payload = result_payload(sampled)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, ensure_ascii=False, separators=(",", ":")),
        encoding="utf-8",
    )
    save_figures(payload, args.figure_dir)


if __name__ == "__main__":
    main()
