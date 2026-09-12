"""Regressions found while reviewing the four-question delivery."""

import importlib.util
from pathlib import Path
import sys

import numpy as np
import pytest

from problem1 import ChamberHistory, Problem1Solution, result_payload
from problem3 import Problem3Result, problem3_payload, truncate_solution_at_threshold


def test_fractional_second_is_not_silently_rounded_at_large_times():
    solution = Problem1Solution(np.array([10800.01]), np.array([0.0, 0.02]),
                                np.full((1, 2), 50.0), np.full((1, 2), 0.2))
    with pytest.raises(ValueError, match="whole seconds"):
        result_payload(solution)


@pytest.mark.parametrize("column", [0, 1, 2])
def test_history_rejects_nonfinite_observation_columns(column):
    columns = [np.array([0.0, 60.0]), np.array([28.0, 30.0]), np.array([0.02, 0.03])]
    columns[column][1] = np.nan
    with pytest.raises(ValueError, match="finite"):
        ChamberHistory(*columns)


def test_threshold_interpolation_when_wettest_node_changes():
    solution = Problem1Solution(
        np.array([60.0, 120.0, 180.0]), np.array([0.0, 0.02]), np.full((3, 2), 50.0),
        np.array([[0.20, 0.19], [0.14, 0.149], [0.13, 0.14]]),
    )
    actual = truncate_solution_at_threshold(solution, threshold=0.15, output_interval_s=60.0)
    expected = 60.0 + 60.0 * (0.19 - 0.15) / (0.19 - 0.149)
    assert actual.time_s[-2] == pytest.approx(expected, abs=1e-10)
    assert np.max(actual.moisture_concentration[-2]) == pytest.approx(0.15, abs=1e-14)
    assert actual.time_s[-1] == 120.0


def test_critical_moisture_validation_uses_absolute_tolerance_only():
    solution = Problem1Solution(
        np.array([60.0, 90.0, 120.0]), np.array([0.0, 0.02]), np.full((3, 2), 50.0),
        np.array([[0.16, 0.10], [0.150001, 0.09], [0.14999, 0.08]]),
    )
    result = Problem3Result(solution, 90.0, 120.0, 0.150001, 0.0, 50.0, 0.05, 0.15, 60.0)
    with pytest.raises(ValueError, match="unrounded threshold"):
        problem3_payload(result)


@pytest.mark.parametrize("question", [1, 3])
def test_runner_rejects_invalid_richardson_grid_ratio(question, monkeypatch):
    path = Path(__file__).resolve().parent / f"run_problem{question}.py"
    spec = importlib.util.spec_from_file_location(f"run_problem{question}", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(sys, "argv", [str(path), "--coarse-radial-intervals", "8", "--fine-radial-intervals", "24"])
    with pytest.raises(ValueError, match="twice"):
        module.main()
