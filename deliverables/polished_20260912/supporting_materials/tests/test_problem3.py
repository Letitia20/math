import numpy as np
import drying_model.problem3 as problem3

from drying_model.problem1 import ChamberHistory, Problem1Parameters, Problem1Solution
from drying_model.problem2 import solve_problem2
from drying_model.problem3 import (
    Problem3Result,
    extend_chamber_history_to_plateau,
    problem3_payload,
    richardson_extrapolate_drying_time,
    solve_problem3,
    truncate_solution_at_threshold,
)


def test_plateau_extension_preserves_observations_and_uses_tail_mean():
    history = ChamberHistory(
        time_s=np.array([0.0, 60.0, 120.0, 180.0]),
        temperature_c=np.array([20.0, 30.0, 40.0, 50.0]),
        moisture_concentration=np.array([0.4, 0.3, 0.2, 0.1]),
    )

    extended = extend_chamber_history_to_plateau(
        history,
        end_time_s=600.0,
        averaging_window_s=120.0,
        transition_s=60.0,
    )

    np.testing.assert_array_equal(extended.time_s[:4], history.time_s)
    np.testing.assert_array_equal(extended.temperature_c[:4], history.temperature_c)
    np.testing.assert_array_equal(
        extended.moisture_concentration[:4],
        history.moisture_concentration,
    )
    np.testing.assert_array_equal(extended.time_s[-2:], np.array([240.0, 600.0]))
    np.testing.assert_allclose(extended.temperature_c[-2:], 40.0)
    np.testing.assert_allclose(extended.moisture_concentration[-2:], 0.2)


def test_problem3_solver_stops_on_full_field_threshold():
    history = ChamberHistory(
        time_s=np.array([0.0, 100.0]),
        temperature_c=np.array([50.0, 50.0]),
        moisture_concentration=np.array([0.05, 0.05]),
    )
    parameters = Problem1Parameters(radius_m=1.0e-4)

    result = solve_problem3(
        history,
        radial_intervals=8,
        horizon_s=1000.0,
        output_interval_s=60.0,
        threshold=0.2,
        averaging_window_s=100.0,
        parameters=parameters,
        relative_tolerance=1.0e-9,
        max_step_s=1.0,
    )

    assert result.drying_time_s == result.solution.time_s[-1]
    assert result.strict_completion_time_s is None
    assert result.drying_time_s % 60.0 != 0.0
    np.testing.assert_allclose(result.maximum_moisture, 0.2, atol=1.0e-9)
    assert result.maximum_radius_m in result.solution.radius_m
    assert np.all(
        np.max(result.solution.moisture_concentration[:-1], axis=1) > 0.2
    )


def test_problem3_payload_keeps_exact_endpoint_and_rounds_only_outputs():
    solution = Problem1Solution(
        time_s=np.array([60.0, 120.0, 135.567891, 180.0]),
        radius_m=np.array([0.0, 0.01, 0.02]),
        temperature_c=np.full((4, 3), 50.0),
        moisture_concentration=np.array(
            [
                [0.20006, 0.19004, 0.18003],
                [0.16006, 0.15504, 0.15103],
                [0.15, 0.14994, 0.14991],
                [0.1499, 0.1497, 0.1495],
            ]
        ),
    )
    result = Problem3Result(
        solution=solution,
        drying_time_s=135.567891,
        strict_completion_time_s=180.0,
        maximum_moisture=0.15,
        maximum_radius_m=0.0,
        plateau_temperature_c=50.0,
        plateau_moisture_concentration=0.05,
        threshold=0.15,
        output_interval_s=60.0,
    )

    payload = problem3_payload(result)

    assert payload["time_s"] == [60, 120, 135.567891, 180]
    assert payload["radius_cm"] == [0.0, 1.0, 2.0]
    assert payload["moisture_concentration"] == [
        [0.2001, 0.19, 0.18],
        [0.1601, 0.155, 0.151],
        [0.15, 0.1499, 0.1499],
        [0.1499, 0.1497, 0.1495],
    ]
    assert payload["drying_time_s"] == 135.567891
    assert payload["drying_time_h"] == round(135.567891 / 3600.0, 10)
    assert payload["strict_completion_time_s"] == 180
    assert payload["strict_completion_time_h"] == 0.05


def test_problem3_payload_distinguishes_nearby_critical_and_regular_times():
    solution = Problem1Solution(
        time_s=np.array([206880.0, 206939.5, 206940.0]),
        radius_m=np.array([0.0, 0.02]),
        temperature_c=np.full((3, 2), 50.0),
        moisture_concentration=np.array(
            [[0.1501, 0.05], [0.15, 0.05], [0.1499, 0.05]]
        ),
    )
    result = Problem3Result(
        solution=solution,
        drying_time_s=206939.5,
        strict_completion_time_s=206940.0,
        maximum_moisture=0.15,
        maximum_radius_m=0.0,
        plateau_temperature_c=50.0,
        plateau_moisture_concentration=0.05,
        threshold=0.15,
        output_interval_s=60.0,
    )

    payload = problem3_payload(result)

    assert payload["time_s"][-2:] == [206939.5, 206940]


def test_richardson_extrapolates_second_order_drying_time_error():
    assert richardson_extrapolate_drying_time(104.0, 101.0, order=2) == 100.0

    with np.testing.assert_raises(ValueError):
        richardson_extrapolate_drying_time(104.0, 101.0, order=0)


def test_radial_moisture_validation_rejects_an_outward_increase():
    valid = Problem1Solution(
        time_s=np.array([60.0, 120.0]),
        radius_m=np.array([0.0, 0.01, 0.02]),
        temperature_c=np.full((2, 3), 50.0),
        moisture_concentration=np.array(
            [[0.30, 0.20, 0.10], [0.15, 0.14, 0.05]]
        ),
    )
    invalid = Problem1Solution(
        time_s=valid.time_s,
        radius_m=valid.radius_m,
        temperature_c=valid.temperature_c,
        moisture_concentration=np.array(
            [[0.30, 0.10, 0.20], [0.15, 0.14, 0.05]]
        ),
    )

    assert problem3.validate_radially_nonincreasing_moisture(valid) <= 0.0
    with np.testing.assert_raises_regex(ValueError, "outward increase"):
        problem3.validate_radially_nonincreasing_moisture(invalid)


def test_full_solution_moisture_change_matches_integrated_surface_flux():
    history = ChamberHistory(
        time_s=np.array([0.0, 600.0]),
        temperature_c=np.array([50.0, 50.0]),
        moisture_concentration=np.array([0.05, 0.05]),
    )
    parameters = Problem1Parameters()
    solution = solve_problem2(
        history,
        radial_intervals=16,
        output_times_s=np.arange(0.0, 601.0),
        parameters=parameters,
        relative_tolerance=1.0e-9,
        max_step_s=0.5,
    )

    diagnostics = problem3.moisture_balance_diagnostics(
        solution,
        history,
        parameters=parameters,
    )

    assert diagnostics["relative_max_residual"] < 1.0e-6
    assert abs(diagnostics["final_residual"]) < 1.0e-8


def test_threshold_interpolation_keeps_crossing_and_first_strict_regular_row():
    times = np.array([60.0, 120.0, 130.0, 140.0, 180.0])
    moisture = np.array(
        [
            [0.20, 0.18],
            [0.16, 0.14],
            [0.155, 0.135],
            [0.145, 0.13],
            [0.10, 0.09],
        ]
    )
    solution = Problem1Solution(
        time_s=times,
        radius_m=np.array([0.0, 0.02]),
        temperature_c=np.column_stack([times, times + 1.0]),
        moisture_concentration=moisture,
    )

    truncated = truncate_solution_at_threshold(
        solution,
        threshold=0.15,
        output_interval_s=60.0,
    )

    np.testing.assert_allclose(
        truncated.time_s,
        np.array([60.0, 120.0, 135.0, 180.0]),
    )
    np.testing.assert_allclose(truncated.moisture_concentration[-2, 0], 0.15)
    assert np.max(truncated.moisture_concentration[-1]) < 0.15
    np.testing.assert_allclose(truncated.temperature_c[-2], np.array([135.0, 136.0]))
