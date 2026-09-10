import math
from pathlib import Path

import numpy as np

from drying_model.problem1 import (
    ChamberHistory,
    diffusivity_q1,
    nodal_control_volumes,
    radial_flux_divergence,
    interpolate_history,
    load_chamber_history_csv,
    result_payload,
    sample_solution,
    solve_problem1,
)


def test_q1_diffusivity_uses_entire_minus_089_over_c_exponent():
    concentration = np.array([2.55, 1.0, 0.5])

    actual = diffusivity_q1(concentration)

    expected = 7.0e-9 * np.exp(-0.89 / concentration)
    np.testing.assert_allclose(actual, expected, rtol=1e-14, atol=0.0)


def test_chamber_history_is_linearly_interpolated_in_seconds():
    time = np.array([0.0, 60.0, 120.0])
    values = np.array([28.0, 30.0, 29.0])

    assert math.isclose(interpolate_history(30.0, time, values), 29.0)
    assert math.isclose(interpolate_history(90.0, time, values), 29.5)


def test_radial_flux_divergence_preserves_uniform_equilibrium():
    state = np.full(9, 2.55)

    rhs = radial_flux_divergence(
        state,
        coefficient=np.full_like(state, 4.0e-9),
        radius_m=0.02,
        exchange_coefficient=8.0e-7,
        ambient_value=2.55,
    )

    np.testing.assert_allclose(rhs, 0.0, atol=1e-18)


def test_closed_radial_domain_conserves_volume_weighted_state():
    state = np.array([2.6, 2.5, 2.4, 2.2, 2.0], dtype=float)
    radius_m = 0.02

    rhs = radial_flux_divergence(
        state,
        coefficient=np.full_like(state, 5.0e-9),
        radius_m=radius_m,
        exchange_coefficient=0.0,
        ambient_value=0.0,
    )
    volumes = nodal_control_volumes(state.size, radius_m)

    assert abs(float(np.dot(rhs, volumes))) < 1e-18


def test_surface_exchange_removes_exact_integrated_amount():
    state = np.full(5, 2.55)
    radius_m = 0.02
    transfer = 8.0e-7
    ambient = 0.02

    rhs = radial_flux_divergence(
        state,
        coefficient=np.full_like(state, 5.0e-9),
        radius_m=radius_m,
        exchange_coefficient=transfer,
        ambient_value=ambient,
    )
    volumes = nodal_control_volumes(state.size, radius_m)
    expected_per_length = -2.0 * math.pi * radius_m * transfer * (2.55 - ambient)

    assert math.isclose(float(np.dot(rhs, volumes)), expected_per_length, rel_tol=1e-14)


def test_solver_preserves_uniform_state_at_matching_ambient_conditions():
    history = ChamberHistory(
        time_s=np.array([0.0, 60.0]),
        temperature_c=np.array([28.0, 28.0]),
        moisture_concentration=np.array([2.55, 2.55]),
    )

    solution = solve_problem1(
        history,
        radial_intervals=8,
        output_times_s=np.array([0.0, 30.0, 60.0]),
        max_step_s=5.0,
    )

    np.testing.assert_allclose(solution.temperature_c, 28.0, atol=1e-11)
    np.testing.assert_allclose(solution.moisture_concentration, 2.55, atol=1e-11)


def test_solver_creates_expected_surface_to_center_gradients():
    history = ChamberHistory(
        time_s=np.array([0.0, 120.0]),
        temperature_c=np.array([50.0, 50.0]),
        moisture_concentration=np.array([0.02, 0.02]),
    )

    solution = solve_problem1(
        history,
        radial_intervals=12,
        output_times_s=np.array([0.0, 120.0]),
        max_step_s=2.0,
    )

    final_temperature = solution.temperature_c[-1]
    final_moisture = solution.moisture_concentration[-1]
    assert 28.0 <= final_temperature[0] < final_temperature[-1] < 50.0
    assert 0.02 < final_moisture[-1] < final_moisture[0] <= 2.55


def test_output_can_start_after_initial_time_without_resetting_initial_state():
    history = ChamberHistory(
        time_s=np.array([0.0, 10.0]),
        temperature_c=np.array([50.0, 50.0]),
        moisture_concentration=np.array([0.02, 0.02]),
    )

    with_initial = solve_problem1(
        history,
        radial_intervals=6,
        output_times_s=np.array([0.0, 1.0]),
        max_step_s=0.2,
    )
    without_initial = solve_problem1(
        history,
        radial_intervals=6,
        output_times_s=np.array([1.0]),
        max_step_s=0.2,
    )

    np.testing.assert_allclose(without_initial.temperature_c[0], with_initial.temperature_c[1])
    np.testing.assert_allclose(
        without_initial.moisture_concentration[0], with_initial.moisture_concentration[1]
    )


def test_load_chamber_history_csv_reads_chinese_headers(tmp_path: Path):
    source = tmp_path / "attachment1.csv"
    source.write_text(
        "时间,温度,水分浓度\n0,28,0.01963\n60,28.528,0.02002\n",
        encoding="utf-8-sig",
    )

    history = load_chamber_history_csv(source)

    np.testing.assert_array_equal(history.time_s, [0.0, 60.0])
    np.testing.assert_allclose(history.temperature_c, [28.0, 28.528])
    np.testing.assert_allclose(history.moisture_concentration, [0.01963, 0.02002])


def test_sample_solution_interpolates_to_requested_physical_radii():
    history = ChamberHistory(
        time_s=np.array([0.0, 10.0]),
        temperature_c=np.array([40.0, 40.0]),
        moisture_concentration=np.array([0.02, 0.02]),
    )
    solution = solve_problem1(
        history,
        radial_intervals=8,
        output_times_s=np.array([0.0, 10.0]),
        max_step_s=1.0,
    )

    sampled = sample_solution(solution, radius_cm=np.array([0.0, 0.5, 1.0, 1.5, 2.0]))

    assert sampled.temperature_c.shape == (2, 5)
    assert sampled.moisture_concentration.shape == (2, 5)
    np.testing.assert_allclose(sampled.radius_m, np.array([0.0, 0.005, 0.01, 0.015, 0.02]))


def test_result_payload_uses_seconds_centimetres_and_four_decimal_values():
    solution = solve_problem1(
        ChamberHistory(
            time_s=np.array([0.0, 1.0]),
            temperature_c=np.array([28.0, 28.0]),
            moisture_concentration=np.array([2.55, 2.55]),
        ),
        radial_intervals=2,
        output_times_s=np.array([1.0]),
        max_step_s=0.2,
    )

    payload = result_payload(solution)

    assert payload["time_s"] == [1]
    assert payload["radius_cm"] == [0.0, 1.0, 2.0]
    assert payload["temperature_c"] == [[28.0, 28.0, 28.0]]
    assert payload["moisture_concentration"] == [[2.55, 2.55, 2.55]]
