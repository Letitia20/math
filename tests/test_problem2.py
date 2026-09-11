import numpy as np

from drying_model.problem1 import (
    ChamberHistory,
    Problem1Parameters,
    diffusivity_q1,
    solve_problem1,
)
from drying_model.problem2 import (
    density_q2,
    diffusivity_q2,
    heat_capacity_q2,
    solve_problem2,
    solve_problem2_sampled_in_chunks,
    solve_variable_property_fixed_cylinder,
    thermal_conductivity_q2,
)


def test_q2_material_laws_match_appendix_3_exactly():
    concentration = np.array([0.15, 1.0, 2.55])
    temperature_c = np.array([28.0, 40.0, 50.0])

    np.testing.assert_allclose(density_q2(concentration), 650.0 + 128.0 * concentration)
    np.testing.assert_allclose(
        heat_capacity_q2(concentration),
        1450.0 + 2736.0 * concentration / (concentration + 1.0),
    )
    np.testing.assert_allclose(
        thermal_conductivity_q2(concentration),
        0.21 + 0.38 * concentration / (concentration + 1.0),
    )
    np.testing.assert_allclose(
        diffusivity_q2(concentration, temperature_c),
        2.4e-3
        * np.exp(-0.45 / concentration)
        * np.exp(-3850.0 / (temperature_c + 273.15)),
    )


def test_q2_material_laws_reject_nonphysical_inputs():
    for law in (density_q2, heat_capacity_q2, thermal_conductivity_q2):
        with np.testing.assert_raises(ValueError):
            law(np.array([-1.0]))
    with np.testing.assert_raises(ValueError):
        diffusivity_q2(np.array([0.0]), np.array([28.0]))
    with np.testing.assert_raises(ValueError):
        diffusivity_q2(np.array([1.0]), np.array([-273.15]))


def test_q2_uniform_equilibrium_remains_constant():
    history = ChamberHistory(
        time_s=np.array([0.0, 20.0]),
        temperature_c=np.array([28.0, 28.0]),
        moisture_concentration=np.array([2.55, 2.55]),
    )

    solution = solve_problem2(
        history,
        radial_intervals=8,
        output_times_s=np.array([0.0, 10.0, 20.0]),
        max_step_s=1.0,
    )

    np.testing.assert_allclose(solution.temperature_c, 28.0, atol=1.0e-11)
    np.testing.assert_allclose(solution.moisture_concentration, 2.55, atol=1.0e-11)


def test_generic_coupled_solver_degenerates_to_problem1_laws():
    history = ChamberHistory(
        time_s=np.array([0.0, 30.0]),
        temperature_c=np.array([45.0, 45.0]),
        moisture_concentration=np.array([0.02, 0.02]),
    )
    parameters = Problem1Parameters()
    output_times = np.array([0.0, 10.0, 30.0])
    expected = solve_problem1(
        history,
        radial_intervals=10,
        output_times_s=output_times,
        parameters=parameters,
        relative_tolerance=1.0e-9,
        max_step_s=0.5,
    )

    actual = solve_variable_property_fixed_cylinder(
        history,
        radial_intervals=10,
        output_times_s=output_times,
        density_law=lambda concentration: np.full_like(concentration, parameters.density_kg_m3),
        heat_capacity_law=lambda concentration: np.full_like(
            concentration, parameters.heat_capacity_j_kg_k
        ),
        conductivity_law=lambda concentration: np.full_like(
            concentration, parameters.thermal_conductivity_w_m_k
        ),
        diffusivity_law=lambda concentration, temperature_c: diffusivity_q1(concentration),
        parameters=parameters,
        relative_tolerance=1.0e-9,
        max_step_s=0.5,
    )

    np.testing.assert_allclose(actual.temperature_c, expected.temperature_c, atol=2.0e-8)
    np.testing.assert_allclose(
        actual.moisture_concentration,
        expected.moisture_concentration,
        atol=2.0e-8,
    )


def test_q2_hot_dry_environment_creates_expected_surface_gradients():
    history = ChamberHistory(
        time_s=np.array([0.0, 120.0]),
        temperature_c=np.array([50.0, 50.0]),
        moisture_concentration=np.array([0.02, 0.02]),
    )

    solution = solve_problem2(
        history,
        radial_intervals=16,
        output_times_s=np.array([120.0]),
        max_step_s=1.0,
    )

    assert solution.temperature_c[0, -1] > solution.temperature_c[0, 0]
    assert solution.moisture_concentration[0, -1] < solution.moisture_concentration[0, 0]
    assert np.all(solution.moisture_concentration > 0.0)
def test_q2_solver_can_continue_from_a_previous_chunk_without_resetting_state():
    history = ChamberHistory(
        time_s=np.array([0.0, 120.0]),
        temperature_c=np.array([50.0, 50.0]),
        moisture_concentration=np.array([0.05, 0.05]),
    )
    monolithic = solve_problem2(
        history,
        radial_intervals=16,
        output_times_s=np.array([60.0, 120.0]),
        relative_tolerance=1.0e-9,
        max_step_s=0.5,
    )
    first_chunk = solve_problem2(
        history,
        radial_intervals=16,
        output_times_s=np.array([60.0]),
        relative_tolerance=1.0e-9,
        max_step_s=0.5,
    )

    second_chunk = solve_problem2(
        history,
        radial_intervals=16,
        output_times_s=np.array([120.0]),
        initial_time_s=60.0,
        initial_temperature_c=first_chunk.temperature_c[-1],
        initial_moisture_concentration=first_chunk.moisture_concentration[-1],
        relative_tolerance=1.0e-9,
        max_step_s=0.5,
    )

    np.testing.assert_allclose(
        second_chunk.temperature_c[-1],
        monolithic.temperature_c[-1],
        atol=2.0e-7,
    )
    np.testing.assert_allclose(
        second_chunk.moisture_concentration[-1],
        monolithic.moisture_concentration[-1],
        atol=2.0e-9,
    )


def test_chunked_q2_sampling_matches_monolithic_solution():
    history = ChamberHistory(
        time_s=np.array([0.0, 120.0]),
        temperature_c=np.array([50.0, 50.0]),
        moisture_concentration=np.array([0.05, 0.05]),
    )
    output_times = np.arange(1.0, 121.0)
    sample_radius_cm = np.array([0.0, 1.0, 2.0])
    monolithic = solve_problem2(
        history,
        radial_intervals=16,
        output_times_s=output_times,
        relative_tolerance=1.0e-9,
        max_step_s=0.5,
    )

    actual = solve_problem2_sampled_in_chunks(
        history,
        radial_intervals=16,
        output_times_s=output_times,
        sample_radius_cm=sample_radius_cm,
        chunk_duration_s=30.0,
        relative_tolerance=1.0e-9,
        max_step_s=0.5,
    )

    np.testing.assert_array_equal(actual.time_s, output_times)
    np.testing.assert_allclose(actual.radius_m, sample_radius_cm * 0.01)
    np.testing.assert_allclose(
        actual.temperature_c,
        monolithic.temperature_c[:, [0, 8, 16]],
        atol=3.0e-7,
    )
    np.testing.assert_allclose(
        actual.moisture_concentration,
        monolithic.moisture_concentration[:, [0, 8, 16]],
        atol=3.0e-9,
    )
