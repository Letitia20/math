import numpy as np

from drying_model.problem1 import ChamberHistory, Problem1Parameters
from drying_model.problem4 import (
    RadiusHistory,
    density_q4,
    diffusivity_q4,
    heat_capacity_q4,
    interpolate_radius,
    moisture_balance_diagnostics,
    sample_moving_solution,
    solve_problem4,
    thermal_conductivity_q4,
    Problem4Solution,
)


def test_q4_material_laws_match_appendix_4():
    c = np.array([0.15, 1.0, 2.55])
    t = np.array([28.0, 40.0, 50.0])
    np.testing.assert_allclose(density_q4(c), 760.0 + 90.0 * c)
    np.testing.assert_allclose(heat_capacity_q4(c), 1850.0 + 2150.0 * c / (c + 1.0))
    np.testing.assert_allclose(thermal_conductivity_q4(c), 0.12 + 0.20 * c / (c + 1.0))
    np.testing.assert_allclose(
        diffusivity_q4(c, t), 4.2e-4 * np.exp(-0.30 / c) * np.exp(-3850.0 / (t + 273.15))
    )


def test_radius_history_interpolation_is_monotone_and_holds_tail():
    history = RadiusHistory(np.array([0.0, 10.0, 20.0]), np.array([0.02, 0.018, 0.018]))
    values = interpolate_radius(np.array([0.0, 5.0, 30.0]), history)
    np.testing.assert_allclose(values, [0.02, 0.019, 0.018])


def test_uniform_equilibrium_remains_constant_with_shrinking_radius():
    chamber = ChamberHistory(np.array([0.0, 100.0]), np.array([28.0, 28.0]), np.array([2.55, 2.55]))
    radius = RadiusHistory(np.array([0.0, 50.0, 100.0]), np.array([0.02, 0.018, 0.016]))
    solution = solve_problem4(chamber, radius, 8, np.array([0.0, 50.0, 100.0]), max_step_s=2.0)
    np.testing.assert_allclose(solution.temperature_c, 28.0, atol=2.0e-10)
    np.testing.assert_allclose(solution.moisture_concentration, 2.55, atol=2.0e-10)


def test_hot_dry_boundary_creates_expected_gradients():
    chamber = ChamberHistory(np.array([0.0, 120.0]), np.array([50.0, 50.0]), np.array([0.02, 0.02]))
    radius = RadiusHistory(np.array([0.0, 120.0]), np.array([0.02, 0.019]))
    solution = solve_problem4(chamber, radius, 16, np.array([120.0]), max_step_s=1.0)
    assert solution.temperature_c[0, -1] > solution.temperature_c[0, 0]
    assert solution.moisture_concentration[0, -1] < solution.moisture_concentration[0, 0]
    assert np.all(solution.moisture_concentration > 0.0)


def test_sampling_uses_nan_for_fixed_radii_outside_shrunken_body():
    solution = Problem4Solution(
        time_s=np.array([0.0, 1.0]), xi=np.array([0.0, 0.5, 1.0]),
        temperature_c=np.full((2, 3), 30.0),
        moisture_concentration=np.array([[1.0, 0.8, 0.6], [1.0, 0.8, 0.6]]),
    )
    radius = RadiusHistory(np.array([0.0, 1.0]), np.array([0.02, 0.01]))
    sampled = sample_moving_solution(solution, radius, np.array([0.0, 1.0, 1.5, 2.0]))
    assert np.isfinite(sampled[0, :]).all()
    assert np.isfinite(sampled[1, :2]).all()
    assert np.isnan(sampled[1, 2:]).all()


def test_threshold_event_reports_exact_crossing():
    chamber = ChamberHistory(np.array([0.0, 200.0]), np.array([50.0, 50.0]), np.array([0.02, 0.02]))
    radius = RadiusHistory(np.array([0.0, 200.0]), np.array([0.0001, 0.0001]))
    result = solve_problem4(
        chamber, radius, 8, np.arange(20.0, 201.0, 20.0), maximum_moisture_threshold=0.2,
        parameters=Problem1Parameters(radius_m=0.0001), max_step_s=0.5,
    )
    assert result.time_s[-1] % 20.0 != 0.0
    np.testing.assert_allclose(np.max(result.moisture_concentration[-1]), 0.2, atol=1.0e-8)


def test_moving_domain_moisture_change_matches_surface_flux():
    chamber = ChamberHistory(np.array([0.0, 600.0]), np.array([50.0, 50.0]), np.array([0.05, 0.05]))
    radius = RadiusHistory(np.array([0.0, 600.0]), np.array([0.02, 0.018]))
    solution = solve_problem4(
        chamber, radius, 24, np.arange(0.0, 601.0), relative_tolerance=1.0e-9, max_step_s=0.5
    )
    diagnostics = moisture_balance_diagnostics(solution, chamber, radius)
    assert diagnostics["relative_max_residual"] < 2.0e-6
