import numpy as np
import pytest

from drying_model.problem1 import ChamberHistory, Problem1Parameters
from drying_model.problem4 import (
    RadiusHistory,
    density_q4,
    density_shrinkage_compatibility,
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


def test_constant_radius_reduces_to_fixed_cylinder_equations():
    from drying_model.problem2 import solve_variable_property_fixed_cylinder

    chamber = ChamberHistory(np.array([0.0, 600.0]), np.array([45.0, 50.0]), np.array([0.05, 0.03]))
    radius = RadiusHistory(np.array([0.0, 600.0]), np.array([0.02, 0.02]))
    times = np.array([0.0, 100.0, 600.0])
    moving = solve_problem4(chamber, radius, 24, times, relative_tolerance=1e-9, max_step_s=1.0)
    fixed = solve_variable_property_fixed_cylinder(
        chamber, 24, times, density_law=density_q4, heat_capacity_law=heat_capacity_q4,
        conductivity_law=thermal_conductivity_q4, diffusivity_law=diffusivity_q4,
        relative_tolerance=1e-9, max_step_s=1.0,
    )
    np.testing.assert_allclose(moving.temperature_c, fixed.temperature_c, rtol=1e-8, atol=1e-7)
    np.testing.assert_allclose(moving.moisture_concentration, fixed.moisture_concentration, rtol=1e-8, atol=1e-9)


@pytest.mark.parametrize("method", ["linear", "pchip"])
def test_solver_holds_radius_tail_but_requires_explicit_chamber_extension(method):
    chamber = ChamberHistory(np.array([0.0, 100.0]), np.array([50.0, 50.0]), np.array([0.02, 0.02]))
    radius = RadiusHistory(np.array([0.0, 20.0, 40.0]), np.array([0.02, 0.019, 0.018]), interpolation=method)
    result = solve_problem4(chamber, radius, 8, [40.0, 100.0])
    assert result.time_s[-1] == 100.0
    assert np.isfinite(result.moisture_concentration).all()
    np.testing.assert_allclose(interpolate_radius([40.0, 100.0], radius), 0.018)
    with pytest.raises(ValueError, match="chamber history"):
        solve_problem4(chamber, radius, 8, [101.0])


def test_pchip_preserves_radius_observations_and_monotonicity():
    radius = RadiusHistory(np.array([0.0, 10.0, 20.0, 30.0]), np.array([0.02, 0.016, 0.015, 0.015]), interpolation="pchip")
    np.testing.assert_allclose(interpolate_radius(radius.time_s, radius), radius.radius_m, rtol=0, atol=1e-15)
    samples = interpolate_radius(np.linspace(0.0, 60.0, 601), radius)
    assert np.max(np.diff(samples)) <= 1e-14
    assert np.min(samples) >= 0.015 - 1e-14


def test_density_compatibility_detects_mass_loss_despite_uniform_moisture():
    solution = Problem4Solution(
        time_s=np.array([0.0, 100.0]), xi=np.linspace(0.0, 1.0, 9),
        temperature_c=np.full((2, 9), 28.0), moisture_concentration=np.full((2, 9), 2.55),
    )
    radius = RadiusHistory(solution.time_s, np.array([0.02, 0.01]))
    result = density_shrinkage_compatibility(solution, radius)
    np.testing.assert_allclose(result["mass_ratio_to_initial"], [1.0, 0.25])
    expected_initial = np.pi * 0.02**2 * (760.0 + 90.0 * 2.55) / 3.55
    np.testing.assert_allclose(result["implied_dry_mass_per_length_kg_m"][0], expected_initial)
