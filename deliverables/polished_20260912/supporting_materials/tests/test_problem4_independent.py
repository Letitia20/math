"""Independent finite-volume checks against analytic geometry and balances."""

from dataclasses import replace

import numpy as np

from drying_model.problem1 import ChamberHistory, Problem1Parameters
from drying_model.problem4 import RadiusHistory
from drying_model.problem4_independent import (
    cell_flux_divergence,
    reconstruct_centre,
    solve_independent_problem4,
    surface_from_series_resistance,
)


def test_quadratic_cylindrical_laplacian_and_axis_reconstruction():
    n = 40
    radius = 0.017
    xi = (np.arange(n) + 0.5) / n
    values = 2.0 + 0.3 * xi**2
    diffusion = 2.0e-9
    divergence, _ = cell_flux_divergence(values, np.full(n, diffusion), radius, 0.0, 0.0)
    # Exclude the imposed zero-flux outer face, which this polynomial does not satisfy.
    np.testing.assert_allclose(divergence[:-1], 4 * 0.3 * diffusion / radius**2, rtol=2e-11)
    np.testing.assert_allclose(reconstruct_centre(values), 2.0, atol=1e-15)


def test_arbitrary_field_has_exact_discrete_surface_flux_balance():
    n = 35
    radius = 0.016
    values = 0.2 + np.exp(-np.linspace(0.0, 2.0, n))
    diffusion = np.linspace(1e-9, 2e-9, n)
    exchange = 8e-7
    ambient = 0.05
    divergence, surface = cell_flux_divergence(values, diffusion, radius, exchange, ambient)
    faces = np.arange(n + 1) / n
    physical_volumes_per_two_pi = 0.5 * radius**2 * np.diff(faces**2)
    np.testing.assert_allclose(
        np.dot(divergence, physical_volumes_per_two_pi),
        -radius * exchange * (surface - ambient), rtol=1e-13,
    )
    # The reported surface satisfies both pieces of the series resistance.
    inner_flux = diffusion[-1] * (values[-1] - surface) / (radius / (2 * n))
    np.testing.assert_allclose(inner_flux, exchange * (surface - ambient), rtol=1e-13)


def test_shrinking_domain_preserves_uniform_equilibrium():
    times = np.array([0.0, 100.0, 200.0])
    chamber = ChamberHistory(times, np.full(3, 28.0), np.full(3, 2.55))
    radius = RadiusHistory(times, np.array([0.02, 0.015, 0.012]))
    result = solve_independent_problem4(chamber, radius, 12, times, threshold=None)
    np.testing.assert_allclose(result.moisture, 2.55, atol=1e-12)
    np.testing.assert_allclose(result.temperature_c, 28.0, atol=1e-12)


def test_zero_exchange_keeps_initial_fields_during_shrinkage():
    times = np.array([0.0, 100.0, 200.0])
    chamber = ChamberHistory(times, np.full(3, 70.0), np.full(3, 0.05))
    radius = RadiusHistory(times, np.array([0.02, 0.015, 0.012]))
    parameters = replace(Problem1Parameters(), heat_transfer_coefficient_w_m2_k=0.0,
                         mass_transfer_coefficient_m_s=0.0)
    result = solve_independent_problem4(chamber, radius, 12, times, threshold=None, parameters=parameters)
    np.testing.assert_allclose(result.moisture, 2.55, atol=1e-12)
    np.testing.assert_allclose(result.temperature_c, 28.0, atol=1e-12)
