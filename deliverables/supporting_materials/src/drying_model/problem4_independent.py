"""Independent cell-centred finite-volume check of the Problem 4 model.

This implements the same continuum assumptions with different spatial unknowns
and boundary closure. It does not call the production solver, its material-law
functions, or its flux-divergence routine. Unknowns approximate point values at
the midpoints of N complete annular cells; neither the axis nor the surface is
a state node. Centre values are reconstructed by even quadratic extrapolation.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp
from scipy.sparse import bmat, diags

from drying_model.problem1 import ChamberHistory, Problem1Parameters
from drying_model.problem4 import RadiusHistory


@dataclass(frozen=True)
class IndependentSolution:
    time_s: NDArray[np.float64]
    xi_centres: NDArray[np.float64]
    temperature_c: NDArray[np.float64]
    moisture: NDArray[np.float64]
    surface_temperature_c: NDArray[np.float64]
    surface_moisture: NDArray[np.float64]
    centre_moisture: NDArray[np.float64]
    drying_time_s: float | None
    endpoint_cell_maximum: float | None
    endpoint_reconstructed_maximum: float | None


def appendix4_coefficients(moisture: ArrayLike, temperature_c: ArrayLike):
    """Transcribe Appendix 4 separately from the production implementation."""
    concentration = np.asarray(moisture, dtype=float)
    kelvin = np.asarray(temperature_c, dtype=float) + 273.15
    if np.any(concentration <= 0.0) or np.any(kelvin <= 0.0):
        raise ValueError("Positive concentration and absolute temperature required")
    rho_cp = (760.0 + 90.0 * concentration) * (
        1850.0 + 2150.0 * concentration / (concentration + 1.0)
    )
    conductivity = 0.12 + 0.20 * concentration / (concentration + 1.0)
    diffusion = 4.2e-4 * np.exp(-0.30 / concentration - 3850.0 / kelvin)
    return rho_cp, conductivity, diffusion


def reconstruct_centre(values: ArrayLike) -> NDArray[np.float64]:
    """Even quadratic continuation from xi=dx/2 and 3dx/2 to xi=0.

    These are point-centred unknowns, not volume averages. For a smooth even
    exact profile this formula's reconstruction error is O(dx**4); the solved
    point values still carry the spatial scheme's discretization error.
    """
    field = np.asarray(values, dtype=float)
    if field.shape[-1] < 2:
        raise ValueError("At least two cells are needed for axis reconstruction")
    return (9.0 * field[..., 0] - field[..., 1]) / 8.0


def surface_from_series_resistance(
    last_value: ArrayLike,
    last_coefficient: ArrayLike,
    radius_m: ArrayLike,
    cell_count: int,
    exchange: float,
    ambient: ArrayLike,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return actual-surface value and outward flux using a half-cell resistor.

    J=(u_last-u_inf)/(R*dx/(2*a_last)+1/h), u_surface=u_inf+J/h.
    The coefficient is frozen at the last cell centre. That local closure is
    a truncation approximation, checked by spatial convergence, not an exact
    nonlinear surface solve.
    """
    last = np.asarray(last_value, dtype=float)
    coefficient = np.asarray(last_coefficient, dtype=float)
    ambient = np.asarray(ambient, dtype=float)
    if cell_count < 2 or exchange < 0.0 or np.any(coefficient <= 0.0):
        raise ValueError("Invalid boundary conductance or cell count")
    if exchange == 0.0:
        return last.copy(), np.zeros_like(last)
    resistance_ratio = exchange * np.asarray(radius_m) / (2.0 * cell_count * coefficient)
    surface = (last + resistance_ratio * ambient) / (1.0 + resistance_ratio)
    return surface, exchange * (surface - ambient)


def cell_flux_divergence(
    values: ArrayLike,
    coefficient: ArrayLike,
    radius_m: float,
    exchange: float,
    ambient: float,
) -> tuple[NDArray[np.float64], float]:
    """Annular flux balance in xi with arithmetic internal-face coefficients."""
    field = np.asarray(values, dtype=float)
    transport = np.asarray(coefficient, dtype=float)
    if field.ndim != 1 or field.shape != transport.shape or radius_m <= 0.0:
        raise ValueError("Mismatched fields or nonpositive radius")
    cell_count = field.size
    dx = 1.0 / cell_count
    faces = np.arange(cell_count + 1, dtype=float) * dx
    volumes = 0.5 * (faces[1:] ** 2 - faces[:-1] ** 2)
    flux = np.zeros(cell_count + 1)
    flux[1:-1] = (
        faces[1:-1]
        * 0.5 * (transport[:-1] + transport[1:])
        * np.diff(field) / dx
    )
    surface, outward = surface_from_series_resistance(
        field[-1], transport[-1], radius_m, cell_count, exchange, ambient
    )
    # flux stores xi*a*u_xi. At xi=1 this equals -R*J.
    flux[-1] = -radius_m * float(outward)
    return np.diff(flux) / (radius_m**2 * volumes), float(surface)


def solve_independent_problem4(
    history: ChamberHistory,
    radius_history: RadiusHistory,
    cell_count: int,
    output_times_s: ArrayLike,
    *,
    threshold: float | None = 0.15,
    parameters: Problem1Parameters = Problem1Parameters(),
    rtol: float = 2.0e-9,
    max_step_s: float = 60.0,
) -> IndependentSolution:
    """Solve directly for T,C and terminate on the reconstructed global maximum.

    The checked set includes all cell-centred values, the even-reconstructed
    axis and the surface. In this monotone drying case its maximum is the
    reconstructed axis. This avoids mistaking xi=dx/2 for the actual centre.
    """
    times = np.asarray(output_times_s, dtype=float)
    if cell_count < 2 or times.ndim != 1 or times.size == 0:
        raise ValueError("At least two cells and nonempty output times required")
    if np.any(np.diff(times) <= 0.0) or times[0] < history.time_s[0]:
        raise ValueError("Invalid output-time ordering")
    if times[-1] > min(history.time_s[-1], radius_history.time_s[-1]):
        raise ValueError("Output horizon exceeds supplied history")
    if threshold is not None and not (0.0 < threshold < parameters.initial_moisture_concentration):
        raise ValueError("Invalid moisture threshold")
    n = cell_count
    xi = (np.arange(n, dtype=float) + 0.5) / n
    tri = diags([np.ones(n - 1), np.ones(n), np.ones(n - 1)], [-1, 0, 1], format="csr")
    sparsity = bmat([[tri, tri], [tri, tri]], format="csr")

    def surroundings(time_s):
        radius = float(np.interp(time_s, radius_history.time_s, radius_history.radius_m))
        ambient_t = float(np.interp(time_s, history.time_s, history.temperature_c))
        ambient_c = float(np.interp(time_s, history.time_s, history.moisture_concentration))
        return radius, ambient_t, ambient_c

    def rhs(time_s, state):
        temperature, moisture = state[:n], state[n:]
        capacity, conductivity, diffusivity = appendix4_coefficients(moisture, temperature)
        radius, ambient_t, ambient_c = surroundings(time_s)
        heat_div, _ = cell_flux_divergence(
            temperature, conductivity, radius,
            parameters.heat_transfer_coefficient_w_m2_k, ambient_t,
        )
        mass_div, _ = cell_flux_divergence(
            moisture, diffusivity, radius,
            parameters.mass_transfer_coefficient_m_s, ambient_c,
        )
        return np.concatenate([heat_div / capacity, mass_div])

    def reconstructed_maximum(time_s, state):
        temperature, moisture = state[:n], state[n:]
        _, _, diffusion = appendix4_coefficients(moisture[-1], temperature[-1])
        radius, _, ambient_c = surroundings(time_s)
        surface, _ = surface_from_series_resistance(
            moisture[-1], diffusion, radius, n,
            parameters.mass_transfer_coefficient_m_s, ambient_c,
        )
        return max(float(np.max(moisture)), float(reconstruct_centre(moisture)), float(surface))

    def event(time_s, state):
        return reconstructed_maximum(time_s, state) - threshold

    event.terminal = True
    event.direction = -1.0
    initial_state = np.concatenate([
        np.full(n, parameters.initial_temperature_c),
        np.full(n, parameters.initial_moisture_concentration),
    ])
    result = solve_ivp(
        rhs, (float(history.time_s[0]), float(times[-1])), initial_state,
        t_eval=times, method="BDF", rtol=rtol,
        atol=np.concatenate([np.full(n, 1.0e-9), np.full(n, 1.0e-11)]),
        max_step=max_step_s, jac_sparsity=sparsity,
        events=event if threshold is not None else None,
    )
    if not result.success:
        raise RuntimeError(result.message)
    endpoint = cell_maximum = reconstructed_max = None
    result_times, states = result.t, result.y
    if threshold is not None:
        if not result.t_events[0].size:
            raise RuntimeError("Threshold not reached within supplied horizon")
        endpoint = float(result.t_events[0][0])
        endpoint_state = result.y_events[0][0]
        cell_maximum = float(np.max(endpoint_state[n:]))
        reconstructed_max = reconstructed_maximum(endpoint, endpoint_state)
        result_times = np.append(result_times, endpoint)
        states = np.column_stack([states, endpoint_state])
    temperatures, moistures = states[:n].T, states[n:].T
    radii = np.interp(result_times, radius_history.time_s, radius_history.radius_m)
    ambient_ts = np.interp(result_times, history.time_s, history.temperature_c)
    ambient_cs = np.interp(result_times, history.time_s, history.moisture_concentration)
    _, last_k, last_d = appendix4_coefficients(moistures[:, -1], temperatures[:, -1])
    surface_t, _ = surface_from_series_resistance(
        temperatures[:, -1], last_k, radii, n,
        parameters.heat_transfer_coefficient_w_m2_k, ambient_ts,
    )
    surface_c, _ = surface_from_series_resistance(
        moistures[:, -1], last_d, radii, n,
        parameters.mass_transfer_coefficient_m_s, ambient_cs,
    )
    return IndependentSolution(
        result_times, xi, temperatures, moistures, surface_t, surface_c,
        reconstruct_centre(moistures), endpoint, cell_maximum, reconstructed_max,
    )


def sample_independent_moisture(
    solution: IndependentSolution,
    radius_history: RadiusHistory,
    physical_radii_cm: ArrayLike,
) -> NDArray[np.float64]:
    """Piecewise-linear interpolation including separately reconstructed ends."""
    radii_m = np.asarray(physical_radii_cm, dtype=float) * 0.01
    if radii_m.ndim != 1 or np.any(radii_m < 0.0):
        raise ValueError("Nonnegative one-dimensional physical radii required")
    result = np.full((solution.time_s.size, radii_m.size), np.nan)
    xi = np.concatenate([[0.0], solution.xi_centres, [1.0]])
    for row, time_s in enumerate(solution.time_s):
        radius = float(np.interp(time_s, radius_history.time_s, radius_history.radius_m))
        field = np.concatenate([
            [solution.centre_moisture[row]], solution.moisture[row],
            [solution.surface_moisture[row]],
        ])
        valid = radii_m <= radius + 1.0e-12
        result[row, valid] = np.interp(radii_m[valid] / radius, xi, field)
    return result
