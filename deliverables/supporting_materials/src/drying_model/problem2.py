"""Coupled variable-property model for Problem 2."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import solve_ivp
from scipy.sparse import bmat, diags

from drying_model.problem1 import (
    ChamberHistory,
    Problem1Parameters,
    Problem1Solution,
    interpolate_history,
    radial_flux_divergence,
    sample_solution,
)

ScalarLaw = Callable[[NDArray[np.float64]], NDArray[np.float64]]
CoupledLaw = Callable[
    [NDArray[np.float64], NDArray[np.float64]],
    NDArray[np.float64],
]


def _concentration_array(concentration: ArrayLike) -> NDArray[np.float64]:
    values = np.asarray(concentration, dtype=float)
    if np.any(~np.isfinite(values)) or np.any(values < 0.0):
        raise ValueError("Moisture concentration must be finite and non-negative")
    return values


def density_q2(concentration: ArrayLike) -> NDArray[np.float64]:
    """Return Appendix 3 density in kilograms per cubic metre."""
    values = _concentration_array(concentration)
    return 650.0 + 128.0 * values


def heat_capacity_q2(concentration: ArrayLike) -> NDArray[np.float64]:
    """Return Appendix 3 specific heat capacity in joules per kilogram-kelvin."""
    values = _concentration_array(concentration)
    return 1450.0 + 2736.0 * values / (values + 1.0)


def thermal_conductivity_q2(concentration: ArrayLike) -> NDArray[np.float64]:
    """Return Appendix 3 thermal conductivity in watts per metre-kelvin."""
    values = _concentration_array(concentration)
    return 0.21 + 0.38 * values / (values + 1.0)


def diffusivity_q2(
    concentration: ArrayLike,
    temperature_c: ArrayLike,
) -> NDArray[np.float64]:
    """Return Appendix 3 moisture diffusivity using absolute temperature."""
    moisture = _concentration_array(concentration)
    temperature_k = np.asarray(temperature_c, dtype=float) + 273.15
    if np.any(~np.isfinite(temperature_k)) or np.any(temperature_k <= 0.0):
        raise ValueError("Temperature must be finite and above absolute zero")
    if np.any(moisture <= 0.0):
        raise ValueError("Moisture concentration must be positive for diffusivity")
    moisture, temperature_k = np.broadcast_arrays(moisture, temperature_k)
    return 2.4e-3 * np.exp(-0.45 / moisture) * np.exp(-3850.0 / temperature_k)


def solve_variable_property_fixed_cylinder(
    history: ChamberHistory,
    radial_intervals: int,
    output_times_s: ArrayLike,
    *,
    density_law: ScalarLaw,
    heat_capacity_law: ScalarLaw,
    conductivity_law: ScalarLaw,
    diffusivity_law: CoupledLaw,
    parameters: Problem1Parameters = Problem1Parameters(),
    relative_tolerance: float = 1.0e-8,
    max_step_s: float = 5.0,
    maximum_moisture_threshold: float | None = None,
    initial_time_s: float | None = None,
    initial_temperature_c: ArrayLike | None = None,
    initial_moisture_concentration: ArrayLike | None = None,
) -> Problem1Solution:
    """Solve coupled heat and moisture transport on a fixed cylindrical radius."""
    if radial_intervals < 2:
        raise ValueError("At least two radial intervals are required")
    output_times = np.asarray(output_times_s, dtype=float)
    if output_times.ndim != 1 or output_times.size == 0:
        raise ValueError("Output times must be a non-empty one-dimensional array")
    if np.any(np.diff(output_times) <= 0.0):
        raise ValueError("Output times must be strictly increasing")
    start_time = history.time_s[0] if initial_time_s is None else float(initial_time_s)
    if (
        not np.isfinite(start_time)
        or start_time < history.time_s[0]
        or start_time > output_times[0]
        or output_times[-1] > history.time_s[-1]
    ):
        raise ValueError("Output times must lie within the chamber history")
    if maximum_moisture_threshold is not None and (
        not np.isfinite(maximum_moisture_threshold)
        or maximum_moisture_threshold <= 0.0
        or maximum_moisture_threshold >= parameters.initial_moisture_concentration
    ):
        raise ValueError("Maximum moisture threshold must lie between zero and the initial value")

    node_count = radial_intervals + 1
    radii = np.linspace(0.0, parameters.radius_m, node_count)
    tridiagonal = diags(
        [np.ones(node_count - 1), np.ones(node_count), np.ones(node_count - 1)],
        offsets=[-1, 0, 1],
        shape=(node_count, node_count),
        format="csr",
    )
    sparsity = bmat(
        [[tridiagonal, tridiagonal], [tridiagonal, tridiagonal]],
        format="csr",
    )

    def coupled_rhs(time_s: float, state: NDArray[np.float64]) -> NDArray[np.float64]:
        temperature_c = state[:node_count]
        concentration = np.exp(state[node_count:])
        density = np.asarray(density_law(concentration), dtype=float)
        heat_capacity = np.asarray(heat_capacity_law(concentration), dtype=float)
        conductivity = np.asarray(conductivity_law(concentration), dtype=float)
        diffusivity = np.asarray(
            diffusivity_law(concentration, temperature_c),
            dtype=float,
        )
        expected_shape = concentration.shape
        for name, values in (
            ("density", density),
            ("heat capacity", heat_capacity),
            ("thermal conductivity", conductivity),
            ("diffusivity", diffusivity),
        ):
            if values.shape != expected_shape or np.any(~np.isfinite(values)) or np.any(values <= 0.0):
                raise ValueError(f"{name} law must return finite positive nodal values")

        ambient_temperature = interpolate_history(
            time_s,
            history.time_s,
            history.temperature_c,
        )
        temperature_rate = radial_flux_divergence(
            temperature_c,
            conductivity,
            parameters.radius_m,
            parameters.heat_transfer_coefficient_w_m2_k,
            ambient_temperature,
        ) / (density * heat_capacity)

        ambient_moisture = interpolate_history(
            time_s,
            history.time_s,
            history.moisture_concentration,
        )
        concentration_rate = radial_flux_divergence(
            concentration,
            diffusivity,
            parameters.radius_m,
            parameters.mass_transfer_coefficient_m_s,
            ambient_moisture,
        )
        log_concentration_rate = concentration_rate / concentration
        return np.concatenate([temperature_rate, log_concentration_rate])

    supplied_initial_fields = (
        initial_temperature_c is not None,
        initial_moisture_concentration is not None,
    )
    if any(supplied_initial_fields) and not all(supplied_initial_fields):
        raise ValueError("Both initial fields must be supplied together")
    if all(supplied_initial_fields):
        temperature_initial = np.asarray(initial_temperature_c, dtype=float)
        moisture_initial = np.asarray(initial_moisture_concentration, dtype=float)
        if (
            temperature_initial.shape != (node_count,)
            or moisture_initial.shape != (node_count,)
            or np.any(~np.isfinite(temperature_initial))
            or np.any(~np.isfinite(moisture_initial))
            or np.any(moisture_initial <= 0.0)
        ):
            raise ValueError("Initial fields must be finite positive nodal arrays")
    else:
        temperature_initial = np.full(node_count, parameters.initial_temperature_c)
        moisture_initial = np.full(
            node_count,
            parameters.initial_moisture_concentration,
        )
    initial_state = np.concatenate(
        [temperature_initial, np.log(moisture_initial)]
    )
    absolute_tolerance = np.concatenate(
        [np.full(node_count, 1.0e-9), np.full(node_count, 1.0e-10)]
    )

    threshold_event = None
    if maximum_moisture_threshold is not None:

        def threshold_event(time_s: float, state: NDArray[np.float64]) -> float:
            del time_s
            return float(
                np.max(np.exp(state[node_count:])) - maximum_moisture_threshold
            )

        threshold_event.terminal = True
        threshold_event.direction = -1.0

    result = solve_ivp(
        coupled_rhs,
        (start_time, float(output_times[-1])),
        initial_state,
        method="BDF",
        t_eval=output_times,
        rtol=relative_tolerance,
        atol=absolute_tolerance,
        max_step=max_step_s,
        jac_sparsity=sparsity,
        events=threshold_event,
    )
    if not result.success:
        raise RuntimeError(f"Coupled solver failed: {result.message}")
    result_time = result.t
    result_state = result.y
    if maximum_moisture_threshold is not None:
        if not result.t_events or result.t_events[0].size == 0:
            raise RuntimeError("Maximum moisture threshold was not reached")
        event_time = float(result.t_events[0][0])
        event_state = result.y_events[0][0]
        if result_time.size and np.isclose(
            result_time[-1],
            event_time,
            rtol=0.0,
            atol=1.0e-9,
        ):
            result_time = result_time.copy()
            result_state = result_state.copy()
            result_time[-1] = event_time
            result_state[:, -1] = event_state
        else:
            result_time = np.append(result_time, event_time)
            result_state = np.column_stack([result_state, event_state])
    return Problem1Solution(
        time_s=result_time,
        radius_m=radii,
        temperature_c=result_state[:node_count].T,
        moisture_concentration=np.exp(result_state[node_count:].T),
    )


def solve_problem2(
    history: ChamberHistory,
    radial_intervals: int,
    output_times_s: ArrayLike,
    *,
    parameters: Problem1Parameters = Problem1Parameters(),
    relative_tolerance: float = 1.0e-8,
    max_step_s: float = 5.0,
    maximum_moisture_threshold: float | None = None,
    initial_time_s: float | None = None,
    initial_temperature_c: ArrayLike | None = None,
    initial_moisture_concentration: ArrayLike | None = None,
) -> Problem1Solution:
    """Solve Problem 2 with the four Appendix 3 material-property laws."""
    return solve_variable_property_fixed_cylinder(
        history,
        radial_intervals,
        output_times_s,
        density_law=density_q2,
        heat_capacity_law=heat_capacity_q2,
        conductivity_law=thermal_conductivity_q2,
        diffusivity_law=diffusivity_q2,
        parameters=parameters,
        relative_tolerance=relative_tolerance,
        max_step_s=max_step_s,
        maximum_moisture_threshold=maximum_moisture_threshold,
        initial_time_s=initial_time_s,
        initial_temperature_c=initial_temperature_c,
        initial_moisture_concentration=initial_moisture_concentration,
    )


def solve_problem2_sampled_in_chunks(
    history: ChamberHistory,
    radial_intervals: int,
    output_times_s: ArrayLike,
    sample_radius_cm: ArrayLike,
    *,
    chunk_duration_s: float = 600.0,
    parameters: Problem1Parameters = Problem1Parameters(),
    relative_tolerance: float = 1.0e-8,
    max_step_s: float = 5.0,
) -> Problem1Solution:
    """Solve Problem 2 in bounded-memory chunks and retain sampled radii only."""
    output_times = np.asarray(output_times_s, dtype=float)
    if output_times.ndim != 1 or output_times.size == 0:
        raise ValueError("Output times must be a non-empty one-dimensional array")
    if np.any(~np.isfinite(output_times)) or np.any(np.diff(output_times) <= 0.0):
        raise ValueError("Output times must be finite and strictly increasing")
    if (
        output_times[0] < history.time_s[0]
        or output_times[-1] > history.time_s[-1]
    ):
        raise ValueError("Output times must lie within the chamber history")
    if not np.isfinite(chunk_duration_s) or chunk_duration_s <= 0.0:
        raise ValueError("Chunk duration must be finite and positive")

    requested_radius_cm = np.asarray(sample_radius_cm, dtype=float)
    node_count = radial_intervals + 1
    nodal_radius = np.linspace(0.0, parameters.radius_m, node_count)
    initial_solution = Problem1Solution(
        time_s=np.array([history.time_s[0]]),
        radius_m=nodal_radius,
        temperature_c=np.full((1, node_count), parameters.initial_temperature_c),
        moisture_concentration=np.full(
            (1, node_count),
            parameters.initial_moisture_concentration,
        ),
    )
    sampled_initial = sample_solution(initial_solution, requested_radius_cm)
    sampled_temperature = np.empty(
        (output_times.size, sampled_initial.radius_m.size),
        dtype=float,
    )
    sampled_moisture = np.empty_like(sampled_temperature)

    next_output_index = 0
    current_time = float(history.time_s[0])
    current_temperature: NDArray[np.float64] | None = None
    current_moisture: NDArray[np.float64] | None = None
    if np.isclose(output_times[0], current_time, rtol=0.0, atol=1.0e-12):
        sampled_temperature[0] = sampled_initial.temperature_c[0]
        sampled_moisture[0] = sampled_initial.moisture_concentration[0]
        next_output_index = 1

    final_time = float(output_times[-1])
    while current_time < final_time:
        chunk_end = min(current_time + chunk_duration_s, final_time)
        chunk_stop = int(np.searchsorted(output_times, chunk_end, side="right"))
        requested_times = output_times[next_output_index:chunk_stop]
        if requested_times.size and requested_times[0] <= current_time:
            raise RuntimeError("Chunked output indexing did not advance")
        append_chunk_end = (
            requested_times.size == 0
            or not np.isclose(requested_times[-1], chunk_end, rtol=0.0, atol=1.0e-12)
        )
        integration_times = (
            np.append(requested_times, chunk_end)
            if append_chunk_end
            else requested_times
        )
        chunk_solution = solve_problem2(
            history,
            radial_intervals=radial_intervals,
            output_times_s=integration_times,
            parameters=parameters,
            relative_tolerance=relative_tolerance,
            max_step_s=max_step_s,
            initial_time_s=None if current_temperature is None else current_time,
            initial_temperature_c=current_temperature,
            initial_moisture_concentration=current_moisture,
        )
        sampled_chunk = sample_solution(chunk_solution, requested_radius_cm)
        retained_count = requested_times.size
        if retained_count:
            sampled_temperature[next_output_index:chunk_stop] = (
                sampled_chunk.temperature_c[:retained_count]
            )
            sampled_moisture[next_output_index:chunk_stop] = (
                sampled_chunk.moisture_concentration[:retained_count]
            )
        current_temperature = chunk_solution.temperature_c[-1].copy()
        current_moisture = chunk_solution.moisture_concentration[-1].copy()
        current_time = chunk_end
        next_output_index = chunk_stop

    if next_output_index != output_times.size:
        raise RuntimeError("Chunked solver did not populate every requested output time")
    return Problem1Solution(
        time_s=output_times.copy(),
        radius_m=sampled_initial.radius_m,
        temperature_c=sampled_temperature,
        moisture_concentration=sampled_moisture,
    )
