"""Reproduce a separately implemented cell-centred discretization of Problem 4."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import platform
import time

import numpy as np
import scipy

from problem1 import load_chamber_history_csv
from problem3 import extend_chamber_history_to_plateau
from problem4 import load_radius_history_csv
from problem4_independent import (
    sample_independent_moisture,
    solve_independent_problem4,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", nargs="+", type=int, default=[100, 200, 400, 800, 1600, 3200])
    parser.add_argument("--output", type=Path, default=Path("problem4_independent_validation.json"))
    args = parser.parse_args()
    raw_chamber = load_chamber_history_csv("attachment1.csv")
    chamber = extend_chamber_history_to_plateau(raw_chamber, 72 * 3600.0)
    radius = load_radius_history_csv("attachment2.csv")
    table_hours = np.arange(6.0, 48.1, 6.0)
    requested_times = np.concatenate([[0.0], table_hours * 3600.0, [72 * 3600.0]])
    reference_table = np.asarray([
        [1.7196, 1.5377, 1.0226, 0.4207],
        [0.7377, 0.6546, 0.4083, 0.1669],
        [0.4088, 0.3687, 0.2397, 0.0892],
        [0.2851, 0.2611, 0.1790, 0.0673],
        [0.2263, 0.2094, 0.1493, 0.0595],
        [0.1928, 0.1797, 0.1317, 0.0560],
        [0.1712, 0.1604, 0.1200, 0.0541],
        [0.1561, 0.1468, 0.1116, 0.0530],
    ])
    grids = {}
    for cells in args.cells:
        started = time.perf_counter()
        solution = solve_independent_problem4(chamber, radius, cells, requested_times)
        sampled = sample_independent_moisture(solution, radius, [0.0, 0.5, 1.0])
        table = np.column_stack([sampled[1:9], solution.surface_moisture[1:9]])
        full_profiles = np.column_stack([
            solution.centre_moisture, solution.moisture, solution.surface_moisture,
        ])
        grids[str(cells)] = {
            "cells": cells,
            "drying_time_s": solution.drying_time_s,
            "drying_time_h": solution.drying_time_s / 3600.0,
            "endpoint_cell_maximum": solution.endpoint_cell_maximum,
            "endpoint_reconstructed_maximum": solution.endpoint_reconstructed_maximum,
            "endpoint_centre_moisture": float(solution.centre_moisture[-1]),
            "endpoint_surface_moisture": float(solution.surface_moisture[-1]),
            "table_unrounded": table.tolist(),
            "table_rounded_4dp": np.round(table, 4).tolist(),
            "maximum_outward_increase_at_stored_times": float(np.max(np.diff(full_profiles, axis=1))),
            "elapsed_s": time.perf_counter() - started,
        }
        print(f"cells={cells}: endpoint {solution.drying_time_s:.9f} s = "
              f"{solution.drying_time_s / 3600:.10f} h", flush=True)

    output = {
        "purpose": "Numerical cross-check of the same assumed continuum model; not independent experimental validation.",
        "method": {
            "geometry": "N complete annular cells; point unknowns at xi=(j+1/2)/N",
            "internal_faces": "Arithmetic mean transport coefficient with centred gradient",
            "outer_boundary": "Last-centre half-cell resistance R/(2*N*a_last) plus convection resistance 1/h",
            "centre_reconstruction": "(9*C[0]-C[1])/8 from an even quadratic",
            "event": "Maximum of reconstructed centre, all cell-centre values and reconstructed surface equals 0.15",
            "time_integrator": "SciPy BDF for T and C directly; rtol=2e-9, atol(T)=1e-9, atol(C)=1e-11, max_step=60s",
            "material_laws": "Independent transcription of Appendix 4 with Kelvin in diffusivity",
            "shared_inputs": "Same Attachment 1 and 2 CSV inputs, linear radius and chamber interpolation, final-hour chamber means with 60s transition",
            "limitation": "Arithmetic-face and frozen half-cell coefficient errors, and interior linear sampling, require convergence checks. The reconstruction is not a rigorous bound on continuous subcell maxima.",
        },
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__},
        "chamber_plateau": {"temperature_c": float(chamber.temperature_c[-1]), "moisture": float(chamber.moisture_concentration[-1])},
        "table_time_h": table_hours.tolist(),
        "table_columns": ["0 cm", "0.5 cm", "1.0 cm", "actual surface"],
        "baseline_report": {"drying_time_h": 51.091955, "table_6_rounded_4dp": reference_table.tolist()},
        "unverified_docx_reference_time_h": 51.091999,
        "grids": grids,
    }
    extrapolations = {}
    for coarse in args.cells:
        fine = 2 * coarse
        if str(fine) not in grids:
            continue
        coarse_value, fine_value = grids[str(coarse)], grids[str(fine)]
        extrap_time = fine_value["drying_time_s"] + (
            fine_value["drying_time_s"] - coarse_value["drying_time_s"]
        ) / 3.0
        extrap_table = np.asarray(fine_value["table_unrounded"]) + (
            np.asarray(fine_value["table_unrounded"]) - np.asarray(coarse_value["table_unrounded"])
        ) / 3.0
        rounded = np.round(extrap_table, 4)
        mismatches = np.argwhere(rounded != reference_table)
        extrapolations[f"{coarse}/{fine}"] = {
            "assumed_order": 2,
            "drying_time_s": extrap_time,
            "drying_time_h": extrap_time / 3600.0,
            "delta_from_baseline_s": extrap_time - 51.091955 * 3600.0,
            "delta_from_unverified_docx_reference_s": extrap_time - 51.091999 * 3600.0,
            "table_unrounded": extrap_table.tolist(),
            "table_rounded_4dp": rounded.tolist(),
            "maximum_absolute_difference_from_rounded_table_6": float(np.max(np.abs(extrap_table - reference_table))),
            "rounded_table_mismatch_count": int(len(mismatches)),
            "rounded_table_mismatches": [
                {"time_h": float(table_hours[row]), "column": output["table_columns"][col],
                 "baseline": float(reference_table[row, col]), "independent": float(rounded[row, col])}
                for row, col in mismatches
            ],
        }
    output["second_order_richardson"] = extrapolations
    orders = {}
    for coarse in args.cells:
        if str(2 * coarse) in grids and str(4 * coarse) in grids:
            a, b, c = [grids[str(k)]["drying_time_s"] for k in (coarse, 2 * coarse, 4 * coarse)]
            orders[f"{coarse}/{2 * coarse}/{4 * coarse}"] = float(np.log2(abs((a - b) / (b - c))))
    output["observed_endpoint_orders"] = orders
    output["interpretation"] = (
        "The endpoint differences decrease with refinement, but observed orders remain below two. "
        "The second-order Richardson values are diagnostic extrapolations, not certified error bounds. "
        "Agreement with a production-solver extrapolation does not validate the shared physical assumptions. "
        "Four-decimal table agreement is assessed against the rounded published table; its raw differences "
        "therefore include reference rounding and are not discretization error estimates."
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"extrapolated_endpoints_h": {key: row["drying_time_h"] for key, row in extrapolations.items()},
                      "observed_orders": orders}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
