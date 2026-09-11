"""Fill the official Problem 4 workbook template with numeric results."""

from __future__ import annotations

import argparse
import json
import math
from copy import copy
from pathlib import Path

import openpyxl
from openpyxl.comments import Comment


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("template", type=Path)
    parser.add_argument("payload", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    payload = json.loads(args.payload.read_text(encoding="utf-8"))
    times = payload["time_s"]
    radii = payload["radius_cm"]
    matrix = payload["moisture_concentration"]
    surface = payload["surface_moisture_concentration"]
    if len(radii) != 21 or len(times) != len(matrix) or len(times) != len(surface):
        raise ValueError("Payload dimensions do not match the result4 template")
    if any(len(row) != len(radii) for row in matrix):
        raise ValueError("Moisture matrix has an invalid shape")

    workbook = openpyxl.load_workbook(args.template)
    if "Sheet1" not in workbook.sheetnames:
        raise ValueError("Official result4 template must contain Sheet1")
    sheet = workbook["Sheet1"]
    # The template has a header followed by placeholder rows.  Preserve its
    # sheet name and styling, replacing only the data region with values.
    sheet.cell(1, 1).value = "时间\\到药材中心的距离"
    for col, radius in enumerate(radii, start=2):
        sheet.cell(1, col).value = radius
    sheet.cell(1, len(radii) + 2).value = "药材表面"
    sheet.cell(1, 1).comment = Comment(
        "时间单位：s；固定物理半径和药材表面含水率单位：kg/kg。\n"
        "固定半径超过该时刻收缩半径的单元留空。", "OpenAI"
    )
    required_rows = len(times) + 1
    if sheet.max_row > required_rows:
        sheet.delete_rows(required_rows + 1, sheet.max_row - required_rows)
    for row_index, time in enumerate(times, start=2):
        sheet.cell(row_index, 1).value = time
        for col_index, value in enumerate(matrix[row_index - 2], start=2):
            sheet.cell(row_index, col_index).value = None if value is None or (isinstance(value, float) and math.isnan(value)) else value
        value = surface[row_index - 2]
        sheet.cell(row_index, len(radii) + 2).value = None if value is None or (isinstance(value, float) and math.isnan(value)) else value

    # Match the supplied template's centered Chinese layout and numeric format.
    sheet.column_dimensions["A"].width = 24
    for col in range(2, len(radii) + 3):
        sheet.column_dimensions[openpyxl.utils.get_column_letter(col)].width = 11
    for row in sheet.iter_rows(min_row=1, max_row=required_rows, min_col=1, max_col=len(radii) + 2):
        for cell in row:
            alignment = copy(cell.alignment)
            alignment.horizontal = "center"
            alignment.vertical = "center"
            cell.alignment = alignment
            if cell.row == 1:
                font = copy(cell.font)
                font.name = "宋体"
                font.size = 10
                cell.font = font
            elif cell.column == 1:
                cell.number_format = "0.######"
            else:
                cell.number_format = "0.0000"
    sheet.freeze_panes = "B2"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    workbook.save(args.output)

    # A result workbook should contain values only; formulas would violate the
    # official template's numeric-output convention.
    check = openpyxl.load_workbook(args.output, data_only=False)
    formulas = [cell.coordinate for row in check.active.iter_rows() for cell in row if isinstance(cell.value, str) and cell.value.startswith("=")]
    if formulas:
        raise RuntimeError(f"Unexpected formulas in result workbook: {formulas[:5]}")
    print(f"saved {args.output} ({required_rows} rows x {len(radii) + 2} columns)")


if __name__ == "__main__":
    main()
