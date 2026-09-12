"""Restore the complete, rounded, per-second result2.xlsx from the lossless NPZ."""
from pathlib import Path
import argparse
import numpy as np
from openpyxl import Workbook
from openpyxl.cell import WriteOnlyCell

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,default=Path('result2_full_process.npz'))
    parser.add_argument('--output',type=Path,default=Path('result2_full_process.xlsx'))
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f'Refusing to overwrite {args.output}; choose another --output')
    wb = Workbook(write_only=True)
    with np.load(args.input,allow_pickle=False) as archive:
        for name in ('温度','水分浓度'):
            ws = wb.create_sheet(name)
            ws.append(['时间\\到药材中心的距离']+[i/10 for i in range(21)])
            for row in archive[name]:
                cells = [int(row[0])]
                for value in row[1:]:
                    cell = WriteOnlyCell(ws,float(value)); cell.number_format='0.0000'; cells.append(cell)
                ws.append(cells)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    wb.save(args.output)
    print(args.output)

if __name__=='__main__': main()
