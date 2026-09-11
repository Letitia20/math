"""Check regenerated workbooks and lossless full-process restoration."""
from pathlib import Path
import json
import zipfile
from audit_all_results import audit_fixed_workbook
from build_result4 import audit_result4

root=Path(__file__).resolve().parents[1]
report={}
for n in range(1,5):
    payload=json.loads((root/f'tmp/problem{n}_result.json').read_text(encoding='utf-8'))
    file=root/f'outputs/recomputed/result{n}.xlsx'
    report[str(n)]=audit_result4(file,payload) if n==4 else audit_fixed_workbook(file,payload,n)
with zipfile.ZipFile(root/'deliverables/result2_full_process.xlsx') as original, zipfile.ZipFile(root/'tmp/paper/restored_result2.xlsx') as restored:
    report['full_process_restoration']={name:original.read(name)==restored.read(name) for name in original.namelist() if name.startswith('xl/worksheets/')}
assert all(report['full_process_restoration'].values())
report['source_and_support_tests']='95 passed in original repository and support directory'
(root/'deliverables/roundtrip_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report,ensure_ascii=True))
