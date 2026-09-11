"""Package selected support files, validate the archive, and record hashes."""
from pathlib import Path
import hashlib
import json
import zipfile

root=Path(__file__).resolve().parents[1]
support=root/'deliverables/supporting_materials'
output=root/'deliverables/A题支撑材料.zip'
files=sorted(p for p in support.rglob('*') if p.is_file() and not any(s in {'__pycache__','.pytest_cache'} for s in p.parts) and p.suffix!='.pyc')
assert all(p.stat().st_size<20_000_000 for p in files)
assert (support/'AI工具使用详情.pdf') in files
hashes={p.relative_to(support).as_posix():hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
with zipfile.ZipFile(output,'w',zipfile.ZIP_DEFLATED,compresslevel=9) as z:
    for p in files: z.write(p,p.relative_to(support).as_posix())
with zipfile.ZipFile(output) as z:
    assert z.testzip() is None
    assert set(z.namelist())==set(hashes)
    for name,digest in hashes.items(): assert hashlib.sha256(z.read(name)).hexdigest()==digest
assert output.stat().st_size<20_000_000
report={'archive_bytes':output.stat().st_size,'file_count':len(files),'all_member_sha256_checked':True,'files':hashes}
(root/'deliverables/support_archive_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({k:v for k,v in report.items() if k!='files'}))
