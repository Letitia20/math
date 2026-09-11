"""Inspect final Word/PDF structure and render contact sheets for visual review."""
from pathlib import Path
import json
import zipfile
from lxml import etree
import pymupdf
from PIL import Image, ImageOps, ImageDraw

root=Path(__file__).resolve().parents[1]
out=root/'deliverables'
pdf=pymupdf.open(out/'A题论文_提交预览.pdf')
texts=[p.get_text() for p in pdf]
appendix=next(i+1 for i,t in enumerate(texts) if '附录 A' in t)
assert appendix <= 32, f'Body too long: appendix starts on {appendix}'
assert '摘要' in texts[0].replace(' ','')
assert '关键词' in texts[0], 'Abstract spills past page 1'
assert '1 问题重述' in texts[1], 'Main text should start on page 2'
assert out.joinpath('A题论文_提交预览.pdf').stat().st_size < 20_000_000
with zipfile.ZipFile(out/'A题论文_可编辑版.docx') as z:
    xml=etree.fromstring(z.read('word/document.xml'))
    ns={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','m':'http://schemas.openxmlformats.org/officeDocument/2006/math'}
    equations=len(xml.findall('.//m:oMath',ns))
    tables=xml.findall('.//w:tbl',ns)
    assert len(tables)==10, f'Expected symbols + 9 tables, got {len(tables)}'
    assert equations>=70
    assert b'MATHPLACEHOLDER' not in z.read('word/document.xml')
    assert b'{{TABLE' not in z.read('word/document.xml')
    for m in xml.findall('.//w:pgMar',ns):
        assert all(int(m.get('{'+ns['w']+'}'+k))>=1418 for k in ('top','bottom','left','right'))
report={'total_pages':len(pdf),'appendix_first_page':appendix,'abstract_pages':1,'body_pages':appendix-2,'native_editable_equations':equations,'editable_tables':len(tables),'pdf_bytes':out.joinpath('A题论文_提交预览.pdf').stat().st_size}
review=root/'tmp/paper/layout';review.mkdir(parents=True,exist_ok=True)
for batch in range(0,appendix,6):
    sheet=Image.new('RGB',(1260,2*615),'#d0d0d0');draw=ImageDraw.Draw(sheet)
    for j,idx in enumerate(range(batch,min(batch+6,appendix))):
        pix=pdf[idx].get_pixmap(matrix=pymupdf.Matrix(0.7,0.7))
        im=Image.frombytes('RGB',[pix.width,pix.height],pix.samples)
        im.thumbnail((400,570));x=(j%3)*420+10;y=(j//3)*615+30
        sheet.paste(im,(x,y));draw.text((x,y-20),f'Page {idx+1}',fill='black')
    sheet.save(review/f'pages_{batch+1}.png')
for idx in (0,2,appendix-2):
    pdf[idx].get_pixmap(matrix=pymupdf.Matrix(1.5,1.5)).save(review/f'detail_{idx+1}.png')
(out/'layout_audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(report))
