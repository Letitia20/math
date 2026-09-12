"""Finalize editable Word math, then audit and package the revised document.

Run `python finalize_and_audit.py math` after render_paper.mjs;
export PDFs with ../../tools/export_paper_pdf.ps1 -Directory <this directory>;
then run `python finalize_and_audit.py audit`.
"""
from pathlib import Path
import copy
import hashlib
import json
import re
import shutil
import sys
import zipfile

from lxml import etree

OUT = Path(__file__).resolve().parent
ROOT = OUT.parents[1]
NS = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main',
      'm': 'http://schemas.openxmlformats.org/officeDocument/2006/math'}


def finalize_math():
    from latex2mathml.converter import convert
    transform = etree.XSLT(etree.parse('C:/Program Files/Microsoft Office/root/Office16/MML2OMML.XSL'))
    equations = json.loads((OUT / 'equations.json').read_text(encoding='utf-8'))
    mapping = {}
    for item in equations:
        # Keep the evaluation bar/subscript attached to the material derivative.
        latex = item['latex'].replace(r'\left.', '').replace(r'\right|', r'\vert')
        if latex.startswith(r'\begin{aligned}'):
            # Office's MathML XSLT flattens the converter's aligned mtable.
            # Build native equation-array rows explicitly, preserving each equation.
            rows = latex.removeprefix(r'\begin{aligned}').removesuffix(r'\end{aligned}').split(r'\\')
            node = etree.Element('{' + NS['m'] + '}oMath', nsmap={'m': NS['m']})
            arr = etree.SubElement(node, '{' + NS['m'] + '}eqArr')
            props = etree.SubElement(arr, '{' + NS['m'] + '}eqArrPr')
            jc = etree.SubElement(props, '{' + NS['m'] + '}baseJc')
            jc.set('{' + NS['m'] + '}val', 'center')
            for row in rows:
                entry = etree.SubElement(arr, '{' + NS['m'] + '}e')
                rendered = transform(etree.fromstring(convert(row.strip()).encode())).getroot()
                for child in rendered:
                    entry.append(copy.deepcopy(child))
        else:
            node = transform(etree.fromstring(convert(latex).encode())).getroot()
        for run in node.xpath('//*[local-name()="r"]'):
            props = etree.Element('{' + NS['w'] + '}rPr')
            size = etree.SubElement(props, '{' + NS['w'] + '}sz')
            size.set('{' + NS['w'] + '}val', str(item['size']))
            # OMML requires mathematical properties before Word run properties.
            run.insert(1 if run.find('./m:rPr', NS) is not None else 0, props)
        mapping[item['key']] = node
    for filename in ('A题论文_可编辑版.docx', 'AI工具使用详情_可编辑版.docx'):
        dest = OUT / filename
        with zipfile.ZipFile(dest) as z:
            content = {name: z.read(name) for name in z.namelist()}
        doc = etree.fromstring(content['word/document.xml'])
        # Word defaults may hang full-width Chinese punctuation past the margin.
        # docx-js currently drops an explicit false value for this property.
        later_props = {'topLinePunct', 'autoSpaceDE', 'autoSpaceDN', 'bidi', 'adjustRightInd',
                       'snapToGrid', 'spacing', 'ind', 'contextualSpacing', 'mirrorIndents',
                       'suppressOverlap', 'jc', 'textDirection', 'textAlignment',
                       'textboxTightWrap', 'outlineLvl', 'divId', 'cnfStyle', 'rPr', 'sectPr'}
        for paragraph in doc.findall('.//w:p', NS):
            props = paragraph.find('./w:pPr', NS)
            if props is None:
                props = etree.Element('{' + NS['w'] + '}pPr')
                paragraph.insert(0, props)
            for existing in props.findall('./w:overflowPunct', NS):
                props.remove(existing)
            overflow = etree.Element('{' + NS['w'] + '}overflowPunct')
            overflow.set('{' + NS['w'] + '}val', '0')
            index = next((i for i, child in enumerate(props)
                          if etree.QName(child).localname in later_props), len(props))
            props.insert(index, overflow)
        count = 0
        for t in doc.xpath('//w:t', namespaces=NS):
            if t.text in mapping:
                run = t.getparent()
                run.getparent().replace(run, copy.deepcopy(mapping[t.text]))
                count += 1
        content['word/document.xml'] = etree.tostring(doc, xml_declaration=True, encoding='UTF-8', standalone=True)
        assert b'MATHPLACEHOLDER' not in content['word/document.xml']
        with zipfile.ZipFile(dest, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
            for name, data in content.items():
                z.writestr(name, data)
        print(filename, count, 'editable equations')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit():
    import pymupdf
    from PIL import Image, ImageDraw

    old = (ROOT / 'deliverables/paper_complete.md').read_text(encoding='utf-8')
    new = (OUT / 'paper_complete.md').read_text(encoding='utf-8')
    old_source = (ROOT / 'deliverables/paper_source.md').read_text(encoding='utf-8')
    new_source = (OUT / 'paper_source.md').read_text(encoding='utf-8')
    # Compare math in source, before the renderer adds equation layout only.
    old_display = re.findall(r'\$\$(.*?)\$\$', old_source, re.S)
    new_display = re.findall(r'\$\$(.*?)\$\$', new_source, re.S)
    assert old_display == new_display
    old_tables = re.findall(r'(?m)^\|.*\|$', old)
    new_tables = re.findall(r'(?m)^\|.*\|$', new)
    assert old_tables == new_tables, 'A numeric/symbol table changed'
    new_codes = re.findall(r'```python\n(.*?)```', new, re.S)
    assert len(new_codes) == 4, 'Appendix must contain only four core model files'
    references = new_source.split('## 参考文献', 1)[1].split('<!-- APPENDIX -->', 1)[0]
    assert not re.search(r'(?i)\bAI\b|Scientific Agent Skills|Nature Skills|Codex', references)
    assert re.findall(r'(?m)^\[(\d+)\]', references) == ['1', '2', '3']
    for n in ('33.5753', '36.7856', '2.5500', '1.5102', '49.8495', '49.9664',
              '1.7662', '1.0081', '57.4740', '57.4833', '51.0920', '51.1000',
              '129.8481', '0.890289', '10.97%'):
        assert n in new
    with zipfile.ZipFile(OUT / 'A题论文_可编辑版.docx') as z:
        doc = etree.fromstring(z.read('word/document.xml'))
        equations = len(doc.findall('.//m:oMath', NS))
        arrays = doc.findall('.//m:eqArr', NS)
        assert len(arrays) == 4 and all(len(x.findall('./m:e', NS)) == 2 for x in arrays)
        tables = doc.findall('.//w:tbl', NS)
        assert len(tables) == 10
        assert equations >= 79
        assert b'MATHPLACEHOLDER' not in z.read('word/document.xml')
        for margin in doc.findall('.//w:pgMar', NS):
            assert all(int(margin.get('{' + NS['w'] + '}' + x)) >= 1418
                       for x in ('top', 'bottom', 'left', 'right'))
        for table in tables:
            vertical = table.find('./w:tblPr/w:tblBorders/w:insideV', NS)
            assert vertical is not None and vertical.get('{' + NS['w'] + '}val') in ('none', 'nil')
    pdf = pymupdf.open(OUT / 'A题论文_提交预览.pdf')
    texts = [p.get_text() for p in pdf]
    appendix = next(i + 1 for i, t in enumerate(texts) if re.search(r'(?m)^附录 A 支撑材料与复现说明', t))
    assert '关键词' in texts[0] and '1 问题重述' in texts[1]
    assert appendix - 2 <= 30
    assert all(abs(p.rect.width-595.3) < 1 and abs(p.rect.height-841.9) < 1 for p in pdf)
    for i, t in enumerate(texts):
        assert t.strip(), f'Empty page {i+1}'
    # Right-aligned numbering must survive the Word/PDF export for all displays.
    main_text = '\n'.join(texts[:appendix-1])
    for n in range(1, 23):
        assert re.search(r'\(' + str(n) + r'\)', main_text), f'Missing equation number {n}'
    review = OUT / 'layout_review'
    if review.exists():
        shutil.rmtree(review)
    review.mkdir()
    for batch in range(0, len(pdf), 12):
        sheet = Image.new('RGB', (1600, 3 * 585), '#dedede')
        draw = ImageDraw.Draw(sheet)
        for j, idx in enumerate(range(batch, min(batch+12, len(pdf)))):
            pix = pdf[idx].get_pixmap(matrix=pymupdf.Matrix(.65, .65))
            im = Image.frombytes('RGB', [pix.width, pix.height], pix.samples)
            im.thumbnail((380, 550))
            x, y = (j % 4) * 400 + 10, (j // 4) * 585 + 25
            sheet.paste(im, (x, y))
            draw.text((x, y-20), f'Page {idx+1}', fill='black')
        sheet.save(review / f'pages_{batch+1:02}.png')
    for idx in range(min(appendix, len(pdf))):
        pdf[idx].get_pixmap(matrix=pymupdf.Matrix(1.2, 1.2)).save(review / f'detail_{idx+1:02}.png')
    page_bounds = []
    for i, page in enumerate(pdf):
        blocks = [b for b in page.get_text('blocks') if len(b) > 4 and b[4].strip() != str(i+1)]
        if blocks:
            bounds = [min(b[0] for b in blocks), min(b[1] for b in blocks),
                      max(b[2] for b in blocks), max(b[3] for b in blocks)]
            page_bounds.append({'page': i+1, 'bounds': bounds})
    # Preserve the single-level package and refresh only the required disclosure PDF.
    support = OUT / 'supporting_materials'
    shutil.copy2(OUT / 'AI工具使用详情.pdf', support / 'AI工具使用详情.pdf')
    assert not any(p.is_dir() for p in support.iterdir()), 'Supporting package is not single-level'
    files = sorted(p for p in support.iterdir() if p.is_file())
    expected_results = {f'result{index}.xlsx' for index in range(1, 5)}
    assert expected_results <= {p.name for p in files}
    canonical_results = ROOT / 'deliverables' / 'supporting_materials' / 'outputs'
    for filename in expected_results:
        assert sha(support / filename) == sha(canonical_results / filename)
    flat_audit = json.loads((support / 'flat_support_audit.json').read_text(encoding='utf-8'))
    assert flat_audit['structure'] == 'single-level'
    assert flat_audit['workbooks_resaved'] is False
    assert set(flat_audit['result_filenames']) == expected_results
    assert flat_audit['python_source_files'] == 25
    assert flat_audit['flat_package_tests'] == '95 passed'
    assert len(list(support.glob('*.py'))) == 25
    with zipfile.ZipFile(OUT / 'A题支撑材料.zip', 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for p in files:
            z.write(p, p.name)
    with zipfile.ZipFile(OUT / 'A题支撑材料.zip') as z:
        assert z.testzip() is None
        assert all('/' not in name and '\\' not in name for name in z.namelist())
        assert set(z.namelist()) == {p.name for p in files}
        for p in files:
            assert hashlib.sha256(z.read(p.name)).hexdigest() == sha(p)
    for filename in ('A题论文_提交预览.pdf', 'A题支撑材料.zip'):
        assert (OUT / filename).stat().st_size < 20_000_000
    def count(text):
        return len(re.sub(r'\s', '', text))
    before_main, after_main = old_source.split('<!-- APPENDIX -->')[0], new_source.split('<!-- APPENDIX -->')[0]
    section_counts = []
    for heading in re.findall(r'(?m)^### .+$', before_main):
        if heading in after_main:
            def body(text):
                return re.split(r'\n##?#? ', text.split(heading, 1)[1], maxsplit=1)[0]
            section_counts.append({'section': heading[4:], 'before_chars': count(body(before_main)),
                                   'after_chars': count(body(after_main))})
    report = {'total_pages': len(pdf), 'abstract_pages': 1, 'body_pages': appendix-2,
              'appendix_first_page': appendix, 'native_editable_equations': equations,
              'editable_three_line_tables': len(tables), 'display_equations_unchanged': len(old_display),
              'table_rows_unchanged': len(old_tables), 'python_listings': len(new_codes),
              'support_python_sources': 25, 'flat_package_tests': '95 passed',
              'support_files': len(files), 'support_subdirectories': 0, 'reference_count': 3,
              'result_workbooks_resaved': False, 'source_main_chars_before': count(before_main),
              'source_main_chars_after': count(after_main), 'section_counts': section_counts,
              'archive_bytes': (OUT/'A题支撑材料.zip').stat().st_size,
              'pdf_bytes': (OUT/'A题论文_提交预览.pdf').stat().st_size,
              'page_text_bounds': page_bounds,
              'sha256': {n: sha(OUT/n) for n in ('A题论文_可编辑版.docx','A题论文_提交预览.pdf','A题支撑材料.zip')}}
    (OUT / 'revision_audit.json').write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ('page_text_bounds','sha256','section_counts')}, ensure_ascii=False))


if __name__ == '__main__':
    {'math': finalize_math, 'audit': audit}[sys.argv[1]]()
