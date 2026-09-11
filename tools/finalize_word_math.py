"""Replace generated math markers with editable native Word OMML equations."""
from pathlib import Path
import json
import zipfile
import copy
from lxml import etree
from latex2mathml.converter import convert

root = Path(__file__).resolve().parents[1]
transform = etree.XSLT(etree.parse('C:/Program Files/Microsoft Office/root/Office16/MML2OMML.XSL'))
ns = {'w': 'http://schemas.openxmlformats.org/wordprocessingml/2006/main'}
equations = json.loads((root/'deliverables/equations.json').read_text(encoding='utf-8'))
mapping = {}
for item in equations:
    latex = item['latex'].replace(r'\left.','').replace(r'\right|',r'\vert')
    mathml = etree.fromstring(convert(latex).encode())
    node = transform(mathml).getroot()
    for run in node.xpath('//*[local-name()="r"]'):
        props = etree.Element('{'+ns['w']+'}rPr')
        size = etree.SubElement(props, '{'+ns['w']+'}sz')
        size.set('{'+ns['w']+'}val', '22')
        run.insert(0, props)
    mapping[item['key']] = node
for filename in ('A题论文_可编辑版.docx', 'AI工具使用详情_可编辑版.docx'):
    dest = root/'deliverables'/filename
    with zipfile.ZipFile(dest) as z:
        content = {name:z.read(name) for name in z.namelist()}
    document = etree.fromstring(content['word/document.xml'])
    count = 0
    for t in document.xpath('//w:t',namespaces=ns):
        if t.text in mapping:
            run = t.getparent()
            run.getparent().replace(run,copy.deepcopy(mapping[t.text]))
            count += 1
    content['word/document.xml'] = etree.tostring(document,xml_declaration=True,encoding='UTF-8',standalone=True)
    assert b'MATHPLACEHOLDER' not in content['word/document.xml']
    with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for name,data in content.items(): z.writestr(name,data)
    print(filename, 'native equations:',count)
