import fs from 'node:fs';
import path from 'node:path';
import {createRequire} from 'node:module';
const require=createRequire(import.meta.url);
const {AlignmentType:A,Document,Footer,HeadingLevel,ImageRun,PageBreak,Packer,Paragraph,Table,TableCell,TableRow,TextRun,WidthType,PageNumber}=require('docx');
const math=[],font='宋体';
function inline(s,size=22){return s.split(/(\$[^$]+\$)/g).filter(Boolean).map(x=>{if(x.startsWith('$')&&x.endsWith('$')){const key=`MATHPLACEHOLDER${math.length}END`;math.push({key,latex:x.slice(1,-1),display:false});return new TextRun({text:key,size,font});}return new TextRun({text:x.replace(/\*\*/g,'').replace(/`/g,''),size,font});});}
function p(s,opts={}){return new Paragraph({alignment:A.JUSTIFIED,spacing:{line:320,after:80},children:inline(s),...opts});}
function h(s,level=1){return new Paragraph({heading:level===1?HeadingLevel.HEADING_1:HeadingLevel.HEADING_2,keepNext:true,spacing:{before:160,after:90},children:[new TextRun({text:s,font:'黑体',color:'000000',bold:true,size:level===1?27:23})]});}
function table(lines){const rows=lines.filter(l=>!/^\|[\s:|-]+\|$/.test(l)).map(l=>l.trim().slice(1,-1).split('|').map(x=>x.trim()));const width=9000,widths=rows[0].map((_,i)=>Math.floor(width/rows[0].length)+(i===0?width%rows[0].length:0));return new Table({width:{size:width,type:WidthType.DXA},columnWidths:widths,rows:rows.map((row,i)=>new TableRow({tableHeader:i===0,cantSplit:true,children:row.map((c,j)=>new TableCell({width:{size:widths[j],type:WidthType.DXA},margins:{top:60,bottom:60,left:60,right:60},children:[p(c,{alignment:A.CENTER,spacing:{line:260,after:0},children:inline(c,19)})]}))}))});}
function walk(dir){return fs.readdirSync(dir,{withFileTypes:true}).flatMap(e=>e.isDirectory()?(['__pycache__','.pytest_cache'].includes(e.name)?[]:walk(path.join(dir,e.name))):(e.name.endsWith('.pyc')?[]:[path.join(dir,e.name)]));}
let source=fs.readFileSync('deliverables/paper_source.md','utf8');
for(let n=1;n<=6;n++){const md=fs.readFileSync(`reports/problem${[1,1,2,2,3,4][n-1]}.md`,'utf8').replaceAll('\r','');const block=md.slice(md.indexOf(`**表 ${n}`)).split('\n\n');if(!block[1]?.trim().startsWith('|'))throw Error(`Missing table ${n}`);source=source.replace(`{{TABLE${n}}}`,block[0].replace(/\*\*/g,'')+'\n\n'+block[1]);}
source=source.replace('{{ENDPOINT3}}','第三问临界时刻为 206906.2447 s（57.4740 h），首个严格达标整分钟为 206940 s（57.4833 h）。').replace('{{ENDPOINT4}}','第四问临界时刻为 183931.0374 s（51.0920 h），首个严格达标整分钟为 183960 s（51.1000 h）。');
const support='deliverables/supporting_materials';
const files=walk(support).map(f=>path.relative(support,f).replaceAll('\\','/')).sort();
source=source.replace('{{FILELIST}}',files.map(f=>'`'+f+'`').join('\n\n'));
const codes=files.filter(f=>/^(src|scripts|tools|tests)\//.test(f)&&f.endsWith('.py'));
source=source.replace('{{SOURCECODE}}',codes.map(f=>'### '+f+'\n\n```python\n'+fs.readFileSync(path.join(support,f),'utf8')+'\n```').join('\n\n'));
fs.writeFileSync('deliverables/paper_complete.md',source);
function parse(text){const children=[],lines=text.split('\n');let code=false;for(let i=0;i<lines.length;i++){const ln=lines[i],s=ln.trim();
 if(s.startsWith('```')){code=!code;continue;}
 if(code){children.push(new Paragraph({spacing:{line:200,lineRule:'exact',after:0},children:[new TextRun({text:ln||' ',font:'Consolas',size:16})]}));continue;}
 if(!s)continue;
 if(s==='<!-- PAGEBREAK -->'||s==='<!-- APPENDIX -->'){children.push(new Paragraph({children:[new PageBreak()],spacing:{after:0}}));continue;}
 if(s.startsWith('<!--'))continue;
 if(s==='$$'){const formula=[];while(++i<lines.length&&lines[i].trim()!=='$$')formula.push(lines[i]);const key=`MATHPLACEHOLDER${math.length}END`;math.push({key,latex:formula.join(' '),display:true});children.push(p(key,{alignment:A.CENTER,spacing:{line:360,before:100,after:100},children:[new TextRun({text:key,font:'Cambria Math',size:22})]}));continue;}
 if(s.startsWith('|')){const block=[s];while(i+1<lines.length&&lines[i+1].trim().startsWith('|'))block.push(lines[++i].trim());children.push(table(block));continue;}
 if(s.startsWith('# ')){children.push(p(s.slice(2),{alignment:A.CENTER,keepNext:true,spacing:{before:0,after:140},children:[new TextRun({text:s.slice(2),font:'黑体',bold:true,size:30})]}));continue;}
 if(s.startsWith('## ')){children.push(s==='## 摘要'?p('摘  要',{alignment:A.CENTER,keepNext:true,children:[new TextRun({text:'摘  要',size:26,font:'黑体',bold:true})]}):h(s.slice(3)));continue;}
 if(s.startsWith('### ')){children.push(h(s.slice(4),2));continue;}
 const img=s.match(/^!\[([^\]]+)\]\(([^)]+)\)/);
 if(img){const f=path.join('reports/figures',path.basename(img[2]));if(!fs.existsSync(f))throw Error(f);const b=fs.readFileSync(f),w=b.readUInt32BE(16),ht=b.readUInt32BE(20);children.push(new Paragraph({alignment:A.CENTER,keepNext:true,children:[new ImageRun({type:'png',data:b,transformation:{width:550,height:550*ht/w}})]}));children.push(p(img[1],{alignment:A.CENTER,children:inline(img[1],19)}));continue;}
 if(/^\[\d+\]/.test(s)){children.push(p(s,{alignment:A.LEFT,spacing:{line:250,after:60},keepLines:true,children:inline(s,20)}));continue;}
 children.push(p(s,{keepNext:/^表 [0-9]/.test(s)}));}return children;}
async function build(text,file){const doc=new Document({creator:'',lastModifiedBy:'',title:'药材烘干模型',styles:{default:{document:{run:{font,size:22},paragraph:{spacing:{line:320,after:80}}}}},sections:[{properties:{page:{size:{width:11906,height:16838},margin:{top:1418,bottom:1418,left:1418,right:1418}}},footers:{default:new Footer({children:[new Paragraph({alignment:A.CENTER,children:[new TextRun({children:[PageNumber.CURRENT],size:18,font})]})]})},children:parse(text)}]});fs.writeFileSync(file,await Packer.toBuffer(doc));}
await build(source,'deliverables/A题论文_可编辑版.docx');
await build(fs.readFileSync('deliverables/AI工具使用详情.md','utf8'),'deliverables/AI工具使用详情_可编辑版.docx');
fs.writeFileSync('deliverables/equations.json',JSON.stringify(math,null,2));
console.log(JSON.stringify({supportFiles:files.length,sourceFiles:codes.length,equations:math.length}));
