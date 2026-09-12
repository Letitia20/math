import fs from 'node:fs';
import path from 'node:path';
import {createRequire} from 'node:module';
import {fileURLToPath} from 'node:url';
const out=path.dirname(fileURLToPath(import.meta.url));
const root=path.resolve(out,'../..');
process.chdir(root);
const require=createRequire(import.meta.url);
const {AlignmentType:A,BorderStyle,Document,Footer,HeadingLevel,ImageRun,PageBreak,Packer,Paragraph,Table,TableCell,TableRow,TextRun,WidthType,PageNumber,TabStopType}=require('docx');
const math=[],font={ascii:'Times New Roman',hAnsi:'Times New Roman',eastAsia:'宋体',cs:'Times New Roman'};
function inline(s,size=24){return s.split(/(\$[^$]+\$)/g).filter(Boolean).map(x=>{if(x.startsWith('$')&&x.endsWith('$')){const key=`MATHPLACEHOLDER${math.length}END`;math.push({key,latex:x.slice(1,-1),display:false,size});return new TextRun({text:key,size,font});}return new TextRun({text:x.replace(/\*\*/g,'').replace(/`/g,''),size,font});});}
function p(s,opts={}){return new Paragraph({alignment:A.JUSTIFIED,indent:{firstLine:480},widowControl:true,spacing:{line:300,after:70},children:inline(s),...opts});}
function h(s,level=1){return new Paragraph({heading:level===1?HeadingLevel.HEADING_1:HeadingLevel.HEADING_2,keepNext:true,keepLines:true,indent:{firstLine:0},spacing:{before:200,after:100},children:[new TextRun({text:s,font:{...font,eastAsia:'黑体'},color:'000000',bold:true,size:level===1?28:24})]});}
const none={style:BorderStyle.NONE,size:0,color:'FFFFFF'};
const rule=(size)=>({style:BorderStyle.SINGLE,size,color:'000000'});
function table(lines){
 const rows=lines.filter(l=>!/^\|[\s:|-]+\|$/.test(l)).map(l=>l.trim().slice(1,-1).split('|').map(x=>x.trim()));
 const width=9000,n=rows[0].length;
 let widths=Array(n).fill(Math.floor(width/n)); widths[0]+=width-widths.reduce((a,b)=>a+b,0);
 if(rows[0][0]==='符号')widths=[1950,4200,2850];
 else if(n===2)widths=[6750,2250];
 else if(n===3)widths=[2000,4200,2800];
 return new Table({width:{size:width,type:WidthType.DXA},columnWidths:widths,borders:{top:none,bottom:none,left:none,right:none,insideHorizontal:none,insideVertical:none},rows:rows.map((row,i)=>new TableRow({tableHeader:i===0,cantSplit:true,children:row.map((c,j)=>new TableCell({width:{size:widths[j],type:WidthType.DXA},borders:{top:i===0?rule(12):none,bottom:i===0?rule(6):i===rows.length-1?rule(12):none,left:none,right:none},margins:{top:65,bottom:65,left:75,right:75},children:[p(c,{indent:{firstLine:0},keepNext:i<rows.length-1,alignment:n===2&&j===0?A.LEFT:A.CENTER,spacing:{line:250,after:0},children:inline(c,21)})]}))}))});
}
function walk(dir){return fs.readdirSync(dir,{withFileTypes:true}).flatMap(e=>e.isDirectory()?(['__pycache__','.pytest_cache'].includes(e.name)?[]:walk(path.join(dir,e.name))):(e.name.endsWith('.pyc')?[]:[path.join(dir,e.name)]));}
let source=fs.readFileSync(path.join(out,'paper_source.md'),'utf8');
for(let n=1;n<=6;n++){const md=fs.readFileSync(`reports/problem${[1,1,2,2,3,4][n-1]}.md`,'utf8').replaceAll('\r','');const block=md.slice(md.indexOf(`**表 ${n}`)).split('\n\n');if(!block[1]?.trim().startsWith('|'))throw Error(`Missing table ${n}`);source=source.replace(`{{TABLE${n}}}`,block[0].replace(/\*\*/g,'')+'\n\n'+block[1]);}
source=source.replace('{{ENDPOINT3}}','第三问临界时刻为 206906.2447 s（57.4740 h），首个严格达标整分钟为 206940 s（57.4833 h）。').replace('{{ENDPOINT4}}','第四问临界时刻为 183931.0374 s（51.0920 h），首个严格达标整分钟为 183960 s（51.1000 h）。');
const support=path.join(out,'supporting_materials');
const files=walk(support).map(f=>path.relative(support,f).replaceAll('\\','/')).sort();
source=source.replace('{{FILELIST}}',files.map(f=>'`'+f+'`').join('\n\n'));
const codes=files.filter(f=>/^problem[1-4]\.py$/.test(f));
source=source.replace('{{SOURCECODE}}',codes.map(f=>'### '+f+'\n\n```python\n'+fs.readFileSync(path.join(support,f),'utf8')+'\n```').join('\n\n'));
fs.writeFileSync(path.join(out,'paper_complete.md'),source);
function parse(text){const children=[],lines=text.split('\n');let code=false;for(let i=0;i<lines.length;i++){const ln=lines[i],s=ln.trim();
 if(s.startsWith('```')){code=!code;continue;}
 if(code){children.push(new Paragraph({widowControl:false,spacing:{line:180,lineRule:'exact',after:0},children:[new TextRun({text:ln||' ',font:{ascii:'Consolas',hAnsi:'Consolas',eastAsia:'宋体'},size:14})]}));continue;}
 if(!s)continue;
 if(s==='<!-- PAGEBREAK -->'||s==='<!-- APPENDIX -->'){children.push(new Paragraph({children:[new PageBreak()],spacing:{after:0}}));continue;}
 if(s.startsWith('<!--'))continue;
 if(s==='$$'){const formula=[];while(++i<lines.length&&lines[i].trim()!=='$$')formula.push(lines[i]);let latex=formula.join(' ');const number=latex.match(/\\qquad\s*\((\d+)\)\s*$/);if(number)latex=latex.slice(0,number.index).trim().replace(/\.$/,'');if(number&&['5','7','12','17'].includes(number[1]))latex='\\begin{aligned}'+latex.replace(/,?\\qquad\s*/g,' \\\\ ')+'\\end{aligned}';const key=`MATHPLACEHOLDER${math.length}END`;math.push({key,latex,display:true,size:23,number:number?.[1]});children.push(p(key,{indent:{firstLine:0},alignment:A.LEFT,keepLines:true,tabStops:[{type:TabStopType.CENTER,position:4450},{type:TabStopType.RIGHT,position:9000}],spacing:{line:340,before:100,after:100},children:[new TextRun({text:'\t'}),new TextRun({text:key,font:'Cambria Math',size:23}),new TextRun({text:'\t'+(number?'('+number[1]+')':''),font,size:22})]}));continue;}
 if(s.startsWith('|')){const block=[s];while(i+1<lines.length&&lines[i+1].trim().startsWith('|'))block.push(lines[++i].trim());children.push(table(block));continue;}
 if(s.startsWith('# ')){children.push(p(s.slice(2),{indent:{firstLine:0},alignment:A.CENTER,keepNext:true,spacing:{before:0,after:180},children:[new TextRun({text:s.slice(2),font:'黑体',bold:true,size:32})]}));continue;}
 if(s.startsWith('## ')){children.push(s==='## 摘要'?p('摘  要',{indent:{firstLine:0},alignment:A.CENTER,keepNext:true,children:[new TextRun({text:'摘  要',size:26,font:'黑体',bold:true})]}):h(s.slice(3)));continue;}
 if(s.startsWith('### ')){children.push(h(s.slice(4),2));continue;}
 const img=s.match(/^!\[([^\]]+)\]\(([^)]+)\)/);
 if(img){const f=path.join(support,path.basename(img[2]));if(!fs.existsSync(f))throw Error(f);const b=fs.readFileSync(f),w=b.readUInt32BE(16),ht=b.readUInt32BE(20);const dw=Math.min(565,310*w/ht);children.push(new Paragraph({alignment:A.CENTER,keepNext:true,spacing:{before:60,after:50},children:[new ImageRun({type:'png',data:b,transformation:{width:dw,height:dw*ht/w}})]}));children.push(p(img[1],{indent:{firstLine:0},alignment:A.CENTER,keepLines:true,spacing:{line:260,after:100},children:inline(img[1],21)}));continue;}
 if(/^\[\d+\]/.test(s)){children.push(p(s,{indent:{left:300,hanging:300},alignment:A.LEFT,spacing:{line:260,after:60},keepLines:true,children:inline(s,21)}));continue;}
 if(/^表 [0-9]/.test(s)&&lines.slice(i+1).find(x=>x.trim())?.trim().startsWith('|')){children.push(p(s,{indent:{firstLine:0},alignment:A.CENTER,keepNext:true,spacing:{before:100,after:70},children:inline(s,21)}));continue;}
 if(s.startsWith('`')&&s.endsWith('`')){children.push(p(s,{indent:{firstLine:0},alignment:A.LEFT,spacing:{line:240,after:0},children:inline(s,19)}));continue;}
 if(s.startsWith('关键词：')){children.push(p(s,{indent:{firstLine:0},children:inline(s,22)}));continue;}
 children.push(p(s));}return children;}
async function build(text,file){const doc=new Document({creator:'',lastModifiedBy:'',title:'药材烘干的传热传质模型',styles:{default:{document:{run:{font,size:24},paragraph:{spacing:{line:300,after:70}}}}},sections:[{properties:{page:{size:{width:11906,height:16838},margin:{top:1418,bottom:1418,left:1418,right:1500}}},footers:{default:new Footer({children:[new Paragraph({alignment:A.CENTER,children:[new TextRun({children:[PageNumber.CURRENT],size:20,font})]})]})},children:parse(text)}]});fs.writeFileSync(file,await Packer.toBuffer(doc));}
await build(source,path.join(out,'A题论文_可编辑版.docx'));
await build(fs.readFileSync(path.join(out,'AI工具使用详情.md'),'utf8'),path.join(out,'AI工具使用详情_可编辑版.docx'));
fs.writeFileSync(path.join(out,'equations.json'),JSON.stringify(math,null,2));
console.log(JSON.stringify({supportFiles:files.length,sourceFiles:codes.length,equations:math.length}));
