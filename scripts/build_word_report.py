"""Render the research README to Word with editable equations and shared content.

Requires Python, python-docx and Pandoc (also distributed with Quarto).
Set PANDOC or pass --pandoc when Pandoc is not on PATH.
"""
from pathlib import Path
import argparse, hashlib, json, os, re, shutil, subprocess, tempfile
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn

ROOT=Path(__file__).resolve().parents[1]

def reference(path):
    d=Document();sec=d.sections[0]
    sec.page_width=Inches(8.27);sec.page_height=Inches(11.69)
    sec.top_margin=sec.bottom_margin=Inches(.7)
    sec.left_margin=sec.right_margin=Inches(.72)
    for style in d.styles:
        if style.type==1:
            style.font.name='Calibri';style.font.size=Pt(11)
            style.font.color.rgb=RGBColor(0,0,0)
            style.paragraph_format.space_after=Pt(5)
            style.paragraph_format.line_spacing=1.04
    for name,size in [('Title',23),('Subtitle',12),('Heading 1',16),('Heading 2',13),('Heading 3',11)]:
        st=d.styles[name];st.font.size=Pt(size)
        st.paragraph_format.keep_with_next=True
        st.paragraph_format.space_before=Pt(10 if name.startswith('Heading') else 0)
    for e in list(d.styles.element.iter(qn('w:pBdr'))):e.getparent().remove(e)
    d.save(path)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--pandoc',default=os.environ.get('PANDOC') or shutil.which('pandoc'))
    args=ap.parse_args()
    if not args.pandoc:ap.error('Pandoc is required: use --pandoc or set PANDOC.')
    source=(ROOT/'README.md').read_text(encoding='utf-8')
    title=source.splitlines()[0].removeprefix('# ')
    text='\n'.join(source.splitlines()[6:])
    text='\n'.join(line for line in text.splitlines() if not line.startswith('[Word manuscript]'))
    # Repository-relative links must also work when Word is downloaded alone.
    text=re.sub(r'(?<!!)\[([^]]+)\]\((?!https?://|#)([^)]+)\)',
        r'[\1](https://github.com/Vusumuzi2Mabasa/Multiancestry_GWAS_Bayesian_Comparison/blob/main/\2)',text)
    # Standard Markdown plus TeX math lets Pandoc emit native Office Math.
    text=re.sub(r'!\[([^]]+)\]\(([^)]+)\)',r'![\1](\2){width=5.9in}',text)
    out=ROOT/'report/Multiancestry_GWAS_Methodology_and_Results.docx'
    out.parent.mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory() as temp:
        temp=Path(temp);ref=temp/'reference.docx';reference(ref)
        md=temp/'manuscript.md';md.write_text(text,encoding='utf-8')
        subprocess.run([args.pandoc,str(md),'--from=markdown-implicit_figures+tex_math_dollars','--to=docx','--shift-heading-level-by=-1',
            '--resource-path='+str(ROOT),'--reference-doc='+str(ref),
            '--metadata=title:'+title,'--metadata=author:Vusumuzi Mabasa',
            '--metadata=subtitle:A paired comparison of Manc-COJO and Bayesian hierarchical variable selection',
            '--output='+str(out)],cwd=ROOT,check=True)
    d=Document(out);d.core_properties.author='Vusumuzi Mabasa';d.core_properties.title=title
    for p in d.paragraphs:
        if p.style.name=='Heading 1' and p.text in ['Methods','Results','Discussion','Reproducibility']:p.paragraph_format.page_break_before=True
        if p.style.name.startswith('Heading'):p.paragraph_format.keep_with_next=True
        if p._p.xpath('.//w:drawing'):p.paragraph_format.keep_with_next=True
        if p.text.startswith(('Figure ','Table ')):
            p.paragraph_format.space_after=Pt(8)
        if p.style.name=='Heading 2' and p.text in [
            'The benefit of borrowing depended on effect similarity']:
            p.paragraph_format.page_break_before=True
        if p.text.startswith('Strong borrowing attenuated the AFR-specific effect'):
            p.paragraph_format.page_break_before=True
    for table in d.tables:
        table.autofit=True
        borders=OxmlElement('w:tblBorders')
        for side in ['top','left','bottom','right','insideH','insideV']:
            e=OxmlElement('w:'+side)
            for k,v in [('val','single'),('sz','4'),('color','D9D9D9')]:e.set(qn('w:'+k),v)
            borders.append(e)
        table._tbl.tblPr.append(borders)
        for i,row in enumerate(table.rows):
            props=row._tr.get_or_add_trPr();props.append(OxmlElement('w:cantSplit'))
            if i==0:props.append(OxmlElement('w:tblHeader'))
            for cell in row.cells:
                cell.vertical_alignment=1
                if i==0:
                    shade=OxmlElement('w:shd');shade.set(qn('w:fill'),'E8EFF4');cell._tc.get_or_add_tcPr().append(shade)
                for p in cell.paragraphs:
                    p.paragraph_format.space_before=Pt(3);p.paragraph_format.space_after=Pt(3)
                    for r in p.runs:r.font.size=Pt(9);r.bold=(i==0)
    for e in list(d.styles.element.iter(qn('w:pBdr'))):e.getparent().remove(e)
    d.save(out)
    import docx
    (ROOT/'results/word_report_software.json').write_text(json.dumps({
        'builder':'scripts/build_word_report.py','source':'README.md',
        'source_sha256':hashlib.sha256(source.encode()).hexdigest(),
        'pandoc':subprocess.check_output([args.pandoc,'--version'],text=True).splitlines()[0],
        'python_docx':docx.__version__, 'author':'Vusumuzi Mabasa',
        'equations':'Native Office Math (OMML)','analyses_rerun':False},indent=2)+'\n',encoding='utf-8')
    print(out)

if __name__=='__main__':main()
