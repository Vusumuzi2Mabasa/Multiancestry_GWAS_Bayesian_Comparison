"""Check research manuscript provenance, numerical tables, equations and links."""
from pathlib import Path
import csv, hashlib, json, re
from docx import Document
root=Path(__file__).resolve().parents[1]
source=(root/'README.md').read_bytes();provenance=json.loads((root/'results/word_report_software.json').read_text())
assert provenance['source_sha256']==hashlib.sha256(source).hexdigest(),'Rebuild Word after changing README'
d=Document(root/'report/Multiancestry_GWAS_Methodology_and_Results.docx')
assert d.core_properties.author=='Vusumuzi Mabasa'
assert len(d.inline_shapes)==4 and len(d.tables)==5
math_count=len(d.element.xpath('.//m:oMath'));assert math_count>=30
texts=[p.text for p in d.paragraphs]
for heading in ['Abstract','Introduction','Methods','Results','Discussion','Conclusion','Reproducibility','References']:
    assert heading in texts,heading
assert any('This is a simulation study' in t for t in texts)
assert any('every phenotype is simulated' in t for t in texts)
for i,p in enumerate(d.paragraphs):
    if p._p.xpath('.//w:drawing'):assert len(d.paragraphs[i-1].text)>100,'Interpretation missing before figure'
with (root/'results/extension/summary.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
lookup={(r['scenario'],r['n_african'],r['reference'],r['method']):r for r in rows}
methods={'Bayes ρ = 0':'Bayes rho 0','Bayes ρ = 0.8':'Bayes rho 0.8','MultiSuSiE':'MultiSuSiE','Manc-COJO':'Manc-COJO'}
assert [len(t.rows) for t in d.tables]==[4,7,17,17,7]
for row in d.tables[2].rows[1:]:
    a=[c.text for c in row.cells];r=lookup[('null',a[0],a[1].lower(),methods[a[2]])]
    assert a[3]==f"{float(r['false_selection_rate']):.3f}"
    assert a[4]==f"{float(r['false_lower']):.3f}–{float(r['false_upper']):.3f}"
for row in d.tables[3].rows[1:]:
    a=[c.text for c in row.cells];r=lookup[(a[0].lower(),a[1],'matched',methods[a[2]])]
    assert a[3:]==[f"{float(r[k]):.3f}" for k in ['mean_recall','rmse_african','mean_fdp']]
links=[r.target_ref for r in d.part.rels.values() if r.reltype.endswith('/hyperlink')]
assert all(x.startswith(('https://','http://','#')) for x in links)
for dest in re.findall(r'\]\(([^)]+)\)',source.decode()):
    if not dest.startswith(('https://','#')):assert (root/dest).is_file(),dest
audit=dict(author='Vusumuzi Mabasa',source_sha256=hashlib.sha256(source).hexdigest(),native_math_objects=math_count,
    figures=4,tables=5,numerical_tables_match=True,interpretations_before_figures=True,simulation_label=True,portable_links=True,
    visual_review='See render audit after page-by-page inspection')
(root/'results/extension/manuscript_validation.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
