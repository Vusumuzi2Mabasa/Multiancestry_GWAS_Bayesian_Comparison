"""Check narrative provenance, editable equations and tables against saved data."""
from pathlib import Path
import csv, hashlib, json, math, re
from docx import Document

root=Path(__file__).resolve().parents[1]
source=(root/'README.md').read_bytes()
provenance=json.loads((root/'results/word_report_software.json').read_text())
assert provenance['source_sha256']==hashlib.sha256(source).hexdigest(),'Rebuild Word after changing README'
d=Document(root/'report/Multiancestry_GWAS_Methodology_and_Results.docx')
assert d.core_properties.author=='Vusumuzi Mabasa'
assert len(d.inline_shapes)==4 and len(d.tables)==5
math_count=len(d.element.xpath('.//m:oMath'));assert math_count>=20
texts=[p.text for p in d.paragraphs]
for heading in ['Abstract','Introduction','Methods','Results','Discussion','Conclusion','Reproducibility','References']:
    assert heading in texts,heading
assert any('All participants, genotypes and phenotypes are synthetic' in t for t in texts)
for i,p in enumerate(d.paragraphs):
    if p._p.xpath('.//w:drawing'):
        assert len(d.paragraphs[i-1].text)>100,'Missing interpretation before figure'
def rows(name):
    with (root/'results'/name).open() as f:return list(csv.DictReader(f,delimiter='\t'))
lookup={(r['scenario'],int(r['n_afr']),r['method']):r for r in rows('summary.tsv')}
methods={'Manc-COJO':'Manc_COJO','Bayesian ρ = 0':'Bayes_rho0_PIP0.95','Bayesian ρ = 0.8':'Bayes_rho0.8_PIP0.95'}
table=d.tables[2];assert len(table.rows)==13
for row in table.rows[1:]:
    a=[c.text for c in row.cells];r=lookup[(a[0].lower(),int(a[1].replace(',','')),methods[a[2]])]
    assert a[3:]==[f'{v:.3f}' for v in [float(r['recall_mean']),math.sqrt(float(r['mse_afr_mean'])),float(r['fdp_mean'])]]
null={(r['method'],int(r['n_afr'])):int(r['any_false_positive']) for r in rows('null_event_counts.tsv')}
nullmethods=['Manc_COJO','Bayes_rho0_PIP0.95','Bayes_rho0.8_PIP0.95','Bayes_rho0_PIP0.5','Bayes_rho0.8_PIP0.5']
for row,method in zip(d.tables[3].rows[1:],nullmethods):
    assert [c.text for c in row.cells[1:]]==[f'{null[method,n]}/20' for n in [250,500,1000,2000]]
links=[r.target_ref for r in d.part.rels.values() if r.reltype.endswith('/hyperlink')]
assert all(x.startswith(('https://','http://','#')) for x in links)
for dest in re.findall(r'\]\(([^)]+)\)',source.decode()):
    if not dest.startswith(('https://','#')):assert (root/dest).is_file(),dest
print(f'PASS: author, source hash, sections, synthetic label, {math_count} native math objects, four interpreted figures, five tables, result values and portable links.')
