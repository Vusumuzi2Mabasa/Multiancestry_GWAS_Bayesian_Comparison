"""Validate completed paired outputs and write a reproducibility manifest."""
from pathlib import Path
import csv,hashlib,json,platform,sys
import numpy as np
from simulate_compare import ROOT,IDS

config=json.loads((ROOT/'results/config.json').read_text())
assert (ROOT/'results/status.txt').read_text().strip()=='COMPLETE_WITH_MANC'
with (ROOT/'results/replicate_metrics.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
expected=config['replicates']*len(config['scenarios'])*len(config['n_afr'])
methods={x['method'] for x in r};assert len(methods)==5
for method in methods:
    z=[x for x in r if x['method']==method]
    assert len(z)==expected,(method,len(z))
    assert len({(x['scenario'],x['n_afr'],x['replicate']) for x in z})==expected
    assert all(np.isfinite(float(x['mse_afr'])) and int(x['selected'])<=len(IDS) for x in z)
assert len(list((ROOT/'data').rglob('manc.command.json')))==expected
assert len(list((ROOT/'data').rglob('manc.stdout.txt')))==expected
for p in (ROOT/'data').rglob('manc.jma.cojo'):
    with p.open() as f:
        rows=list(csv.DictReader(f,delimiter='\t'))
    for x in rows:
        assert x['A1']=='A' and float(x['pJ.ma'])<=config['cojo_p']*1.0001,x
with (ROOT/'results/variant_results.tsv').open() as f:v=list(csv.DictReader(f,delimiter='\t'))
assert len(v)==expected*2*len(IDS)
assert all(0<=float(x['pip'])<=1+1e-12 for x in v)
environment=dict(python=sys.version,numpy=np.__version__,platform=platform.platform(),
    executable=sys.executable,paired_datasets=expected,methods=sorted(methods),
    checks='Complete paired rows, finite MSE, PIP range, actual command logs, allele orientation and selected joint p-values')
(ROOT/'results/environment_and_audit.json').write_text(json.dumps(environment,indent=2))
files=[p for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and '.git' not in p.parts and p.name!='checksums.sha256']
(ROOT/'checksums.sha256').write_text(''.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(ROOT).as_posix()+'\n' for p in sorted(files)))
print(f'PASS: {expected} datasets x {len(methods)} method/threshold combinations; {len(v)} Bayesian variant rows; provenance and checksums recorded.')
