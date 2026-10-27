"""Record software, input and design identities without machine paths or donor IDs."""
from pathlib import Path
import sys, os, json, hashlib, platform, importlib.metadata as meta
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tmp/runtime'))
def main():
    config=json.loads((ROOT/'extension/design.json').read_text())
    data=dict(python=platform.python_version(),platform=platform.system(),
        packages={x:meta.version(x) for x in ['numpy','scipy','numba','tqdm','matplotlib','MultiSuSiE']},
        multisusie_commit=config['multisusie_commit'],
        manc_binary_sha256=hashlib.sha256(Path(os.environ.get('MANC_COJO',str(ROOT/('tools/manc_cojo_win.exe' if os.name=='nt' else 'tools/manc_cojo')))).read_bytes()).hexdigest(),
        file_sha256={str(p.relative_to(ROOT)).replace('\\','/'):hashlib.sha256(p.read_bytes()).hexdigest()
            for p in [ROOT/'extension/design.json',ROOT/'extension/benchmark.py',ROOT/'scripts/simulate_compare.py',ROOT/'extension/prepare_reference.py',ROOT/'extension/reference_downloads.json']},
        source_permission='1000 Genomes Phase 3 no embargo after final publication; cite consortium',
        source_permission_url='https://www.internationalgenome.org/IGSR_disclaimer/',
        references=[dict(method='MultiSuSiE',doi='10.1038/s41588-025-02450-5',source='https://github.com/jordanero/MultiSuSiE',license='GPL-3.0'),
            dict(method='Manc-COJO',doi='10.1038/s41588-026-02799-1',source='https://github.com/light156/multi-ancestry-COJO',version='1.1.0'),
            dict(data='1000 Genomes Phase 3',doi='10.1038/nature15393')])
    (ROOT/'results/extension/software_provenance.json').write_text(json.dumps(data,indent=2)+'\n')
    print(json.dumps(data['packages'],indent=2))
if __name__=='__main__':main()
