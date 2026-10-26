"""Add real Manc-COJO results to an existing simulation; no data regeneration."""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from simulate_compare import ROOT,IDS,run_cojo,metrics,table

ap=argparse.ArgumentParser();ap.add_argument('--manc',type=Path,required=True);a=ap.parse_args()
binary=a.manc.resolve();assert binary.is_file(),binary
config=json.loads((ROOT/'results/config.json').read_text())
config['manc_binary']=str(binary);config['manc_sha256']=hashlib.sha256(binary.read_bytes()).hexdigest()
(ROOT/'results/config.json').write_text(json.dumps(config,indent=2))
with (ROOT/'results/replicate_metrics.tsv').open() as f:
    rows=[r for r in csv.DictReader(f,delimiter='\t') if r['method']!='Manc_COJO']
(ROOT/'results/status.txt').write_text('MANC_RUNNING\n')
for s in config['scenarios']:
    for rep in range(1,config['replicates']+1):
        for n in config['n_afr']:
            d=ROOT/'data'/s/f'rep_{rep:03d}'/f'AFR_{n}'
            truth=[]
            for cohort in ['EUR','AFR']:
                with (d/f'{cohort}.truth.tsv').open() as f:
                    t=list(csv.DictReader(f,delimiter='\t'))
                assert [r['SNP'] for r in t]==IDS
                truth.append([float(r['beta']) for r in t])
            selected,effects,elapsed=run_cojo(binary,d,config['cojo_p'])
            rows.append(metrics('Manc_COJO',selected,effects,np.array(truth),s,n,rep,elapsed))
            table(ROOT/'results/replicate_metrics.tsv',rows)
        print(f'Manc-COJO {s} {rep}/{config["replicates"]} done',flush=True)
(ROOT/'results/status.txt').write_text('COMPLETE_WITH_MANC\n')
print('COMPLETE_WITH_MANC')
