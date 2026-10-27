"""Independent checks of frozen thresholds, decisions, effects and output completeness."""
from pathlib import Path
import sys,csv,json,hashlib,math,argparse
from collections import defaultdict
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tmp/runtime'))
import numpy as np
def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--cell',type=int);args=ap.parse_args()
    cell_ids=[args.cell] if args.cell is not None else list(range(16))
    out=ROOT/'results/extension';design=json.loads((ROOT/'extension/design.json').read_text())
    expected=hashlib.sha256((ROOT/'extension/design.json').read_bytes()+(ROOT/'extension/benchmark.py').read_bytes()).hexdigest()
    counts=defaultdict(int)
    for cell in cell_ids:
        p=out/'cells'/f'cell_{cell:02d}';done=json.loads((p/'complete.json').read_text());t=json.loads((p/'thresholds.json').read_text())
        assert done['design_code_sha256']==expected==t['design_code_sha256']
        assert done['threshold_sha256']==hashlib.sha256((p/'thresholds.json').read_bytes()).hexdigest()
        cal=read(p/'calibration_scores.tsv');metrics=read(p/'metrics.tsv');variants=read(p/'variant_results.tsv')
        by=defaultdict(list)
        for row in cal:by[row['method']].append(float(row['score']))
        assert len(by)==4
        for method,scores in by.items():
            assert len(scores)==300
            h=math.ceil((len(scores)+1)*.95)
            assert sorted(scores)[h-1]==t['thresholds'][method]
        vby=defaultdict(list)
        for v in variants:vby[(v['scenario'],v['replicate'],v['method'])].append(v)
        assert len(metrics)==1600 and len(vby)==1600
        seen=set()
        for row in metrics:
            key=(row['scenario'],row['replicate'],row['method']);assert key not in seen;seen.add(key)
            v=vby[key];assert len(v)==8 and len({x['index'] for x in v})==8
            selected=np.array([int(x['selected']) for x in v]);truth=np.array([float(x['true_african'])!=0 or float(x['true_ceu'])!=0 for x in v])
            assert selected.sum()==int(row['selected'])
            assert np.sum(selected*truth)==int(row['tp'])
            assert np.sum(selected*(~truth))==int(row['fp'])
            err=np.array([float(x['estimate_african'])-float(x['true_african']) for x in v])
            assert math.isclose(np.mean(err**2),float(row['mse_african']),rel_tol=1e-10,abs_tol=1e-12)
            if row['method']!='Manc-COJO':
                pip=np.array([float(x['pip']) for x in v]);assert np.array_equal(selected,pip>t['thresholds'][row['method']])
            elif selected.any():
                assert float(row['meta_entry_score'])>=t['thresholds']['Manc-COJO']-1e-8
            counts[row['scenario']]+=1
        for scenario,n in [('null',200),('shared',100),('heterogeneous',100)]:
            for method in by:
                assert {int(r['replicate']) for r in metrics if r['scenario']==scenario and r['method']==method}==set(range(n))
    # Distinct phase namespaces generate distinct RNG states, including repeated numeric replicate IDs.
    states=[tuple(np.random.SeedSequence([design['seed'],cell,phase,rep]).generate_state(4))
        for cell in cell_ids for phase,n in enumerate([300,200,100,100]) for rep in range(n)]
    assert len(states)==len(set(states))==700*len(cell_ids)
    report=dict(all_passed=True,cells=len(cell_ids),unique_dataset_random_states=len(states),method_fit_counts=dict(counts),
        frozen_threshold_hashes=True,quantiles_recomputed=True,variant_decisions_recomputed=True,
        effect_mse_recomputed=True,manc_entry_implication_checked=True,all_replicates_present=True)
    destination=out/'independent_audit.json' if args.cell is None else ROOT/f'tmp/audit_cell_{args.cell}.json'
    destination.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
