"""Independent audit and execution gate for the subsequent INLA extension."""
from engine import *
import json,csv,hashlib,subprocess,datetime
from manc import table
def main():
    cfg=json.loads((ROOT/'revision/protocol.json').read_text());details=[];identities=set();total=0;cal_nonconverged=0
    source_files=['revision/protocol.json','revision/engine.py','revision/manc.py','revision/benchmark.py','revision/prepare_reference.py']
    expected_hash=hashlib.sha256(b''.join((ROOT/p).read_bytes() for p in source_files)).hexdigest()
    frozen=subprocess.check_output(['git','show','3d414e2:revision/protocol.json'],cwd=ROOT)
    assert frozen==(ROOT/'revision/protocol.json').read_bytes(),'Protocol differs from pre-simulation commit'
    assert json.loads((ROOT/'results/revision/engine_tests.json').read_text())['passed']
    refcheck=json.loads((ROOT/'results/revision/reference/validation.json').read_text());assert refcheck['all_passed']
    for entry in refcheck['arrays']:
        arr=np.load(ROOT/f'data/revision/reference/region_{entry["region"]}.npz')[entry['pool']]
        assert hashlib.sha256(arr.tobytes()).hexdigest()==entry['sha256']
    assert cfg['sensitivity']['seed']!=cfg['validation']['seed']
    for phase in ['sensitivity','validation']:
        regions=cfg['sensitivity_regions'] if phase=='sensitivity' else cfg['fresh_validation_regions']
        for region in regions:
            for pop in ['YRI','LWK']:
                path=ROOT/f'results/revision/{phase}/r{region}_{pop}';meta=json.loads((path/'complete.json').read_text());th=json.loads((path/'thresholds.json').read_text())
                assert meta['threshold_sha256']==hashlib.sha256((path/'thresholds.json').read_bytes()).hexdigest()
                assert meta['code_sha256']==th['code_sha256']==expected_hash
                with (path/'calibration_scores.tsv').open() as f:cal=list(csv.DictReader(f,delimiter='\t'))
                with (path/'metrics.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
                count=sum(cfg[phase][s] for s in ['null','shared','heterogeneous','one','three'])
                assert len(rows)==count*len(th['arms'])==meta['all_metrics']
                for a in th['arms']:
                    scores=[float(r['score']) for r in cal if int(r['arm'])==a['arm']];assert len(scores)==cfg[phase]['calibration']
                    assert np.isfinite(scores).all()
                    cutoff=sorted(scores)[math.ceil((len(scores)+1)*.95)-1];assert cutoff==th['cutoffs'][str(a['arm'])]
                    cal_nonconverged+=sum(r['converged']=='0' for r in cal if int(r['arm'])==a['arm'])
                    rr=[r for r in rows if int(r['arm'])==a['arm']]
                    for scenario in ['null','shared','heterogeneous','one','three']:
                        assert sum(r['scenario']==scenario for r in rr)==cfg[phase][scenario]
                    for r in rr:
                        assert r['status'] in ['ok','nonconverged','actual_manc','analytic_no_entry','no_eligible']
                        assert int(r['tp'])+int(r['fp'])==int(r['selected'])
                        assert math.isfinite(float(r['mse_african'])) and float(r['mse_african'])>=0
                        if a['method']!='Manc':
                            for group in ['causal','noncausal']:
                                assert 0<=int(r['covered_'+group])<=int(r['n_'+group]);assert float(r['width_'+group])>=0
                    details.append(dict(phase=phase,region=region,population=pop,arm=a['arm'],rows=len(rr),nonconverged=sum(r['converged']=='0' for r in rr),actual_manc=sum(r['status']=='actual_manc' for r in rr),screened_manc=sum(r['status']=='analytic_no_entry' for r in rr)))
                b=np.load(ROOT/f'data/revision/banks/{phase}/r{region}_{pop}.npz');assert len(b['q'])==count
                # Regenerate one complete cohort pair per scenario independently from
                # the prepared haplotypes and public seeds, checking input identity.
                arr=np.load(ROOT/f'data/revision/reference/region_{region}.npz')
                for sid,scenario in enumerate(['calibration','null','shared','heterogeneous','one','three']):
                    if scenario=='calibration':continue
                    j=int(np.flatnonzero((b['scenario']==scenario)&(b['replicate']==0))[0])
                    truth=np.zeros((2,8));indices=cfg['causal_indices_zero_based'].get(scenario,[])
                    if indices:truth[:,indices]=cfg['effects'][scenario]
                    np.testing.assert_array_equal(truth,b['truth'][j])
                    rng=np.random.default_rng(np.random.SeedSequence([cfg[phase]['seed'],region,{'YRI':0,'LWK':1}[pop],sid,0]))
                    for k,(p,n) in enumerate([('CEU',2000),(pop,250)]):
                        h=arr[p+'_study'];draw=rng.integers(len(h),size=(n,2));g=h[draw].sum(axis=1);x=g-g.mean(0)
                        y=x@truth[k]+rng.normal(size=n);y-=y.mean();q=x.T@y;cross=x.T@x;d=np.diag(cross)
                        np.testing.assert_allclose(q,b['q'][j,k],atol=1e-10);np.testing.assert_allclose(d,b['d'][j,k],atol=1e-10)
                        np.testing.assert_allclose(cross/np.sqrt(d[:,None]*d[None,:]),b['R'][j,k],atol=1e-12)
                for j in range(count):
                    h=hashlib.sha256(b['q'][j].tobytes()).hexdigest();assert h not in identities;identities.add(h)
                    d=(b['n'][j,:,None]-1)*b['varY'][j,:,None]/((b['n'][j,:,None]-2)*b['se'][j]**2+b['b'][j]**2)
                    np.testing.assert_allclose(d,b['d'][j],rtol=1e-10);np.testing.assert_allclose(d*b['b'][j],b['q'][j],atol=1e-10)
                    assert np.linalg.eigvalsh(b['R'][j]).min()>0
                total+=len(b['q'])
    table(ROOT/'results/revision/baseline_audit_details.tsv',details)
    files=sorted((ROOT/'results/revision').glob('*/*/complete.json'))
    out=dict(passed=True,baseline_completed_before_inla=True,evaluation_datasets=total,sensitivity_evaluation_datasets=880,validation_evaluation_datasets=3600,calibration_datasets=3200,unique_summary_fingerprints=len(identities),nonconverged=sum(d['nonconverged'] for d in details),actual_manc=sum(d['actual_manc'] for d in details),screened_manc=sum(d['screened_manc'] for d in details),completion_hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files})
    out.update(calibration_nonconverged=cal_nonconverged,audit_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),code_sha256=expected_hash)
    gatefile=ROOT/'results/revision/baseline_gate.json'
    if gatefile.exists():
        previous=json.loads(gatefile.read_text())
        if {k:v for k,v in previous.items() if k!='audit_utc'}=={k:v for k,v in out.items() if k!='audit_utc'}:out=previous
    gatefile.write_text(json.dumps(out,indent=2)+'\n',newline='\n');print(out)
if __name__=='__main__':main()
