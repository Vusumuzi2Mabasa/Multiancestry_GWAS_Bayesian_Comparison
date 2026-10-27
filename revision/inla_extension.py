"""INLA execution gate, exact Gaussian controls and independent quadrature."""
from engine import *
import json,argparse,csv,subprocess,hashlib,time
from scipy.linalg import block_diag
from scipy.special import logsumexp
from manc import table
RSCRIPT=os.environ.get('RSCRIPT','C:/Program Files/R/R-4.5.1/bin/Rscript.exe')
def cases():
    cfg=json.loads((ROOT/'revision/protocol.json').read_text());out=[]
    for region in range(6):
        arr=dict(np.load(ROOT/f'data/revision/reference/region_{region}.npz'))
        for pop in ['YRI','LWK']:
            b=np.load(ROOT/f'data/revision/banks/validation/r{region}_{pop}.npz')
            for scenario in ['null','shared','heterogeneous']:
                ids=np.flatnonzero((b['scenario']==scenario)&(b['replicate']<20))
                assert len(ids)==20
                for j in ids:
                    for ld in ['in_sample','matched']:
                        rs=b['R'][j] if ld=='in_sample' else [arr[p+'_reference_R'] for p in ['CEU',pop]]
                        Q=[r*np.sqrt(d[:,None]*d[None,:]) for r,d in zip(rs,b['d'][j])];B=[];y=[]
                        for k in range(2):
                            lower=np.linalg.cholesky(Q[k]);B.append(lower.T/np.sqrt(b['varY'][j,k]));y.append(np.linalg.solve(lower,b['q'][j,k])/np.sqrt(b['varY'][j,k]))
                        out.append(dict(case_id=f'r{region}_{pop}_{scenario}_{b["replicate"][j]}_{ld}',region=region,population=pop,scenario=scenario,replicate=int(b['replicate'][j]),ld=ld,B=block_diag(*B).tolist(),y=np.concatenate(y).tolist(),truth=b['truth'][j].reshape(16).tolist()))
    return out
T=np.block([[np.eye(8),np.eye(8),np.zeros((8,8))],[np.eye(8),np.zeros((8,8)),np.eye(8)]])
def conditional(B,y,delta_sd):
    A=B@T;v=np.r_[np.full(8,.15**2),np.full(16,delta_sd**2)]
    va=v[:,None]*A.T;S=np.eye(16)+A@va
    inv=np.linalg.inv(S);mean=v*(A.T@(inv@y));cov=np.diag(v)-va@inv@va.T
    return T@mean,np.maximum(np.diag(T@cov@T.T),0),-.5*(np.linalg.slogdet(S)[1]+y@inv@y+16*np.log(2*np.pi))
def exact(z,mode='adaptive',points=401,log_range=(-14,1)):
    B=np.array(z['B']);y=np.array(z['y'])
    if mode=='fixed':
        m,v,_=conditional(B,y,.075);sd=np.sqrt(v)
        return dict(mean=m,sd=sd,lower=m-1.95996398454*sd,upper=m+1.95996398454*sd,tail_mass=0.)
    grid=np.linspace(*log_range,points);tau=np.exp(grid);lam=-np.log(.05)/.2
    fits=[conditional(B,y,t) for t in tau];m=np.array([x[0] for x in fits]);v=np.array([x[1] for x in fits]);lw=np.array([x[2] for x in fits])+np.log(lam)+grid-lam*tau
    lw[[0,-1]]+=np.log(.5);w=np.exp(lw-logsumexp(lw));mean=w@m;sd=np.sqrt(w@(v+m*m)-mean*mean)
    # No atoms for this continuous hierarchy; direct mixture CDF is appropriate.
    limits=[]
    for j in range(16):
        ss=np.sqrt(v[:,j]);bound=max(1.,float(np.max(np.abs(m[:,j])+12*ss)))
        limits.append([brentq(lambda q:float(w@ndtr((q-m[:,j])/ss))-p,-bound,bound,xtol=1e-10) for p in [.025,.975]])
    return dict(mean=mean,sd=sd,lower=np.array(limits)[:,0],upper=np.array(limits)[:,1],tail_mass=float(w[:3].sum()+w[-3:].sum()))
def launch(shard):
    name,items=shard;work=ROOT/'tmp/inla_inputs';work.mkdir(parents=True,exist_ok=True);out=ROOT/'results/revision/inla/shards';out.mkdir(parents=True,exist_ok=True)
    inp=work/f'{name}.json';dest=out/f'{name}.json';inp.write_text(json.dumps(items)+'\n')
    stamp=hashlib.sha256(inp.read_bytes()+(ROOT/'revision/inla_fit.R').read_bytes()+(ROOT/'revision/inla_mixture.R').read_bytes()+(ROOT/'results/revision/baseline_gate.json').read_bytes()+(ROOT/'results/revision/software_provenance.json').read_bytes()).hexdigest();meta=out/f'{name}.meta.json'
    if dest.exists() and meta.exists() and json.loads(meta.read_text())['input_code_sha256']==stamp:
        rows=json.loads(dest.read_text())
        if len(rows)==len(items):return name,'cached'
    with (work/f'{name}.log').open('w') as f:
        rtemp=ROOT/'tmp/r_temp';rtemp.mkdir(exist_ok=True)
        renv=dict(os.environ,TMPDIR=str(rtemp),TMP=str(rtemp),TEMP=str(rtemp))
        p=subprocess.run([RSCRIPT,'--vanilla',str(ROOT/'revision/inla_fit.R'),str(inp),str(dest)],cwd=ROOT,env=renv,stdout=f,stderr=subprocess.STDOUT)
    if p.returncode:raise RuntimeError((work/f'{name}.log').read_text()[-4000:])
    meta.write_text(json.dumps(dict(input_code_sha256=stamp,cases=len(items)))+'\n')
    return name,'complete'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--stage',choices=['validate','run','audit'],required=True);ap.add_argument('--workers',type=int,default=2);args=ap.parse_args()
    gate=json.loads((ROOT/'results/revision/baseline_gate.json').read_text());assert gate['passed']
    allcases=cases();out=ROOT/'results/revision/inla';out.mkdir(parents=True,exist_ok=True)
    # Evenly spaced region/population/scenario/LD contexts; deterministic, outcome-blind.
    fixed=[dict(allcases[i],mode='fixed') for i in np.linspace(0,len(allcases)-1,24,dtype=int)]
    adaptive=[dict(allcases[i],mode='adaptive') for i in np.linspace(0,len(allcases)-1,48,dtype=int)]
    if args.stage=='validate':
        launch(('validate_fixed',fixed));launch(('validate_adaptive',adaptive))
        checks=[];ok=True
        for mode,inputs in [('fixed',fixed),('adaptive',adaptive)]:
            rr=json.loads((out/f'shards/validate_{mode}.json').read_text())
            assert len(inputs)==len(rr)
            for z,r in zip(inputs,rr):
                assert z['case_id']==r['case_id'] and mode==r['mode']
                e=exact(z,mode);meanerr=float(np.max(np.abs(np.array(r['mean'])-e['mean'])));sderr=float(np.max(np.abs(np.array(r['sd'])/e['sd']-1)));qerr=float(max(np.max(abs(np.array(r[k])-e[k])) for k in ['lower','upper']))
                refine=exact(z,mode,801);referr=float(max(np.max(abs(e[k]-refine[k])) for k in ['mean','sd','lower','upper']))
                extended=exact(z,mode,1001,(-18,1.5));rangeerr=float(max(np.max(abs(e[k]-extended[k])) for k in ['mean','sd','lower','upper']))
                passed=meanerr<(1e-4 if mode=='fixed' else .005) and sderr<(.01 if mode=='fixed' else .05) and (mode=='fixed' or qerr<.015) and referr<1e-5 and rangeerr<1e-5
                ok &= passed
                checks.append(dict(case_id=z['case_id'],mode=mode,mean_max_abs_error=meanerr,sd_max_relative_error=sderr,quantile_max_abs_error=qerr,quadrature_refinement_max_error=referr,quadrature_extended_range_max_error=rangeerr,endpoint_grid_mass=e['tail_mass'],passed=int(passed)))
        table(out/'numerical_validation.tsv',checks)
        report=dict(passed=bool(ok),fixed_cases=24,adaptive_cases=48,max_mean_error=max(x['mean_max_abs_error'] for x in checks),max_sd_relative_error=max(x['sd_max_relative_error'] for x in checks),baseline_gate_sha256=hashlib.sha256((ROOT/'results/revision/baseline_gate.json').read_bytes()).hexdigest())
        report.update(inla_r_sha256=hashlib.sha256((ROOT/'revision/inla_fit.R').read_bytes()).hexdigest(),mixture_r_sha256=hashlib.sha256((ROOT/'revision/inla_mixture.R').read_bytes()).hexdigest(),extension_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        (out/'validation.json').write_text(json.dumps(report,indent=2)+'\n');print(report);assert ok,'INLA failed prespecified numerical validation; do not run the primary experiment'
    elif args.stage=='run':
        validation=json.loads((out/'validation.json').read_text());assert validation['passed']
        assert validation['inla_r_sha256']==hashlib.sha256((ROOT/'revision/inla_fit.R').read_bytes()).hexdigest()
        assert validation['mixture_r_sha256']==hashlib.sha256((ROOT/'revision/inla_mixture.R').read_bytes()).hexdigest()
        assert validation['extension_code_sha256']==hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
        from concurrent.futures import ThreadPoolExecutor
        items=[dict(z,mode='adaptive') for z in allcases]
        shards=[(f'primary_{i:02}',items[i*60:(i+1)*60]) for i in range(24)]
        with ThreadPoolExecutor(max_workers=args.workers) as ex:
            for result in ex.map(launch,shards):print(result,flush=True)
    else:
        results=[]
        for i in range(24):results.extend(json.loads((out/f'shards/primary_{i:02}.json').read_text()))
        lookup={r['case_id']:r for r in results};assert len(lookup)==len(allcases)==1440
        metrics=[];details=[]
        for z in allcases:
            r=lookup[z['case_id']];assert r['status']=='ok';truth=np.array(z['truth']).reshape(2,8)
            for k in ['mean','sd','lower','upper']:assert len(r[k])==16 and np.isfinite(r[k]).all()
            assert (np.array(r['sd'])>0).all() and (np.array(r['lower'])<=np.array(r['upper'])).all()
            for method,f in [('Adaptive INLA',r),('Fixed Gaussian',exact(z,'fixed'))]:
                mean=np.array(f['mean']).reshape(2,8);lo=np.array(f['lower']).reshape(2,8);hi=np.array(f['upper']).reshape(2,8)
                causal=truth[1]!=0
                row={k:z[k] for k in ['case_id','region','population','scenario','replicate','ld']};row.update(method=method,mse_african=float(np.mean((mean[1]-truth[1])**2)),specific_estimate=float(mean[1,6]) if z['scenario']=='heterogeneous' else '',specific_bias=float(mean[1,6]-truth[1,6]) if z['scenario']=='heterogeneous' else '',seconds=r['seconds'] if method=='Adaptive INLA' else '')
                for group,mask in [('causal',causal),('noncausal',~causal)]:row.update({f'n_{group}':int(mask.sum()),f'covered_{group}':int(((lo[1]<=truth[1])&(hi[1]>=truth[1])&mask).sum()),f'width_{group}':float((hi[1]-lo[1])[mask].sum())})
                metrics.append(row)
                for j in range(8):details.append(dict(case_id=z['case_id'],method=method,index=j,truth_african=truth[1,j],estimate_african=mean[1,j],lower_african=lo[1,j],upper_african=hi[1,j]))
        table(out/'metrics.tsv',metrics);table(out/'variant_results.tsv',details)
        (out/'audit.json').write_text(json.dumps(dict(passed=True,primary_contexts=1440,all_statuses_ok=True,continuous_prior_no_pip=True),indent=2)+'\n');print('INLA audit passed: 1440 primary contexts')
if __name__=='__main__':main()
