"""Paired sensitivity experiment and independent-seed validation; no INLA here."""
from engine import *
import json,csv,hashlib,time,argparse
from manc import table,eligible,entry_score,run as run_manc
CFG=json.loads((ROOT/'revision/protocol.json').read_text())
def digest():
    files=['revision/protocol.json','revision/engine.py','revision/manc.py','revision/benchmark.py','revision/prepare_reference.py']
    return hashlib.sha256(b''.join((ROOT/p).read_bytes() for p in files)).hexdigest()
def arms(phase):
    out=[]
    def add(method,ld,variance='total',rho=0,tau=.15,pi=.2,L=5,freq=.2):
        a=dict(method=method,ld=ld,variance=variance,rho=rho,tau=tau,pi=pi,L=L,freq=freq)
        if a not in out:out.append(a)
    for ld in CFG[phase+'_ld']:
        for v in (CFG['variance_sensitivity'] if phase=='sensitivity' else ['total']):
            for rho in [0,.8]:add('Exact',ld,v,rho=rho)
            add('MultiSuSiE',ld,v,rho=.75)
        for freq in ([.1,.2,1] if phase=='sensitivity' and ld!='in_sample' else [.2]):add('Manc',ld,freq=freq)
    if phase=='sensitivity':
        for L in [1,2,3]:add('MultiSuSiE','matched',rho=.75,L=L)
        for rho in [0,.5,.95]:add('MultiSuSiE','matched',rho=rho)
        for pi in [.05,.4]:add('Exact','matched',rho=.8,pi=pi)
        for tau in [.075,.3]:add('Exact','matched',rho=.8,tau=tau)
    for i,a in enumerate(out):a['arm']=i
    return out
def generate(arr,pop,phase,region,scenario,rep):
    sid=['calibration','null','shared','heterogeneous','one','three'].index(scenario)
    rng=np.random.default_rng(np.random.SeedSequence([CFG[phase]['seed'],region,{'YRI':0,'LWK':1}[pop],sid,rep]))
    truth=np.zeros((2,8));idx=CFG['causal_indices_zero_based'].get(scenario,[])
    if idx:truth[:,idx]=CFG['effects'][scenario]
    stats=[];studyR=[]
    for k,(p,n) in enumerate(zip(['CEU',pop],[2000,250])):
        h=arr[p+'_study'];ix=rng.integers(len(h),size=(n,2));g=h[ix[:,0]]+h[ix[:,1]];x=g-g.mean(0)
        y=x@truth[k]+rng.normal(size=n);y-=y.mean();Q=x.T@x;d=np.diag(Q);q=x.T@y;b=q/d
        se=np.sqrt((y@y-q*q/d)/(n-2)/d);r=Q/np.sqrt(d[:,None]*d[None,:])
        stats.append(dict(b=b,se=se,d=d,q=q,n=n,varY=float(y@y/(n-1)),freq=g.mean(0)/2));studyR.append(r)
    return stats,studyR,truth
def references(a,arr,pop,stats,studyR):
    if a['ld']=='in_sample':return studyR,np.array([s['freq'] for s in stats])
    rp=pop if a['ld']=='matched' else {'YRI':'LWK','LWK':'YRI'}[pop]
    return [arr[p+'_reference_R'] for p in ['CEU',rp]],np.array([arr[p+'_reference'].mean(0) for p in ['CEU',rp]])
def fit(a,stats,rs,intervals):
    Q=[r*np.sqrt(s['d'][:,None]*s['d'][None,:]) for r,s in zip(rs,stats)];q=[s['q'] for s in stats]
    variance=[s['varY'] for s in stats] if a['variance']=='total' else [1.,1.]
    if a['method']=='Exact':return custom(Q,q,variance,a['rho'],a['tau'],a['pi'],intervals)
    return multisusie(Q,q,stats,variance,a['L'],a['rho'],intervals)
def cell(phase,region,pop,smoke=False):
    out=ROOT/('tmp/revision_smoke' if smoke else 'results/revision')/phase/f'r{region}_{pop}';out.mkdir(parents=True,exist_ok=True)
    stamp=digest();complete=out/'complete.json'
    if complete.exists():assert json.loads(complete.read_text())['code_sha256']==stamp;return str(out)+' cached'
    arr=dict(np.load(ROOT/f'data/revision/reference/region_{region}.npz'));vs=json.loads((ROOT/f'data/revision/reference/region_{region}.json').read_text());aa=arms(phase)
    counts={k:(20 if k=='calibration' else 2) if smoke else CFG[phase][k] for k in ['calibration','null','shared','heterogeneous','one','three']}
    work=ROOT/'tmp/revision_manc'/phase/f'r{region}_{pop}';scores={a['arm']:[] for a in aa};cal=[];metrics=[];variants=[];bank=[]
    for scenario,count in counts.items():
        if scenario=='null':
            rank=math.ceil((counts['calibration']+1)*.95)
            cuts={i:float(np.sort(s)[rank-1]) for i,s in scores.items()}
            table(out/'calibration_scores.tsv',cal)
            (out/'thresholds.json').write_text(json.dumps(dict(code_sha256=stamp,seed=CFG[phase]['seed'],rank=rank,cutoffs=cuts,arms=aa),indent=2)+'\n',newline='\n')
            threshold_hash=hashlib.sha256((out/'thresholds.json').read_bytes()).hexdigest()
        for rep in range(count):
            stats,studyR,truth=generate(arr,pop,phase,region,scenario,rep)
            if scenario!='calibration':bank.append(dict(scenario=scenario,rep=rep,truth=truth,stats=stats,R=studyR))
            for a in aa:
                t=time.perf_counter();i=a['arm'];rs,fr=references(a,arr,pop,stats,studyR);keep=eligible(stats,fr,a['freq']) if a['method']=='Manc' else np.ones(8,bool)
                status='ok';converged=True
                if a['method']=='Manc':
                    score=entry_score(stats,keep)
                    if scenario=='calibration':f=None
                    else:
                        selected,mean,status=run_manc(work,stats,vs,rs,fr,keep,cuts[i]);f=dict(mean=mean,lower=None,upper=None)
                else:
                    f=fit(a,stats,rs,scenario!='calibration');score=float(f['pip'].max());converged=f['converged']
                    if not converged:status='nonconverged'
                    assert np.isfinite(f['pip']).all() and np.isfinite(f['mean']).all()
                    if scenario!='calibration':selected=f['pip']>cuts[i]
                if scenario=='calibration':scores[i].append(score);cal.append(dict(arm=i,replicate=rep,score=score,converged=int(converged)));continue
                causal=np.any(truth!=0,axis=0);nc=int(causal.sum());tp=int((selected&causal).sum());fp=int((selected&~causal).sum());err=f['mean']-truth
                row=dict(phase=phase,region=region,population=pop,scenario=scenario,replicate=rep,**a,status=status,converged=int(converged),seconds=time.perf_counter()-t,
                    selected=int(selected.sum()),tp=tp,fp=fp,recall=tp/nc if nc else '',fdp=fp/max(1,int(selected.sum())),any_false_positive=int(fp>0),
                    ld_recall_ceu=match_recovery(selected,causal,studyR[0]) if nc else '',ld_recall_african=match_recovery(selected,causal,studyR[1]) if nc else '',
                    retained_causal=int((keep&causal).sum()),retained_noncausal=int((keep&~causal).sum()),mse_african=float(np.mean(err[1]**2)),
                    specific_estimate=float(f['mean'][1,6]) if scenario in ['heterogeneous','three'] else '',specific_bias=float(err[1,6]) if scenario in ['heterogeneous','three'] else '')
                for group,mask in [('causal',truth[1]!=0),('noncausal',truth[1]==0)]:
                    row['n_'+group]=int(mask.sum());row['covered_'+group]='';row['width_'+group]=''
                    if f['lower'] is not None:
                        row['covered_'+group]=int(((f['lower'][1]<=truth[1])&(f['upper'][1]>=truth[1])&mask).sum())
                        row['width_'+group]=float((f['upper'][1]-f['lower'][1])[mask].sum())
                for label,mask in [('causal',causal),('noncausal',~causal)]:
                    for stat in ['mean','max']:
                        row['pip_'+label+'_'+stat]=float(getattr(np,stat)(f['pip'][mask])) if 'pip' in f and mask.any() else ''
                metrics.append(row)
                if phase=='validation':
                    for j in range(8):variants.append(dict(arm=i,scenario=scenario,replicate=rep,index=j,truth_ceu=truth[0,j],truth_african=truth[1,j],estimate_ceu=f['mean'][0,j],estimate_african=f['mean'][1,j],selected=int(selected[j]),pip=f['pip'][j] if 'pip' in f else '',lower_african=f['lower'][1,j] if f['lower'] is not None else '',upper_african=f['upper'][1,j] if f['upper'] is not None else ''))
        print(f'{phase} r{region} {pop}: {scenario} {count}',flush=True)
    table(out/'metrics.tsv',metrics);table(out/'variant_results.tsv',variants)
    bankdir=ROOT/'data/revision/banks'/phase;bankdir.mkdir(parents=True,exist_ok=True)
    # Publicly reproducible summary-only bank, not committed: all source seeds are in the protocol.
    np.savez_compressed(bankdir/f'r{region}_{pop}.npz',scenario=np.array([b['scenario'] for b in bank]),replicate=np.array([b['rep'] for b in bank]),truth=np.array([b['truth'] for b in bank]),R=np.array([b['R'] for b in bank]),**{key:np.array([[s[key] for s in b['stats']] for b in bank]) for key in ['b','se','d','q','n','varY','freq']})
    complete.write_text(json.dumps(dict(code_sha256=stamp,threshold_sha256=threshold_hash,calibration_datasets=counts['calibration'],evaluation_datasets=sum(counts.values())-counts['calibration'],method_arms=len(aa),all_metrics=len(metrics)),indent=2)+'\n',newline='\n')
    return f'{phase} r{region} {pop} complete'
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--phase',choices=['sensitivity','validation'],required=True);ap.add_argument('--workers',type=int,default=2);ap.add_argument('--smoke',action='store_true');args=ap.parse_args()
    from concurrent.futures import ProcessPoolExecutor,as_completed
    regions=CFG[args.phase+'_regions'] if args.phase=='sensitivity' else CFG['fresh_validation_regions']
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        fs=[ex.submit(cell,args.phase,r,p,args.smoke) for r in regions for p in ['YRI','LWK']]
        for f in as_completed(fs):print(f.result(),flush=True)
if __name__=='__main__':main()
