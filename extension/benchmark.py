"""Reference-informed benchmark; calibration is frozen before held-out evaluation."""
from pathlib import Path
import os
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='1'
os.environ['MKL_NUM_THREADS']='1';os.environ['NUMBA_NUM_THREADS']='1'
import sys, json, csv, math, hashlib, subprocess, argparse, time, contextlib, io
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tmp/runtime'));sys.path.insert(0,str(ROOT/'scripts'))
import numpy as np
from scipy.special import ndtr
from MultiSuSiE import multisusie_rss
from simulate_compare import bayes, write_bed, table
CFG=json.loads((ROOT/'extension/design.json').read_text())
METHODS=['Bayes rho 0','Bayes rho 0.8','MultiSuSiE','Manc-COJO']

def calibration_cutoff(scores,alpha=.05):
    """Strict exceedance of a finite-sample order statistic; ties conservative."""
    k=math.ceil((len(scores)+1)*(1-alpha))
    if k>len(scores):return float('inf')
    return float(np.sort(scores)[k-1])

def summarize(g,y):
    x=g-g.mean(0);yc=y-y.mean();d=np.sum(x*x,axis=0);q=x.T@yc
    b=q/d;se=np.sqrt((yc@yc-q*q/d)/(len(y)-2)/d)
    return dict(b=b,se=se,d=d,q=q,varY=float(np.var(y,ddof=1)),freq=g.mean(0)/2,n=len(y))

def synthetic(h,n,beta,rng):
    ix=rng.integers(0,len(h),size=(n,2));g=h[ix[:,0]]+h[ix[:,1]]
    y=(g-g.mean(0))@beta+rng.normal(size=n)
    return summarize(g,y)

def reference_files(work,arrays,variants,pops,condition):
    refs=[]
    for k,pop in enumerate(pops):
        rp=({'YRI':'LWK','LWK':'YRI'}[pop] if condition=='swapped' and k==1 else pop)
        h=arrays[rp+'_reference'];g=h[::2]+h[1::2];prefix=work/f'ref{k}'
        write_bed(prefix,g)
        prefix.with_suffix('.bim').write_text(''.join(f"22\t{v['id']}\t0\t{v['position']}\t{v['alt']}\t{v['ref']}\n" for v in variants))
        refs.append(prefix)
    return refs

def cojo(work,refs,stats,variants,threshold):
    for k,s in enumerate(stats):
        table(work/f'study{k}.ma',[dict(SNP=v['id'],A1=v['alt'],A2=v['ref'],freq=s['freq'][j],
            b=s['b'][j],se=s['se'][j],p=max(float(2*ndtr(-abs(s['b'][j]/s['se'][j]))),1e-300),N=s['n']) for j,v in enumerate(variants)])
    out=work/'manc'; output=work/'manc.jma.cojo'
    if output.exists():output.unlink() # fixed task-local file, prevents stale empty-fit results
    binary=Path(os.environ.get('MANC_COJO',str(ROOT/('tools/manc_cojo_win.exe' if os.name=='nt' else 'tools/manc_cojo'))))
    cmd=[str(binary),'--bfile',*[str(p) for p in refs],'--cojo-file',str(work/'study0.ma'),str(work/'study1.ma'),
        '--cojo-slct','--cojo-p',str(threshold),'--diff-freq','1','1','--thread-num','1','--out',str(out)]
    run=subprocess.run(cmd,capture_output=True,text=True,timeout=45)
    log=run.stdout+'\n'+run.stderr
    (work/'last_command.json').write_text(json.dumps(cmd,indent=2)+'\n')
    (work/'last_stdout.txt').write_text(log)
    if run.returncode:raise RuntimeError(log[-2500:])
    selected=np.zeros(8,dtype=bool);effects=np.zeros((2,8))
    if output.exists():
        with output.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
        # Manc uses whitespace in some platform releases.
        if rows and 'SNP' not in rows[0]:
            lines=[x.split() for x in output.read_text().splitlines() if x.strip()]
            rows=[dict(zip(lines[0],x)) for x in lines[1:]]
        lookup={v['id']:i for i,v in enumerate(variants)}
        for row in rows:
            j=lookup[row['SNP']]
            if row['A1']!=variants[j]['alt']:raise ValueError('Effect allele mismatch')
            selected[j]=True;effects[:,j]=[float(row['bJ.1']),float(row['bJ.2'])]
    elif not any(x in log.lower() for x in ['no snp','no significant','0 snp','no variant','not significant']):
        raise RuntimeError('Missing Manc output without explicit empty-selection message: '+log[-1200:])
    if not np.isfinite(effects).all():raise ValueError('Nonfinite Manc effects')
    return selected,effects

def fit_bayes(stats,rs):
    fits={};cohorts=[]
    for s,r in zip(stats,rs):
        # d recovered identically from b, se, varY and n; study genotypes not supplied to estimators.
        d=(s['varY']*(s['n']-1))/((s['n']-2)*s['se']**2+s['b']**2)
        cohorts.append((r*np.sqrt(d[:,None]*d[None,:]),d*s['b'],s['varY']))
    for rho in CFG['bayesian_correlations']:
        pip,mean,sd=bayes(cohorts,rho,CFG['bayesian_tau'],CFG['bayesian_inclusion_prior'])
        fits[f'Bayes rho {rho:g}']=(pip,mean,True,0)
    with contextlib.redirect_stdout(io.StringIO()):
        fit=multisusie_rss(R_list=rs,population_sizes=[s['n'] for s in stats],
            b_list=[s['b'] for s in stats],s_list=[s['se'] for s in stats],varY_list=[s['varY'] for s in stats],
            L=CFG['multisusie_L'],rho=CFG['multisusie_rho'],estimate_residual_variance=False,
            single_population_mac_thresh=0,max_iter=CFG['multisusie_max_iter'],tol=CFG['multisusie_tolerance'],
            float_type=np.float64,min_abs_corr=.5)
    fits['MultiSuSiE']=(fit.pip,fit.coef,bool(fit.converged),int(fit.niter))
    for pip,mean,_,_ in fits.values():
        if not np.isfinite(pip).all() or not np.isfinite(mean).all():raise ValueError('Nonfinite Bayesian fit')
    return fits

def meta_score(stats):
    weights=np.stack([1/s['se']**2 for s in stats]);bs=np.stack([s['b'] for s in stats])
    z=np.sum(weights*bs,axis=0)/np.sqrt(weights.sum(0))
    return -math.log10(max(float(np.min(2*ndtr(-np.abs(z)))),1e-300))

def cells():
    return [(l,p,n,c) for l in range(2) for p in ['YRI','LWK'] for n in CFG['n_african'] for c in CFG['reference_conditions']]

def run_cell(cell_id,smoke=False):
    locus,pop,n,condition=cells()[cell_id]
    out=ROOT/('tmp/smoke' if smoke else 'results/extension/cells')/f'cell_{cell_id:02d}';out.mkdir(parents=True,exist_ok=True)
    stamp=hashlib.sha256((ROOT/'extension/design.json').read_bytes()+(ROOT/'extension/benchmark.py').read_bytes()).hexdigest()
    done=out/'complete.json'
    if done.exists():
        assert json.loads(done.read_text())['design_code_sha256']==stamp,'Checkpoint differs from code/design; use a new output location'
        return cell_id,'cached'
    work=ROOT/'tmp/benchmark'/f'cell_{cell_id:02d}';work.mkdir(parents=True,exist_ok=True)
    arrays=dict(np.load(ROOT/f'data/prepared/locus_{locus}.npz'))
    variants=json.loads((ROOT/f'data/prepared/locus_{locus}.json').read_text())
    pops=['CEU',pop];refpop=pop if condition=='matched' else {'YRI':'LWK','LWK':'YRI'}[pop]
    rs=[arrays['CEU_reference_R'],arrays[refpop+'_reference_R']]
    refs=reference_files(work,arrays,variants,pops,condition)
    counts=[30,3,3,3] if smoke else [CFG['calibration_replicates'],CFG['null_evaluation_replicates'],CFG['alternative_replicates'],CFG['alternative_replicates']]
    scores={m:[] for m in METHODS};calrows=[];metrics=[];variantrows=[];thresholds=None
    for phase,(scenario,count) in enumerate(zip(['calibration','null','shared','heterogeneous'],counts)):
        if phase==1:
            thresholds={m:calibration_cutoff(scores[m],CFG['target_locus_false_selection_rate']) for m in METHODS}
            table(out/'calibration_scores.tsv',calrows)
            (out/'thresholds.json').write_text(json.dumps(dict(cell_id=cell_id,thresholds=thresholds,calibration_replicates=counts[0],
                rule='score strictly exceeds order statistic ceil((B+1)*0.95)',design_code_sha256=stamp),indent=2)+'\n')
        for rep in range(count):
            rng=np.random.default_rng(np.random.SeedSequence([CFG['seed'],cell_id,phase,rep]))
            truth=np.zeros((2,8))
            if scenario in ('shared','heterogeneous'):truth[:,CFG['causal_marker_indices_zero_based']]=CFG[scenario+'_effects']
            stats=[synthetic(arrays[p+'_study'],size,truth[k],rng) for k,(p,size) in enumerate(zip(pops,[CFG['n_ceu'],n]))]
            fits=fit_bayes(stats,rs);ms=meta_score(stats)
            if phase==0:
                for method,(pip,_,converged,_) in fits.items():
                    # Algorithm output after the fixed iteration budget defines the score;
                    # nonconvergence is retained and reported, never removed selectively.
                    score=float(pip.max());scores[method].append(score)
                    calrows.append(dict(cell_id=cell_id,replicate=rep,method=method,score=score,converged=int(converged)))
                scores['Manc-COJO'].append(ms);calrows.append(dict(cell_id=cell_id,replicate=rep,method='Manc-COJO',score=ms,converged=1))
                continue
            causal=np.any(truth!=0,axis=0)
            selections={m:pip>thresholds[m] for m,(pip,_,_,_) in fits.items()}
            # Actual stepwise fit at calibrated entry threshold, not filtering a permissive path.
            sel,eff=cojo(work,refs,stats,variants,10**(-thresholds['Manc-COJO']))
            selections['Manc-COJO']=sel;fits['Manc-COJO']=(np.full(8,np.nan),eff,True,0)
            for method,(pip,effects,converged,niter) in fits.items():
                selected=selections[method];tp=int(np.sum(selected&causal));fp=int(np.sum(selected&~causal))
                metrics.append(dict(cell_id=cell_id,locus=locus,population=pop,n_african=n,reference=condition,scenario=scenario,replicate=rep,
                    method=method,selected=int(selected.sum()),tp=tp,fp=fp,recall=tp/2 if causal.any() else '',any_false_positive=int(fp>0),
                    mse_african=float(np.mean((effects[1]-truth[1])**2)),converged=int(converged),iterations=niter,
                    population_specific_recovered=int(selected[6]) if scenario=='heterogeneous' else '',
                    fixed_pip95_selected=int(np.sum(pip>.95)) if method!='Manc-COJO' else '',meta_entry_score=ms))
                for j in range(8):variantrows.append(dict(cell_id=cell_id,scenario=scenario,replicate=rep,method=method,variant=variants[j]['id'],
                    index=j,pip=float(pip[j]) if np.isfinite(pip[j]) else '',selected=int(selected[j]),true_ceu=truth[0,j],true_african=truth[1,j],
                    estimate_ceu=effects[0,j],estimate_african=effects[1,j]))
        print(f'cell {cell_id} {scenario} completed {count}',flush=True)
    table(out/'metrics.tsv',metrics);table(out/'variant_results.tsv',variantrows)
    done.write_text(json.dumps(dict(cell_id=cell_id,design_code_sha256=stamp,evaluation_datasets=sum(counts[1:]),calibration_datasets=counts[0],
        threshold_sha256=hashlib.sha256((out/'thresholds.json').read_bytes()).hexdigest()),indent=2)+'\n')
    return cell_id,'completed'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--smoke',action='store_true');ap.add_argument('--workers',type=int,default=2);ap.add_argument('--cell',type=int)
    args=ap.parse_args()
    if args.cell is not None:print(run_cell(args.cell,args.smoke));return
    from concurrent.futures import ProcessPoolExecutor,as_completed
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        fs=[pool.submit(run_cell,i,args.smoke) for i in range(len(cells()))]
        for f in as_completed(fs):print(f.result(),flush=True)

if __name__=='__main__':main()
