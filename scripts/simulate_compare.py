"""Small-locus pilot: actual Manc-COJO versus enumerated Bayesian model.

Requires numpy only. No substitute COJO implementation is used if the binary
is absent. Run from the project root; see README.md for design/limitations.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
os.environ.setdefault('OMP_NUM_THREADS', '1')
from pathlib import Path
import argparse, csv, hashlib, itertools, json, math, subprocess, time
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
M=8
IDS=[f'rsSIM{j+1:03d}' for j in range(M)]
N_AFR=[250,500,1000,2000]
N_EUR=2000

def table(path, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def genotypes(n,rho,rng):
    # Two independent haplotypes; two independent four-marker blocks.
    h=np.empty((2*n,M),dtype=np.int8)
    for j in range(M):
        fresh=rng.binomial(1,.3,2*n)
        h[:,j]=fresh if j%4==0 else np.where(rng.random(2*n)<rho,h[:,j-1],fresh)
    return h[::2]+h[1::2]

def write_bed(prefix,g):
    n,p=g.shape
    # PLINK BED SNP-major: 00=A1/A1, 10=heterozygous, 11=A2/A2.
    # Our g is A-allele dosage (A1), so g=2 maps to 00, g=0 to 11.
    code=np.array([3,2,0],dtype=np.uint8)[g]
    with prefix.with_suffix('.bed').open('wb') as f:
        f.write(bytes([0x6c,0x1b,0x01]))
        for j in range(p):
            c=np.pad(code[:,j],(0,(-n)%4))
            f.write((c[::4]|(c[1::4]<<2)|(c[2::4]<<4)|(c[3::4]<<6)).tobytes())
    prefix.with_suffix('.bim').write_text(''.join(f'1\t{IDS[j]}\t0\t{100000+j*1000}\tA\tG\n' for j in range(p)))
    prefix.with_suffix('.fam').write_text(''.join(f'P{i+1} P{i+1} 0 0 0 -9\n' for i in range(n)))

def summarize(prefix,g,y,beta):
    n=len(y);x=g-g.mean(axis=0);yc=y-y.mean()
    q=x.T@yc;Q=x.T@x;d=np.diag(Q);b=q/d
    se=np.sqrt((yc@yc-q*q/d)/(n-2)/d)
    rows=[]
    for j in range(M):
        # Normal approximation, as used for summary-data COJO input.
        p=math.erfc(abs(b[j]/se[j])/math.sqrt(2))
        rows.append(dict(SNP=IDS[j],A1='A',A2='G',freq=g[:,j].mean()/2,
                         b=b[j],se=se[j],p=max(p,1e-300),N=n))
    table(prefix.with_suffix('.ma'),rows);write_bed(prefix,g)
    R=Q/np.sqrt(d[:,None]*d[None,:])
    np.savetxt(prefix.with_suffix('.ld.tsv'),R,delimiter='\t')
    table(prefix.with_suffix('.truth.tsv'),[dict(SNP=IDS[j],beta=beta[j]) for j in range(M)])
    # Bayesian analysis reconstructs all cross-products from summary b/se,
    # LD and reference variances. It never takes y or individual phenotypes.
    table(prefix.with_suffix('.variance.tsv'),[dict(SNP=IDS[j],sum_squares=d[j]) for j in range(M)])
    return read_summary(prefix)

def read_summary(prefix):
    with prefix.with_suffix('.ma').open() as f:r=list(csv.DictReader(f,delimiter='\t'))
    with prefix.with_suffix('.variance.tsv').open() as f:v=list(csv.DictReader(f,delimiter='\t'))
    b=np.array([float(t['b']) for t in r]);se=np.array([float(t['se']) for t in r])
    d=np.array([float(t['sum_squares']) for t in v]);n=int(r[0]['N'])
    R=np.loadtxt(prefix.with_suffix('.ld.tsv'));Q=R*np.sqrt(d[:,None]*d[None,:]);q=d*b
    syy=float(np.median(se**2*d*(n-2)+b*b*d))
    sigma2=(syy-q@np.linalg.solve(Q,q))/(n-M-1)
    assert sigma2>0 and np.linalg.eigvalsh(Q).min()>0
    return Q,q,sigma2

def bayes(cohorts,rho,tau=.15,pi=.2):
    """Exact sum over 2^8 shared inclusion configurations.

    gamma_j~Bernoulli(pi); if gamma_j=0 both effects are exactly zero.
    If active: mu_j~N(0,rho*tau^2), delta_jk~N(0,(1-rho)*tau^2),
    beta_jk=mu_j+delta_jk. Conditional residual variance is plug-in estimated.
    The normal likelihood is integrated analytically; no MCMC approximation.
    """
    configurations=np.array(list(itertools.product([0,1],repeat=M)),dtype=bool)
    lw=[];means=[];seconds=[]
    for mask in configurations:
        idx=np.flatnonzero(mask);s=len(idx);full=np.zeros(2*M);moment=np.zeros(2*M)
        logprior=s*math.log(pi)+(M-s)*math.log(1-pi)
        if s==0:lw.append(logprior);means.append(full);seconds.append(moment);continue
        prior=tau*tau*np.block([[np.eye(s),rho*np.eye(s)],[rho*np.eye(s),np.eye(s)]])
        pinv=np.linalg.inv(prior);H=pinv.copy();u=np.empty(2*s)
        for k,(Q,q,sigma2) in enumerate(cohorts):
            sl=slice(k*s,(k+1)*s);H[sl,sl]+=Q[np.ix_(idx,idx)]/sigma2;u[sl]=q[idx]/sigma2
        cov=np.linalg.inv(H);mean=cov@u
        logbf=-.5*np.linalg.slogdet(prior)[1]-.5*np.linalg.slogdet(H)[1]+.5*u@mean
        for k in range(2):
            sl=slice(k*s,(k+1)*s);target=k*M+idx
            full[target]=mean[sl];moment[target]=np.diag(cov)[sl]+mean[sl]**2
        lw.append(logprior+logbf);means.append(full);seconds.append(moment)
    lw=np.array(lw);weights=np.exp(lw-lw.max());weights/=weights.sum()
    pip=weights@configurations;mean=weights@np.array(means)
    sd=np.sqrt(np.maximum(0,weights@np.array(seconds)-mean*mean))
    assert abs(weights.sum()-1)<1e-10 and np.all((pip>=0)&(pip<=1+1e-12))
    return pip,mean.reshape(2,M),sd.reshape(2,M)

def run_cojo(binary,directory,threshold):
    out=directory/'manc'
    command=[str(binary),'--bfile',str(directory/'EUR'),str(directory/'AFR'),
             '--cojo-file',str(directory/'EUR.ma'),str(directory/'AFR.ma'),
             '--cojo-slct','--cojo-p',str(threshold),'--thread-num','1','--out',str(out)]
    t=time.perf_counter();result=subprocess.run(command,capture_output=True,text=True,timeout=60)
    (directory/'manc.stdout.txt').write_text(result.stdout+'\n'+result.stderr)
    (directory/'manc.command.json').write_text(json.dumps(command,indent=2))
    if result.returncode:raise RuntimeError(f'Manc-COJO failed: {directory}\n{result.stdout[-1500:]} {result.stderr}')
    file=directory/'manc.jma.cojo'
    effects=np.zeros((2,M));selected=np.zeros(M,dtype=bool)
    if file.exists():
        lines=[l.split() for l in file.read_text().splitlines() if l.strip()]
        if lines:
            header=lines[0]
            for values in lines[1:]:
                r=dict(zip(header,values));j=IDS.index(r['SNP']);selected[j]=True
                # Output effects are aligned to cohort-1 BIM A1=A.
                assert r['A1']=='A',r
                effects[0,j]=float(r['bJ.1']);effects[1,j]=float(r['bJ.2'])
    else:
        # A successful empty selection must be explicitly documented in log.
        text=(result.stdout+result.stderr+out.with_suffix('.log').read_text()).lower()
        assert any(s in text for s in ['no snp','no significant','0 snp','no variant','not significant']),text[-2000:]
    return selected,effects,time.perf_counter()-t

def metrics(method,selected,effects,truth,scenario,n,rep,seconds):
    causal=np.any(truth!=0,axis=0);tp=int(np.sum(selected&causal));fp=int(np.sum(selected&~causal))
    afr_only=(truth[0]==0)&(truth[1]!=0)
    causal_blocks={int(j//4) for j in np.flatnonzero(causal)}
    selected_blocks={int(j//4) for j in np.flatnonzero(selected)}
    return dict(scenario=scenario,n_eur=N_EUR,n_afr=n,replicate=rep,method=method,
        true_signals=int(causal.sum()),selected=int(selected.sum()),tp=tp,fp=fp,
        recall=tp/causal.sum() if causal.any() else '',fdp=fp/max(1,int(selected.sum())),
        any_false_positive=int(fp>0),afr_only_recall=float(np.mean(selected[afr_only])) if afr_only.any() else '',
        causal_block_recall=len(causal_blocks&selected_blocks)/len(causal_blocks) if causal_blocks else '',
        mse_afr=float(np.mean((effects[1]-truth[1])**2)),mse_eur=float(np.mean((effects[0]-truth[0])**2)),
        mse_afr_causal=float(np.mean((effects[1,causal]-truth[1,causal])**2)) if causal.any() else '',seconds=seconds)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--replicates',type=int,default=20)
    ap.add_argument('--seed',type=int,default=20261026);ap.add_argument('--manc',type=Path)
    ap.add_argument('--cojo-p',type=float,default=5e-8);args=ap.parse_args()
    binary=args.manc.resolve() if args.manc else None
    assert args.replicates>0
    allmetrics=[];variants=[]
    config=dict(seed=args.seed,replicates=args.replicates,n_eur=N_EUR,n_afr=N_AFR,markers=M,
        scenarios=['shared','heterogeneous','null'],rho_genotype_EUR=.8,rho_genotype_AFR=.35,
        bayes_tau=.15,bayes_inclusion_prior=.2,bayes_effect_correlations=[0,.8],
        pip_thresholds=[.5,.95],cojo_p=args.cojo_p,
        manc_binary=str(binary) if binary else None,
        manc_sha256=hashlib.sha256(binary.read_bytes()).hexdigest() if binary else None)
    (ROOT/'results').mkdir(exist_ok=True);(ROOT/'results/config.json').write_text(json.dumps(config,indent=2))
    for scenario_id,scenario in enumerate(config['scenarios']):
      for rep in range(1,args.replicates+1):
        rng=np.random.default_rng(np.random.SeedSequence([args.seed,scenario_id,rep]))
        gE=genotypes(N_EUR,.8,rng);gA=genotypes(max(N_AFR),.35,rng)
        truth=np.zeros((2,M))
        if scenario=='shared':truth[:,1]=.18;truth[:,5]=-.16
        if scenario=='heterogeneous':truth[:,1]=[.18,.09];truth[:,5]=[0,.28]
        yE=(gE-gE.mean(axis=0))@truth[0]+rng.normal(size=N_EUR)
        yA=(gA-gA.mean(axis=0))@truth[1]+rng.normal(size=max(N_AFR))
        for n in N_AFR:
          d=ROOT/'data'/scenario/f'rep_{rep:03d}'/f'AFR_{n}';d.mkdir(parents=True,exist_ok=True)
          cohorts=[summarize(d/'EUR',gE,yE,truth[0]),summarize(d/'AFR',gA[:n],yA[:n],truth[1])]
          for rho in [0,.8]:
            t=time.perf_counter();pip,effects,sd=bayes(cohorts,rho);elapsed=time.perf_counter()-t
            name=f'Bayes_rho{rho:g}'
            for threshold in [.5,.95]:
              allmetrics.append(metrics(f'{name}_PIP{threshold}',pip>=threshold,effects,truth,scenario,n,rep,elapsed))
            for j in range(M):variants.append(dict(scenario=scenario,replicate=rep,n_afr=n,method=name,
                SNP=IDS[j],pip=pip[j],beta_eur=effects[0,j],beta_afr=effects[1,j],
                sd_eur=sd[0,j],sd_afr=sd[1,j],truth_eur=truth[0,j],truth_afr=truth[1,j]))
          if binary:
            selected,effects,elapsed=run_cojo(binary,d,args.cojo_p)
            allmetrics.append(metrics('Manc_COJO',selected,effects,truth,scenario,n,rep,elapsed))
          # Checkpoints preserve useful work if a later external run fails.
          table(ROOT/'results/replicate_metrics.tsv',allmetrics);table(ROOT/'results/variant_results.tsv',variants)
        print(f'{scenario} replicate {rep}/{args.replicates} done',flush=True)
    status='COMPLETE_WITH_MANC' if binary else 'BAYES_ONLY_MANC_NOT_RUN'
    (ROOT/'results/status.txt').write_text(status+'\n')
    print(status,flush=True)

if __name__=='__main__':main()
