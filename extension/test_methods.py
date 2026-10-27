"""Numerical and file-interface tests on independent toy data (no tuning)."""
from pathlib import Path
import json, time, math
import benchmark as b
import numpy as np
from simulate_compare import genotypes

def main():
    checks={};rng=np.random.default_rng(583029)
    checks['finite_rank_cutoff']=b.calibration_cutoff(list(range(300)))==285
    checks['ties_conservative']=not (1>b.calibration_cutoff([1]*300))
    variants=[dict(id=f'rsTEST{i}',position=20000000+i*1000,ref='G',alt='A') for i in range(8)]
    arrays={};stats=[];rs=[]
    truth=np.zeros(8);truth[1]=.25;truth[6]=.2
    for pop,n in [('CEU',800),('YRI',400)]:
        g=genotypes(n,.5,rng);y=(g-g.mean(0))@truth+rng.normal(size=n);s=b.summarize(g,y);stats.append(s)
        arrays[pop+'_reference']=np.stack([g//2,(g+1)//2],axis=1).reshape(2*n,8)
        r=np.corrcoef(g,rowvar=False);rs.append(r)
        d=s['varY']*(n-1)/((n-2)*s['se']**2+s['b']**2)
        checks[pop+'_summary_reconstruction']=bool(np.allclose(d,s['d'],rtol=1e-10))
    start=time.perf_counter();fits=b.fit_bayes(stats,rs);seconds=time.perf_counter()-start
    checks['bayesian_outputs_finite']=all(np.isfinite(x[0]).all() and np.isfinite(x[1]).all() for x in fits.values())
    work=b.ROOT/'tmp/method_tests';work.mkdir(parents=True,exist_ok=True)
    refs=b.reference_files(work,arrays,variants,['CEU','YRI'],'matched')
    # An entry boundary test exercises the executable on either side of its first-step score.
    score=b.meta_score(stats);p=10**(-score)
    sel0,_=b.cojo(work,refs,stats,variants,p*.999)
    sel1,_=b.cojo(work,refs,stats,variants,min(p*1.001,1))
    checks['manc_first_entry_boundary']=bool(not sel0.any() and sel1.any())
    checks['multisusie_converged']=fits['MultiSuSiE'][2]
    report=dict(checks=checks,bayesian_three_fits_seconds=seconds,
        multisusie_iterations=fits['MultiSuSiE'][3],all_passed=all(checks.values()))
    (b.ROOT/'results/extension').mkdir(parents=True,exist_ok=True)
    (b.ROOT/'results/extension/method_tests.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
    assert all(checks.values())

if __name__=='__main__':main()
