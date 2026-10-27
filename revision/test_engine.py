"""Executable numerical controls; no simulated outcome determines the design."""
from engine import *
import json
import sys
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'extension'))
from simulate_compare import bayes,write_bed
from benchmark import summarize,cojo
from MultiSuSiE import multisusie_rss
import manc
def main():
    rng=np.random.default_rng(916);arr=dict(np.load(ROOT/'data/revision/reference/region_0.npz'))
    vs=json.loads((ROOT/'data/revision/reference/region_0.json').read_text());stats=[];qs=[];Q=[];rs=[];freq=[]
    for p,n in [('CEU',2000),('YRI',250)]:
        h=arr[p+'_study'];ix=rng.integers(len(h),size=(n,2));g=h[ix[:,0]]+h[ix[:,1]];x=g-g.mean(0)
        y=x@np.array([0,.22,0,0,0,0,.28,0])+rng.normal(size=n);s=summarize(g,y);stats.append(s)
        rs.append(arr[p+'_reference_R']);Q.append(rs[-1]*np.sqrt(s['d'][:,None]*s['d'][None,:]));qs.append(s['q'])
        hr=arr[p+'_reference'];freq.append(hr.mean(0));
    v=[s['varY'] for s in stats];f=custom(Q,qs,v,intervals=True);old=bayes(list(zip(Q,qs,v)),.8)
    for a,b in zip([f['pip'],f['mean'],f['sd']],old):np.testing.assert_allclose(a,b,atol=1e-10)
    with contextlib.redirect_stdout(io.StringIO()):
        oldm=multisusie_rss(R_list=rs,population_sizes=[2000,250],b_list=[s['b'] for s in stats],s_list=[s['se'] for s in stats],varY_list=v,L=5,rho=.75,estimate_residual_variance=False,single_population_mac_thresh=0,max_iter=200,tol=.001,float_type=np.float64,min_abs_corr=.5)
    new=multisusie(Q,qs,stats,v,intervals=True)
    np.testing.assert_allclose(new['pip'],oldm.pip,atol=1e-7);np.testing.assert_allclose(new['mean'],oldm.coef,atol=1e-7)
    np.testing.assert_allclose(mixture_quantiles([1],[0],[1]),[-1.9599639845,1.9599639845],atol=1e-8)
    assert mixture_quantiles([.99,.01],[0,1],[0,1])==[0.,0.]
    sign=np.ones(8);sign[1]=-1
    flip=custom([x*sign[:,None]*sign[None,:] for x in Q],[x*sign for x in qs],v)
    np.testing.assert_allclose(flip['pip'],f['pip'],atol=1e-10);np.testing.assert_allclose(flip['mean'],f['mean']*sign,atol=1e-10)
    work=ROOT/'tmp/revision_tests';work.mkdir(exist_ok=True);refs=[]
    for k,p in enumerate(['CEU','YRI']):
        h=arr[p+'_reference'];g=h[::2]+h[1::2];pre=work/f'ref{k}';write_bed(pre,g)
        pre.with_suffix('.bim').write_text(''.join(f"22\t{z['id']}\t0\t{z['position']}\t{z['alt']}\t{z['ref']}\n" for z in vs));refs.append(pre)
    s1,e1=cojo(work,refs,stats,vs,.001)
    s2,e2,_=manc.run(work/'native',stats,vs,rs,np.array(freq),np.ones(8,bool),3,screen=False)
    np.testing.assert_array_equal(s1,s2);np.testing.assert_allclose(e1,e2,atol=1e-5)
    score=manc.entry_score(stats,np.ones(8,bool))
    below=manc.run(work/'below',stats,vs,rs,np.array(freq),np.ones(8,bool),score+.01,screen=False)[0]
    above=manc.run(work/'above',stats,vs,rs,np.array(freq),np.ones(8,bool),score-.01,screen=False)[0]
    assert not below.any() and above.any()
    out=dict(passed=True,exact_enumeration_tolerance=1e-10,multisusie_api_tolerance=1e-7,bed_native_manc_tolerance=1e-5,atom_quantiles=True,allele_flip_invariance=True,manc_entry_boundary=True)
    (ROOT/'results/revision/engine_tests.json').write_text(json.dumps(out,indent=2)+'\n');print(out)
if __name__=='__main__':main()
