"""Summaries and scientific figures from completed, independently evaluated cells."""
import sys, json, csv, math, itertools
from pathlib import Path
from collections import defaultdict
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tmp/runtime'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from benchmark import table, CFG, METHODS
OUT=ROOT/'results/extension';FIG=ROOT/'figures/extension'

def read(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def wilson(x,n):
    z=1.95996398454;p=x/n;d=1+z*z/n;c=(p+z*z/(2*n))/d
    h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return c-h,c+h

def summarize(rows,keys):
    groups=defaultdict(list)
    for r in rows:groups[tuple(r[k] for k in keys)].append(r)
    result=[]
    for key,rs in sorted(groups.items()):
        n=len(rs);false=sum(int(r['any_false_positive']) for r in rs);lo,hi=wilson(false,n)
        rec=[float(r['recall']) for r in rs if r['recall']!=''];mse=[float(r['mse_african']) for r in rs]
        result.append(dict(zip(keys,key),replicates=n,mean_recall=float(np.mean(rec)) if rec else '',
            false_selection_rate=false/n,false_lower=lo,false_upper=hi,rmse_african=math.sqrt(np.mean(mse)),
            mean_fdp=float(np.mean([int(r['fp'])/max(1,int(r['selected'])) for r in rs])),
            nonconverged=sum(1-int(r['converged']) for r in rs),
            population_specific_recall=float(np.mean([float(r['population_specific_recovered']) for r in rs if r['population_specific_recovered']!=''])) if rs[0]['scenario']=='heterogeneous' else ''))
    return result

def main():
    FIG.mkdir(parents=True,exist_ok=True);rows=[];vs=[];cal=[];thresholds=[]
    for i in range(16):
        p=OUT/'cells'/f'cell_{i:02d}';done=json.loads((p/'complete.json').read_text())
        r=read(p/'metrics.tsv');v=read(p/'variant_results.tsv');c=read(p/'calibration_scores.tsv')
        assert len(r)==1600 and len(v)==12800 and len(c)==1200,(i,len(r),len(v),len(c))
        assert done['evaluation_datasets']==400 and done['calibration_datasets']==300
        rows+=r;vs+=v;cal+=c
        t=json.loads((p/'thresholds.json').read_text())
        for method,score in t['thresholds'].items():thresholds.append(dict(cell_id=i,method=method,score=score,
            entry_p=10**(-score) if method=='Manc-COJO' else '',calibration_replicates=300))
    table(OUT/'replicate_metrics.tsv',rows);table(OUT/'thresholds.tsv',thresholds)
    summary=summarize(rows,['scenario','n_african','reference','method']);table(OUT/'summary.tsv',summary)
    table(OUT/'cell_summary.tsv',summarize(rows,['cell_id','locus','population','scenario','n_african','reference','method']))
    table(OUT/'calibration_convergence.tsv',[dict(method=m,nonconverged=sum(1-int(r['converged']) for r in cal if r['method']==m),
        fits=sum(r['method']==m for r in cal)) for m in METHODS])
    # Fixed-cell bootstrap: resample paired replicates within each locus/population.
    paired=[];rng=np.random.default_rng(673299)
    by=defaultdict(dict)
    for r in rows:by[(r['cell_id'],r['scenario'],r['replicate'])][r['method']]=r
    for scenario,n,reference in itertools.product(['shared','heterogeneous'],['250','1000'],['matched','swapped']):
        for method in METHODS[:-1]:
            pairs=[rs for rs in by.values() if rs[method]['scenario']==scenario and rs[method]['n_african']==n and rs[method]['reference']==reference]
            strata=defaultdict(list)
            for rs in pairs:strata[rs[method]['cell_id']].append(rs)
            for metric in ['recall','mse_african']:
                arrays=[np.array([float(rs[method][metric])-float(rs['Manc-COJO'][metric]) for rs in group]) for group in strata.values()]
                boot=np.mean([a[rng.integers(len(a),size=(2000,len(a)))].mean(1) for a in arrays],axis=0)
                paired.append(dict(scenario=scenario,n_african=n,reference=reference,method=method,comparator='Manc-COJO',metric=metric,
                    difference=np.mean([a.mean() for a in arrays]),lower=np.quantile(boot,.025),upper=np.quantile(boot,.975),paired_datasets=len(pairs)))
    table(OUT/'paired_differences.tsv',paired)
    meta={r['cell_id']:r for r in rows};effectgroups=defaultdict(lambda:defaultdict(list))
    for v in vs:
        if v['scenario']=='heterogeneous' and v['index']=='6':
            r=meta[v['cell_id']]
            effectgroups[(r['n_african'],r['reference'],v['method'])][v['cell_id']].append(v)
    effectsummary=[]
    for (n,reference,method),groups in sorted(effectgroups.items()):
        a=[np.array([float(v['estimate_african']) for v in group]) for group in groups.values()]
        ceu=[float(v['estimate_ceu']) for group in groups.values() for v in group]
        mean=float(np.mean([x.mean() for x in a]));se=float(np.sqrt(sum(x.var(ddof=1)/len(x) for x in a)/len(a)**2))
        effectsummary.append(dict(n_african=n,reference=reference,method=method,true_african=.28,mean_african=mean,bias_african=mean-.28,
            mean_lower=mean-1.95996398454*se,mean_upper=mean+1.95996398454*se,
            rmse_african=float(np.sqrt(np.mean((np.concatenate(a)-.28)**2))),
            true_ceu=0,mean_ceu=float(np.mean(ceu)),replicates=sum(map(len,a))))
    table(OUT/'population_specific_effects.tsv',effectsummary)
    plt.rcParams.update({'font.size':13,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})
    colors=['#2166ac','#67a9cf','#1b7837','#b2182b']
    from plot_reference import plot_ld
    plot_ld()
    plt.rcParams.update({'font.size':13})
    fig,axs=plt.subplots(1,2,figsize=(9,4.3),constrained_layout=True)
    for ax,n in zip(axs,['250','1000']):
        for j,m in enumerate(METHODS):
            rs=[next(r for r in summary if r['scenario']=='null' and r['n_african']==n and r['reference']==ref and r['method']==m) for ref in ['matched','swapped']]
            x=np.arange(2)+(j-1.5)*.14;y=np.array([r['false_selection_rate'] for r in rs]);lo=np.array([r['false_lower'] for r in rs]);hi=np.array([r['false_upper'] for r in rs])
            ax.errorbar(x,y,yerr=[y-lo,hi-y],fmt='o',color=colors[j],label=m,capsize=3)
        ax.axhline(.05,color='black',ls='--',lw=1);ax.set(xticks=[0,1],xticklabels=['Matched','Swapped'],ylabel='Null locus false selection rate',title=f'African cohort n = {n}');ax.set_ylim(bottom=0)
    axs[1].legend(fontsize=11);fig.savefig(FIG/'null_calibration.png');plt.close(fig)
    fig,axs=plt.subplots(2,2,figsize=(9,6.4),constrained_layout=True)
    for i,scenario in enumerate(['shared','heterogeneous']):
        for j,reference in enumerate(['matched','swapped']):
            ax=axs[i,j]
            for k,m in enumerate(METHODS):
                ys=[next(r['mean_recall'] for r in summary if r['scenario']==scenario and r['n_african']==n and r['reference']==reference and r['method']==m) for n in ['250','1000']]
                ax.plot([250,1000],ys,'o-',color=colors[k],label=m)
            ax.set(title=f'{scenario.capitalize()} · {reference} LD',xlabel='African cohort size',ylabel='Exact causal variant recall',ylim=(-.02,1.02),xticks=[250,1000])
    axs[0,1].legend(fontsize=11);fig.savefig(FIG/'signal_recovery.png');plt.close(fig)
    from plot_tradeoff import plot_tradeoff
    plot_tradeoff()
    fig,axs=plt.subplots(2,2,figsize=(9,7),constrained_layout=True)
    for ax,scenario in zip(axs[0],['shared','heterogeneous']):
        for k,m in enumerate(METHODS):
            ys=[next(r['rmse_african'] for r in summary if r['scenario']==scenario and r['n_african']==n and r['reference']=='matched' and r['method']==m) for n in ['250','1000']]
            ax.plot([250,1000],ys,'o-',color=colors[k],label=m)
        ax.set(title=scenario.capitalize()+' effects · matched LD',xlabel='African cohort size',ylabel='African effect RMSE',xticks=[250,1000])
    for ax,reference in zip(axs[1],['matched','swapped']):
        for k,m in enumerate(METHODS):
            es=[next(r for r in effectsummary if r['n_african']==n and r['reference']==reference and r['method']==m) for n in ['250','1000']]
            y=np.array([r['mean_african'] for r in es]);lo=np.array([r['mean_lower'] for r in es]);hi=np.array([r['mean_upper'] for r in es])
            ax.errorbar([250,1000],y,yerr=[y-lo,hi-y],fmt='o-',color=colors[k],label=m,capsize=3)
        ax.axhline(.28,color='black',ls='--',lw=1)
        ax.set(title='Population-specific · '+reference+' LD',xlabel='African cohort size',ylabel='Mean African effect estimate',xticks=[250,1000])
    axs[0,1].legend(fontsize=11);fig.savefig(FIG/'effect_error.png');plt.close(fig)
    validation=dict(cells=16,calibration_datasets=4800,held_out_null_datasets=3200,alternative_datasets=3200,
        evaluation_method_fits=len(rows),variant_results=len(vs),all_counts_pass=True,
        nonconverged_evaluation_fits=sum(1-int(r['converged']) for r in rows),calibration_and_evaluation_seed_namespaces_disjoint=True)
    (OUT/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation,indent=2))

if __name__=='__main__':main()
