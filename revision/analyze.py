"""Reproducible summaries, paired Monte Carlo contrasts and manuscript figures."""
from engine import *
import json,csv,collections
from manc import table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
COLORS=['#254e70','#bd6b29','#2d8266','#8b4d78','#c34343','#6c6c6c']
def read(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def num(rows,k):return np.array([float(r[k]) for r in rows if r.get(k,'')!=''])
def avg(rows,k):
    a=num(rows,k);return float(a.mean()) if len(a) else float('nan')
def mcse(rows,k):
    groups=collections.defaultdict(list)
    for r in rows:
        if r.get(k,'')!='':groups[(r['region'],r['population'])].append(float(r[k]))
    return float(np.sqrt(sum(np.var(a,ddof=1)/len(a) for a in groups.values()))/len(groups)) if groups else float('nan')
def primary(r):return r['variance']=='total' and ((r['method']=='Exact' and r['tau']=='0.15' and r['pi']=='0.2' and r['rho'] in ['0','0.8']) or (r['method']=='MultiSuSiE' and r['L']=='5' and r['rho']=='0.75') or (r['method']=='Manc' and r['freq']=='0.2'))
def label(r):return 'Exact rho '+r['rho'] if r['method']=='Exact' else ('MultiSuSiE L5' if r['method']=='MultiSuSiE' else 'Manc-COJO')
def summarize(rows,keys):
    groups=collections.defaultdict(list)
    for r in rows:groups[tuple(r[k] for k in keys)].append(r)
    output=[]
    for key,rr in groups.items():
        out=dict(zip(keys,key));out['datasets']=len(rr)
        for k in ['any_false_positive','recall','fdp','ld_recall_ceu','ld_recall_african','mse_african','specific_estimate','specific_bias','retained_causal','retained_noncausal','pip_causal_mean','pip_noncausal_mean','pip_noncausal_max','seconds']:
            out[k]=avg(rr,k);out[k+'_mcse']=mcse(rr,k)
        out['rmse_african']=math.sqrt(out['mse_african'])
        out['rmse_african_mcse']=out['mse_african_mcse']/(2*out['rmse_african']) if out['rmse_african'] else 0.
        for group in ['causal','noncausal']:
            denom=num(rr,'n_'+group).sum();covered=num(rr,'covered_'+group);width=num(rr,'width_'+group)
            out['coverage_'+group]=float(covered.sum()/denom) if len(covered) and denom else float('nan')
            out['width_'+group]=float(width.sum()/denom) if len(width) and denom else float('nan')
            clustered=[dict(region=r['region'],population=r['population'],coverage=float(r['covered_'+group])/float(r['n_'+group])) for r in rr if r.get('covered_'+group,'')!='' and float(r['n_'+group])>0]
            out['coverage_'+group+'_mcse']=mcse(clustered,'coverage')
        out['nonconverged']=sum(r.get('converged')=='0' for r in rr);output.append(out)
    return output
def save(fig,name):
    fig.tight_layout();fig.savefig(ROOT/f'figures/revision/{name}.png',dpi=190,bbox_inches='tight');plt.close(fig)
def main(baseline_only=False):
    out=ROOT/'results/revision';figdir=ROOT/'figures/revision';figdir.mkdir(parents=True,exist_ok=True)
    sensitivity=[r for p in sorted((out/'sensitivity').glob('*/metrics.tsv')) for r in read(p)]
    validation=[r for p in sorted((out/'validation').glob('*/metrics.tsv')) for r in read(p)]
    for r in validation:r['name']=label(r)
    vs=summarize(validation,['name','ld','scenario']);table(out/'validation_summary.tsv',vs)
    ss=summarize(sensitivity,['arm','method','ld','variance','rho','tau','pi','L','freq','scenario']);table(out/'sensitivity_summary.tsv',ss)
    # Null calibration data and PIP distributions are distinct from alternative FDP.
    thresholds=[]
    for phase in ['sensitivity','validation']:
        for p in sorted((out/phase).glob('*/thresholds.json')):
            d=json.loads(p.read_text());
            for a in d['arms']:thresholds.append(dict(phase=phase,stratum=p.parent.name,**a,cutoff=d['cutoffs'][str(a['arm'])]))
    table(out/'calibrated_thresholds.tsv',thresholds)
    # Paired changes relative to the primary prior, preserving reference/population strata.
    base={}
    for r in sensitivity:
        if primary(r) and r['ld']=='matched':base[(r['region'],r['population'],r['scenario'],r['replicate'],label(r))]=r
    paired=[]
    for arm in sorted(set(r['arm'] for r in sensitivity),key=int):
        rr=[r for r in sensitivity if r['arm']==arm];a=rr[0]
        if a['method']=='Manc':name='Manc-COJO'
        elif a['method']=='MultiSuSiE':name='MultiSuSiE L5'
        else:name='Exact rho '+a['rho']
        for scenario in ['null','shared','heterogeneous','one','three']:
            changes=[]
            for r in rr:
                if r['scenario']!=scenario:continue
                b=base[(r['region'],r['population'],scenario,r['replicate'],name)]
                changes.append(dict(region=r['region'],population=r['population'],**{k:float(r[k])-float(b[k]) for k in ['mse_african','fdp','specific_estimate'] if r[k]!=''}))
            for k in ['mse_african','fdp','specific_estimate']:
                if any(k in c for c in changes):paired.append(dict(arm=arm,scenario=scenario,metric=k,paired_change=avg(changes,k),mcse=mcse(changes,k)))
    table(out/'paired_sensitivity_contrasts.tsv',paired)
    names=['Exact rho 0','Exact rho 0.8','MultiSuSiE L5','Manc-COJO']
    plt.rcParams.update({'font.size':13.5,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(11,7))
    ld=read(out/'reference/ld.tsv')
    for region,ax in enumerate(axes.flat):
        rr=[r for r in ld if r['region']==str(region) and r['pool']=='YRI_reference'];R=np.array([float(r['r']) for r in rr]).reshape(8,8)
        im=ax.imshow(R,vmin=-1,vmax=1,cmap='RdBu_r');ax.set_title(f'Region {region+1}: {20+5*region} Mb');ax.set_xticks([0,3,7],[1,4,8]);ax.set_yticks([0,3,7],[1,4,8])
    fig.subplots_adjust(right=.89,hspace=.4);cax=fig.add_axes([.92,.2,.02,.6]);fig.colorbar(im,cax=cax,label='Signed correlation');fig.savefig(figdir/'reference_ld.png',dpi=190,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(12,4))
    for i,name in enumerate(names):
        rr=[r for r in validation if r['name']==name and r['ld']=='matched' and r['scenario']=='null'];p=avg(rr,'any_false_positive');se=mcse(rr,'any_false_positive')
        axes[0].errorbar(i,p,yerr=1.96*se,fmt='o',color=COLORS[i],capsize=3)
        for ax,k in zip(axes[1:],['recall','fdp']):
            y=[next(r[k] for r in vs if r['name']==name and r['ld']=='matched' and r['scenario']==s) for s in ['one','shared','heterogeneous','three']]
            ax.plot(range(4),y,'o-',label=name,color=COLORS[i])
            if k=='fdp':
                control=[next(r[k] for r in vs if r['name']==name and r['ld']=='in_sample' and r['scenario']==s) for s in ['one','shared','heterogeneous','three']]
                ax.plot(range(4),control,'--',color=COLORS[i],alpha=.55)
    axes[0].axhline(.05,ls='--',color='gray');axes[0].set_xticks(range(4),['Exact 0','Exact .8','Multi L5','Manc'],rotation=30);axes[0].set_ylabel('Held-out null false-selection rate')
    for ax,k in zip(axes[1:],['Exact causal recall','Mean false discovery proportion']):ax.set_xticks(range(4),['One','Shared','Hetero.','Three'],rotation=25);ax.set_ylabel(k);ax.set_ylim(-.02,1.02)
    axes[2].legend(fontsize=11);save(fig,'calibration_recovery')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for scenario,col in zip(['one','shared','heterogeneous','three'],COLORS):
        vals=[]
        for L in [1,2,3,5]:
            rr=[r for r in sensitivity if r['method']=='MultiSuSiE' and r['ld']=='matched' and r['variance']=='total' and r['rho']=='0.75' and int(r['L'])==L and r['scenario']==scenario]
            vals.append(avg(rr,'fdp'))
        axes[0].plot([1,2,3,5],vals,'o-',label=scenario,color=col)
        vals=[]
        for rho in [0,.5,.75,.95]:
            rr=[r for r in sensitivity if r['method']=='MultiSuSiE' and r['ld']=='matched' and r['variance']=='total' and float(r['rho'])==rho and r['L']=='5' and r['scenario']==scenario];vals.append(math.sqrt(avg(rr,'mse_african')))
        axes[1].plot([0,.5,.75,.95],vals,'o-',label=scenario,color=col)
    axes[0].set(xlabel='Maximum single-effect components L',ylabel='Mean false discovery proportion');axes[1].set(xlabel='Prior effect correlation',ylabel='African all-marker RMSE');axes[0].legend(fontsize=11);save(fig,'multisusie_sensitivity')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    pr=[r for r in sensitivity if r['method']=='Exact' and r['ld']=='matched' and r['variance']=='total' and r['rho']=='0.8' and r['scenario']=='shared']
    specs=[(.2,.15),(.05,.15),(.4,.15),(.2,.075),(.2,.3)]
    for i,(pi,tau) in enumerate(specs):
        rr=[r for r in pr if float(r['pi'])==pi and float(r['tau'])==tau];axes[0].bar(i,avg(rr,'fdp'),color=COLORS[i])
    axes[0].set_xticks(range(5),['.2/.15','.05/.15','.4/.15','.2/.075','.2/.3'],rotation=30);axes[0].set(xlabel='Inclusion prior / slab SD',ylabel='Shared-effect mean FDP')
    for i,name in enumerate(names[:3]):
        rr=[r for r in validation if r['name']==name and r['ld']=='matched' and r['scenario']=='shared'];axes[1].boxplot(num(rr,'pip_noncausal_max'),positions=[i],widths=.5,showfliers=False)
    axes[1].set_xticks(range(3),['Exact 0','Exact .8','Multi L5']);axes[1].set_ylabel('Maximum noncausal PIP per dataset');save(fig,'prior_pip_sensitivity')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for i,name in enumerate(names[:3]):
        vals=[]
        for variance in ['total','known_noise']:
            rr=[r for r in sensitivity if primary(dict(r,variance='total')) and label(r)==name and r['ld']=='matched' and r['variance']==variance and r['scenario']=='heterogeneous'];vals.append(avg(rr,'specific_estimate'))
        axes[0].plot([0,1],vals,'o-',color=COLORS[i],label=name)
    axes[0].axhline(.28,ls='--',color='gray');axes[0].set_xticks([0,1],['Total phenotype variance','Known noise variance']);axes[0].set_ylabel('Mean population-specific effect');axes[0].legend(fontsize=11)
    for ld,style in [('matched','o-'),('swapped','s--')]:
        yy=[];ret=[]
        for freq in [.1,.2,1]:
            rr=[r for r in sensitivity if r['method']=='Manc' and r['ld']==ld and float(r['freq'])==freq and r['scenario']=='heterogeneous'];yy.append(avg(rr,'recall'));ret.append(avg(rr,'retained_causal')/2)
        axes[1].plot([.1,.2,1],yy,style,label=ld+' recall');axes[1].plot([.1,.2,1],ret,style,alpha=.45,label=ld+' causal retention')
    axes[1].set(xlabel='Frequency difference limit',ylabel='Proportion',ylim=(-.02,1.02));axes[1].legend(fontsize=11);save(fig,'variance_frequency_sensitivity')
    if baseline_only:
        print('Baseline summaries and five figures written; INLA results not read')
        return
    # INLA comparison is restricted to exactly the prespecified common subset.
    ir=read(out/'inla/metrics.tsv');pairedrows=[]
    for r in validation:
        if int(r['replicate'])<20 and r['scenario'] in ['null','shared','heterogeneous']:
            r=dict(r,method=r['name']);pairedrows.append(r)
    pairedrows+=ir;isum=summarize(pairedrows,['method','ld','scenario']);table(out/'inla/comparison_summary.tsv',isum)
    methods=names+['Fixed Gaussian','Adaptive INLA']
    # Paired primary comparisons use the same dataset, population and LD source.
    paired_lookup={(r['region'],r['population'],r['scenario'],r['replicate'],r['ld'],r['method']):r for r in pairedrows}
    contrasts=[]
    for comparator in methods[:-1]:
        for ld in ['in_sample','matched']:
            for scenario in ['null','shared','heterogeneous']:
                differences=[]
                for r in pairedrows:
                    if r['method']!='Adaptive INLA' or r['ld']!=ld or r['scenario']!=scenario:continue
                    b=paired_lookup[(r['region'],r['population'],scenario,r['replicate'],ld,comparator)]
                    differences.append(dict(region=r['region'],population=r['population'],mse=float(r['mse_african'])-float(b['mse_african'])))
                delta=avg(differences,'mse');se=mcse(differences,'mse')
                contrasts.append(dict(comparator=comparator,ld=ld,scenario=scenario,datasets=len(differences),mse_difference=delta,mcse=se,mc_lower=delta-1.96*se,mc_upper=delta+1.96*se))
    table(out/'inla/paired_contrasts.tsv',contrasts)
    fig,axes=plt.subplots(1,2,figsize=(11,4))
    for i,name in enumerate(methods):
        yy=[next(r['rmse_african'] for r in isum if r['method']==name and r['ld']=='matched' and r['scenario']==s) for s in ['null','shared','heterogeneous']]
        se=[next(r['rmse_african_mcse'] for r in isum if r['method']==name and r['ld']=='matched' and r['scenario']==s) for s in ['null','shared','heterogeneous']]
        axes[0].errorbar(range(3),yy,yerr=1.96*np.array(se),fmt='o-',label=name,color=COLORS[i],capsize=2)
        row=next(r for r in isum if r['method']==name and r['ld']=='matched' and r['scenario']=='heterogeneous');axes[1].bar(i,row['specific_estimate'],color=COLORS[i])
    axes[0].set_xticks(range(3),['Null','Shared','Heterogeneous']);axes[0].set_ylabel('African all-marker RMSE');axes[0].legend(fontsize=11)
    axes[1].axhline(.28,ls='--',color='gray');axes[1].set_xticks(range(6),['Exact 0','Exact .8','Multi L5','Manc','Fixed','INLA'],rotation=35);axes[1].set_ylabel('Mean population-specific effect');save(fig,'inla_effect_estimation')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    for i,name in enumerate([m for m in methods if m!='Manc-COJO']):
        rr=[r for r in isum if r['method']==name and r['ld']=='matched' and r['scenario']=='heterogeneous'][0]
        y=rr['coverage_causal'];delta=1.96*rr['coverage_causal_mcse']
        axes[0].bar(i,y,yerr=np.array([[min(y,delta)],[min(1-y,delta)]]),capsize=3,color=COLORS[methods.index(name)]);axes[1].bar(i,rr['width_causal'],color=COLORS[methods.index(name)])
    for ax in axes:ax.set_xticks(range(5),['Exact 0','Exact .8','Multi L5','Fixed','INLA'],rotation=30)
    axes[0].axhline(.95,ls='--',color='gray');axes[0].set_ylabel('Causal 95% interval coverage');axes[0].set_ylim(0,1.02);axes[1].set_ylabel('Mean causal interval width');save(fig,'interval_coverage')
    print('Summaries and seven manuscript figures written')
if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--baseline-only',action='store_true');args=ap.parse_args();main(args.baseline_only)
