"""Show recovery together with incorrect assignment under matched external LD."""
from pathlib import Path
import sys,csv
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tmp/runtime'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def plot_tradeoff():
    with (ROOT/'results/extension/summary.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    methods=['Bayes rho 0','Bayes rho 0.8','MultiSuSiE','Manc-COJO'];colors=['#2166ac','#67a9cf','#1b7837','#b2182b']
    plt.rcParams.update({'font.size':13,'axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':190})
    fig,axs=plt.subplots(2,2,figsize=(9,6.4),constrained_layout=True)
    for i,scenario in enumerate(['shared','heterogeneous']):
        for j,metric in enumerate(['mean_recall','mean_fdp']):
            ax=axs[i,j]
            for k,m in enumerate(methods):
                y=[float(next(r[metric] for r in rows if r['scenario']==scenario and r['reference']=='matched' and r['n_african']==n and r['method']==m)) for n in ['250','1000']]
                ax.plot([250,1000],y,'o-',color=colors[k],label=m)
            label='Exact causal variant recall' if metric=='mean_recall' else 'Mean false discovery proportion'
            ax.set(title=scenario.capitalize()+' effects',xlabel='African cohort size',ylabel=label,xticks=[250,1000])
            ax.set_ylim((-.02,1.02) if metric=='mean_recall' else (0,.38))
    axs[0,1].legend(fontsize=11);fig.savefig(ROOT/'figures/extension/decision_tradeoff.png');plt.close(fig)
if __name__=='__main__':plot_tradeoff()
