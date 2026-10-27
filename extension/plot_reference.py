"""Plot the actual external LD panels and record descriptive matrix checks."""
from pathlib import Path
import sys,csv
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tmp/runtime'))
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
def plot_ld():
    out=ROOT/'figures/extension';out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':14,'savefig.dpi':190})
    fig,axs=plt.subplots(2,3,figsize=(10,6.5),constrained_layout=True);rows=[]
    for locus in range(2):
        a=np.load(ROOT/f'data/prepared/locus_{locus}.npz')
        for j,pop in enumerate(['CEU','YRI','LWK']):
            r=a[pop+'_reference_R'];ax=axs[locus,j]
            im=ax.imshow(r,vmin=-1,vmax=1,cmap='RdBu_r')
            ax.set(title=f'Locus {locus+1} · {pop}',xlabel='Variant index',ylabel='Variant index')
            ax.set_xticks(range(8),range(1,9));ax.set_yticks(range(8),range(1,9))
            r2=r[np.triu_indices(8,1)]**2
            rows.append(dict(locus=locus,population=pop,reference_donors=40,maximum_pairwise_r2=float(r2.max()),
                mean_pairwise_r2=float(r2.mean()),minimum_eigenvalue=float(np.linalg.eigvalsh(r).min())))
    fig.colorbar(im,ax=axs,shrink=.7,label='Signed LD correlation r')
    fig.savefig(out/'reference_ld.png');plt.close(fig)
    with (ROOT/'results/extension/ld_summary.tsv').open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    return rows
if __name__=='__main__':print(plot_ld())
