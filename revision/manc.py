"""Native Manc-COJO signed-LD interface, with explicit common-panel filtering."""
import csv,subprocess,math
from pathlib import Path
from engine import np,ROOT
from scipy.special import ndtr
def table(path,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    if not rows:return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def eligible(stats,freq,limit):return np.all(np.abs(np.stack([s['freq'] for s in stats])-freq)<=limit+1e-12,axis=0)
def entry_score(stats,keep):
    if not keep.any():return 0.
    w=np.array([1/s['se']**2 for s in stats]);b=np.array([s['b'] for s in stats])
    z=np.sum(w*b,axis=0)/np.sqrt(w.sum(0))
    return -math.log10(max(float(np.min(2*ndtr(-abs(z[keep])))),1e-300))
def run(work,stats,variants,rs,freq,keep,cutoff,screen=True):
    work=Path(work);work.mkdir(parents=True,exist_ok=True)
    sel=np.zeros(8,bool);eff=np.zeros((2,8))
    if not keep.any():return sel,eff,'no_eligible'
    if screen and entry_score(stats,keep)<=cutoff:return sel,eff,'analytic_no_entry'
    refs=[];mas=[]
    for k in range(2):
        pre=work/f'ld{k}';refs.append(str(pre));ma=work/f's{k}.ma';mas.append(str(ma))
        pre.with_suffix('.bim').write_text(''.join(f"22\t{v['id']}\t0\t{v['position']}\t{v['alt']}\t{v['ref']}\n" for v in variants),newline='\n')
        pre.with_suffix('.ld').write_text(''.join(f"{a['id']} {b['id']} {rs[k][i,j]:.17g}\n" for i,a in enumerate(variants) for j,b in enumerate(variants) if j>=i),newline='\n')
        table(pre.with_suffix('.frq'),[dict(CHR=22,SNP=v['id'],A1=v['alt'],A2=v['ref'],MAF=freq[k,j],NCHROBS=80) for j,v in enumerate(variants)])
        s=stats[k];table(ma,[dict(SNP=v['id'],A1=v['alt'],A2=v['ref'],freq=s['freq'][j],b=s['b'][j],se=s['se'][j],p=max(float(2*ndtr(-abs(s['b'][j]/s['se'][j]))),1e-300),N=s['n']) for j,v in enumerate(variants) if keep[j]])
    output=work/'fit.jma.cojo'
    if output.exists():output.unlink()
    cmd=[str(ROOT/'tools/manc_cojo_win.exe'),'--ld',*refs,'--cojo-file',*mas,'--cojo-slct','--cojo-p',str(10**(-cutoff)),'--diff-freq','1','1','--thread-num','1','--out',str(work/'fit')]
    r=subprocess.run(cmd,capture_output=True,text=True,timeout=45);log=r.stdout+'\n'+r.stderr
    (work/'last_stdout.txt').write_text(log)
    if r.returncode:raise RuntimeError(log[-2500:])
    if output.exists():
        lines=[x.split() for x in output.read_text().splitlines() if x.strip()];lookup={v['id']:i for i,v in enumerate(variants)}
        for row in [dict(zip(lines[0],x)) for x in lines[1:]]:
            j=lookup[row['SNP']];assert row['A1']==variants[j]['alt']
            sel[j]=True;eff[:,j]=[float(row['bJ.1']),float(row['bJ.2'])]
    elif not any(x in log.lower() for x in ['no snp','no significant','0 snp','no variant','not significant']):raise RuntimeError(log[-1500:])
    assert np.isfinite(eff).all()
    return sel,eff,'actual_manc'
