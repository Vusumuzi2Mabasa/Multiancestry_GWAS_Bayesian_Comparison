"""Six prespecified public-reference regions; preserve the original benchmark files."""
from pathlib import Path
import sys,json,gzip,csv,hashlib
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'tmp/runtime'))
import numpy as np
def table(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def main():
    cfg=json.loads((ROOT/'revision/protocol.json').read_text());windows=cfg['reference_windows_grch37']
    raw=ROOT/'data/reference';out=ROOT/'data/revision/reference';out.mkdir(parents=True,exist_ok=True)
    pub=ROOT/'results/revision/reference';pub.mkdir(parents=True,exist_ok=True)
    splits=json.loads((ROOT/'data/prepared/donor_splits.json').read_text())
    names=[x for p in ['CEU','YRI','LWK'] for role in ['reference','study'] for x in splits[p][role]]
    assert len(names)==len(set(names))
    records=[[] for w in windows];hs=[[] for w in windows]
    with gzip.open(raw/'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz','rt') as f:
        for line in f:
            if line.startswith('##'):continue
            if line.startswith('#CHROM'):
                samples=line.rstrip().split('\t')[9:];ix=[9+samples.index(n) for n in names];continue
            pos=int(line.split('\t',2)[1])
            if pos>windows[-1][2]:break
            wi=next((i for i,(_,start,end) in enumerate(windows) if start<=pos<=end),None)
            if wi is None:continue
            t=line.rstrip().split('\t');ref,alt=t[3:5]
            if t[6]!='PASS' or len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT' or {ref,alt} in ({'A','T'},{'C','G'}):continue
            gt=[t[j].split(':')[0] for j in ix]
            if any(g not in ('0|0','0|1','1|0','1|1') for g in gt):continue
            records[wi].append(dict(id=f'22:{pos}:{ref}:{alt}',position=pos,ref=ref,alt=alt))
            hs[wi].append(np.array([[int(g[0]),int(g[2])] for g in gt],dtype=np.int8))
    allvariants=[];qc=[];ld=[];identities=[]
    for wi in range(len(windows)):
        h=np.stack(hs[wi],axis=2);blocks={}
        for pop in splits:
            for role,ns in splits[pop].items():
                blocks[f'{pop}_{role}']=h[[names.index(x) for x in ns]].reshape(-1,h.shape[2])
        good=np.ones(h.shape[2],dtype=bool)
        for b in blocks.values():good&=np.minimum(b.mean(0),1-b.mean(0))>=.05
        candidates=np.flatnonzero(good);pos=np.array([r['position'] for r in records[wi]])
        targets=np.linspace(windows[wi][1],windows[wi][2],8)
        ranked=[sorted(candidates,key=lambda j:abs(pos[j]-t)) for t in targets]
        order=[ranked[t][rank] for rank in range(min(20,len(candidates))) for t in range(8)]+candidates.tolist();chosen=[]
        for j in dict.fromkeys(order):
            use=chosen+[j];accept=True
            if chosen:
                for b in blocks.values():
                    if np.max(np.corrcoef(b[:,use],rowvar=False)[-1,:-1]**2)>=.95:accept=False;break
                    g=b[::2,use]+b[1::2,use]
                    if np.linalg.eigvalsh(np.corrcoef(g,rowvar=False)).min()<1e-5:accept=False;break
            if accept:chosen.append(j)
            if len(chosen)==8:break
        assert len(chosen)==8,(wi,len(chosen))
        chosen=sorted(chosen);variants=[records[wi][j] for j in chosen];arrays={}
        for j,v in enumerate(variants):allvariants.append(dict(region=wi,index=j,**v))
        for key,b in blocks.items():
            x=b[:,chosen];g=x[::2]+x[1::2];r=np.corrcoef(g,rowvar=False)
            arrays[key]=x;arrays[key+'_R']=r
            qc.append(dict(region=wi,pool=key,donors=len(x)//2,min_maf=float(np.minimum(x.mean(0),1-x.mean(0)).min()),min_eigenvalue=float(np.linalg.eigvalsh(r).min())))
            identities.append(dict(region=wi,pool=key,sha256=hashlib.sha256(x.tobytes()).hexdigest()))
            for a in range(8):
                for b1 in range(8):ld.append(dict(region=wi,pool=key,i=a,j=b1,r=r[a,b1]))
        np.savez_compressed(out/f'region_{wi}.npz',**arrays)
        (out/f'region_{wi}.json').write_text(json.dumps(variants,indent=2)+'\n',newline='\n')
        print(f'region {wi}: {len(candidates)} common candidates, 8 retained',flush=True)
    table(pub/'variants.tsv',allvariants);table(pub/'qc.tsv',qc);table(pub/'ld.tsv',ld)
    (pub/'validation.json').write_text(json.dumps(dict(all_passed=True,regions=6,markers=48,donor_overlap=0,arrays=identities),indent=2)+'\n',newline='\n')
if __name__=='__main__':main()
