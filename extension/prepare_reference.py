"""Extract prespecified Phase 3 windows and split donors before phenotypes exist."""
from pathlib import Path
import sys, json, gzip, hashlib, csv,argparse
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tmp/runtime'))
import numpy as np

def table(path,rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--vcf',type=Path);args=ap.parse_args()
    cfg=json.loads((ROOT/'extension/design.json').read_text())
    raw=ROOT/'data/reference'; out=ROOT/'data/prepared';out.mkdir(exist_ok=True)
    pub=ROOT/'results/extension';pub.mkdir(parents=True,exist_ok=True)
    panel={}
    for line in (raw/'integrated_call_samples_v3.20130502.ALL.panel').read_text().splitlines()[1:]:
        t=line.split()
        if t[1] in cfg['populations']:panel[t[0]]=t[1]
    rng=np.random.default_rng(cfg['reference_seed']);splits={}; donors=[]
    for pop in cfg['populations']:
        names=np.array(sorted(k for k,v in panel.items() if v==pop))
        names=rng.permutation(names)
        splits[pop]={'reference':names[:40].tolist(),'study':names[40:].tolist()}
        assert not set(splits[pop]['reference']) & set(splits[pop]['study'])
        for role,ids in splits[pop].items():
            donors.append(dict(population=pop,role=role,donors=len(ids),
                donor_list_sha256=hashlib.sha256(('\n'.join(ids)+'\n').encode()).hexdigest()))
    (out/'donor_splits.json').write_text(json.dumps(splits,indent=2)+'\n')
    table(pub/'donor_counts.tsv',donors)
    windows=cfg['windows_grch37']; variants=[[] for _ in windows];haplotypes=[[] for _ in windows]
    names=list(panel); indices=None
    with gzip.open(args.vcf or raw/'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz','rt') as f:
        for line in f:
            if line.startswith('##'):continue
            if line.startswith('#CHROM'):
                fields=line.rstrip().split('\t')
                samples=fields[9:];indices=[9+samples.index(x) for x in names];continue
            pos=int(line.split('\t',2)[1])
            if pos>windows[-1][2]:break
            wi=next((i for i,(_,start,end) in enumerate(windows) if start<=pos<=end),None)
            if wi is None:continue
            fields=line.rstrip().split('\t')
            ref,alt=fields[3:5]
            if len(ref)!=1 or len(alt)!=1 or ref not in 'ACGT' or alt not in 'ACGT':continue
            if {ref,alt} in ({'A','T'},{'C','G'}) or fields[6]!='PASS':continue
            gt=[fields[j].split(':')[0] for j in indices]
            if any(x not in ('0|0','0|1','1|0','1|1') for x in gt):continue
            h=np.array([[int(x[0]),int(x[2])] for x in gt],dtype=np.int8)
            variants[wi].append(dict(chromosome=22,position=pos,id=f'22:{pos}:{ref}:{alt}',ref=ref,alt=alt,rsid=fields[2]))
            haplotypes[wi].append(h)
    qc=[]; selectedrows=[];ldrows=[]
    for wi in range(len(windows)):
        if not haplotypes[wi]:raise ValueError(f'No eligible variants in window {wi}')
        h=np.stack(haplotypes[wi],axis=2); blocks={}
        for pop in cfg['populations']:
            for role,ids in splits[pop].items():
                blocks[f'{pop}_{role}']=h[[names.index(x) for x in ids]].reshape(-1,h.shape[2])
        eligible=np.ones(h.shape[2],dtype=bool)
        for b in blocks.values():eligible&=np.minimum(b.mean(0),1-b.mean(0))>=cfg['minimum_maf_each_split']
        candidates=np.flatnonzero(eligible)
        # Coordinate-spaced traversal, then deterministic fill; no phenotype-informed selection.
        targets=np.linspace(windows[wi][1],windows[wi][2],cfg['markers_per_window'])
        positions=np.array([x['position'] for x in variants[wi]])
        ranked=[sorted(candidates,key=lambda j:abs(positions[j]-target)) for target in targets]
        order=[ranked[t][rank] for rank in range(min(20,len(candidates))) for t in range(len(targets))]
        order.extend(candidates.tolist()); chosen=[]
        for j in dict.fromkeys(order):
            if all(not chosen or np.max(np.corrcoef(b[:,chosen+[j]],rowvar=False)[-1,:-1]**2)<cfg['maximum_pairwise_r2_each_split'] for b in blocks.values()):
                chosen.append(j)
            if len(chosen)==cfg['markers_per_window']:break
        if len(chosen)!=cfg['markers_per_window']:raise ValueError('Too few QC-passing noncollinear markers')
        chosen=sorted(chosen);vs=[variants[wi][j] for j in chosen]
        for j,v in enumerate(vs):selectedrows.append(dict(locus=wi,index=j,**v))
        arrays={}
        for key,b in blocks.items():
            b=b[:,chosen];arrays[key]=b
            # LD from diploid reference donors, not artificially paired haplotypes.
            g=b[::2]+b[1::2];r=np.corrcoef(g,rowvar=False)
            if np.linalg.eigvalsh(r).min()<1e-6:raise ValueError(f'Singular reference {key}')
            arrays[key+'_R']=r
            for j,v in enumerate(vs):qc.append(dict(locus=wi,pool=key,variant=v['id'],maf=min(b[:,j].mean(),1-b[:,j].mean())))
            for j in range(len(vs)):
                for k in range(len(vs)):ldrows.append(dict(locus=wi,pool=key,i=j,j=k,r=r[j,k]))
        np.savez_compressed(out/f'locus_{wi}.npz',**arrays)
        (out/f'locus_{wi}.json').write_text(json.dumps(vs,indent=2)+'\n')
        print('locus',wi,'eligible',len(candidates),'selected',len(chosen),flush=True)
    table(pub/'variants.tsv',selectedrows);table(pub/'reference_qc.tsv',qc);table(pub/'ld_matrices.tsv',ldrows)
    (pub/'reference_validation.json').write_text(json.dumps(dict(donor_overlap=0,windows=len(windows),markers_per_window=8,
        build='GRCh37',allele_dosage='ALT',phase_complete=True,source='1000 Genomes Phase 3 v5b'),indent=2)+'\n')

if __name__=='__main__':main()
