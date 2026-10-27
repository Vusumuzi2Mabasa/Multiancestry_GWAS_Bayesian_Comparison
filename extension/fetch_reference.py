"""Download public Phase 3 chr22 once; do not publish raw reference genotypes."""
from pathlib import Path
import hashlib,json,urllib.request,time,gzip
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/reference';OUT.mkdir(parents=True,exist_ok=True)
BASE='https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/'
FILES=['README_phase3_callset_20150220','README_known_issues_20200731',
 'README_integrated_call_samples_v3.20250704.ALL','integrated_call_samples_v3.20130502.ALL.panel',
 'ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz']
manifest=[]
for name in FILES:
    path=OUT/name
    if not path.exists():
        part=path.with_suffix(path.suffix+'.part')
        if name.endswith('.vcf.gz'):
            total=205612353
            # Bounded HTTP ranges survive intermittent connections without restarting.
            while not part.exists() or part.stat().st_size<total:
                start=part.stat().st_size if part.exists() else 0
                from concurrent.futures import ThreadPoolExecutor
                def chunk(bounds):
                    lo,hi=bounds
                    for attempt in range(8):
                        try:
                            req=urllib.request.Request(BASE+name,headers={'Range':f'bytes={lo}-{hi}'})
                            with urllib.request.urlopen(req,timeout=30) as response:
                                if response.status!=206 or not response.headers.get('Content-Range','').startswith(f'bytes {lo}-'):
                                    raise ValueError('Server did not honour resume range')
                                block=response.read()
                            if len(block)!=hi-lo+1:raise ValueError('Truncated range')
                            return block
                        except Exception as exc:
                            print(f'retry {attempt+1} at {lo}: {exc}',flush=True)
                            if attempt==7:raise
                            time.sleep(2)
                bounds=[(lo,min(lo+512*1024-1,total-1)) for lo in range(start,min(start+8*512*1024,total),512*1024)]
                with ThreadPoolExecutor(max_workers=8) as pool:
                    for bounds_item,block in zip(bounds,pool.map(chunk,bounds)):
                        with part.open('ab') as f:f.write(block)
                        print(f'{name}: {bounds_item[1]+1}/{total}',flush=True)
            if part.stat().st_size!=total:raise ValueError('Unexpected VCF size')
        else:
            with urllib.request.urlopen(BASE+name,timeout=120) as response,part.open('wb') as f:
                while block:=response.read(1024*1024):f.write(block)
        part.replace(path)
    h=hashlib.sha256()
    with path.open('rb') as f:
        while block:=f.read(1024*1024):h.update(block)
    if name.endswith('.vcf.gz'):
        with gzip.open(path,'rb') as f:
            while f.read(8*1024*1024):pass
    manifest.append(dict(file=name,url=BASE+name,bytes=path.stat().st_size,sha256=h.hexdigest()))
    print(name,path.stat().st_size,flush=True)
(ROOT/'extension/reference_downloads.json').write_text(json.dumps(manifest,indent=2)+'\n')
