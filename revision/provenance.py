"""Public software identities without personal environment variables or paths."""
from engine import *
import json,hashlib,subprocess,platform,importlib.metadata
def main():
    rp=os.environ.get('RSCRIPT','C:/Program Files/R/R-4.5.1/bin/Rscript.exe')
    expr="suppressPackageStartupMessages(library(INLA)); cat(as.character(getRversion()),as.character(packageVersion('INLA')),as.character(packageVersion('Matrix')),as.character(packageVersion('jsonlite')),sep='|')"
    r=subprocess.run([rp,'--vanilla','-e',expr],capture_output=True,text=True,check=True)
    vals=r.stdout.strip().split('|');assert len(vals)==4
    out=dict(python=platform.python_version(),R=vals[0],INLA=vals[1],Matrix=vals[2],jsonlite=vals[3],platform='Windows',R_startup='--vanilla',
      packages={p:importlib.metadata.version(p) for p in ['numpy','scipy','numba','MultiSuSiE','matplotlib']},
      multisusie_source_commit='a89e091d45297d983118438f116c3eb5bc1219ed',manc_sha256=hashlib.sha256((ROOT/'tools/manc_cojo_win.exe').read_bytes()).hexdigest(),
      reference_vcf_sha256=hashlib.sha256((ROOT/'data/reference/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz').read_bytes()).hexdigest(),
      protocol_sha256=hashlib.sha256((ROOT/'revision/protocol.json').read_bytes()).hexdigest(),
      software_citations=['10.1038/s41588-026-02799-1','10.1038/s41588-025-02450-5','10.1111/j.1467-9868.2008.00700.x'])
    (ROOT/'results/revision/software_provenance.json').write_text(json.dumps(out,indent=2)+'\n',newline='\n');print(out)
if __name__=='__main__':main()
