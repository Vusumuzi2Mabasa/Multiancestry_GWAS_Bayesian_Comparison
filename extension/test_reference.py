"""Audit real-panel allele encoding, donor separation, frequencies and LD."""
import json, hashlib
import benchmark as b
import numpy as np

def main():
    prepared=b.ROOT/'data/prepared';splits=json.loads((prepared/'donor_splits.json').read_text())
    assert all(not set(x['reference'])&set(x['study']) for x in splits.values())
    checks=[];work=b.ROOT/'tmp/reference_tests';work.mkdir(exist_ok=True)
    for locus in range(2):
        arrays=dict(np.load(prepared/f'locus_{locus}.npz'))
        variants=json.loads((prepared/f'locus_{locus}.json').read_text())
        assert len({v['id'] for v in variants})==8
        assert all({v['ref'],v['alt']} not in ({'A','T'},{'C','G'}) for v in variants)
        # Recreate separately so the reused task-local ref1 prefix is checked immediately.
        for pop in ['CEU','YRI','LWK']:
            prefix=b.reference_files(work,arrays,variants,['CEU',pop] if pop!='CEU' else ['CEU','CEU'],'matched')[1]
            raw=prefix.with_suffix('.bed').read_bytes();assert raw[:3]==b'\x6c\x1b\x01'
            codes=np.frombuffer(raw[3:],dtype=np.uint8).reshape(8,10)
            expanded=np.stack([(codes>>shift)&3 for shift in [0,2,4,6]],axis=2).reshape(8,40).T
            assert not np.any(expanded==1),'Unexpected BED missing value'
            decoded=np.array([2,-99,1,0])[expanded]
            h=arrays[pop+'_reference'];expected=h[::2]+h[1::2]
            assert np.array_equal(decoded,expected)
            assert np.allclose(np.corrcoef(decoded,rowvar=False),arrays[pop+'_reference_R'],atol=1e-12)
            a1=[line.split()[4] for line in prefix.with_suffix('.bim').read_text().splitlines()]
            assert a1==[v['alt'] for v in variants]
        for key,h in arrays.items():
            if key.endswith('_R'):continue
            maf=np.minimum(h.mean(0),1-h.mean(0));assert np.min(maf)>=.05
            r=arrays[key+'_R'];mineig=float(np.linalg.eigvalsh(r).min());assert mineig>1e-6
            checks.append(dict(locus=locus,pool=key,haplotypes=len(h),minimum_maf=float(maf.min()),minimum_eigenvalue=mineig,
                haplotype_array_sha256=hashlib.sha256(h.tobytes()).hexdigest()))
    report=dict(all_passed=True,donor_overlap=0,bed_alt_dosage_roundtrip=True,reference_ld_roundtrip=True,pools=checks)
    (b.ROOT/'results/extension/reference_tests.json').write_text(json.dumps(report,indent=2)+'\n')
    print('PASS: donor separation, allele alignment, BED dosage roundtrip, LD, MAF and positive definiteness')
if __name__=='__main__':main()
