"""Checks that target statistical calculations and binary encoding."""
import tempfile
from pathlib import Path
import numpy as np
from simulate_compare import bayes,write_bed,genotypes,summarize,IDS,M

rng=np.random.default_rng(901)
with tempfile.TemporaryDirectory() as td:
    p=Path(td)/'check';g=genotypes(103,.7,rng);write_bed(p,g)
    raw=p.with_suffix('.bed').read_bytes();assert raw[:3]==bytes([108,27,1])
    back=np.empty_like(g);decode={0:2,2:1,3:0};width=(len(g)+3)//4
    for j in range(M):
        b=raw[3+j*width:3+(j+1)*width]
        back[:,j]=[decode[(b[i//4]>>(2*(i%4)))&3] for i in range(len(g))]
    assert np.array_equal(g,back),'BED dosage/allele orientation mismatch'
    y=rng.normal(size=len(g));Q,q,s2=summarize(p,g,y,np.zeros(M))
    x=g-g.mean(axis=0);yc=y-y.mean()
    assert np.allclose(Q,x.T@x) and np.allclose(q,x.T@yc)
    residual=yc-x@np.linalg.solve(Q,q)
    assert np.allclose(s2,residual@residual/(len(g)-M-1))

# With orthogonal markers, rho=0 and zero scores, each marker's BF is
# 1/(1+tau^2*Qjj); enumeration must give its closed-form inclusion probability.
Q=np.eye(M)*100;q=np.zeros(M);pi=.2;tau=.15
pip,mean,sd=bayes([(Q,q,1),(Q,q,1)],rho=0,tau=tau,pi=pi)
bf=1/(1+tau*tau*100);expected=pi*bf/(1-pi+pi*bf)
assert np.allclose(pip,expected) and np.allclose(mean,0)
# Swapping populations must swap posterior effects and leave shared PIPs intact.
q1=rng.normal(size=M)*8;q2=rng.normal(size=M)*8
p1,m1,_=bayes([(Q,q1,1),(Q,q2,1)],rho=.8)
p2,m2,_=bayes([(Q,q2,1),(Q,q1,1)],rho=.8)
assert np.allclose(p1,p2) and np.allclose(m1,m2[::-1])
print('PASS: BED round-trip; summary likelihood reconstruction; analytic null; ancestry-swap invariance.')
