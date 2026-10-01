"""Direct complex-state check, independent of invariant moment matrices.
This is a numerical cross-check, not the exact proof supplied separately.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from functools import lru_cache
import json,sys,numpy as np,sympy as s
from sympy.physics.wigner import clebsch_gordan
root=Path(__file__).resolve().parent;p=root/'global_purity2';sys.path.insert(0,str(root))
from angular_momentum import tensors
T=tensors(5);j=s.Rational(5,2);data=json.loads((p/'features.json').read_text());cert=json.loads((p/'certificate.json').read_text());C={};Q={}
for bid,d in data.items():
 C[bid]=np.diag([float(s.sympify(z)) for z in d['sigma']])@np.array(s.Matrix(d['K']),dtype=float)@np.array(s.Matrix(cert[bid]['W']),dtype=float);Q[bid]=np.array(s.Matrix(cert[bid]['Q']),dtype=float)
@lru_cache(None)
def cg(a,b,c,m,n):return float(clebsch_gordan(a,b,c,m,n,m+n))
@lru_cache(None)
def pairs(L1,L2,K,M):return [(m,M-m,cg(L1,L2,K,m,M-m)) for m in range(-L1,L1+1) if abs(M-m)<=L2]
def rhs(rho,x):
 @lru_cache(None)
 def f(lab,M):
  L1,L2,K=lab
  if L2==0:return 1.
  if L1==0:return x[L2,M]
  return sum(c*x[L1,m]*x[L2,n] for m,n,c in pairs(L1,L2,K,M))
 total=0.
 for bid,d in data.items():
  twoJ,_,fs=d['meta'];J=s.Rational(twoJ,2);G=np.zeros(Q[bid].shape,complex)
  for twoM in range(-twoJ,twoJ+1,2):
   M=s.Rational(twoM,2);F=np.zeros((len(fs),6),complex)
   for i,lab in enumerate(fs):
    K=lab[2]
    for a in range(6):
     mu=j-a;m=int(M-mu)
     if abs(m)<=K:F[i,a]=cg(j,K,J,mu,m)*f(tuple(lab),m)
   V=C[bid].T@F;G+=V@rho@V.conj().T/(twoJ+1)
  total+=np.trace(Q[bid]@G).real
 return total
rng=np.random.default_rng(187625);worst=0.
for trial in range(24):
 x={}
 for L in [3,4,5]:
  x[L,0]=rng.normal()
  for m in range(1,L+1):x[L,m]=rng.normal()+1j*rng.normal();x[L,-m]=(-1)**m*x[L,m].conjugate()
 X=sum(v*T[k] for k,v in x.items());scale=(.14 if trial<12 else .65)/np.linalg.norm(X,2);x={k:v*scale for k,v in x.items()};rho=np.eye(6)/6+X*scale
 expected=13/18-np.trace(rho@rho).real;actual=rhs(rho,x);err=abs(actual-expected);worst=max(worst,err);assert err<2e-9,(trial,actual,expected)
 print('complex test',trial,'error',err,flush=True)
# The stated rank-two optimizer must saturate the identity.
rho=np.zeros((6,6));rho[np.ix_([0,5],[0,5])]=1/12;rho[np.ix_([1,4],[1,4])]=5/12;x={k:np.trace(v.conj().T@rho) for k,v in T.items() if k[0] in [3,4,5]};assert abs(rhs(rho,x))<2e-10
print('Independent complex-state identity checks passed; maximum error',worst,flush=True)
