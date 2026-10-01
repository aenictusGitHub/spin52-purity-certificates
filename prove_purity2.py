"""Exact global certificate: Tr(rho^2) <= 13/18 for 2-AC spin-5/2 states.

Default: independently regenerate all rational localizing matrices from
Clebsch--Gordan coefficients and the supplied exact invariant bases, then
check the polynomial identity and rational positive-definiteness proofs.
--quick skips regeneration of the cached localizing coefficient matrices.
No optimization solver or floating-point positivity test is used.
"""
from pathlib import Path
from fractions import Fraction as F
from collections import Counter,defaultdict
from itertools import combinations_with_replacement,product
from functools import lru_cache
from math import factorial,prod
import ast,json,sys,numpy as np,sympy as s
from sympy.physics.wigner import clebsch_gordan
root=Path(__file__).resolve().parent/'global_purity2';j=s.Rational(5,2)
if hasattr(sys,'set_int_max_str_digits'):sys.set_int_max_str_digits(0)
def frac(z):
 z=s.sympify(z);assert z.is_Rational
 return F(int(z.p),int(z.q))
@lru_cache(None)
def cg(a,b,c,m,n):return clebsch_gordan(a,b,c,m,n,m+n)
@lru_cache(None)
def tensor(L,a,b):return s.simplify(s.sqrt(s.Rational(2*L+1,6))*cg(j,L,j,j-b,b-a))
cL={L:s.simplify(tensor(L,0,0)*s.sqrt(s.binomial(2*L,L))) for L in [3,4,5]}
assert cL=={3:s.Rational(5,3),4:s.sqrt(10)/2,5:s.S.One}
D=json.loads((root/'invariants_degree2to4.json').read_text());D.update({'(5,'+k[1:]:v for k,v in json.loads((root/'invariants_degree5.json').read_text()).items()})
assert {ast.literal_eval(k) for k in D}=={(d,n3,n4,d-n3-n4) for d in range(2,6) for n3 in range(d+1) for n4 in range(d-n3+1) if (d-n4)%2==0}
mom={():[(0,F(1))]};offset=1
for typ,d in D.items():
 tuples=[tuple(map(tuple,t)) for t in d['tuples']];deg=len(tuples[0]);counts=tuple(sum(L==l for L,m in tuples[0]) for l in [3,4,5]);ts={0:[],1:[]}
 for parts in product(*[list(combinations_with_replacement(range(-L,L+1),n)) for L,n in zip([3,4,5],counts)]):
  M=sum(map(sum,parts))
  if M in ts:ts[M].append(tuple((L,m) for L,ms in zip([3,4,5],parts) for m in ms))
 assert tuples==ts[0];V=d['basis'];assert len(V)==len(ts[0])-len(ts[1])
 for v in V:
  raised=defaultdict(int)
  for t,z in zip(tuples,v):
   if not z:continue
   for (L,m),n in Counter(t).items():
    if m==L:continue
    q=list(t);q.remove((L,m));q.append((L,m+1));raised[tuple(sorted(q))]+=n*(L+m+1)*z
  assert all(z==0 for z in raised.values())
 for k,v in enumerate(V):assert any(z and all(w[i]==0 for l,w in enumerate(V) if l!=k) for i,z in enumerate(v))
 for i,t in enumerate(tuples):
  mult=factorial(deg)//prod(factorial(n) for n in Counter(t).values());bins=prod(int(s.binomial(2*L,L+m)) for L,m in t);assert mult==d['multiplicity'][i] and bins==d['binomial_product'][i]
  mom[t]=[(offset+k,F(v[i],mult)) for k,v in enumerate(V) if v[i]]
 offset+=len(V)
assert offset==99
print('Exact invariant bases: complete, independent, and annihilated by raising.',flush=True)
@lru_cache(None)
def feature(lab,M):
 L1,L2,K=lab
 if L2==0:return {():s.S.One} if M==0 else {}
 if L1==0:return {((L2,M),):s.sqrt(s.binomial(2*L2,L2+M))/cL[L2]} if abs(M)<=L2 else {}
 out=defaultdict(lambda:s.S.Zero)
 for m in range(-L1,L1+1):
  n=M-m
  if abs(n)<=L2:out[tuple(sorted(((L1,m),(L2,n))))]+=cg(L1,L2,K,m,n)*s.sqrt(s.binomial(2*L1,L1+m)*s.binomial(2*L2,L2+n))/(cL[L1]*cL[L2])
 return {k:s.simplify(v) for k,v in out.items() if v}
features=json.loads((root/'features.json').read_text());polys={}
for bid,d in features.items():
 twoJ,_,fs=d['meta'];J=s.Rational(twoJ,2);original=[];K=np.array([[F(z) for z in row] for row in d['K']],dtype=object)
 for lab,sigma in zip(fs,d['sigma']):
  sigma=s.sympify(sigma);rank=lab[2];vv=[]
  for a in range(6):
   mu=j-a;m=int(J-mu);vv.append({key:frac(s.simplify(sigma*cg(j,rank,J,mu,m)*v*s.sqrt(s.binomial(5,a)))) for key,v in feature(tuple(lab),m).items()} if abs(m)<=rank else {})
  original.append(vv)
 P=[]
 for u in range(K.shape[1]):
  vv=[]
  for a in range(6):
   pp=defaultdict(F)
   for i in range(K.shape[0]):
    for key,c in original[i][a].items():pp[key]+=K[i,u]*c
   vv.append({key:c for key,c in pp.items() if c})
  P.append(vv)
 stored=[[{tuple(map(tuple,key)):F(z) for key,z in row} for row in pp] for pp in d['P']];assert P==stored;polys[bid]=P
print('All coupled feature polynomials verified from exact Clebsch--Gordan coefficients.',flush=True)
RR=[[[] for b in range(6)] for a in range(6)]
for a in range(6):
 for b in range(6):
  for L in [3,4,5]:
   if abs(b-a)>L:continue
   val=frac(s.simplify(tensor(L,a,b)*s.sqrt(s.binomial(2*L,L+b-a))/(cL[L]*s.sqrt(s.binomial(5,a)*s.binomial(5,b)))))
   if val:RR[a][b].append((((L,b-a),),val))
Hstored=json.loads((root/'moment_matrices.json').read_text());cert=json.loads((root/'certificate.json').read_text());sumcoeff=np.full(99,F(0),dtype=object)
for bid,P in polys.items():
 H=np.array([[[F(z) for z in row] for row in h] for h in Hstored[bid]],dtype=object);n=len(P)
 if '--quick' not in sys.argv:
  for u in range(n):
   for v in range(u,n):
    poly=defaultdict(F)
    for a in range(6):
     for b in range(6):
      for mi,ci in P[u][a].items():
       for mj,cj in P[v][b].items():
        ms=mi+tuple((L,-m) for L,m in mj);value=ci*cj*(-1)**abs(sum(m for L,m in mj))
        if a==b and sum(L for L,m in ms)%2==0:poly[tuple(sorted(ms))]+=value/F(6*int(s.binomial(5,a)))
        for key,fac in RR[a][b]:
         if sum(L for L,m in ms+key)%2==0:poly[tuple(sorted(ms+key))]+=value*fac
    co=[F(0)]*99
    for ms,c in poly.items():
     if c:
      for k,q in mom.get(ms,[]):co[k]+=c*q
    assert all(H[k,u,v]==H[k,v,u]==c for k,c in enumerate(co))
 d=cert[bid];W=np.array([[F(z) for z in row] for row in d['W']],dtype=object);Q=np.array([[F(z) for z in row] for row in d['Q']],dtype=object);R=np.array([[F(z) for z in row] for row in d['R']],dtype=object);assert np.all(Q==Q.T)
 G=R@Q@R.T;assert all(G[i,i]>sum(abs(G[i,k]) for k in range(G.shape[0]) if k!=i) for i in range(G.shape[0]));Y=W@Q@W.T
 sumcoeff+=np.array([sum((Y*hk).flat,F(0)) for hk in H],dtype=object)
 print('Block',bid,': exact coefficients and rational positive definiteness passed.',flush=True)
S=[F(0)]*99
for L in [3,4,5]:
 for m in range(-L,L+1):
  fac=frac(s.simplify((-1)**abs(m)*s.binomial(2*L,L+m)/cL[L]**2))
  for k,q in mom[tuple(sorted(((L,m),(L,-m))))]:S[k]+=fac*q
stored=[F(z) for z in json.loads((root/'objective.json').read_text())['objective']];assert S==stored;target=-np.array(S,dtype=object);target[0]=F(5,9);assert np.all(sumcoeff==target)
# Exact feasible optimizer in the descending-m basis.
rho=s.zeros(6)
for a,b,w in [(0,5,s.Rational(1,12)),(1,4,s.Rational(5,12))]:
 for i in [a,b]:
  for k in [a,b]:rho[i,k]=w
assert s.trace(rho)==1 and s.trace(rho*rho)==s.Rational(13,18)
assert rho.eigenvals()=={s.S.Zero:4,s.Rational(1,6):1,s.Rational(5,6):1}
for L in [1,2]:
 for m in range(-L,L+1):assert s.simplify(sum(tensor(L,a,a+m)*rho[a,a+m] for a in range(6) if 0<=a+m<6))==0
print('EXACT IDENTITY: 13/18 - Tr(rho^2) = sum of nine nonnegative terms.')
print('The rank-two state attains equality. Hence P_max(5,2) = 13/18 globally.')
