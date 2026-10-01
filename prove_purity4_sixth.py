"""Exact coefficient verification of S^2(1/8-S)>=0 for rank-five X.
No optimization, numerical eigensolver, invariant projection, or stored
moment matrix is used: expand all seven nonnegative terms directly.
"""
from pathlib import Path
from functools import lru_cache
from collections import defaultdict
from fractions import Fraction as Q
import sympy as s
from sympy.physics.wigner import clebsch_gordan
root=Path(__file__).resolve().parent
j=s.Rational(5,2);B5=[int(s.binomial(5,a)) for a in range(6)]
@lru_cache(None)
def cg(a,b,c,m,n):return clebsch_gordan(a,b,c,m,n,m+n)
@lru_cache(None)
def feature(deg,K,M,pair=-1):
 out=defaultdict(lambda:s.S.Zero)
 for a in range(-5,6):
  for b in range(-5,6):
   if deg==2:
    if a+b!=M:continue
    key=tuple(sorted((a,b)));val=cg(5,5,K,a,b)*s.sqrt(s.binomial(10,5+a)*s.binomial(10,5+b))
   else:
    c=M-a-b
    if abs(c)>5 or abs(a+b)>pair:continue
    key=tuple(sorted((a,b,c)));val=cg(5,5,pair,a,b)*cg(pair,5,K,a+b,c)*s.sqrt(s.binomial(10,5+a)*s.binomial(10,5+b)*s.binomial(10,5+c))
   if val:out[key]+=val
 return {k:s.simplify(v) for k,v in out.items() if v}
def rationalize(rows):
 lead=next((v for r in rows for v in r.values() if v),None)
 if lead is None:return Q(0),[{} for r in rows]
 gamma=s.simplify(lead**2);assert gamma.is_Rational
 pp=[]
 for r in rows:
  rr={}
  for k,v in r.items():
   val=s.simplify(v/lead);assert val.is_Rational
   if val:rr[k]=Q(int(val.p),int(val.q))
  pp.append(rr)
 return Q(int(gamma.p),int(gamma.q)),pp
def conj(p):return {tuple(sorted(-m for m in key)):v*(-1)**abs(sum(key)) for key,v in p.items()}
def mul(a,b):
 out=defaultdict(Q)
 for ka,va in a.items():
  for kb,vb in b.items():out[tuple(sorted(ka+kb))]+=va*vb
 return {k:v for k,v in out.items() if v}
def add(out,p,fac=Q(1)):
 for k,v in p.items():out[k]+=fac*v
# X_ab = sqrt(binomial(5,a) binomial(5,b)) (-1)^b y_(b-a).
R=[[{(b-a,):Q((-1)**b)} for b in range(6)] for a in range(6)]
X2=[[defaultdict(Q) for b in range(6)] for a in range(6)]
for a in range(6):
 for b in range(6):
  for c in range(6):add(X2[a][b],mul(R[a][c],R[c][b]),Q(B5[c]))
terms=[
 ('local',5,[(0,s.Rational(45,2)),(2,9*s.sqrt(1365)),(4,2*s.sqrt(273))],Q(1094426,4929645591)),
 ('local',7,[(2,16*s.sqrt(33)/11),(4,s.S.One)],Q(131329055,252802338)),
 ('local',9,[(2,-2*s.sqrt(33)/11),(4,s.S.One)],Q(3843346489,252802338)),
 ('local',11,[(4,43*s.sqrt(149226)/244188),(6,s.sqrt(16302)/114),(8,s.S.One)],Q(607184206272,6025122389)),
 ('scalar',3,[(2,18*s.sqrt(26)/13),(4,s.S.One)],Q(2002382525,14156930928)),
 ('scalar',4,[(2,s.S.One)],Q(12329190117,337069784)),
 ('scalar',6,[(2,-6*s.sqrt(3)),(4,s.S.One)],Q(162694181,758407014)),
]
total=defaultdict(Q)
for typ,spin,components,weight in terms:
 assert weight>0;poly=defaultdict(Q)
 if typ=='local':
  J=s.Rational(spin,2)
  for twoM in range(-spin,spin+1,2):
   M=s.Rational(twoM,2);rows=[]
   for a in range(6):
    mu=j-a;m=M-mu;r=defaultdict(lambda:s.S.Zero)
    for K,w in components:
     if abs(m)>K:continue
     for key,c in feature(2,K,int(m)).items():r[key]+=w*cg(j,K,J,mu,m)*c*s.sqrt(B5[a])
    rows.append({key:s.simplify(c) for key,c in r.items() if c})
   gamma,P=rationalize(rows)
   for a in range(6):
    for b in range(6):
     pp=mul(P[a],conj(P[b]));fac=gamma/Q(spin+1)
     if a==b:add(poly,pp,fac/Q(36*B5[a]))
     add(poly,mul(pp,X2[a][b]),-fac)
 else:
  K=spin
  for M in range(-K,K+1):
   r=defaultdict(lambda:s.S.Zero)
   for pair,w in components:
    for key,c in feature(3,K,M,pair).items():r[key]+=w*c
   gamma,rows=rationalize([{key:s.simplify(c) for key,c in r.items() if c}]);add(poly,mul(rows[0],conj(rows[0])),gamma/Q(2*K+1))
 add(total,poly,weight);print('Verified positive term',typ,spin,flush=True)
S=defaultdict(Q)
for m in range(-5,6):S[tuple(sorted((m,-m)))]+=Q((-1)**abs(m)*int(s.binomial(10,5+m)))
target=defaultdict(Q);S2=mul(S,S);add(target,S2,Q(1,8));add(target,mul(S2,S),Q(-1))
residual={k:total.get(k,Q(0))-target.get(k,Q(0)) for k in set(total)|set(target)};residual={k:v for k,v in residual.items() if v}
assert not residual,residual
print('EXACT FULL POLYNOMIAL IDENTITY VERIFIED: S^2(1/8-S) = seven nonnegative terms.')
print('Therefore every 4-AC spin-5/2 density matrix has purity <= 7/24.')

# Exact attainment by the Table-I rank-four state.
rho=s.Matrix([[9,0,0,0,0,s.sqrt(99)],[0,15,0,0,0,0],[0,0,0,0,0,0],[0,0,0,20,0,0],[0,0,0,0,5,0],[s.sqrt(99),0,0,0,0,11]])/60
assert rho.eigenvals()=={s.S.Zero:2,s.Rational(1,3):2,s.Rational(1,4):1,s.Rational(1,12):1}
assert s.trace(rho*rho)==s.Rational(7,24)
for L in range(1,5):
 for m in range(-L,L+1):
  value=sum(s.sqrt(s.Rational(2*L+1,6))*cg(j,L,j,j-(a+m),m)*rho[a,a+m] for a in range(6) if 0<=a+m<6)
  assert s.simplify(value)==0
print('The stated rank-four optimizer is positive, 4-AC, and has purity exactly 7/24.')
