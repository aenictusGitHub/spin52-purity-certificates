"""Exact three-copy positivity certificate: Pmax(5,3)=1/2."""
from functools import lru_cache
from itertools import permutations
import sympy as s
from sympy.physics.wigner import clebsch_gordan,wigner_3j
j=s.Rational(5,2);sq=s.sqrt
@lru_cache(None)
def T(L,m,a):
 # Matrix element from Dicke index a to a-m.
 if not 0<=a-m<6:return s.S(0)
 return sq(s.Rational(2*L+1,6))*clebsch_gordan(j,L,j,j-a,m,j-a+m)
def ket(J,vs):
 ks=[K for K in range(6) if abs(K-j)<=J<=K+j]
 out={}
 for a in range(6):
  for b in range(6):
   c=3*j-a-b-J
   if c.is_integer and 0<=c<6:
    c=int(c);m1=j-a;m2=j-b;m3=j-c
    v=sum(vk*clebsch_gordan(j,j,K,m1,m2,m1+m2)*clebsch_gordan(K,j,J,m1+m2,m3,J) for K,vk in zip(ks,vs))
    v=s.simplify(v)
    if v:out[a,b,c]=v
 assert s.simplify(sum(v*v for v in out.values()))==1
 return out
def expectation(v,Ls,ms):
 terms=[]
 for src,x in v.items():
  dst=tuple(a-m for a,m in zip(src,ms))
  if dst not in v:continue
  terms.append(x*v[dst]*s.prod(T(L,m,a) for L,m,a in zip(Ls,ms,src)))
 return s.simplify(sum(terms))
def row(v,flip):
 rr=[s.Rational(1,216)]
 for L in [4,5]:
  val=0
  for ij in [(0,1),(0,2),(1,2)]:
   for m in range(-L,L+1):
    LL=[0,0,0];ms=[0,0,0]
    LL[ij[0]]=LL[ij[1]]=L;ms[ij[0]]=m;ms[ij[1]]=-m
    # T00=I/sqrt6, so insert sqrt6 to recover the spectator identity.
    val+=s.S(-1)**m*sq(6)*expectation(v,tuple(LL),tuple(ms))*((-1)**LL[0] if flip else 1)
  rr.append(s.simplify(val/(6*(2*L+1))))
 for triple in [(4,4,4),(4,5,5)]:
  val=0
  for LL in set(permutations(triple)):
   for m1 in range(-LL[0],LL[0]+1):
    for m2 in range(-LL[1],LL[1]+1):
     m3=-m1-m2
     if abs(m3)>LL[2]:continue
     c=wigner_3j(*LL,m1,m2,m3)*((-1)**LL[0] if flip else 1)
     if c:val+=c*expectation(v,LL,(m1,m2,m3))
  rr.append(s.simplify(val))
 return s.Matrix([rr])
# Coupled basis |((j,j)K,j)J,M=J>, with increasing allowed K.
cases=[(s.Rational(3,2),[0,sq(s.Rational(5,7)),0,-sq(s.Rational(2,7))],False),
       (s.Rational(5,2),[sq(s.Rational(2,9)),0,-sq(s.Rational(5,18)),0,-sq(s.Rational(1,2)),0],False),
       (s.Rational(3,2),[sq(s.Rational(24,65)),sq(s.Rational(7,13)),sq(s.Rational(6,65)),0],True)]
rows=[]
for J,vs,flip in cases:
 r=row(ket(J,vs),flip);rows.append(r);print('row',r,flush=True)
result=s.simplify(sum((w*r for w,r in zip([s.Rational(84,25),s.Rational(624,25),s.Rational(1092,25)],rows)),s.zeros(1,5)))
assert result==s.Matrix([[s.Rational(1,3),-1,-1,0,0]])
print('EXACT IDENTITY VERIFIED:',result,flush=True)
