"""Verify the four-variable discriminant specification of the exact q2 roof.
All identities and interval inequalities use rational arithmetic. Polynomial
translation prevents cancellation in interval evaluation of the discriminant.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb
import sympy as s,mpmath as mp,json
p=Path(__file__).resolve().parent
raw=json.loads((p/'q2_reduced_polynomial.json').read_text())
A,C,u,v=s.symbols('A C u v')
P=(7920-4400*u)*A*A+81312-96800*v
Q=(243-400*u)*A*A+12870-8800*v
L=1584*A+8712;M=1584*A*A
E=(288-400*u)*A**4+(4752-8800*u-8800*v)*A*A+156816-193600*v
K=s.expand((P+825*C*C)*(Q*C*C-2*M*C+E)-L*L*C*C)
coeffs={tuple(ex):F(co) for ex,co in raw['Phi']}
Phi=sum(s.Rational(co)*A**ex[0]*u**ex[1]*v**ex[2] for ex,co in coeffs.items())
assert s.expand(s.discriminant(K,C)-s.Rational(raw['discriminant_scale'])*P*Phi)==0
assert s.degree(Phi,A)==22
print('Exact quartic-discriminant identity verified; degree in A is 22.',flush=True)
# Completing the square gives the exact elimination of the B amplitude.
B=s.symbols('B')
H=s.sympify(json.loads((p/'q2_algebraic_system.json').read_text())['H'],locals={str(t):t for t in [A,B,C,u,v]})
assert s.expand((P+825*C*C)*H-((P+825*C*C)*B-L*C)**2-K)==0
print('Exact amplitude elimination verified.',flush=True)

def add(x,y):return x[0]+y[0],x[1]+y[1]
def sub(x,y):return x[0]-y[1],x[1]-y[0]
def mul(x,y):
 z=[a*b for a in x for b in y];return min(z),max(z)
def I(x):return F(x),F(x)
def amp(x):return max(abs(t) for t in x)
def translate(pol,center):
 out=pol
 for j,x in enumerate(center):
  new={}
  for ex,co in out.items():
   for k in range(ex[j]+1):
    ee=list(ex);ee[j]=k;ee=tuple(ee)
    new[ee]=new.get(ee,F(0))+co*comb(ex[j],k)*x**(ex[j]-k)
  out={ex:co for ex,co in new.items() if co}
 return out

def derivative(pol,j):
 out={}
 for ex,co in pol.items():
  if ex[j]:
   ee=list(ex);ee[j]-=1;out[tuple(ee)]=co*ex[j]
 return out

def eval_box(pol,rad):
 # Every coordinate is in [-rad,rad]. An even monomial is nonnegative.
 zero=(0,)*3;out=I(pol.get(zero,F(0)))
 for ex,co in pol.items():
  if sum(ex)==0:continue
  top=rad**sum(ex);term=(F(0),top) if all(k%2==0 for k in ex) else (-top,top)
  out=add(out,mul(I(co),term))
 return out
old=json.loads((p/'q2_algebraic_verification.json').read_text())
coarse=old['coarse'];indices=[0,3,6,7]
center=[F(coarse['center'][j]) for j in indices];r=F(coarse['radius'])
blocks=[translate(coeffs,[center[i],center[2],center[3]]) for i in [0,1]]
polys=[];Jbox=[];fc=[];Jc=[]
for k,poly in enumerate(blocks):
 for ff in [poly,derivative(poly,0)]:
  fc.append(ff.get((0,0,0),F(0)))
  row=[I(0)]*4;cr=[F(0)]*4
  for i,j in enumerate([k,2,3]):
   der=derivative(ff,i);row[j]=eval_box(der,r);cr[j]=der.get((0,0,0),F(0))
  Jbox.append(row);Jc.append(cr)
mp.mp.dps=160
Jm=mp.matrix([[mp.mpf(t.numerator)/t.denominator for t in row] for row in Jc]);Cm=Jm**-1
inv=[[F(int(mp.nint(Cm[i,j]*10**110)),10**110) for j in range(4)] for i in range(4)]
beta=F(0);disp=F(0)
for i in range(4):
 rowsum=F(0)
 for j in range(4):
  z=I(0)
  for k in range(4):z=add(z,mul(I(inv[i][k]),Jbox[k][j]))
  rowsum+=amp(sub(I(int(i==j)),z))
 beta=max(beta,rowsum)
 disp=max(disp,abs(sum(inv[i][k]*fc[k] for k in range(4))))
assert beta<1 and disp+r*beta<r
# The fine box for the previously certified roof projects strictly inside X.
fcenter=[F(old['fine']['center'][j]) for j in indices];fr=F(old['fine']['radius'])
assert all(c-r<z-fr and z+fr<c+r for c,z in zip(center,fcenter))
# A generic quadratic common root of the quartic and its derivative has
# nonzero second derivative; this licenses the discriminant-derivative step.
# The original root is an interior stationary zero, so K=K_C=K_A=0.
subres=s.subresultants(K,s.diff(K,C),C)
S1=next(t for t in subres if s.degree(t,C)==1)
# Rational interval enclosure for low-degree regularity polynomials.
def eval_general(expr,cc,rad):
 pp=s.Poly(expr,A,C,u,v);out=I(0)
 for ex,co in pp.terms():
  t=I(str(co))
  for x,k in zip(cc,ex):
   for _ in range(k):t=mul(t,(x-rad,x+rad))
  out=add(out,t)
 return out
regs=[]
# Use the original fine certificate, so natural intervals remain narrow.
for ai,ci in [(0,2),(3,5)]:
 cc=[F(old['fine']['center'][i]) for i in [ai,ci,6,7]]
 regular=[eval_general(t,cc,fr) for t in [P,s.diff(K,C,2),s.Poly(S1,C).coeff_monomial(C),Q]]
 assert regular[0][0]>0 and regular[1][0]>0 and regular[2][0]*regular[2][1]>0 and regular[3][0]>0
 regs.append([[str(x) for x in t] for t in regular])
report={'variables':['A_minus','A_plus','u','v'],'center':[str(x) for x in center],'radius':str(r),'contraction_bound':str(beta),'selfmap_bound':str(disp+r*beta),'rational_inverse':[[str(x) for x in row] for row in inv],'regularity_intervals':regs}
(p/'q2_reduced_verification.json').write_text(json.dumps(report,indent=2))
print('Unique reduced root certified: contraction <=',float(beta),'selfmap bound <=',float(disp+r*beta),flush=True)
print('Its box contains the projection of the certified original root.',flush=True)
print('Quartic-root regularity verified: same globally proved q2, with four unknowns.',flush=True)
