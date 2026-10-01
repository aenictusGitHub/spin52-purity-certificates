"""Exact global certificate for the analytic q1 value (no numerical PSD tests).
All polynomial identities use an algebraic number field; positivity is checked
with rational interval arithmetic and integer square roots.
"""
from pathlib import Path
from fractions import Fraction as F
from math import isqrt
import json
import sympy as s
from sympy.polys.rings import ring
root=Path(__file__).resolve().parent
S=s.sqrt(6071+448*s.sqrt(22))
K=s.QQ.algebraic_field(s.sqrt(2),s.sqrt(11),S)
cv=K.from_sympy
R,a,b,c,d=ring('a,b,c,d',K)
x=[a,b,c,d];n=sum(t*t for t in x)
y=cv((-286-8*s.sqrt(22)+4*S)/(6*(s.sqrt(11)+4*s.sqrt(2))))
q=cv((694-4*s.sqrt(22)-2*S)/675)
meta=json.loads((root/'q1_field.json').read_text())
def coeff(v):return K([s.QQ(t) for t in v])
ell=[coeff(v) for v in meta['ell']]
Gs=[[[coeff(t) for t in row] for row in json.loads((root/f'q1_matrix{k}.json').read_text())] for k in [0,1]]
# Compact analytic Gram matrices; their identity is verified independently below.
la,lb,lc,ld=ell
p=K.one/2;r=K.one/4
b22=cv(s.Rational(28,25))-lc;b33=cv(s.Rational(8,5))-ld
a11=cv(s.Rational(52,25))-la
a22=cv(s.Rational(99,100))-lc;a44=cv(s.Rational(24,25))-ld
a12=-cv(3*s.sqrt(11)/25)+r+y*b22
a13=-cv(12*s.sqrt(2)/25)-r+y*b33
a24=(la+9*ld-cv(s.Rational(383,50)))/8
a33=cv(s.Rational(49,25))-lc-ld-2*a24
Z=K.zero
analyticGs=[[[a11,a12,a13,Z],[a12,a22,Z,a24],[a13,Z,a33,Z],[Z,a24,Z,a44]],[[p,r,Z],[r,b22,Z],[Z,Z,b33]]]
assert analyticGs==Gs
Gs=analyticGs
M=-cv(s.Rational(3,2))*a*a+cv(s.Rational(3,2))*b*b-c*c*cv(s.Rational(1,4))-d*d*cv(s.Rational(1,2))
L=cv(s.Rational(3,2))*b*c+cv(s.sqrt(11)/2)*a*c+cv(s.sqrt(8))*a*d
f=n*n-cv(s.Rational(4,25))*(M*M+L*L)
hs=[[a*b-y*a*a,c*c-4*a*a,c*d-4*a*a,d*d-4*a*a], [a*d-a*c,b*c-y*a*c,b*d-y*a*c]]
represented=sum(G[i][j]*h[i]*h[j] for G,h in zip(Gs,hs) for i in range(len(h)) for j in range(len(h)))
assert represented==f-n*sum(e*t*t for e,t in zip(ell,x))
assert sum(e*cv(w) for e,w in zip(ell,[s.Rational(1,12),s.Rational(1,4),s.Rational(1,3),s.Rational(1,3)]))==q
assert ell[1]==cv(s.Rational(16,25))
# The ensemble objective is exactly the same algebraic number.
Q=(93*y*y-cv(6*(s.sqrt(11)+4*s.sqrt(2)))*y+cv(551-8*s.sqrt(22)))/(75*(y*y+9))
assert Q==q
# Rational intervals.
def I(v):return (F(v),F(v))
def add(x,y):return (x[0]+y[0],x[1]+y[1])
def neg(x):return (-x[1],-x[0])
def sub(x,y):return add(x,neg(y))
def mul(x,y):
 t=[a*b for a in x for b in y];return min(t),max(t)
def div(x,y):
 assert y[0]*y[1]>0
 return mul(x,(1/y[1],1/y[0]))
def sqrt(x):
 assert x[0]>=0
 den=10**100
 def lower(t):return F(isqrt(t.numerator*den*den//t.denominator),den)
 return lower(x[0]),lower(x[1])+F(1,den)
primitive=add(add(sqrt(I(2)),sqrt(I(11))),sqrt(add(I(6071),mul(I(448),sqrt(I(22))))))
def enclosure(t):
 val=I(0)
 for co in t.to_list():val=add(mul(val,primitive),I(str(co)))
 return val
assert enclosure(y)[0]>0 and enclosure(y*y)[1]<3
report={'q1_interval':[str(v) for v in enclosure(q)],'ldl_pivot_intervals':[]}
for G in Gs:
 B=[[enclosure(t) for t in row] for row in G];pivots=[]
 for j in range(len(B)):
  piv=B[j][j];assert piv[0]>0
  pivots.append([str(v) for v in piv])
  for i in range(j+1,len(B)):
   for k in range(j+1,len(B)):B[i][k]=sub(B[i][k],div(mul(B[i][j],B[j][k]),piv))
 report['ldl_pivot_intervals'].append(pivots)
 print('Certified positive LDL pivots:',[float(F(t[0])) for t in pivots])
(root/'q1_verification.json').write_text(json.dumps(report,indent=2))
print('Exact polynomial identity, ensemble equality, and global positivity verified.')
print('q1 = (694 - 4 sqrt(22) - 2 sqrt(6071 + 448 sqrt(22)))/675')
print('q1 approximately',float(enclosure(q)[0]))
