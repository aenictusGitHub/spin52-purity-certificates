"""Certify the algebraic q2 roof: root existence/uniqueness, SOS, and ensemble.
Verification uses exact rational polynomial identities and rational intervals.
Decimal approximations select the root, but no numerical residual or numerical
positive-semidefinite test is accepted as a certificate.
"""
from fractions import Fraction as F
from pathlib import Path
from math import isqrt
import json,sympy as s,mpmath as mp
p=Path(__file__).resolve().parent
data=json.loads((p/'q2_algebraic_system.json').read_text())
xx=s.symbols('a b c e f g u v');a,b,c,e,f,g,u,v=xx
A,B,C,D=s.symbols('A B C D');ss=s.symbols('s:4')

def expr(terms):return sum(s.Rational(co)*s.prod(x**k for x,k in zip(xx,ex)) for ex,co in terms)
polys=[expr(poly) for poly in data['contact_polynomials']]
Q=s.Matrix([[expr(poly) for poly in row] for row in data['gram']])
G=s.sympify(data['G'],locals={str(x):x for x in [A,B,C,D,u,v]})
H=G.subs(D,1)
expected=[]
for tr in [(a,b,c),(e,f,g)]:expected += [H.subs(dict(zip([A,B,C],tr)))]+[s.diff(H,x).subs(dict(zip([A,B,C],tr))) for x in [A,B,C]]
assert all(s.expand(x-y)==0 for x,y in zip(polys,expected))
T=s.Matrix([[0,0,a,e],[1,0,b,f],[0,1,c,g],[0,0,1,1]])
trans=G.subs(dict(zip([A,B,C,D],T*s.Matrix(ss))),simultaneous=True)
mon=s.Matrix([ss[i]*ss[j] for i in range(4) for j in range(i+1,4)])
h1,da1,db1,dc1,h2,da2,db2,dc2=polys
remainder=h1*ss[2]**4+db1*ss[0]*ss[2]**3+dc1*ss[1]*ss[2]**3+(4*h1+(e-a)*da1+(f-b)*db1+(g-c)*dc1)*ss[2]**3*ss[3]
remainder+=h2*ss[3]**4+db2*ss[0]*ss[3]**3+dc2*ss[1]*ss[3]**3+(4*h2+(a-e)*da2+(b-f)*db2+(c-g)*dc2)*ss[2]*ss[3]**3
assert s.expand(trans-(mon.T*Q*mon)[0]-remainder)==0
print('Exact SOS identity modulo the eight contact equations verified.',flush=True)
# Independently check the real-amplitude bound and its rational rescaling.
y=s.symbols('x:8',real=True);zz=s.symbols('z:4');ww=s.symbols('w:4')
f2=s.sympify((p/'poly_t2.txt').read_text(),locals=dict(zip(map(str,y),y)))
fc=s.expand(f2.subs({y[k]:(zz[k]+ww[k])/2 for k in range(4)}).subs({y[k+4]:(zz[k]-ww[k])/(2*s.I) for k in range(4)}))
for powers,co in s.Poly(fc,*(zz+ww)).terms():
 if powers[:4]!=powers[4:]:assert co.is_negative
real=f2.subs(dict.fromkeys(y[4:],0));scaled=real.subs(dict(zip(y[:4],[A,s.sqrt(11)*B,C,s.sqrt(22)*D])),simultaneous=True)
n=A*A+11*B*B+C*C+22*D*D
assert s.expand(G-400*(scaled-n*(u*A*A+s.Rational(18,25)*11*B*B+s.Rational(297,400)*C*C+22*v*D*D)))==0
print('Complex phases are bounded by the real-amplitude polynomial exactly.',flush=True)
# Interval arithmetic, all endpoints Fraction.
def I(x):return (F(x),F(x))
def add(x,y):return x[0]+y[0],x[1]+y[1]
def neg(x):return -x[1],-x[0]
def sub(x,y):return add(x,neg(y))
def mul(x,y):
 t=[a*b for a in x for b in y];return min(t),max(t)
def div(x,y):
 assert y[0]*y[1]>0
 return mul(x,(1/y[1],1/y[0]))
def abmax(x):return max(abs(x[0]),abs(x[1]))
def serial(poly):return [(ex,F(str(co))) for ex,co in s.Poly(poly,*xx).terms()]
P=[serial(poly) for poly in polys]
J=[[serial(s.diff(poly,x)) for x in xx] for poly in polys]
QQ=[[serial(Q[i,j]) for j in range(6)] for i in range(6)]
def evalpoly(poly,box):
 out=I(0)
 for ex,co in poly:
  term=I(co)
  for k,degree in enumerate(ex):
   for _ in range(degree):term=mul(term,box[k])
  out=add(out,term)
 return out
# Use the high-precision solve only to specify rational centers and inverse guesses.
mp.mp.dps=160
initial=json.loads((p/'q2_600digits.json').read_text())['solutions']
vs=list(map(mp.mpf,initial));vals=[mp.sqrt(22)*vs[0],mp.sqrt(2)*vs[1],mp.sqrt(22)*vs[2],mp.sqrt(22)*vs[3],mp.sqrt(2)*vs[4],mp.sqrt(22)*vs[5],vs[6],vs[7]]
report={}
for title,digits,radius in [('coarse',10,F(1,10**9)),('fine',65,F(1,10**60))]:
 center=[F(int(mp.nint(t*10**digits)),10**digits) for t in vals]
 box=[(t-radius,t+radius) for t in center];point=[I(t) for t in center]
 if title=='coarse':coarse_box=box
 else:assert all(a<=c and d<=b for (a,b),(c,d) in zip(coarse_box,box))
 Fc=[evalpoly(poly,point)[0] for poly in P]
 Jc=[[evalpoly(poly,point)[0] for poly in row] for row in J]
 Jbox=[[evalpoly(poly,box) for poly in row] for row in J]
 Jnum=mp.matrix([[mp.mpf(t.numerator)/t.denominator for t in row] for row in Jc]);Cnum=Jnum**-1
 inv=[[F(int(mp.nint(Cnum[i,j]*10**90)),10**90) for j in range(8)] for i in range(8)]
 # Infinity-norm contraction test for x -> x - inv*f(x), not a residual-only test.
 contraction=F(0);error=F(0)
 for i in range(8):
  rowsum=F(0)
  for j in range(8):
   z=I(0)
   for k in range(8):z=add(z,mul(I(inv[i][k]),Jbox[k][j]))
   rowsum+=abmax(sub(I(int(i==j)),z))
  contraction=max(contraction,rowsum)
  displacement=abs(sum(inv[i][k]*Fc[k] for k in range(8)))
  error=max(error,displacement+radius*rowsum)
 assert contraction<1 and error<radius
 # C is nonsingular: ||I-C J(center)||<1 proves this as J is square.
 # Banach's theorem now gives exactly one root in the box.
 Gbox=[[evalpoly(poly,box) for poly in row] for row in QQ]
 pivots=[]
 for j in range(6):
  pivot=Gbox[j][j];assert pivot[0]>0,pivot;pivots.append(pivot)
  for i in range(j+1,6):
   for k in range(j+1,6):Gbox[i][k]=sub(Gbox[i][k],div(mul(Gbox[i][j],Gbox[j][k]),pivot))
 assert sub(box[0],box[3])[1]<0 # invertibility of the four-zero change of coordinates
 # Positive exact orbit weights: solve a/d population reconstruction first.
 aa,bb,cc,ee,ff,gg,uu,vv=box
 a2=mul(aa,aa);e2=mul(ee,ee)
 t1=div(sub(I(F(1,12)),mul(I(F(1,66)),e2)),sub(a2,e2))
 t2=sub(I(F(1,66)),t1)
 bweight=sub(I(F(1,4)),mul(I(11),add(mul(t1,mul(bb,bb)),mul(t2,mul(ff,ff)))))
 cweight=sub(I(F(1,3)),add(mul(t1,mul(cc,cc)),mul(t2,mul(gg,gg))))
 n1=add(add(a2,mul(I(11),mul(bb,bb))),add(mul(cc,cc),I(22)))
 n2=add(add(e2,mul(I(11),mul(ff,ff))),add(mul(gg,gg),I(22)))
 weights=[bweight,cweight,mul(t1,n1),mul(t2,n2)]
 assert all(w[0]>0 for w in weights)
 qinterval=add(I(F(171,400)),div(add(uu,mul(I(4),vv)),I(12)))
 report[title]={'center':[str(t) for t in center],'radius':str(radius),'rational_inverse':[[str(t) for t in row] for row in inv],'contraction_bound':str(contraction),'selfmap_error_bound':str(error),'ldl_pivots':[[str(t) for t in z] for z in pivots],'weights':[[str(t) for t in z] for z in weights],'q2_interval':[str(t) for t in qinterval]}
 print(title,'root certified; contraction <=',float(contraction),'selfmap error <=',float(error),'radius',float(radius),flush=True)
 print('LDL lower bounds', [float(t[0]) for t in pivots],flush=True)
 print('Positive ensemble weights', [float(t[0]) for t in weights],flush=True)
(p/'q2_algebraic_verification.json').write_text(json.dumps(report,indent=2))
print('GLOBAL EXACT q2 = 171/400 + (u+4v)/12, at the uniquely isolated algebraic root.',flush=True)
print('q2 =',mp.nstr(mp.mpf(171)/400+(vals[6]+4*vals[7])/12,62),flush=True)
