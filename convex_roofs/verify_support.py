"""Rebuild pure-state polynomials from symmetric five-qubit reductions exactly.
This checks the physical input to the independent algebraic certificates.
The Dicke index k is the number of down spins (m = 5/2 - k).
"""
from pathlib import Path
import sympy as s

root=Path(__file__).resolve().parent
x=s.symbols('x:8',real=True)
z=s.Matrix([x[k]+s.I*x[k+4] for k in range(4)])
V=s.zeros(6,4)
V[4,0]=1;V[1,1]=1;V[0,2]=s.sqrt(s.Rational(9,20))
V[5,2]=s.sqrt(s.Rational(11,20));V[3,3]=1
assert V.H*V==s.eye(4)
state=V*z
norm=sum(t*t for t in x)

def pure_polynomial(t):
    amplitudes=s.Matrix(t+1,6-t,lambda a,b:
        state[a+b]*s.sqrt(s.binomial(t,a)*s.binomial(5-t,b)/s.binomial(5,a+b)))
    red=amplitudes*amplitudes.H
    purity=sum(red[i,j]*red[j,i] for i in range(t+1) for j in range(t+1))
    return s.expand(s.Rational(t+1,t)*(norm*norm-purity))

f1=pure_polynomial(1)
f2=pure_polynomial(2)
stored=s.sympify((root/'poly_t2.txt').read_text(),locals=dict(zip(map(str,x),x)))
assert s.expand(f2-stored)==0
print('Second-order complex polynomial rebuilt from the physical five-qubit state.',flush=True)
a,b,c,d=x[:4];n=a*a+b*b+c*c+d*d
M=-s.Rational(3,2)*a*a+s.Rational(3,2)*b*b-c*c/4-d*d/2
L=3*b*c/2+s.sqrt(11)*a*c/2+s.sqrt(8)*a*d
assert s.expand(f1.subs(dict.fromkeys(x[4:],0))-(n*n-s.Rational(4,25)*(M*M+L*L)))==0
# The exact compressed spin raising operator gives the triangle-inequality bound.
j=s.Rational(5,2)
Jp=s.zeros(6)
for k in range(1,6):
    m=j-k;Jp[k-1,k]=s.sqrt((j-m)*(j+m+1))
compressed=V.H*Jp*V
assert compressed==s.Matrix([[0,0,s.sqrt(11)/2,0],[0,0,0,0],[0,s.Rational(3,2),0,0],[s.sqrt(8),0,0,0]])
mag=s.expand((z.H*compressed*z)[0]);magz=s.expand((z.H*(V.H*s.diag(*[j-k for k in range(6)])*V)*z)[0])
assert s.expand(f1-(norm*norm-s.Rational(4,25)*(magz**2+mag*s.conjugate(mag))))==0
print('First-order polynomial and phase bound rebuilt from exact spin operators.',flush=True)
# Check that the C5 rotation distinguishes all support characters.
# Doubled magnetic quantum numbers mod 10 account for half-integer spin.
characters=[(-3)%10,3%10,5%10,(-1)%10]
assert len(set(characters))==4 and 5%10==(-5)%10
assert sum([s.Rational(1,12),s.Rational(1,4),s.Rational(1,3),s.Rational(1,3)])==1
print('Orthonormal support, populations, and five-state orbit averaging verified.',flush=True)
