"""Exact feasibility and attainment for the three global purity bounds."""
import sympy as s
from math import comb

def project(v):return v*v.conjugate().T

def reduction(rho,t):
    N=rho.rows-1
    return s.Matrix(t+1,t+1,lambda a,c:s.simplify(sum(
        rho[a+b,c+b]*comb(N-t,b)*s.sqrt(s.Rational(comb(t,a)*comb(t,c),comb(N,a+b)*comb(N,c+b)))
        for b in range(N-t+1))))

def main():
    g=s.Matrix([1,0,0,0,0,1])/s.sqrt(2)
    h=s.Matrix([0,1,0,0,1,0])/s.sqrt(2)
    u=s.Matrix([0,s.sqrt(s.Rational(5,6)),0,0,0,1/s.sqrt(6)])
    v=s.Matrix([1/s.sqrt(6),0,0,0,-s.sqrt(s.Rational(5,6)),0])
    r4=s.Matrix([[9,0,0,0,0,s.sqrt(99)],[0,15,0,0,0,0],[0,0,0,0,0,0],
                 [0,0,0,20,0,0],[0,0,0,0,5,0],[s.sqrt(99),0,0,0,0,11]])/60
    for t,rho,P in [(2,project(g)/6+5*project(h)/6,s.Rational(13,18)),
                    (3,(project(u)+project(v))/2,s.Rational(1,2)),
                    (4,r4,s.Rational(7,24))]:
        assert rho==rho.conjugate().T and s.trace(rho)==1
        assert all(e>=0 for e in rho.eigenvals())
        assert reduction(rho,t)==s.eye(t+1)/(t+1)
        assert s.trace(rho*rho)==P
        print(f'EXACT optimizer t={t}: PSD, trace 1, maximally mixed t-reduction, purity {P}.',flush=True)
if __name__=='__main__':main()
