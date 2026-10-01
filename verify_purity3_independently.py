"""Direct numerical three-copy check; supplementary to the exact proof."""
import numpy as np
import sympy as s
from sympy.physics.wigner import clebsch_gordan
from angular_momentum import tensors
rng=np.random.default_rng(47381)
# Rebuild the three witness projectors from angular momentum coupling and
# check the purity identity directly on random 3-AC density matrices.
j=s.Rational(5,2)
cases=[(s.Rational(3,2),[0,np.sqrt(5/7),0,-np.sqrt(2/7)],False,84/25),
 (s.Rational(5,2),[np.sqrt(2/9),0,-np.sqrt(5/18),0,-np.sqrt(.5),0],False,624/25),
 (s.Rational(3,2),[np.sqrt(24/65),np.sqrt(7/13),np.sqrt(6/65),0],True,1092/25)]
projectors=[]
for J,coeff,flip,w in cases:
 ks=[K for K in range(6) if abs(K-j)<=J<=K+j];P=np.zeros((216,216))
 for twoM in range(-int(2*J),int(2*J)+1,2):
  M=s.Rational(twoM,2);v=np.zeros((6,6,6))
  for a in range(6):
   for b in range(6):
    c=3*j-a-b-M
    if c.is_integer and 0<=c<6:
     c=int(c);m1=j-a;m2=j-b;m3=j-c
     v[a,b,c]=sum(q*float(clebsch_gordan(j,j,K,m1,m2,m1+m2)*clebsch_gordan(K,j,J,m1+m2,m3,M)) for K,q in zip(ks,coeff))
  P+=np.outer(v.ravel(),v.ravel())/float(2*J+1)
 projectors.append((P,flip,w))
Ts=tensors(5);H=[]
for L in [4,5]:
 H.append(Ts[L,0])
 for M in range(1,L+1):
  T=Ts[L,M];H.extend([(T+T.conj().T)/np.sqrt(2),1j*(T-T.conj().T)/np.sqrt(2)])
H=np.array(H);U=np.zeros((6,6))
for k in range(6):U[5-k,k]=(-1)**k
worst=0
for _ in range(12):
 X=np.einsum('a,aij->ij',rng.normal(size=len(H)),H);X*=.14/np.linalg.norm(X,2);rho=np.eye(6)/6+X
 flip=U@rho.conj()@U.T
 rr=np.kron(np.kron(rho,rho),rho);tr=np.kron(np.kron(flip,rho),rho)
 val=sum(w*np.trace(P@(tr if fl else rr)).real for P,fl,w in projectors)
 err=abs(val-(.5-np.trace(rho@rho).real));worst=max(worst,err)
assert worst<2e-14
print('Independent three-copy identity tests: PASS; max residual',worst,flush=True)
