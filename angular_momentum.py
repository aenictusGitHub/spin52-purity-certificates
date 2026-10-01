"""Hilbert--Schmidt normalized spin tensors for independent numerical checks."""
import numpy as np
from functools import lru_cache
from sympy import Rational
from sympy.physics.wigner import clebsch_gordan

@lru_cache(None)
def tensors(N):
    j=Rational(N,2); out={}
    for L in range(N+1):
        for M in range(-L,L+1):
            T=np.zeros((N+1,N+1),complex)
            for c in range(N+1):
                m=j-c; mp=m+M; row=j-mp
                if 0<=row<=N:
                    T[int(row),c]=np.sqrt((2*L+1)/(N+1))*float(clebsch_gordan(j,L,j,m,M,mp))
            out[L,M]=T
    return out
