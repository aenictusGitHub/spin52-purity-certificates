Exact second-order global purity certificate for spin 5/2
=========================================================

Statement: every 2-anticoherent spin-5/2 density matrix satisfies
Tr(rho^2) <= 13/18. The rank-two state in Table I attains equality.
This proof has no assumed density-matrix rank, reality, or eigenbasis.

Run from the project directory:
    python3 prove_purity2.py
This regenerates the coefficient matrices, verifies the complete invariant
bases and coupled features, proves positivity by exact rational inequalities,
and checks the attaining state. NumPy and SymPy are required; no numerical
optimization solver is used. --quick skips only coefficient regeneration.

The separate verify_purity2_independently.py performs a numerical cross-check
on complex Hermitian matrices using the physical coupled matrices directly.
It is independent of the invariant moment matrices and is not the exact proof.

Data conventions
----------------
Physical basis: m = 5/2, 3/2, ..., -5/2; row index a=0,...,5.
T_Lm are Hilbert--Schmidt normalized, with standard Clebsch--Gordan phases.
For 2-AC states rho=I/6+sum_{L=3,4,5;m} x_Lm T_Lm, use variables
    x_Lm = sqrt(binomial(2L,L+m)) y_Lm / c_L,
    (c_3,c_4,c_5) = (5/3, sqrt(10)/2, 1).
Hermiticity means y_L,-m = (-1)^m conjugate(y_Lm).

features.json:
  Keys 1,...,9 correspond to J=1/2,3/2,...,17/2.
  meta records 2J and the ordered features (L1,L2,K).
  (0,0,0) is the constant; (0,L,L) is x_L; otherwise the feature is
  (x_L1 tensor x_L2)_K. Only features that couple j=5/2 to J occur.
  sigma is an explicit real square-root scaling vector; K is rational.
  P gives the rational polynomial components at M=J after scaling the
  physical component a by sqrt(binomial(5,a)). The verifier derives them
  independently from the coupling convention.

certificate.json:
  W, Q, and R are rational matrices. In the notation of the manuscript,
      C_J = diag(sigma_J) K_J W_J.
  The exact certificate is
      13/18 - Tr(rho^2) = sum_J Tr(Q_J C_J^T G_J(rho) C_J),
  with the positive rotational/time-reversal averages G_J defined in the
  proof note. Each R_J Q_J R_J^T has positive diagonal and is strictly
  diagonally dominant. Thus every Q_J is positive definite exactly.

invariants_degree2to4.json and invariants_degree5.json:
  Complete invariant bases for each multipole-rank multiset of even total
  rank. The integer basis vectors are null vectors of the raising operator
  in the unnormalized symmetric spin basis. Moving one (L,m) component to
  (L,m+1) has coefficient n_(L,m) (L+m+1). Completeness follows from the
  difference between the weight-zero and weight-one dimensions.
  For a sorted monomial t, its averaged y-moment is sum_k v_k[t] q_k / mult[t],
  where mult[t] is the multinomial multiplicity. There are 98 nonconstant
  invariant coordinates: 3 quadratic, 4 cubic, 26 quartic, and 65 quintic.

moment_matrices.json:
  The 99 rational coefficient matrices for each kernel-restricted Gram
  matrix, in the invariant-coordinate ordering used by the verifier.
  Coordinate zero is the constant 1. These cached matrices are regenerated
  and compared exactly by the default verification, so they are not trusted
  as an unexplained numerical input.

objective.json:
  The exact expansion of S=Tr(rho^2)-1/6 in the same coordinates.
  The first three nonconstant coordinates occur as
      S = 13860 q_1 + 504 q_2 + (756/5) q_3.
  Higher-degree terms cancel exactly in the certificate.

The high-precision semidefinite searches were used only to discover the
certificate. The proof is the supplied rational identity and positivity
verification, independent of solver tolerances or numerical convergence.
