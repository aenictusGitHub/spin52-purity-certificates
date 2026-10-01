# Methods and relation to the literature

## 1. The optimization problem and what counts as a proof

A spin-5/2 density matrix is a positive semidefinite $6\times6$ Hermitian matrix of trace one. In the symmetric five-qubit representation, $t$-anticoherence means that its $t$-qubit reduced state is $I_{t+1}/(t+1)$. Equivalently, its irreducible multipoles of ranks $1,\ldots,t$ vanish. The relationship between anticoherence and symmetric reductions is developed in [Baguette and Martin (2017)](https://doi.org/10.1103/PhysRevA.96.032304); the accompanying manuscript gives the mixed-state setting.

The feasible set is convex, but the objective $P=\operatorname{Tr}\rho^2$ is also convex and is **maximized**. The original problem therefore is not a convex semidefinite minimization problem. Finding a feasible state of purity $P_*$ supplies a lower bound on the maximum, even if many local searches return that same state. A global result requires a second argument: $P\le P_*$ for every feasible state.

This repository supplies both parts. `verify_optimizers.py` checks attaining density matrices exactly. The three proof scripts express the difference between a proposed upper bound and the purity as nonnegative quantities. Their identities apply to arbitrary feasible states, so no search over ranks or eigenbases remains.

## 2. Shared conventions: multipoles, rotations, and time reversal

Use the descending magnetic basis $m=5/2,3/2,\ldots,-5/2$ and Hilbert–Schmidt orthonormal tensors $T_{LM}$. Their matrix elements are constructed from Clebsch–Gordan coefficients with the phases of [Varshalovich, Moskalev, and Khersonskii (1988)](https://doi.org/10.1142/0270). In particular,

$$\rho=I_6/6+\sum_{L=t+1}^5\sum_{M=-L}^L x_{LM}T_{LM},\qquad P=1/6+\sum_{L,M}|x_{LM}|^2.$$

Hermiticity is imposed as $x_{L,-M}=(-1)^M x_{LM}^*$. Complex components are retained; a proof for real coefficient vectors alone would not suffice.

Spin time reversal is the antiunitary $\Theta$ satisfying $\Theta\mathbf J\Theta^{-1}=-\mathbf J$. It preserves positivity and flips a Hermitian rank-$L$ multipole by $(-1)^L$. Averaging a scalar expression over time reversal removes terms of odd total multipole rank. Averaging over rotations, or equivalently summing within coupled angular-momentum sectors, removes orientation dependence. Neither operation requires a density matrix to be invariant under time reversal or to be real.

These symmetry reductions make the certificate search much smaller. Their general mathematical rationale is related to the invariant sum-of-squares and semidefinite reductions of [Gatermann and Parrilo (2004)](https://arxiv.org/abs/math/0211450). The specific continuous spin-rotation couplings and certificate matrices here are stated explicitly in the proof note and checked from their defining coefficients.

## 3. Order two: a matrix-localizing certificate of degree five

For a 2-AC state the remaining multipoles have ranks 3, 4, and 5. Form a feature list from a constant, these linear multipoles, and their coupled quadratic products. Couple each feature to the physical spin $j=5/2$ and group by total spin $J$.

A matrix of the form $F\rho F^\dagger$ is positive semidefinite whenever $\rho\succeq0$. This is the core of the localizing construction: multiply the physical positivity constraint by polynomial features without losing its sign. Rotation and time-reversal averages give positive matrices $G_J(\rho)$. The certificate is

$$\frac{13}{18}-P=\sum_{J=1/2}^{17/2}\operatorname{Tr}\!\left[Q_J C_J^{\mathsf T}G_J(\rho)C_J\right].$$

The nine rational matrices $Q_J$ have dimensions $7,11,14,12,11,6,3,1,1$. Each summand is nonnegative because both factors under the trace are positive semidefinite. For a quadratic feature list, $F\rho F^\dagger$ has polynomial degree at most five. The stored identity cancels all terms above degree two exactly.

The coefficient comparison uses complete invariant bases: 3 quadratic, 4 cubic, 26 quartic, and 65 quintic coordinates, plus the constant. The verifier checks that the proposed basis vectors are annihilated by the raising operator, are independent, and have the dimension obtained from the difference between weight-zero and weight-one multiplicities. It reconstructs the features from exact Clebsch–Gordan coefficients and, in the default full run, regenerates every localizing coefficient matrix.

Positivity is checked without numerical eigenvalue tolerances. For each block, a supplied rational congruence $R_J$ makes $R_JQ_JR_J^{\mathsf T}$ strictly diagonally dominant with positive diagonal. This proves the transformed real symmetric matrix positive definite. It also implies that $R_J$ is nonsingular, so $Q_J$ is positive definite. All comparisons use exact rational arithmetic.

Polynomial positivity certificates and moment/semidefinite relaxations are the general framework developed by [Parrilo (2000)](https://doi.org/10.7907/2K6Y-CH43) and [Lasserre (2001)](https://doi.org/10.1137/S1052623400366802). Here we use one explicit finite matrix-localizing identity. Its validity does not depend on an asymptotic hierarchy converging, nor does this repository claim that every nonnegative polynomial admits a sum-of-squares decomposition of this degree.

## 4. Order three: cancellation in a positive three-copy witness

Only rank-four and rank-five multipoles survive for a 3-AC state. Write their squared norms as $a_4$ and $a_5$, so $P=1/6+a_4+a_5$. Consider rotational averages of $\rho^{\otimes3}$ and $\widetilde\rho\otimes\rho\otimes\rho$, where $\widetilde\rho=\Theta\rho\Theta^{-1}$. Both are positive: each is an average of positive tensor products. This argument makes **no assumption about a partial transpose**.

Three explicitly specified coupled vectors define nonnegative expectations $e_u,e_v,e_w$. Exact angular-momentum algebra gives

$$\frac12-P=\frac{84e_u+624e_v+1092e_w}{25}.$$

Besides the scalar and quadratic coefficients, there are two symmetric cubic invariants, with rank triples $(4,4,4)$ and $(4,5,5)$. `prove_purity3.py` derives their coefficients and verifies that this positive weighted sum cancels both cubic terms. The remaining expression is exactly $1/2-P$.

This is a direct finite witness using positivity, tensor products, and representation theory. The required vectors, phases, and weights appear in the proof note; the exact script uses no stored floating-point witness or optimization result. The separately verified rank-two state attains $P=1/2$.

## 5. Order four: a seven-term degree-six identity

A 4-AC state has $\rho=I_6/6+X$ with only rank-five multipoles. Time reversal sends $X$ to $-X$. Hence positivity of $\rho$ and its time reverse implies $-I_6/6\preceq X\preceq I_6/6$, or

$$B(X)=I_6/36-X^2\succeq0,\qquad S=\operatorname{Tr}X^2=P-1/6.$$

The certificate combines four positive localizing terms built from $B(X)$ and quadratic coupled multipoles with three squared norms of cubic coupled multipoles:

$$S^2(1/8-S)=\sum_{r=1}^4 c_r\mathcal L_{J_r}(w_r)+\sum_{r=1}^3d_r\mathcal C_{K_r}(h_r).$$

All seven weights are explicitly positive rational numbers. Each right-hand term is nonnegative. If $S>0$, divide by $S^2$ to get $S\le1/8$, hence $P\le7/24$. The omitted case $S=0$ is the maximally mixed state and already satisfies the bound.

`prove_purity4_sixth.py` expands all seven terms directly after the substitution $x_m=\sqrt{\binom{10}{5+m}}\,y_m$. Rational square-root scale factors are handled symbolically before the final coefficient comparison, which uses rational arithmetic. It compares every monomial against $S^2(1/8-S)$ and checks the explicit attaining state. Neither cached moment matrices nor a numerical positive-semidefinite test are used for this proof.

The factor $S^2$ is part of this certificate. Dropping it or dividing at the maximally mixed state would not be a valid proof. The bound uses both physical positivity and the rank-five multipole restriction; it is not a statement for all traceless Hermitian $6\times6$ matrices.

## 6. Numerical discovery versus exact verification

The development used numerical searches to suggest candidate optimal states and symmetry-reduced semidefinite certificates. Kernel restrictions and rational reconstruction then produced compact exact data. This numeric-to-exact strategy is closely related to [Peyrl and Parrilo (2008)](https://doi.org/10.1016/j.tcs.2008.09.025), who study rational sum-of-squares recovery from numerical solutions. We cite that general strategy, without claiming that every step of our reconstruction implements their algorithm or meets every hypothesis in that paper.

The published package contains the final certificates and verifiers. It does not require the original optimization trajectory to be reproduced. In the order-two case, the default run rebuilds the cached coefficient matrices; in the order-three and order-four cases, the required coefficients are derived directly. Exact attainment closes the gap between the universal bound and a feasible state.

The optional random complex-state checks exercise a different implementation of the coupled operators. They are useful for detecting transcription and convention errors, but finite sampling cannot establish a global theorem. The exact identities and positivity arguments are the proof.

## References

1. J. Denis, T. Lacaille, J. Martin, and E. Serrano-Ensástiga, *Total, quantum, and classical measures of anticoherence for mixed spin states*, [arXiv:2605.29436](https://arxiv.org/abs/2605.29436) (2026). Accompanying manuscript; its cited repository commit fixes the version of these certificates.
2. D. Baguette and J. Martin, *Anticoherence measures for pure spin states*, Phys. Rev. A **96**, 032304 (2017), [doi:10.1103/PhysRevA.96.032304](https://doi.org/10.1103/PhysRevA.96.032304).
3. D. A. Varshalovich, A. N. Moskalev, and V. K. Khersonskii, *Quantum Theory of Angular Momentum* (World Scientific, 1988), [doi:10.1142/0270](https://doi.org/10.1142/0270).
4. P. A. Parrilo, *Structured Semidefinite Programs and Semialgebraic Geometry Methods in Robustness and Optimization*, Ph.D. thesis, Caltech (2000), [doi:10.7907/2K6Y-CH43](https://doi.org/10.7907/2K6Y-CH43).
5. J. B. Lasserre, *Global optimization with polynomials and the problem of moments*, SIAM J. Optim. **11**, 796–817 (2001), [doi:10.1137/S1052623400366802](https://doi.org/10.1137/S1052623400366802).
6. K. Gatermann and P. A. Parrilo, *Symmetry groups, semidefinite programs, and sums of squares*, J. Pure Appl. Algebra **192**, 95–128 (2004), [arXiv:math/0211450](https://arxiv.org/abs/math/0211450), [doi:10.1016/j.jpaa.2003.12.011](https://doi.org/10.1016/j.jpaa.2003.12.011).
7. H. Peyrl and P. A. Parrilo, *Computing sum of squares decompositions with rational coefficients*, Theor. Comput. Sci. **409**, 269–281 (2008), [doi:10.1016/j.tcs.2008.09.025](https://doi.org/10.1016/j.tcs.2008.09.025).

The last four works are methodological background, not prior claims of the particular spin-5/2 optimum values or identities given here.
