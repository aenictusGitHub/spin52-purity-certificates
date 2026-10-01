# Exact quantum-anticoherence roofs of the rank-four state

The [complete proof note](../docs/convex-roof-certificates.pdf) defines the state and gives all derivations. It preserves the full development removed from the manuscript appendix **Exact roofs for the rank-four state**. These certificates concern the **purity-based** quantum measure.

The results are

$$q_1=\frac{694-4\sqrt{22}-2\sqrt{6071+448\sqrt{22}}}{675}=0.73249884302689881566\ldots,$$

$$q_2=\frac{171}{400}+\frac{u+4v}{12}=0.67058403099168934995\ldots.$$

The proof note specifies $u,v$ through a unique root of a rational polynomial system in four variables. This is an exact algebraic specification with no remaining minimization; neither a compact radical expression nor a minimal polynomial for $q_2$ is claimed. Complementary reductions give $q_3=(8/9)q_2$ and $q_4=(5/8)q_1$.

## Run

From the repository root:

```bash
python verify_roofs.py
```

This also runs as part of `python verify_all.py`. Do not use Python's `-O` option. The checks require SymPy and mpmath, use exact algebraic and rational arithmetic, and do not require an optimization solver.

## Proof obligations and files

1. `verify_support.py` reconstructs the complex pure-state polynomials from the symmetric five-qubit reductions, the first-order phase bound from spin operators, and the rotation characters. This independently checks the physical input to the certificate.
2. `verify_q1.py` checks the radical value, affine-bound polynomial identity, positive Gram matrices, feasible ensemble parameter, and matching ensemble objective. `q1_field.json` and `q1_matrix{0,1}.json` contain coefficients in descending powers of the primitive element $\sqrt2+\sqrt{11}+\sqrt{6071+448\sqrt{22}}$.
3. `verify_q2_algebraic.py` checks the complex-phase bound, rationally rescaled polynomial, Gram identity modulo eight contact equations, unique root by rational interval contraction, strict positive Gram pivots, and positive attaining-ensemble weights. `q2_algebraic_system.json` stores the contact polynomials and Gram data in variables $(A_-,B_-,C_-,A_+,B_+,C_+,u,v)$.
4. `verify_q2_reduced.py` checks the quartic discriminant identity, exact amplitude elimination, unique root of the four-variable specification, and its equivalence to the original certified root. `q2_reduced_polynomial.json` stores $\Phi/202649056051200$, with variables $(A,u,v)$; this nonzero scale leaves the roots unchanged. Degree 22 refers to the auxiliary variable $A$, not to a claimed minimal polynomial degree of $q_2$.

`poly_t2.txt` is the complex-state polynomial in the real amplitudes `x0`–`x3` and imaginary amplitudes `x4`–`x7`. `q2_600digits.json` contains **numerical guesses**, used only to propose rational box centers and inverse matrices. The ensuing rational inequalities, rather than numerical residuals, certify the result. The fine root box encloses $q_2$ in an interval narrower than $10^{-60}$.

The scripts regenerate `q1_verification.json`, `q2_algebraic_verification.json`, and `q2_reduced_verification.json`. These reports record exact rational enclosures, inverse matrices, positive pivots, and ensemble weights. Run the original eight-variable verification before the reduced verifier, which reads its report; the driver enforces this order. A recorded successful run is in [verification/roof-run.txt](../verification/roof-run.txt).

The attaining ensembles have six components for $q_1$ and twelve for $q_2$. The [methods guide](../METHODS.md) explains the connection to supporting affine functionals, symmetry reduction, sum-of-squares certificates, and verified root isolation, with references to the literature.
