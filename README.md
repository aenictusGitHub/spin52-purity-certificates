# Exact purity and convex-roof certificates for spin 5/2

This repository contains the full proofs, exact coefficient data, and verification programs for the global purity optima of mixed spin-5/2 states with anticoherence orders 2, 3, and 4. It accompanies *Total, quantum, and classical measures of anticoherence for mixed spin states* by Jérôme Denis, Tara Lacaille, John Martin, and Eduardo Serrano-Ensástiga ([arXiv:2605.29436](https://arxiv.org/abs/2605.29436)). It collects the material removed from the revised manuscript appendices on purity certificates and exact rank-four convex roofs. Version 1.1.0 adds the full proofs and verification data for the purity-based quantum anticoherence of the rank-four state.

For a spin-5/2 density operator, equivalently a symmetric state of five qubits, define

$$P_{\max}(5,t)=\max\{\mathrm{Tr}(\rho^2):\rho\succeq0,\ \mathrm{Tr}\rho=1,\ \rho_t=I_{t+1}/(t+1)\}.$$

| AC order $t$ | Exact maximum purity | Upper-bound certificate | Rank of an attaining state |
|---|---|---|---|
| 2 | $13/18$ | Nine positive localizing blocks; polynomial degree at most five | 2 |
| 3 | $1/2$ | Positive three-copy witness | 2 |
| 4 | $7/24$ | Seven-term degree-six positivity identity | 4 |

These bounds hold for **all feasible density operators**, without assuming a rank, a real matrix, or a preferred eigenbasis. Explicit states attain each bound. The exact verification does not call a numerical optimization solver.

For the rank-four fourth-order state, the purity-based quantum contributions are also certified globally:

$$q_1=\frac{694-4\sqrt{22}-2\sqrt{6071+448\sqrt{22}}}{675}=0.73249884302689881566\ldots,$$

$$q_2=\frac{171}{400}+\frac{u+4v}{12}=0.67058403099168934995\ldots.$$

Here $u,v$ are specified by a uniquely isolated root of a four-variable polynomial system in the [convex-roof proof note](docs/convex-roof-certificates.pdf). Matching affine lower bounds and attaining six- and twelve-state ensembles prove these values. Complementary orders give $q_3=(8/9)q_2$ and $q_4=(5/8)q_1$.

## Read the proofs

- [Methods and relation to the literature](METHODS.md): why the certificates prove global optima, how rotation symmetry and time reversal are used, and the distinction between numerical discovery and exact verification.
- [Complete purity proof note](docs/purity-certificates.pdf), with [editable LaTeX source](docs/purity-certificates.tex): all three proofs, the attaining states, coefficient conventions, and references.
- [Complete convex-roof proof note](docs/convex-roof-certificates.pdf), with [editable LaTeX source](docs/convex-roof-certificates.tex): state definition, phase minimization, exact formulas, root isolation, positive Gram matrices, and attaining ensembles.
- [Convex-roof data and verification guide](convex_roofs/README.md).
- [Order-two data conventions](global_purity2/README.txt): invariant coordinates, localizing blocks, and the rational positivity checks.
- [Bibliography](docs/references.bib): source references in BibTeX format.

## Verify from a clean checkout

Python 3.10 or later, NumPy, SymPy, and mpmath are sufficient. Use a virtual environment if desired:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python verify_all.py
```

On Windows, activate with `.venv\Scripts\activate` instead. The full command regenerates the order-two coefficient matrices from Clebsch–Gordan coefficients, rather than trusting the stored matrices. Depending on the machine, exact symbolic verification can take several minutes. The progress output identifies the block or certificate being checked. Do not run the verifiers with `python -O`, because the proof checks use assertions.

Optional commands:

```bash
python verify_roofs.py
python verify_all.py --quick
python verify_all.py --numerical-crosschecks
```

`verify_roofs.py` runs only the two convex-roof proofs and their independent support-polynomial check. `--quick` skips only regeneration of the cached order-two localizing coefficient matrices; the invariant-basis, coupling, positivity, polynomial-identity, and optimizer checks still run. Use the default full verification for an independent reproduction of the certificate. The optional numerical cross-checks evaluate coupled operators directly on complex states; they provide an additional implementation check and are **not** the proof.

## Files and proof obligations

| File | What it verifies |
|---|---|
| `verify_roofs.py` | Physical support polynomials, both exact convex roofs, and the four-variable algebraic reduction |
| `convex_roofs/` | Exact roof coefficient data, verifiers, and rational interval reports |
| `verify_optimizers.py` | Exact positivity, trace one, maximally mixed reductions, and attainment for all three states |
| `prove_purity2.py` | Complete invariant bases; regenerated coupled features and coefficient matrices; rational positivity of nine blocks; exact cancellation to $13/18-\mathrm{Tr}\rho^2$ |
| `prove_purity3.py` | Exact three-copy expectations and cancellation to $1/2-\mathrm{Tr}\rho^2$ |
| `prove_purity4_sixth.py` | Direct expansion of the seven nonnegative terms and exact attainment of $7/24$ |
| `verify_purity2_independently.py` | Supplementary numerical check of the order-two identity on complex Hermitian matrices |
| `verify_purity3_independently.py` | Supplementary direct three-copy check on complex 3-AC states |
| `angular_momentum.py` | Tensor construction used by the numerical cross-checks |
| `global_purity2/*.json` | Exact integer/rational data for the order-two certificate |
| `verification/` | Recorded full-run output and runtime versions for this distribution |

The proof note and repository data form a complete statement of the certificate. A small eigenvalue returned by a floating-point solver, or a small numerical residual, is never used as the final proof of positivity. Earlier exploratory optimization runs and the manuscript's unrelated data are omitted.

## Cite and reproduce

Cite the accompanying manuscript and the repository **at the commit referenced in that manuscript**. A commit-specific URL identifies the exact scripts and coefficient data even if the default branch changes. `CITATION.cff` provides bibliographic metadata. The methods references explain the mathematical background; the particular spin-5/2 certificates and attaining states are the results provided here.

`SHA256SUMS` records the hashes of the distributed sources, data, proof note, and verification record. On systems with `shasum`, check them with `shasum -a 256 -c SHA256SUMS` from the repository root.
