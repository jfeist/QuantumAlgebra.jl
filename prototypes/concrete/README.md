# Prototype: compact representation for concrete modes

`ConcreteProto.jl` converts expressions whose operators and parameters all have integer indices (no
sums, δs, expectation values or correlators) to a compact representation with a fixed set of modes,
computes products, sums and commutators there, and converts the result back to a `QuExpr`:

```julia
CH, CA = concrete(H, A)        # common basis (modes and parameters of both expressions)
R = QuExpr(comm(CH, CA))       # equal to comm(H, A)
```

Each monomial is a fixed number of 64-bit words:
- bosons: the exponents (m, n) of a†ᵐaⁿ for each mode, products by Wick's theorem;
- two-level systems with Pauli matrices: a Pauli matrix for each site, products are a single Pauli
  matrix times a phase (bit operations);
- all other operators (fermions, two-level systems with σ±, qudits, spins, so(n) generators): a label
  for a product of operators on each site, in QuantumAlgebra's normal form. Products of two such local
  monomials are computed once with QuantumAlgebra and cached, so the results agree with QuantumAlgebra
  by construction (e.g., QuantumAlgebra does not reduce Jᶻ² = 1/4 for spins with j = 1/2, which a
  matrix representation would). Fermions with the same name give signs when they are reordered;
- parameters: the exponent of each parameter.

Pairs of monomials that commute trivially are skipped in commutators.

`test_concrete.jl` compares products, sums and commutators of random expressions with QuantumAlgebra,
for all kinds of operators and both `use_σpm` settings. `bench_concrete.jl` gives (Julia 1.13,
M-series Mac; all results equal to QuantumAlgebra's; powers are normal-ordered after each product):

| case | terms | QuantumAlgebra | prototype | speedup |
|---|---|---|---|---|
| `use_σpm(false)` (Pauli matrices) | | | | |
| Tavis–Cummings (N=4), [H,a], depth 8 | 22008 | 0.579 s | 0.078 s | 7× |
| Tavis–Cummings (N=4), H⁵ | 34457 | 1.203 s | 0.101 s | 12× |
| XXZ chain (L=10), [H,τˣ₅], depth 8 | 101701 | 3.959 s | 0.141 s | 28× |
| `use_σpm(true)` (σ±) | | | | |
| Tavis–Cummings (N=4), [H,a], depth 8 | 15267 | 0.215 s | 0.020 s | 11× |
| Tavis–Cummings (N=4), H⁵ | 9432 | 0.131 s | 0.0045 s | 29× |
| XXZ chain (L=10), [H,τˣ₅], depth 6 | 37642 | 0.944 s | 0.068 s | 14× |
| other operators | | | | |
| 3 driven Kerr modes, [H,b₁], depth 7 | 1165 | 0.020 s | 0.0014 s | 15× |
| Fermi–Hubbard chain (L=4), [H,c₁↑], depth 6 | 320 | 0.0031 s | 0.0011 s | 3× |
| Λ-type qudits in a cavity (N=3), [H,q¹²₁], depth 6 | 2329 | 0.046 s | 0.0028 s | 16× |

Which basis is efficient depends on the expressions: Tavis–Cummings has fewer terms with σ± (and
integer coefficients, while σ± = (σˣ ± iσʸ)/2 gives rational ones with Pauli matrices), the XXZ chain
has fewer terms with Pauli matrices (σᶻ = 2σ⁺σ⁻ - 1 has two terms). Converting back to `QuExpr` takes
a small fraction of the time.

Possible improvements: bit operations for σ± (as for Pauli matrices) instead of the cached products,
fermion signs with bit masks, choosing the basis of two-level systems independently of `use_σpm`,
expectation values and correlators as further commuting variables.
