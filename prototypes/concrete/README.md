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
  by construction (also where a matrix representation does not exist, e.g., for spins without fixed or
  with symbolic j, or for so(n) generators). Fermions with the same name give signs when they are
  reordered;
- parameters: the exponent of each parameter.

Pairs of monomials that commute trivially are skipped in commutators.

`test_concrete.jl` compares products, sums and commutators of random expressions with QuantumAlgebra,
for all kinds of operators and both `use_σpm` settings. `bench_concrete.jl` gives (Julia 1.13,
M-series Mac; all results equal to QuantumAlgebra's; powers are normal-ordered after each product):

| case | terms | QuantumAlgebra | prototype | speedup |
|---|---|---|---|---|
| `use_σpm(false)` (Pauli matrices) | | | | |
| Tavis–Cummings (N=4), [H,a], depth 8 | 22008 | 0.534 s | 0.075 s | 7× |
| Tavis–Cummings (N=4), H⁴ | 8972 | 0.155 s | 0.0089 s | 17× |
| XXZ chain (L=10), [H,τˣ₅], depth 8 | 101701 | 3.771 s | 0.135 s | 28× |
| `use_σpm(true)` (σ±) | | | | |
| Tavis–Cummings (N=4), [H,a], depth 8 | 15267 | 0.183 s | 0.019 s | 9.5× |
| Tavis–Cummings (N=4), H⁴ | 2649 | 0.011 s | 0.0010 s | 11× |
| XXZ chain (L=10), [H,τˣ₅], depth 6 | 37642 | 0.897 s | 0.072 s | 12.5× |
| other operators | | | | |
| 3 driven Kerr modes, [H,b₁], depth 7 | 1165 | 0.021 s | 0.0014 s | 15× |
| Fermi–Hubbard chain (L=4), [H,c₁↑], depth 6 | 320 | 0.0032 s | 0.0010 s | 3× |
| spin-1/2 XXZ chain (L=8), [H,Sᶻ₄], depth 6 | 6052 | 0.143 s | 0.023 s | 6× |
| spin-1 XXZ chain (L=6), [H,Kᶻ₃], depth 5 | 6840 | 0.104 s | 0.0105 s | 10× |
| collective spin j=5/2, [χJ_z²+Ω(J₊+J₋), J₋], depth 12 | 134 | 0.020 s | 0.0035 s | 6× |
| Λ-type qudits in a cavity (N=3), [H,q¹²₁], depth 6 | 2329 | 0.036 s | 0.0026 s | 14× |

Which basis is efficient depends on the expressions: Tavis–Cummings has fewer terms with σ± (and
integer coefficients, while σ± = (σˣ ± iσʸ)/2 gives rational ones with Pauli matrices), the XXZ chain
has fewer terms with Pauli matrices (σᶻ = 2σ⁺σ⁻ - 1 has two terms). For spins with j = 1/2 (and other
half-integer j), the coefficients are rational (from j(j+1)). Converting back to `QuExpr` takes a small
fraction of the time.

Possible improvements: bit operations for σ± (as for Pauli matrices) instead of the cached products,
fermion signs with bit masks, choosing the basis of two-level systems independently of `use_σpm`,
expectation values and correlators as further commuting variables.
