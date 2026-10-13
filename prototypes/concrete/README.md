# Prototype: compact representation for concrete modes

`ConcreteProto.jl` converts expressions whose operators and parameters all have integer indices (no
sums, δs, expectation values or correlators) to a compact representation with a fixed set of modes,
computes products, sums and commutators there, and converts the result back to a `QuExpr`:

```julia
CH, CA = concrete(H, A)        # common basis (modes and parameters of both expressions)
R = QuExpr(comm(CH, CA))       # equal to comm(H, A)
```

Each monomial is a fixed number of 64-bit words: the exponents (m, n) of a†ᵐaⁿ for each boson mode
(products by Wick's theorem), a Pauli matrix for each two-level system (products are a single Pauli
matrix times a phase), and the exponent of each parameter. Pairs of monomials on different modes are
skipped in commutators.

`test_concrete.jl` compares products, sums and commutators of random expressions with QuantumAlgebra
(for both `use_σpm` settings). `bench_concrete.jl` gives (Julia 1.13, M-series Mac; all results
equal to QuantumAlgebra's):

| case | terms | QuantumAlgebra | prototype | speedup |
|---|---|---|---|---|
| Tavis–Cummings (N=4), [H,a], depth 8 | 22008 | 0.544 s | 0.068 s | 8× |
| Tavis–Cummings (N=4), H⁴ | 1185921 | 0.694 s | 0.154 s | 4.5× |
| XXZ chain (L=10), [H,τˣ₅], depth 8 | 101701 | 3.878 s | 0.082 s | 48× |
| 3 driven Kerr modes, [H,b₁], depth 7 | 1165 | 0.021 s | 0.001 s | 20× |

Converting back to `QuExpr` takes a small fraction of the time. For Tavis–Cummings, the coefficients
are `Complex{Rational{Int}}` (from σ± = (σˣ ± iσʸ)/2), whose arithmetic costs a factor 2.2 compared
to `ComplexF64`.

Not supported yet: fermions (signs from the order of the modes), qudits, spins and so(n) generators
(with a local basis such as matrix units or generalized Pauli matrices), expectation values and
correlators (as further commuting variables).
