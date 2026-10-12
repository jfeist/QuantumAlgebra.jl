```@meta
CurrentModule = QuantumAlgebra
```
```@meta
DocTestSetup = quote
    using QuantumAlgebra
    QuantumAlgebra.auto_normal_form(false)
    QuantumAlgebra.use_σxyz()
    @static if !QuantumAlgebra._DEFINE_DEFAULT_OPS
        @boson_ops a
        @fermion_ops f
        @tlsxyz_ops σ
        @tlspm_ops σ
    end
end
```

# Legacy interface and migration

Before version 1.7, operators were created with constructor functions, e.g.,
`a(:i)` instead of `a[:i]` (after `@boson a[:]`). This interface still works in
all 1.x versions, but it is deprecated and will be removed in version 2.0. Both
interfaces create the same expressions, so they can be mixed while converting
existing code, and all other functions (`normal_form`, `comm`, `heisenberg_eom`,
...) are not affected.

## The legacy interface

- When the package is loaded, the "default" operators `a`, `f`, `σx`, `σy`,
  `σz`, `σp`, `σm` (and the deprecated `adag`, `fdag`) are defined as
  constructor functions, which take an arbitrary number of indices:
  `a(inds...)` and `a'(inds...)` give ``a`` and ``a^{†}`` for a bosonic mode,
  `f(inds...)` and `f'(inds...)` the corresponding fermionic operators,
  `σx(inds...)`, `σy(inds...)`, `σz(inds...)` the Pauli matrices and
  `σp(inds...)`, `σm(inds...)` the excitation and deexcitation operators of a
  two-level system. This can be disabled with the preference
  `"define_default_ops"`.
  ```jldoctest
  julia> a'(:i)*a()
  a†(i) a()

  julia> f'(1,2,:i_9)
  f†(12i₉)
  ```
- Constructor functions for other operators are defined by the macros
  `@boson_ops name` (defining `name` and the deprecated `namedag`),
  `@fermion_ops name`, `@tlsxyz_ops name` (defining `namex`, `namey`, `namez`)
  and `@tlspm_ops name` (defining `namep` and `namem`), and
  `@anticommuting_fermion_group name1 name2 ...` defines fermionic operators that
  mutually anticommute. The functions `boson_ops`, `fermion_ops`, `tlspm_ops`
  and `tlsxyz_ops` return the constructors.
  ```jldoctest
  julia> @boson_ops b
  (b (QuExpr constructor), b† (QuExpr constructor))

  julia> b'(:k)*b(:i)
  b†(k) b(i)
  ```
- Parameters are created with `param(name::Symbol,state='n',inds...)`, where
  `state` is one of `'r'`, `'n'`, or `'c'` for purely real, non-conjugated
  complex, and conjugated complex parameters, or with the string macros
  `Pr"name_inds..."` and `Pc"name_inds..."` for real and complex parameters.
  ```jldoctest
  julia> Pr"g_i,j_2,k"
  g(ij₂k)

  julia> Pr"g_i,j_2,k" == param(:g,'r',:i,:j_2,:k)
  true

  julia> Pc"α_3" == param(:α,3)
  true
  ```

## Migrating to the new interface

The main difference is that in the new interface, operators are variables
instead of functions, and the number of indices of each operator is fixed when
it is defined. Operators that were used with different numbers of indices (which
are treated as different modes) therefore need different names.

| Legacy interface | New interface |
|:-----------------|:--------------|
| `@boson_ops b` | `@boson b` (without indices) or `@boson b[:,:]` (with two indices) |
| `@fermion_ops f` | `@fermion f` or `@fermion f[:]` |
| `@tlspm_ops σ` (defining `σp`, `σm`) | `@tls σ` (with `σ` = ``σ^{-}`` and `σ'` = ``σ^{+}``) |
| `@tlsxyz_ops σ` | `@tlsxyz σ` or `@tlsxyz σ[:]` |
| `@anticommuting_fermion_group c d` | `@anticommuting_fermions c d` (or `c[:] d[:]`) |
| `b()`, `b'()`, `bdag()` | `b`, `b'`, `b'` |
| `b(:i,2)`, `b'(:i,2)` | `b[:i,2]`, `b'[:i,2]` (or `b[i,2]` after `@qindex i`) |
| `b(:j_2)` | `b[:j_2]` (or `b[j[2]]` after `@qindex j[:]`) |
| `Pr"g_i"`, `param(:g,'r',:i)` | `g[:i]` after `@param g[:]` |
| `Pc"α"`, `Pc"α"'` | `α`, `α'` after `@param α::Complex` |
| `boson_ops(:b)` etc. | `boson(:b; ndims=2)` etc. |
| default operators `a`, `f`, `σx`, ... | define them, e.g., `@boson a[:]` |

In addition, QuantumAlgebra currently prints operators with parentheses (e.g.,
`b†(i2)` and `b()`). With `QuantumAlgebra.use_bracket_display(true)`, they are
printed in the syntax of the new interface (e.g., `b†[i,2]` and `b`), which will
become the default in version 2.0.
