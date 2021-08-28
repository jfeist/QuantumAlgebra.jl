# QuantumAlgebra.jl - quantum operator algebra in Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfeist.github.io/QuantumAlgebra.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfeist.github.io/QuantumAlgebra.jl/dev)
[![Build Status](https://github.com/jfeist/QuantumAlgebra.jl/workflows/CI/badge.svg)](https://github.com/jfeist/QuantumAlgebra.jl/actions)
[![Coverage](https://codecov.io/gh/jfeist/QuantumAlgebra.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jfeist/QuantumAlgebra.jl)
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/jfeist/QuantumAlgebra.jl/main?filepath=examples)
[![DOI](https://zenodo.org/badge/211471154.svg)](https://zenodo.org/badge/latestdoi/211471154)

This package does quantum operator algebra (i.e., algebra with non-commuting
operators) in Julia, supporting bosonic, fermionic, and two-level system
operators, with arbitrary names and indices, as well as sums over any of the
indices. It defines an opinionated canonical form (normal ordering plus some
additional rules) to automatically simplify expressions. It is recommended to
use an interface that can display LaTeX formulas (e.g., Jupyter notebooks) for
convenient output formatting. While there is some documentation, it is not
always kept fully up to date, and it is recommended to look at the latest commit
messages to get an idea about new features etc. You can also check out the
notebooks in the `examples` folder, which can be viewed online with
[nbviewer](https://nbviewer.jupyter.org/github/jfeist/QuantumAlgebra.jl/blob/main/examples/)
and tried out interactively with
[Binder](https://mybinder.org/v2/gh/jfeist/QuantumAlgebra.jl/main?filepath=examples).

# Updates in v1.0.0
This is a major revision with some breaking changes. The backend has been almost
completely rewritten to make the code more efficient when dealing with large
expressions, and the interface has been cleaned up in several places.

### Important changes:
- Canonical normal form is **not** automatically enforced by default. In order
  to transform expressions to normal form, use `normal_form(x)`. Since automatic
  conversion to normal form can be convenient for interactive work, it can be
  enabled with `QuantumAlgebra.auto_normal_form(true)`, or alternatively by
  setting the environment variable `QUANTUMALGEBRA_AUTO_NORMAL_FORM` to `"true"`
  (or any value that `parse(Bool,value)` parses as `true`) before `using
  QuantumAlgebra`.
- The function to obtain expectation values is now `expval(A)` (instead of
  `ExpVal`), and `expval_as_corrs(A)` to express an expectation value through a
  correlator / cumulant expansion, e.g., ⟨_AB_⟩ = ⟨_AB_⟩<sub>c</sub> -
  ⟨_A_⟩<sub>c</sub> ⟨_B_⟩<sub>c</sub>, with corresponding extensions for
  products of more operators. Note that for a single operator, ⟨_A_⟩<sub>c</sub>
  = ⟨_A_⟩, but we distinguish the two formally for clarity.
- There is a new function `julia_expression(A)` that converts a QuantumAlgebra
  object to a julia expression, which helps in using QuantumAlgebra to
  programatically derive codes for numerical implementation. The object `A`
  cannot contain any "bare" operators, but only expectation values or
  correlators. See the documentation for more details.
- QuantumAlgebra expressions are now printed in pretty format in the terminal
  etc.

# Overview

The basic functions to create QuantumAlgebra expressions (which are of type
`QuExpr`) are
- `a(inds...)` and `adag(inds...)` for _a_ and _a<sup>†</sup>_, the annihilation
  and creation operators for a bosonic mode.
- `f(inds...)` and `fdag(inds...)` for _f_ and _f<sup>†</sup>_, the annihilation
  and creation operators for a fermionic mode.
- `σx(inds...)`, `σy(inds...)`, `σz(inds...)` for the Pauli matrices
  _σ<sup>x,y,z</sup>_ for a two-level system (TLS).
- `σp(inds...)`, `σm(inds...)` for excitation and deexcitation operators
  _σ<sup>±</sub>_ for a two-level system (TLS).

- **Indices**: All of these functions take an arbitrary number of indices as
  arguments, which can be either integers (1,2,...) or symbolic, where symbolic
  indices must be a single unicode character, with possibly an integer subindex:
  ```julia
  julia> using QuantumAlgebra

  julia> a()
  a()

  julia> adag(:i)
  a†(i)

  julia> fdag(1,2,:i_9)
  f†(12i₉)

  julia> σx(:i_1, 1, :j, :k_2, :μ_2, :◔_1, :😄_121)
  σˣ(i₁1jk₂μ₂◔₁😄₁₂₁)
  ```

- You can define your own bosonic/fermionic/two-level system operators with a
  set of macros:
  - `@boson_ops name` defines new functions `$name()` and `$(name)dag()` for
    bosonic species `name`.
  - `@fermion_ops name` defines new functions `$name()` and `$(name)dag()` for
    fermionic species `name`.
  - `@tlsxyz_ops name` defines new functions `$(name)x()`, `$(name)y()` and
    `$(name)z()` for the Pauli matrices for two-level system species `name`.
  - `@tlspm_ops name` defines new functions `$(name)p()` and `$(name)m()` for
    the two-level system excitation and deexcitation operators for species
    `name`.
  ```julia
  julia> @boson_ops b
  (QuantumAlgebra.OpConstructors.b, QuantumAlgebra.OpConstructors.bdag)

  julia> bdag(:k)*b(:i)
  b†(k) b(i)
  ```
  Operators with different names are assumed to belong to different "species"
  and always commute. For fermions, this is not always desired, since you might
  want to use different named operators to refer to different kinds of states
  for the same species (e.g., localized and itinerant electrons). This can be
  achieved with the macro `@anticommuting_fermion_group`, which creates several
  fermionic operators that mutually anticommute:
  ```julia
  julia> @anticommuting_fermion_group c d

  julia> normal_form(c()*d() + d()*c())
  0
  ```

- `param(name::Symbol,state='n',inds...)` to create a named parameter. `state` must be
  one of `'r'`, `'n'`, or `'c'` for purely real, non-conjugated complex, and
  conjugated complex parameters. More conveniently, parameters can be entered
  with string macros `Pr"name_inds..."` and `Pc"name_inds..."` for real and
  complex parameters:
  ```julia
  julia> Pr"g_i,j_2,k"
  g(ij₂k)

  julia> Pr"g_i,j_2,k" == param(:g,'r',:i,:j_2,:k)
  true

  julia> Pc"α_3" == param(:α,3)
  true
  ```

- Arithmetic operations (`*`, `+`, `-`, `^`, `adjoint`=`'`) are supported
  (exponents must be nonnegative integers), with any `Number` types integrating
  automatically.
  ```julia
  julia> 5*adag(:k)*f(3)*σx(3)
  5 a†(k) f(3) σˣ(3)

  julia> (5//3+4im) * adag(:k)*f(3)*σx(3) + 9.4
  9.4 + (5//3+4i) a†(k) f(3) σˣ(3)

  julia> (a(:i)*f(:k))'
  f†(k) a†(i)
  ```
  If you need a bare number as a QuantumAlgebra expression, you can use
  `x*one(QuExpr)` (or `one(A)`, where `A` is any `QuExpr`).

- `∑(ind,A::QuExpr)` to represent an analytic sum over index `ind`. Since summed
  indices have no semantic meaning, the index within the expression gets
  replaced by a special numbered sum index `#ᵢ`, with `i=1,2,...`.
  ```julia
  julia> ∑(:i,a(:i))
  ∑₁ a(#₁)
  ```

- `normal_form(A::QuExpr)` converts an expression to a well-defined "canonical"
  order. To achieve this canonical form, relevant commutators etc are used, so
  an expression written as a single product can turn into a sum of expressions.
  The order is essentially normal ordering (creation before annihilation
  operators, with σˣʸᶻ in the middle), with some additional conventions to make
  the normal form (hopefully) unique. In some contexts (e.g., interactive work),
  it can be convenient to automatically transform all expressions to normal
  form. This can be enabled by calling `QuantumAlgebra.auto_normal_form(true)`,
  or alternatively by setting the environment variable
  `QUANTUMALGEBRA_AUTO_NORMAL_FORM` to `"true"` (or any value that
  `parse(Bool,value)` parses as `true`) before `using QuantumAlgebra`.
  ```julia
  julia> normal_form(a(:i)*adag(:j))
  δ(ij)  + a†(j) a(i)
  ```

- `expval(A::QuExpr)` to represent an expectation value.
  ```julia
  julia> expval(adag(:j)*a(:i))
  ⟨a†(j) a(i)⟩
  ```

- `expval_as_corrs(A::QuExpr)` to represent an expectation value through its
  correlators, i.e., a cumulant expansion.
  ```julia
  julia> expval_as_corrs(adag(:j)*a(:i))
  ⟨a†(j)⟩c ⟨a(i)⟩c  + ⟨a†(j) a(i)⟩c
  ```

- `comm(A::QuExpr,B::QuExpr)` to calculate the commutator [_A,B_] = _AB - BA_.
  ```julia
  julia> comm(a(),adag())
  -a†() a() + a() a†()

  julia> normal_form(comm(a(),adag()))
  1
  ```

- `Avac(A)` and `vacA(A)` simplify operators by assuming they are applied to the
  vacuum from the left or right, respectively. To be precise, `Avac(A)` returns
  _A'_ such that _A_|0⟩ = _A'_|0⟩, while `vacA(A)` does the same for ⟨0|_A_.
  These functions automatically apply `normal_form` to assure that the operators
  are simplified as much as possible. Note that "vacuum" for two-level systems
  is interpreted as the lower state, `σᶻ|0⟩ = -|0⟩`.
  ```julia
  julia> Avac(a())
  0

  julia> Avac(a(:i)*adag(:j))
  δ(ij)

  julia> Avac(a()*adag()*adag())
  2 a†()

  julia> vacA(a()*adag()*adag())
  0

  julia> Avac(σx())
  σˣ()

  julia> Avac(σz())
  -1
  ```

- `vacExpVal(A,S=1)` calculates the vacuum expectation value
  ⟨0|_S<sup>†</sup>AS_|0⟩, i.e., the expectation value ⟨ψ|_A_|ψ⟩ for the state
  defined by |ψ⟩=_S_|0⟩. The result is guaranteed to not contain any operators.
  ```julia
  julia> vacExpVal(adag()*a())
  0

  julia> vacExpVal(adag()*a(), adag()^4/sqrt(factorial(4)))
  4.000000000000001

  julia> vacExpVal(adag()*a(), adag()^4/sqrt(factorial(big(4))))
  4

  julia> vacExpVal(σx())
  0
  ```

- `julia_expression(A)` to obtain a julia expression that can be used to
  automatically build codes implementing equations derived with QuantumAlgebra.
  Every expectation value or correlator is treated as a separate array. Daggers
  are represented as `ᴴ`, which are valid identifiers that can appear in the
  array names. Note that expectation values and correlators are not
  distinguished, so it is best to have all expressions use the same kind.
  ```julia
  julia> julia_expression(expval_as_corrs(adag(:j)*a(:i)))
  :(aᴴ[j] * a[i] + aᴴa[j, i])
  ```
  Also note that expressions are always treated as arrays, even if they have no
  indices (which gives zero-dimensional arrays). If you are working with scalar
  quantities exclusively, it might be useful to clean up the resulting
  expression (e.g., use `MacroTools` to remove the `[]`).
  ```julia
  julia> julia_expression(expval(adag()*a()*σx()))
  :(aᴴaσˣ[])
  ```

- By default, two-level system operators are represented by the Pauli
  matrices `σˣʸᶻ`, and calling `σp()` and `σm()` will give results expressed through them:
  ```julia
  julia> σp()
  1//2 σˣ() + 1//2i σʸ()

  julia> σm()
  1//2 σˣ() - 1//2i σʸ()
  ```
  This can be changed by calling `QuantumAlgebra.use_σpm(true)`. In this mode,
  `σ⁺` and `σ⁻` are the "fundamental" operators, and all expressions are written in terms of them. Note that mixing conventions within the same expression is not supported, so it is suggested to set this flag once at the beginning of any calculation.
  ```julia
  julia> QuantumAlgebra.use_σpm(true)

  julia> σp()
  σ⁺()

  julia> σx()
  σ⁺() + σ⁻()

  julia> σz()
  -1 + 2 σ⁺() σ⁻()
  ```

## Citing

If you use QuantumAlgebra in academic work, we would appreciate a citation. See
[`CITATION.bib`](CITATION.bib) for the relevant references.
