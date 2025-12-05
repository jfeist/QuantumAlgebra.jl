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

# [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl) - quantum operator algebra in Julia

This package does quantum operator algebra (i.e., algebra with non-commuting
operators) in Julia, supporting bosonic, fermionic, and two-level system
operators, with arbitrary names and indices, as well as sums over any of the
indices. It defines an opinionated canonical form (normal ordering plus some
additional rules) to automatically simplify expressions. It is recommended to
use an interface that can display LaTeX formulas (e.g., Jupyter notebooks) for
convenient output formatting. 

Starting from v1.4, QuantumAlgebra also interoperates with computer algebra
systems (CAS) such as
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) or
[SymPy.jl](https://github.com/JuliaPy/SymPy.jl) /
[SymPyPythonCall.jl](https://github.com/jverzani/SymPyPythonCall.jl), as the
"scalar" prefactors of each quantum term can be arbitrary expressions provided
by these systems. While such expressions do not support symbolic indices in the
same way as QuantumAlgebra, they provide much more flexibility in terms of the
mathematical operations and powerful manipulation functions possible on the
parameters.

Example jupyter notebooks are available in the `examples` folder and can be
viewed online with
[nbviewer](https://nbviewer.jupyter.org/github/jfeist/QuantumAlgebra.jl/blob/main/examples/)
and tried out interactively with
[Binder](https://mybinder.org/v2/gh/jfeist/QuantumAlgebra.jl/main?filepath=examples).

## Release notes / changelog
Please see the [release notes](release_notes.md) for a summary of changes in each version.

## Overview

The basic functions to create QuantumAlgebra expressions (which are of type
`QuExpr`) are
- `a(inds...)` and `a'(inds...)` for ``a`` and ``a^{†}``, the annihilation
  and creation operators for a bosonic mode.
- `f(inds...)` and `f'(inds...)` for ``f`` and ``f^{†}``, the annihilation
  and creation operators for a fermionic mode.
- `σx(inds...)`, `σy(inds...)`, `σz(inds...)` for the Pauli matrices
  ``σ^{x,y,z}`` for a two-level system (TLS).
- `σp(inds...)`, `σm(inds...)` for excitation and deexcitation operators
  ``σ^{±}`` for a two-level system (TLS).

- **Indices**: All of these functions take an arbitrary number of indices as
  arguments, which can be either integers (1,2,...) or symbolic, where symbolic
  indices must be a single unicode character, with possibly an integer subindex:
  ```jldoctest
  julia> using QuantumAlgebra

  julia> a()
  a()

  julia> a'(:i)
  a†(i)

  julia> f'(1,2,:i_9)
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

  Note that for `@boson_ops` and `@fermion_ops`, deprecated `$(name)dag()`
  functions are defined for backward compatibility. These will be removed in a
  future version, as `$(name)'()` is now the preferred syntax for creating an
  adjoint.
  ```jldoctest
  julia> @boson_ops b
  (b (QuExpr constructor), b† (QuExpr constructor))

  julia> b'(:k)*b(:i)
  b†(k) b(i)
  ```
  Operators with different names are assumed to belong to different "species"
  and always commute. For fermions, this is not always desired, since you might
  want to use different named operators to refer to different kinds of states
  for the same species (e.g., localized and itinerant electrons). This can be
  achieved with the macro `@anticommuting_fermion_group`, which creates several
  fermionic operators that mutually anticommute:
  ```jldoctest
  julia> @anticommuting_fermion_group c d

  julia> normal_form(c()*d() + d()*c())
  0
  ```

- `param(name::Symbol,state='n',inds...)` to create a named parameter. `state` must be
  one of `'r'`, `'n'`, or `'c'` for purely real, non-conjugated complex, and
  conjugated complex parameters. More conveniently, parameters can be entered
  with string macros `Pr"name_inds..."` and `Pc"name_inds..."` for real and
  complex parameters:
  ```jldoctest
  julia> Pr"g_i,j_2,k"
  g(ij₂k)

  julia> Pr"g_i,j_2,k" == param(:g,'r',:i,:j_2,:k)
  true

  julia> Pc"α_3" == param(:α,3)
  true
  ```

- Arithmetic operations (`*`, `+`, `-`, `^`, `adjoint`=`'`) are supported
  (exponents must be nonnegative integers), with any `Number` types integrating
  automatically. Division by numbers is also supported.
  ```jldoctest
  julia> 5*a'(:k)*f(3)*σx(3)
  5 a†(k) f(3) σˣ(3)

  julia> (5//3+4im) * a'(:k)*f(3)*σx(3) + 9.4
  9.4 + (5//3+4i) a†(k) f(3) σˣ(3)

  julia> (a(:i)*f(:k))'
  f†(k) a†(i)
  ```
  If you need a bare number as a QuantumAlgebra expression, you can use
  `x*one(QuExpr)` (or `one(A)`, where `A` is any `QuExpr`).

- `∑(ind,A::QuExpr)` to represent an analytic sum over index `ind`. Since summed
  indices have no semantic meaning, the index within the expression gets
  replaced by a special numbered sum index `#ᵢ`, with `i=1,2,...`.
  ```jldoctest
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
  form. This can be enabled by calling `QuantumAlgebra.auto_normal_form(true)`.
  To make the setting permanent, call `QuantumAlgebra.auto_normal_form(true; set_preference=true)`
  or alternatively use [Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl) directly,
  i.e., call `Preferences.set_preferences!(QuantumAlgebra,"auto_normal_form"=>true/false)`.
  ```jldoctest
  julia> normal_form(a(:i)*a'(:j))
  δ(ij)  + a†(j) a(i)
  ```

- `expval(A::QuExpr)` to represent an expectation value.
  ```jldoctest
  julia> expval(a'(:j)*a(:i))
  ⟨a†(j) a(i)⟩
  ```

- `expval_as_corrs(A::QuExpr)` to represent an expectation value through its
  correlators, i.e., a cumulant expansion.
  ```jldoctest
  julia> expval_as_corrs(a'(:j)*a(:i))
  ⟨a†(j)⟩c ⟨a(i)⟩c  + ⟨a†(j) a(i)⟩c
  ```

- `comm(A::QuExpr,B::QuExpr)` to calculate the commutator ``[A,B] = AB - BA``.
  ```jldoctest
  julia> comm(a(),a'())
  -a†() a() + a() a†()

  julia> normal_form(comm(a(),a'()))
  1
  ```

- `Avac(A)` and `vacA(A)` simplify operators by assuming they are applied to the
  vacuum from the left or right, respectively. To be precise, `Avac(A)` returns
  ``A'`` such that ``A|0⟩ = A'|0⟩``, while `vacA(A)` does the same for ``⟨0|A``.
  These functions automatically apply `normal_form` to assure that the operators
  are simplified as much as possible. Note that "vacuum" for two-level systems
  is interpreted as the lower state, ``σ^{z}|0⟩ = -|0⟩``.
  ```jldoctest
  julia> Avac(a())
  0

  julia> Avac(a(:i)*a'(:j))
  δ(ij)

  julia> Avac(a()*a'()*a'())
  2 a†()

  julia> vacA(a()*a'()*a'())
  0

  julia> Avac(σx())
  σˣ()

  julia> Avac(σz())
  -1
  ```
  Both functions can also be called with an optional second argument,
  `Avac(A,modes_in_vacuum)` or `vacA(A,modes_in_vacuum)`, which is an iterable
  over operators (or a single operator) that will be assumed to be in the vacuum
  state, while all others are not. Note that the operators in `modes_in_vacuum`
  do not distinguish by index, i.e., if the modes have indices, all modes with
  the same name are assumed to be in the vacuum state. To avoid confusion, the
  `modes_in_vacuum` argument thus does not accept operators with indices.
    ```jldoctest
  julia> Avac(a(),a())
  0

  julia> Avac(a(),f())
  a()

  julia> Avac(a(:i)*a'(:j),f())
  δ(ij)  + a†(j) a(i)

  julia> Avac(a'()*a()*f()*f'(),f())
  a†() a()

  julia> @boson_ops b;

  julia> Avac(a'()*a()*b()*b'()^2*f()*f'(),(f(),b()))
  2 a†() b†() a()
  ```

- `vacExpVal(A,S=1)` calculates the vacuum expectation value ``⟨0|S^{†}AS|0⟩``,
  i.e., the expectation value ``⟨ψ|A|ψ⟩`` for the state defined by ``|ψ⟩=S|0⟩``.
  The result is guaranteed to not contain any operators.
  ```jldoctest
  julia> vacExpVal(a'()*a())
  0

  julia> vacExpVal(a'()*a(), a'()^4/sqrt(factorial(4)))
  4

  julia> vacExpVal(a'()*a(), a'()^4/sqrt(factorial(big(4))))
  4

  julia> vacExpVal(σx())
  0
  ```
  Like `vacA` and `Avac`, `vacExpVal` also takes an optional `modes_in_vacuum`
  argument, `vacExpVal(A,S,modes_in_vacuum)` (since all arguments are
  positional, `S` has to be given explicitly in this case even if it is just the
  identity operator, i.e., `vacExpVal(A,1,a())`):
  ```jldoctest
  julia> @boson_ops b;

  julia> vacExpVal(a'()*a()*b()^2*b'()^2*f()*f'(), 1, (f(),b()))
  2 a†() a()
  ```

- `heisenberg_eom(A,H,Ls=())` calculates the [Heisenberg equation of
  motion](https://en.wikipedia.org/wiki/Lindbladian#Heisenberg_picture) for
  operator `A` under the action of Hamiltonian `H` and potential Lindblad decay
  terms `Ls`, given by ``\frac{dA}{dt} = i[H,A] + ∑_i γ_i (L_i^† A L_i - \frac{1}{2} \{L_i^†
  L_i,A\})``. The Lindblad decay operators are passed as a tuple (not an array)
  of tuples, where each inner tuple describes one decay operator. The possible
  forms are `(L,)` for decay operator `L`, `(γ,L)` for decay operator `L` with
  rate `γ`, and `(inds,γ,L)` for decay operators summed over the given indices
  (note that this is different from the operator itself being a sum, seen in the
  example below). Finally, `L` can (in all three cases above) be just a single
  operator or a tuple of two operators `L=(X,Y)` to represent off-diagonal
  Lindblad terms ``L_{X,Y}[ρ] = X ρ Y^† - \frac{1}{2} \{Y^† X,ρ\}``.
  ```jldoctest
  julia> H = Pr"ω"*a'()a();

  julia> Ls = ((Pr"γ",a()),);

  julia> normal_form(heisenberg_eom(a(),H,Ls))
  -1//2 γ a() - 1i ω a()

  julia> H = QuExpr();

  julia> Ls = ((:i,a(:i)),);

  julia> normal_form(heisenberg_eom(a(:i),H,Ls))
  -1//2 a(i)

  julia> Ls = ((∑(:i,a(:i)),),);

  julia> normal_form(heisenberg_eom(a(:i),H,Ls))
  -1//2 ∑₁ a(#₁)

  julia> Ls = (((:i,:j),(a(:i),a(:j))),);

  julia> normal_form(heisenberg_eom(a(:i),H,Ls))
  -1//2 ∑₁ a(#₁)
  ```

- `heisenberg_eom_system(H,rhsfilt,Ls=(),ops=nothing)` calculates the system of
  equations of motion for the expectation values of operators appearing in `H`
  and `Ls` (same conventions as for `heisenberg_eom` above). Typically, these
  equation systems are not closed without approximations as equations for
  products of ``n`` operators involve products of ``m>n`` operators, so the
  system has to be truncated. This is achieved with a filter function that
  removes higher-order terms or rewrites them (approximately) in terms of
  lower-order expressions. The function `rhsfilt` is applied to the right-hand
  side of the equations to filter them as desired. If
  `rhsfilt(A::QuExpr)::QuExpr` is a function, it will be applied to the
  calculated right-hand side of the equations. `QuantumAlgebra` comes with two
  predefined constructors for filter functions, `droplen(maxorder::Int)`, which
  leads to all terms of order higher than `maxorder` being neglected, and
  `dropcorr(maxorder::Int)`, where all terms of order higher than `maxorder` are
  rewritten in terms of lower-order expressions up to order `maxorder` and
  higher-order correlators, with those correlations being neglected (i.e.,
  `dropcorr(1)` will replace ``⟨a^† a⟩ = ⟨a^† a⟩_c + ⟨a^†⟩ ⟨a⟩ ≈ ⟨a^†⟩ ⟨a⟩``).
  If `rhsfilt` is a number, it will be interpreted as `droplen(rhsfilt)`.
  Finally, the `ops` argument can be used to specify the operators that should
  be used to "seed" the system of equations, otherwise all operators appearing
  in `H` are used.
  ```jldoctest
  julia> H = Pr"ω"*a'()*a() + Pr"χ"*a'()*(a'()+a())*a();

  julia> Ls = ((Pr"γ",a()),);

  julia> heisenberg_eom_system(H,2,Ls,a())
  dₜ⟨a()⟩ = -1//2 γ ⟨a()⟩  - 1i ω ⟨a()⟩  - 2i χ ⟨a†() a()⟩  - 1i χ ⟨a()²⟩ 
  dₜ⟨a†() a()⟩ = -γ ⟨a†() a()⟩ 
  dₜ⟨a()²⟩ = -2i χ ⟨a()⟩  - γ ⟨a()²⟩  - 2i ω ⟨a()²⟩  
  ```
  The `heisenberg_eom_system` function can also be passed either `ExpVal` or
  `Corr` as a first argument, which will give the equations of motion of the
  expectation values (the default) or correlators (corresponding to a cumulant
  expansion) of the operators.
    ```jldoctest
  julia> H = Pr"ω"*a'()*a() + Pr"χ"*a'()*(a'()+a())*a();

  julia> Ls = ((Pr"γ",a()),);

  julia> heisenberg_eom_system(Corr,H,1,Ls,a())
  dₜ⟨a()⟩c = -1//2 γ ⟨a()⟩c  - 1i ω ⟨a()⟩c  - 2i χ ⟨a†()⟩c ⟨a()⟩c  - 1i χ ⟨a()⟩c² 
  ```

- `julia_expression(A)` to obtain a julia expression that can be used to
  automatically build codes implementing equations derived with QuantumAlgebra.
  Every expectation value or correlator is treated as a separate array. Daggers
  are represented as `ᴴ`, which are valid identifiers that can appear in the
  array names. Note that expectation values and correlators are not
  distinguished, so it is best to have all expressions use the same kind.
  ```jldoctest
  julia> julia_expression(expval_as_corrs(a'(:j)*a(:i)))
  :(aᴴ[j] * a[i] + aᴴa[j, i])
  ```
  Also note that expressions are always treated as arrays, even if they have no
  indices (which gives zero-dimensional arrays). If you are working with scalar
  quantities exclusively, it might be useful to clean up the resulting
  expression (e.g., use `MacroTools` to remove the `[]`).
  ```jldoctest
  julia> julia_expression(expval(a'()*a()*σx()))
  :(aᴴaσˣ[])
  ```

- By default, two-level system operators are represented by the Pauli
  matrices ``σ^{xyz}``, and calling `σp()` and `σm()` will give results expressed through them:
  ```jldoctest
  julia> σp()
  1//2 σˣ() + 1//2i σʸ()

  julia> σm()
  1//2 σˣ() - 1//2i σʸ()
  ```
  This can be changed by calling `QuantumAlgebra.use_σpm(true; set_preference=true/false)`
  (where the value of `set_preference` determines whether this is stored
  permanently using Preferences.jl). In this mode, ``σ^{+}`` and ``σ^{-}`` are the
  "fundamental" operators, and all expressions are written in terms of them.
  Note that mixing conventions within the same expression is not supported, so
  it is suggested to set this flag once at the beginning of any calculation.
  ```jldoctest
  julia> QuantumAlgebra.use_σpm(true)

  julia> σp()
  σ⁺()

  julia> σx()
  σ⁺() + σ⁻()

  julia> σz()
  -1 + 2 σ⁺() σ⁻()
  ```

### Preferences
Several preferences changing the behavior of QuantumAlgebra can be set
permanently (this uses [Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl)):
  - `"define_default_ops"`: if this is set to `false` (default is `true`), the
    "default" operators `a, adag, f, fdag, σx, σy, σz, σp, σm` are not defined
    upon import. Note that changing this value requires restarting the Julia
    session to take effect. The setting can be changed with
    `QuantumAlgebra.set_define_default_ops(true/false)` (which will inform you
    whether a restart is required) or with
    `Preferences.set_preferences!(QuantumAlgebra,"define_default_ops"=>true/false).`
  - `"auto_normal_form"`: Choose whether all expressions are automatically
    converted to normal form upon creation. The default is `false`. It can be
    changed for a single session with
    `QuantumAlgebra.auto_normal_form(true/false)`, and can be made permanent
    with `QuantumAlgebra.auto_normal_form(true/false; set_preference=true)` or
    with
    `Preferences.set_preferences!(QuantumAlgebra,"auto_normal_form"=>true/false)`.
    Note that this could previously be set by defining an environment variable
    `"QUANTUMALGEBRA_AUTO_NORMAL_FORM"`, but this usage has been deprecated and
    will be removed in a future version.
  - `"use_σpm"`: Choose whether for two-level systems, the "basic" operators are
    excitation/deexcitation operators ``σ^{±}`` or the Pauli matrices
    ```σ^{xyz}``. This can be changed in a single session by calling
    `QuantumAlgebra.use_σpm(true/false)`, and can be made permanent with
    `QuantumAlgebra.use_σpm(true/false; set_preference=true)` or with
    `Preferences.set_preferences!(QuantumAlgebra,"use_σpm"=>true/false)`.
  - `"quindices_type"`: Choose the underlying type for storing indices. The
    default is `"Vector"`, which allows an arbitrary number of indices but is
    slower (allocates memory). The alternative is `"NTuple{N}"` (e.g.,
    `"NTuple{5}"`), which uses a fixed-size tuple to store up to `N` indices.
    This is faster (stack-allocated) but limits the number of indices per
    operator to a maximum of `N`.
    This setting can be changed with `QuantumAlgebra.set_quindices_type("Vector"
    / "NTuple{N}")`. Note that changing this value requires restarting the Julia
    session to take effect.

## Citing

If you use QuantumAlgebra in academic work, we would appreciate a citation. See
[`CITATION.bib`](https://github.com/jfeist/QuantumAlgebra.jl/blob/main/CITATION.bib) for the relevant references.
