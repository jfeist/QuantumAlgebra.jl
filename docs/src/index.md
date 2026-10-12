```@meta
CurrentModule = QuantumAlgebra
```
```@meta
DocTestSetup = quote
    using QuantumAlgebra
    QuantumAlgebra.auto_normal_form(false)
    QuantumAlgebra.use_σxyz()
end
```

# [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl) - quantum operator algebra in Julia

This package does quantum operator algebra (i.e., algebra with non-commuting
operators) in Julia, supporting bosonic, fermionic, and two-level system
operators, qudits (d-level systems), spins and so(n) generators, with arbitrary
names and indices, as well as sums over any of the indices. It defines an opinionated canonical form (normal ordering plus some
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

!!! note "New interface in v1.7"
    Version 1.7 introduced a new interface, where operators are variables
    defined with macros (e.g., `@boson a[:]`, used as `a[:i]`) instead of
    constructor functions (e.g., `@boson_ops a`, used as `a(:i)`). The previous
    interface still works, but is deprecated and will be removed in version 2.0.
    See [Legacy interface and migration](@ref) for how to update existing code.

### Defining operators

Operators are defined with macros, which create variables for them:

- `@boson a` defines the annihilation operator `a` of a bosonic mode, with the
  creation operator `a'`. With `@boson a[:]` (or `a[:,:]`, with any number of
  `:`), the operator has indices, given as `a[i]` (or `a[i,j]`). The number of
  indices is fixed by the declaration.
- `@fermion f` defines the annihilation operator `f` of a fermionic mode
  (with `f'` for the creation operator), and `@anticommuting_fermions c d`
  defines several fermionic operators that mutually anticommute.
- `@tlsxyz σ` defines the Pauli matrices `σx`, `σy`, `σz` of a two-level
  system (TLS), and `@tls σ` defines `σ` as the deexcitation operator ``σ^{-}``
  of a TLS, with `σ'` the excitation operator ``σ^{+}``.
- `@qudit d q`, `@spin J`/`@spin j J` and `@so n L` define qudits (``d``-level
  systems), spins, and generators of so(n), see below.
- `@qparam g` defines a real parameter `g`, and `@qparam α::Complex` a complex
  parameter `α` (with `α'` its complex conjugate). Several parameters can be
  defined at once, e.g., `@qparam ω g[:] α::Complex`.

**Indices** can be integers (1,2,...) or symbolic, where symbolic indices must
be a single unicode character, with possibly an integer subindex. They can be
given as symbols (`:i`, `:j_2`) or with index variables defined by `@qindex`,
where `@qindex j[:]` defines an index `j` with a subindex, written `j[2]`:
```jldoctest
julia> @boson a
a()

julia> a'*a
a†() a()

julia> @fermion f[:,:,:]
f[:,:,:]

julia> f'[1,2,:i_9]
f†(12i₉)

julia> @qindex i j[:];

julia> f[i,j[2],:μ_2]
f(ij₂μ₂)

julia> @tlsxyz σ[:]
(σx[:], σy[:], σz[:])

julia> σx[:◔_1]*σy[:😄_121]
σˣ(◔₁) σʸ(😄₁₂₁)

julia> @qparam g[:,:] α::Complex;

julia> g[i,1]*α'
g(i1) α*
```

Operators with different names are assumed to belong to different "species"
and always commute. For fermions, this is not always desired, since you might
want to use different named operators to refer to different kinds of states for
the same species (e.g., localized and itinerant electrons). This can be
achieved with `@anticommuting_fermions`:
```jldoctest
julia> @anticommuting_fermions c d
(c(), d())

julia> normal_form(c*d + d*c)
0
```

Every macro has a function with the same name that returns the operators instead
of assigning them to variables, e.g., `a = boson(:a; ndims=2)` is equivalent to
`@boson a[:,:]` (see [`boson`](@ref), [`fermion`](@ref),
[`anticommuting_fermions`](@ref), [`tls`](@ref), [`tlsxyz`](@ref),
[`qudit`](@ref), [`spin`](@ref), [`so`](@ref), [`qparam`](@ref) and
[`qindex`](@ref)). This is
useful to create operators programmatically, e.g., with names that are generated
in code.

Expressions are currently printed with parentheses for the indices (e.g.,
`f(ij₂μ₂)` and `a()` above). Calling `QuantumAlgebra.use_bracket_display(true)`
prints them in the same syntax as they are entered (e.g., `f[i,j₂,μ₂]` and `a`),
which will become the default in version 2.0 (see [Preferences](@ref)).

For backward compatibility, QuantumAlgebra currently also defines the default
operators `a`, `f`, `σx`, `σy`, `σz`, `σp`, `σm` of the legacy interface when it
is loaded. Defining operators with the same names (as above) replaces them. On
Julia versions before 1.12, this is only possible if the default operator has
not been used yet; the default operators can be disabled with the preference
`"define_default_ops"` (see [Preferences](@ref)).

### Qudits, spins and so(n)

- `@qudit d q` (or `q = qudit(:q, d)`) defines the operators of a `d`-level
  system (qudit), where `q(i,j)` is the transition operator ``|i⟩⟨j|``
  and level 1 is the ground state (the vacuum). With `@qudit d q[:]` (with any
  number of `:`), the qudit has indices, given as `q(i,j)[k]`. Internally, all
  operators are expressed through ``|i⟩⟨1|`` and ``|1⟩⟨i|`` (with ``i ≥ 2``), which
  for `d = 2` are equivalent to ``σ^{+}`` and ``σ^{-}``. `gellmann(q)` returns
  the generalized Gell-Mann matrices of the qudit, normalized as
  ``\mathrm{Tr}(λ_a λ_b) = 2δ_{ab}`` (where the normalization of the diagonal
  matrices uses `Float64` square roots by default, see the docstring for how to keep
  them exact).
  ```jldoctest
  julia> @qudit 3 q[:]
  q[:] (3-level qudit)

  julia> q(2,3)[:i]
  q²³(i)

  julia> normal_form(q(1,2)[:i]*q(2,3)[:i])
  q¹³(i)

  julia> @qudit 3 r;

  julia> r(1,1)
  1 - r²²() - r³³()
  ```

- `@spin J` (or `Jp, Jm, Jz = spin(:J)`) defines the operators `Jp`, `Jm`, `Jz`,
  i.e., ``J_+``, ``J_-``, ``J_z``, of a spin (the Lie algebra su(2) ≅ so(3)), with
  commutation relations ``[J_z,J_±] = ±J_±`` and ``[J_+,J_-] = 2J_z``. Without
  fixed spin `j`, only these are used, so that results hold in any representation
  (e.g., for collective spins). With fixed spin `j` (`@spin j J`, where `j` can
  also be symbolic), products ``J_+ J_-`` on the same site are replaced using the
  Casimir operator, ``J_+ J_- = j(j+1) - J_z^2 + J_z``, and the vacuum is the
  lowest-weight state.
  ```jldoctest
  julia> @spin J;

  julia> normal_form(comm(Jp, Jm))
  2 Jᶻ()

  julia> @spin 1//2 S;

  julia> normal_form(Sp*Sm)
  3//4 + Sᶻ() - Sᶻ()²
  ```

- `@so n L` (or `L = so(:L, n)`) defines the generators ``L_{ab} = -L_{ba}``
  of the Lie algebra so(n) of rotations in ``n`` dimensions, given as `L(a,b)`
  (or `L(a,b)[k]` with `@so n L[:]`). The generators are Hermitian, with
  commutation relations
  ``[L_{ab}, L_{cd}] = i(δ_{ac} L_{bd} - δ_{ad} L_{bc} - δ_{bc} L_{ad} + δ_{bd} L_{ac})``
  (e.g., ``L_{ab} = x_a p_b - x_b p_a``). As for spins without fixed `j`, only
  these relations are used, so that results hold in any representation.
  ```jldoctest
  julia> @so 4 L;

  julia> L(2,1)
  -L¹²()

  julia> normal_form(comm(L(1,2), L(2,3)))
  -1i L¹³()
  ```

### Working with expressions

- Arithmetic operations (`*`, `+`, `-`, `^`, `adjoint`=`'`) are supported
  (exponents must be nonnegative integers), with any `Number` types integrating
  automatically. Division by numbers is also supported.
  ```jldoctest
  julia> @boson a[:]; @fermion f[:]; @tlsxyz σ[:];

  julia> 5*a'[:k]*f[3]*σx[3]
  5 a†(k) f(3) σˣ(3)

  julia> (5//3+4im) * a'[:k]*f[3]*σx[3] + 9.4
  9.4 + (5//3+4i) a†(k) f(3) σˣ(3)

  julia> (a[:i]*f[:k])'
  f†(k) a†(i)
  ```
  If you need a bare number as a QuantumAlgebra expression, you can use
  `x*one(QuExpr)` (or `one(A)`, where `A` is any `QuExpr`).

- `∑(ind,A::QuExpr)` to represent an analytic sum over index `ind`. Since summed
  indices have no semantic meaning, the index within the expression gets
  replaced by a special numbered sum index `#ᵢ`, with `i=1,2,...`.
  ```jldoctest
  julia> @boson a[:]; @qindex i;

  julia> ∑(i,a[i])
  ∑₁ a(#₁)
  ```

- `normal_form(A::QuExpr)` converts an expression to a well-defined "canonical"
  order. To achieve this canonical form, relevant commutators etc are used, so
  an expression written as a single product can turn into a sum of expressions.
  The order is essentially normal ordering (creation before annihilation
  operators, with σˣʸᶻ in the middle), with some additional conventions to make
  the normal form (hopefully) unique. In particular, sum indices are numbered
  canonically, so that terms that only differ in the naming of summed indices
  are identified (e.g., `normal_form(∑(:i,∑(:j,a[:i]*a'[:j])) - ∑(:j,∑(:i,a[:j]*a'[:i])))`
  is zero). In some contexts (e.g., interactive work),
  it can be convenient to automatically transform all expressions to normal
  form. This can be enabled by calling `QuantumAlgebra.auto_normal_form(true)`.
  To make the setting permanent, call `QuantumAlgebra.auto_normal_form(true; set_preference=true)`
  or alternatively use [Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl) directly,
  i.e., call `Preferences.set_preferences!(QuantumAlgebra,"auto_normal_form"=>true/false)`.
  ```jldoctest
  julia> @boson a[:]
  a[:]

  julia> normal_form(a[:i]*a'[:j])
  δ(ij) + a†(j) a(i)
  ```

- `expval(A::QuExpr)` to represent an expectation value.
  ```jldoctest
  julia> @boson a[:];

  julia> expval(a'[:j]*a[:i])
  ⟨a†(j) a(i)⟩
  ```

- `expval_as_corrs(A::QuExpr)` to represent an expectation value through its
  correlators, i.e., a cumulant expansion.
  ```jldoctest
  julia> @boson a[:];

  julia> expval_as_corrs(a'[:j]*a[:i])
  ⟨a†(j)⟩c ⟨a(i)⟩c + ⟨a†(j) a(i)⟩c
  ```

- `comm(A::QuExpr,B::QuExpr)` to calculate the commutator ``[A,B] = AB - BA``
  (in normal form). The commutator is computed directly from the commutators of
  pairs of operators, so that terms that would cancel in ``AB - BA`` are never
  generated.
  ```jldoctest
  julia> @boson a; @boson b[:];

  julia> comm(a,a')
  1

  julia> comm(b'[:i]*b[:j], b'[:k]*b[:k])
  -δ(ik) b†(i) b(j) + δ(jk) b†(i) b(j)
  ```

- `Avac(A)` and `vacA(A)` simplify operators by assuming they are applied to the
  vacuum from the left or right, respectively. To be precise, `Avac(A)` returns
  ``A'`` such that ``A|0⟩ = A'|0⟩``, while `vacA(A)` does the same for ``⟨0|A``.
  These functions automatically apply `normal_form` to assure that the operators
  are simplified as much as possible. Note that "vacuum" for two-level systems
  is interpreted as the lower state, ``σ^{z}|0⟩ = -|0⟩``.
  ```jldoctest
  julia> @boson a; @boson b[:]; @tlsxyz σ;

  julia> Avac(a)
  0

  julia> Avac(b[:i]*b'[:j])
  δ(ij)

  julia> Avac(a*a'*a')
  2 a†()

  julia> vacA(a*a'*a')
  0

  julia> Avac(σx)
  σˣ()

  julia> Avac(σz)
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
  julia> @boson a; @boson b; @fermion f;

  julia> Avac(a,a)
  0

  julia> Avac(a,f)
  a()

  julia> Avac(a'*a*f*f',f)
  a†() a()

  julia> Avac(a'*a*b*b'^2*f*f',(f,b))
  2 a†() b†() a()
  ```

- `vacExpVal(A,S=1)` calculates the vacuum expectation value ``⟨0|S^{†}AS|0⟩``,
  i.e., the expectation value ``⟨ψ|A|ψ⟩`` for the state defined by ``|ψ⟩=S|0⟩``.
  The result is guaranteed to not contain any operators.
  ```jldoctest
  julia> @boson a; @tlsxyz σ;

  julia> vacExpVal(a'*a)
  0

  julia> vacExpVal(a'*a, a'^4/sqrt(factorial(4)))
  4

  julia> vacExpVal(a'*a, a'^4/sqrt(factorial(big(4))))
  4

  julia> vacExpVal(σx)
  0
  ```
  Like `vacA` and `Avac`, `vacExpVal` also takes an optional `modes_in_vacuum`
  argument, `vacExpVal(A,S,modes_in_vacuum)` (since all arguments are
  positional, `S` has to be given explicitly in this case even if it is just the
  identity operator, i.e., `vacExpVal(A,1,a)`):
  ```jldoctest
  julia> @boson a; @boson b; @fermion f;

  julia> vacExpVal(a'*a*b^2*b'^2*f*f', 1, (f,b))
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
  julia> @boson a; @boson b[:]; @qparam ω γ;

  julia> H = ω*a'*a;

  julia> Ls = ((γ,a),);

  julia> normal_form(heisenberg_eom(a,H,Ls))
  -1//2 γ a() - 1i ω a()

  julia> H = QuExpr();

  julia> Ls = ((:i,b[:i]),);

  julia> normal_form(heisenberg_eom(b[:i],H,Ls))
  -1//2 b(i)

  julia> Ls = ((∑(:i,b[:i]),),);

  julia> normal_form(heisenberg_eom(b[:i],H,Ls))
  -1//2 ∑₁ b(#₁)

  julia> Ls = (((:i,:j),(b[:i],b[:j])),);

  julia> normal_form(heisenberg_eom(b[:i],H,Ls))
  -1//2 ∑₁ b(#₁)
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
  julia> @boson a; @qparam ω χ γ;

  julia> H = ω*a'*a + χ*a'*(a'+a)*a;

  julia> Ls = ((γ,a),);

  julia> heisenberg_eom_system(H,2,Ls,a)
  dₜ⟨a()⟩ = -1//2 γ ⟨a()⟩ - 1i ω ⟨a()⟩ - 2i χ ⟨a†() a()⟩ - 1i χ ⟨a()²⟩
  dₜ⟨a†() a()⟩ = -γ ⟨a†() a()⟩
  dₜ⟨a()²⟩ = -2i χ ⟨a()⟩ - γ ⟨a()²⟩ - 2i ω ⟨a()²⟩  
  ```
  The `heisenberg_eom_system` function can also be passed either `ExpVal` or
  `Corr` as a first argument, which will give the equations of motion of the
  expectation values (the default) or correlators (corresponding to a cumulant
  expansion) of the operators.
  ```jldoctest
  julia> @boson a; @qparam ω χ γ;

  julia> H = ω*a'*a + χ*a'*(a'+a)*a;

  julia> Ls = ((γ,a),);

  julia> heisenberg_eom_system(Corr,H,1,Ls,a)
  dₜ⟨a()⟩c = -1//2 γ ⟨a()⟩c - 1i ω ⟨a()⟩c - 2i χ ⟨a†()⟩c ⟨a()⟩c - 1i χ ⟨a()⟩c² 
  ```

- `julia_expression(A)` to obtain a julia expression that can be used to
  automatically build codes implementing equations derived with QuantumAlgebra.
  Every expectation value or correlator is treated as a separate array. Daggers
  are represented as `ᴴ`, which are valid identifiers that can appear in the
  array names. Note that expectation values and correlators are not
  distinguished, so it is best to have all expressions use the same kind.
  ```jldoctest
  julia> @boson a[:];

  julia> julia_expression(expval_as_corrs(a'[:j]*a[:i]))
  :(aᴴ[j] * a[i] + aᴴa[j, i])
  ```
  Also note that expressions are always treated as arrays, even if they have no
  indices (which gives zero-dimensional arrays). If you are working with scalar
  quantities exclusively, it might be useful to clean up the resulting
  expression (e.g., use `MacroTools` to remove the `[]`).
  ```jldoctest
  julia> @boson a; @tlsxyz σ;

  julia> julia_expression(expval(a'*a*σx))
  :(aᴴaσˣ[])
  ```

- By default, two-level system operators are represented by the Pauli
  matrices ``σ^{xyz}``, and the excitation and deexcitation operators defined
  by `@tls` are expressed through them:
  ```jldoctest
  julia> @tls σ
  1//2 σˣ() - 1//2i σʸ()

  julia> σ'
  1//2 σˣ() + 1//2i σʸ()
  ```
  This can be changed by calling `QuantumAlgebra.use_σpm(true; set_preference=true/false)`
  (where the value of `set_preference` determines whether this is stored
  permanently using Preferences.jl). In this mode, ``σ^{+}`` and ``σ^{-}`` are the
  "fundamental" operators, and all expressions are written in terms of them.
  Note that mixing conventions within the same expression is not supported, so
  it is suggested to set this flag once at the beginning of any calculation,
  before defining the operators (operators without indices are created when they
  are defined).
  ```jldoctest
  julia> QuantumAlgebra.use_σpm(true)

  julia> @tls σ; @tlsxyz σ;

  julia> σ'
  σ⁺()

  julia> σx
  σ⁺() + σ⁻()

  julia> σz
  -1 + 2 σ⁺() σ⁻()
  ```

### Preferences
Several preferences changing the behavior of QuantumAlgebra can be set
permanently (this uses [Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl)):
  - `"define_default_ops"`: if this is set to `false` (default is `true`), the
    "default" operators `a, adag, f, fdag, σx, σy, σz, σp, σm` of the [legacy
    interface](@ref "Legacy interface and migration") are not defined upon
    import (they will be removed in version 2.0). Note that changing this value requires restarting the Julia
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
  - `"bracket_display"`: Choose whether operators and parameters are printed in
    the syntax of the interface with operator variables, e.g., `a†[i,j]`, `g[i]`
    and `a` (for `true`), or with parentheses, e.g., `a†(ij)`, `g(i)` and `a()`
    (for `false`, the default). This can be changed in a single session by calling
    `QuantumAlgebra.use_bracket_display(true/false)`, and can be made permanent
    with `QuantumAlgebra.use_bracket_display(true/false; set_preference=true)`.
    Printing with brackets will become the default in version 2.0.
  - `"quindices_type"`: Deprecated and ignored. Terms are now always stored in a
    compact form that supports an arbitrary number of indices per operator.

## Citing

If you use QuantumAlgebra in academic work, we would appreciate a citation. See
[`CITATION.bib`](https://github.com/jfeist/QuantumAlgebra.jl/blob/main/CITATION.bib) for the relevant references.
