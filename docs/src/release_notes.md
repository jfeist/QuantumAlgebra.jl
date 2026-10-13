# Release notes

## (unreleased)
- Operators of qudits with a given matrix: `q(M)` is ``\sum_{ij} M_{ij} |i⟩⟨j|``
  for a ``d×d`` matrix `M` (with numbers or parameters as entries), e.g., for a
  Hamiltonian or for the generators of any Lie algebra in a ``d``-dimensional
  representation.
- Pauli matrices of a two-level system with an empty name are printed as `X`,
  `Y` and `Z`, as in quantum computing (e.g., `X, Y, Z = tlsxyz(""; ndims=1)`,
  or `@tlsxyz P = ""` for variables `Px`, `Py`, `Pz`).

## v1.7.1 (2026-10-13)
A patch release that makes normal forms unique in two more cases (Pauli
matrices with symbolic indices, and spins with a fixed numerical `j`) and speeds
up the handling of coefficients. Results for such Pauli matrices and spins can
change form, since they are now simplified further.

- Faster handling of the coefficients of terms: calculating equations of motion,
  correlators and vacuum expectation values is up to 30% faster.
- The normal form of products of Pauli matrices (``σ^{x,y,z}``) with symbolic
  indices is now unique: two Pauli matrices on the same site are now always
  multiplied, also if there are Pauli matrices between them that might be on the
  same site, e.g., ``σ^x_i σ^y_j σ^z_i = -i σ^y_i σ^y_j + 2i δ_{ij}``.
  Previously, such products stayed as they were, so that equal expressions could
  have different normal forms (and did not cancel).
- Spins with a fixed numerical `j` now use all relations of the ``(2j+1)``-dimensional
  representation, ``J_+^a ∏_{m=-j}^{j-a} (J_z - m) = 0`` and
  ``J_-^b ∏_{m=b-j}^{j} (J_z - m) = 0`` (e.g., ``J_+^2 = 0`` and ``J_z^2 = 1/4``
  for ``j = 1/2``, or ``J_z^3 = J_z`` for ``j = 1``). Previously, only the Casimir
  operator was used, so that equal operators could have different normal forms,
  and powers of ``J_z`` were not reduced. For a symbolic `j`, only the Casimir
  operator is used, as before.
- Qudits, spins and so(n) generators (new in version 1.7) are now marked as
  experimental in the documentation: the rules that are used for them (and
  therefore the form of results) can still change.

## v1.7.0 (2026-10-12)
This is a minor release because it has no breaking changes (existing code keeps
working, although some results are presented differently, see below), but it
contains many new features and prepares the transition to version 2.0. Its main
changes are a new interface for defining operators (the previous interface is
deprecated and will be removed in version 2.0), new kinds of operators (qudits,
spins and so(n) generators), large performance improvements, and support for
using QuantumAlgebra from several threads at the same time. The minimum supported
Julia version is now 1.10 (the current LTS).

### New interface for defining operators
- Operators are now variables that are defined with macros, e.g., `@boson a`
  (used as `a` and `a'`) or `@boson a[:,:]` (used as `a[i,j]` and `a'[i,j]`),
  instead of constructor functions (`@boson_ops a`, used as `a()` and
  `a(:i,:j)`). The number of indices is fixed when an operator is defined. The
  macros are `@boson`, `@fermion`, `@anticommuting_fermions`, `@tls` (for
  ``σ^{-}``, with ``σ^{+}`` given by `σ'`), `@tlsxyz`, and `@qudit`, `@spin`,
  `@so` for the new kinds of operators (see below). Parameters are defined with
  `@qparam`, e.g., `@qparam ω g[:] α::Complex`, and index variables with
  `@qindex`, e.g., `@qindex i j[:]` for `a[i,j[2]]` (indices can also still be
  given as symbols, e.g., `a[:i,:j_2]`).
- Each macro has a function with the same name to create operators
  programmatically, e.g., `boson(:a; ndims=2)` or `qparam(:α, Complex)`.
  Several operators can be defined at once, e.g., `@boson a b[:]` or
  `@qudit 3 s p[:]`.
- The name of an operator or parameter can differ from the name of its
  variable, e.g., `@fermion cu="c↑"` or `@qparam ft="f(t)"` (for names that are
  not valid variable names), and the functions also accept names as strings,
  e.g., `qparam("f(t)")`.
- The modes that are assumed to be in the vacuum state in `Avac`, `vacA` and
  `vacExpVal` are given by the operator variables, where, e.g., `a` after
  `@boson a[:]` stands for all modes of `a`.
- `QuantumAlgebra.use_bracket_display(true)` (or the preference
  `"bracket_display"`) prints expressions in the syntax of the new interface,
  e.g., `a†[i,j]` and `a` instead of `a†(ij)` and `a()`. This will become the
  default in version 2.0.
- The previous interface (constructor functions, `@boson_ops` etc., the default
  operators `a`, `f`, `σx`, …, and parameters created with `param`, `Pr"…"` and
  `Pc"…"`) still works and can be mixed with the new one, but it is deprecated
  and will be removed in version 2.0. The documentation describes how to update
  existing code, and the script `tools/migration/migrate.jl` in the repository
  converts Julia files and Jupyter notebooks automatically where possible (it
  reports the code that has to be converted by hand, e.g., operators that are
  used with different numbers of indices).

### New kinds of operators
- Qudits (d-level systems): `@qudit d q` defines the transition operators
  `q(i,j)` = ``|i⟩⟨j|`` (or `q(i,j)[k]` for qudits with indices), which are
  expressed through ``|i⟩⟨1|`` and ``|1⟩⟨i|`` to obtain a unique normal form
  (products ``|i⟩⟨1|1⟩⟨j|`` are printed as ``|i⟩⟨j|``, e.g., `q²³(k)`).
  `gellmann(q)` returns the generalized Gell-Mann matrices, normalized as
  ``\mathrm{Tr}(λ_a λ_b) = 2δ_{ab}`` (the normalization factors of the diagonal
  matrices are `Float64` by default, and can be kept exact by passing a symbolic
  square root, e.g., `gellmann(q; sqrt = x -> sqrt(SymPyPythonCall.Sym(x)))`).
- Spins: `@spin J` defines `Jp`, `Jm`, `Jz` with the commutation relations of
  su(2), which hold in any representation (e.g., for collective spins). With a
  fixed spin `j` (`@spin j J`, where `j` can also be symbolic), products
  ``J_+ J_-`` are simplified with the Casimir operator, and the vacuum is the
  lowest-weight state.
- so(n): `@so n L` defines the generators `L(a,b)` = ``L_{ab} = -L_{ba}`` of
  so(n), using only their commutation relations.
- Pauli matrices and jump operators (``σ^{±}``) can now be mixed in one
  expression if they act on different two-level systems.

### Changes that affect results
Results are mathematically unchanged, but can look different:
- Summation indices are numbered canonically, so that terms that only differ in
  the naming of summation indices are combined, and sums that vanish by
  symmetry are recognized (e.g., ``∑_{ij} f_i f_j = 0``). This can reduce the
  number of terms considerably, e.g., by a factor 16 for nested commutators
  with sums in a benchmark.
- `comm(A,B)` now returns the commutator in normal form (previously, it returned
  `A*B - B*A` without simplification). It is computed directly from the
  commutators of pairs of operators, so that the terms of ``AB`` and ``BA`` that
  cancel are never generated. `heisenberg_eom` computes the Lindblad terms as
  ``L^† [A,L] + \frac{1}{2} [L^† L, A]`` for the same reason.
- Terms are printed without trailing spaces, and with single spaces between
  terms (e.g., `δ(ij) + a†(j) a(i)` instead of `δ(ij)  + a†(j) a(i)`).
- Names with unicode characters are converted to LaTeX correctly, also with
  combining accents (e.g., `b̄` gives `\bar{b}`), subscripts and superscripts
  (`ω₁`, `α²`), and script or double-struck letters (`ℰ`, `ℝ`). Names that are
  not a single character are enclosed in braces, so that names with subscripts
  (e.g., `ω_c`) can have indices.

### Performance and other improvements
- Terms are now stored compactly as a flat vector of 64-bit words, which
  reduces run time and memory allocations considerably. Together with the
  changes above, compared to v1.6.1, normal ordering of large expressions is
  1.2–1.4 times faster, cumulant expansions (`expval_as_corrs`) 2–4 times,
  commutators 3 times, equations of motion with Lindblad terms 3–4.4 times,
  correlator equations (`corrheis`) 3 times, and nested commutators with sums up
  to 15 times faster in the benchmarks in `bench/runbench.jl`.
- Operators can have an arbitrary number of indices (up to 255) without any
  performance penalty, so the `"quindices_type"` preference introduced in v1.6.0
  has no effect anymore, and `set_quindices_type` is deprecated.
- QuantumAlgebra can be used from several threads at the same time (for
  different expressions).
- PackageExtensionCompat is no longer a dependency (it was only needed for
  Julia versions before 1.9).
- New exported names: the macros listed above, and the functions `boson`,
  `fermion`, `anticommuting_fermions`, `tls`, `tlsxyz`, `qudit`, `gellmann`,
  `spin`, `so`, `qparam` and `qindex`. If these clash with names in your code,
  the variables in your code take precedence (as long as the names have not
  been used before they are defined).

### Bug fixes
- `A ≈ B` threw an error instead of returning `false` when `A` contains a term
  that is not in `B` and whose coefficient cannot be compared with numbers (e.g.,
  a `Symbolics` variable).

## v1.6.1 (2026-10-11)
This is a patch release with bug fixes:
- Fix `expval_as_corrs` and `corr_as_expvals` (and thus `heisenberg_eom_system` for
  correlators) for sums whose operators simplify when the summation index is
  equal to another index, e.g., for products of Pauli matrices. Previously, the
  prefactors of such simplifications were dropped, products that vanish were not
  taken into account, and in some cases an assertion error was thrown. This
  changes (corrects) results, e.g., of the equations of motion for correlators in
  systems with Pauli matrices and sums.
- Fix `normal_form` for products of Pauli matrices where operators that are not
  next to each other were contracted although an operator in between does not
  commute with them (e.g., `σx(:j)*σy(:i)*sz(1)*σz(:j)` lost a term proportional
  to `δ(ij)`).
- Make the result of `heisenberg_eom_system` deterministic. Which of two
  equivalent (mutually conjugate) expectation values or correlators was chosen
  as a variable could depend on the order in which operators were defined and on
  the package version.
- Remove a method of `Base.hash` for vectors of QuantumAlgebra objects (see
  [#32](https://github.com/jfeist/QuantumAlgebra.jl/issues/32)).
- Allow OrderedCollections v2.

## v1.6.0 (2025-12-05)
This is a minor release with internal improvements and new features:
- The new preference`"quindices_type"` can be used to choose the underlying
  type for storing indices. The default is `"Vector"`, which allows an arbitrary
  number of indices but is slower (allocates memory). The alternative is 
  `"NTuple{N}"` (e.g., `"NTuple{5}"`), which uses a fixed-size tuple to store
  up to `N` indices. This is faster (stack-allocated) but limits the number of
  indices per operator to a maximum of `N`. This setting can be changed with
  `QuantumAlgebra.set_quindices_type("Vector" / "NTuple{N}")`. Note that changing
  this value requires restarting the Julia session to take effect. For normal
  usage, `"Vector"` is recommended -- the `"NTuple{N}"` option is mainly intended
  for specific projects where extremely large expressions (100.000+ terms) are
  treated and speed becomes critical. Note that the setting can be changed
  locally within each project using `QuantumAlgebra`.
- Extend documentation and especially docstrings.
- Fix a bug in the prefactors of expressions that are simplified through several contractions in
  `normal_form` (see [#23](https://github.com/jfeist/QuantumAlgebra.jl/issues/23)).

## v1.5.1 (2024-11-20)
This is a patch release with some fixes to improve compatibility with
[Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl).

## v1.5.0 (2024-07-19)
This is a minor release with some bug fixes and new features:
- Add `heisenberg_eom_system` function to generate a system of Heisenberg
  equations of motion for either the expectation values or the cumulants /
  correlators of operator products, starting from a Hamiltonian and Lindblad
  operators describing decoherence. Typically, these equation systems are not
  closed without approximations as equations for products of n operators involve
  products of m>n operators, so the system has to be truncated. This is achieved
  with a filter function that removes higher-order terms or rewrites them
  (approximately) in terms of lower-order expressions. Simple example:
  ```julia
  julia> using QuantumAlgebra
  julia> @boson_ops a;
  julia> H = Pr"ω"*a'()*a() + Pr"χ"*a'()*(a'()+a())*a();
  julia> Ls = ((Pr"γ",a()),);
  julia> EQ = heisenberg_eom_system(H,2,Ls,a())
  dₜ⟨a()⟩ = -1//2 γ ⟨a()⟩  - 1i ω ⟨a()⟩  - 2i χ ⟨a†() a()⟩  - 1i χ ⟨a()²⟩ 
  dₜ⟨a†() a()⟩ = -γ ⟨a†() a()⟩ 
  dₜ⟨a()²⟩ = -2i χ ⟨a()⟩  - γ ⟨a()²⟩  - 2i ω ⟨a()²⟩  
  ```
- Add additional argument `modes_in_vacuum` to `vacA`, `Avac`, and `vacExpVal`,
  to restrict the modes that are considered to be in the vacuum state, with all
  other operators unaffected.
- Fix a LaTeX output error where Latexify misinterpreted a superscript "^+" as a
  sum and replaced "x^+ - y" with "x^- y".

## v1.4.0 (2024-06-29)
This is a minor release with some bug fixes and new features:
- Allow `Number` type inputs for most exported functions so that they do not
  require `QuExpr` wrapping, such that `expval(3)` or `Avac(1)` no longer raise
  errors (closes [#14](https://github.com/jfeist/QuantumAlgebra.jl/issues/14)).
- Support interoperability with symbolic computer algebra systems (CAS) such as
  [Symbolics.jl](https://github.com/JuliaSymbolics/Symbolics.jl) or
  [SymPy.jl](https://github.com/JuliaPy/SymPy.jl) /
  [SymPyPythonCall.jl](https://github.com/jverzani/SymPyPythonCall.jl).
  Expressions provided by these systems can be used as the "scalar" factors of
  each term in a `QuExpr`, and thus provide an alternative to the
  QuantumAlgebra `Param` objects. In contrast to `Param`s, they do not support
  symbolic indices, but provide much more flexibility in terms of the
  mathematical operations and powerful manipulation functions. For example, one
  can use
  ```julia
  julia> using QuantumAlgebra, Symbolics
  julia> @boson_ops a;
  julia> @variables x y;
  julia> H = cos(x)^2 * a() * a'() + sin(x)^2 * a'() * a();
  julia> Symbolics.simplify(normal_form(H))
  cos(x)^2 + a†() a()

  julia> Symbolics.substitute(H, x => y)
  sin(y)^2 a†() a() + cos(y)^2 a() a†()

  julia> Symbolics.substitute(H, x => 2)
  0.826821810431806 a†() a() + 0.17317818956819406 a() a†()

  julia> Symbolics.substitute(H, x => 2; fold=false)
  sin(2)^2 a†() a() + cos(2)^2 a() a†()
  ```

  While overloads for some functions (like `simplify` and `substitute` from
  `Symbolics`) are available, all other functions/manipulations from the CAS
  packages can be applied with the new function `map_scalar_function` which
  applies a function to each scalar factor in a `QuExpr`. This can be used as,
  e.g.,
  ```julia
  julia> map_scalar_function(Symbolics.simplify, normal_form(H))
  cos(x)^2 + a†() a()
  ```

  In order not to add dependencies on several CAS packages to QuantumAlgebra,
  the interoperability helpers are defined in extension modules, using
  [PackageExtensionCompat.jl](https://github.com/cjdoris/PackageExtensionCompat.jl)
  to maintain compatibility with Julia versions before v1.9.

## v1.3.1 (2023-12-22)
This is a patch release with some bug fixes and performance improvements.

## v1.3.0 (2023-09-19)
This is a minor release with some bug fixes in LaTeX output, added unit tests for many functions, and some new features and deprecations:
- The functions that generate operators now can be conjugated directly, e.g., `a'()` is equivalent to `adag()` (and to `a()'`, which was already the case beforehand). The old `_dag` functions for bosons and fermions are deprecated and will be removed in a future version.
- [Preferences.jl](https://github.com/JuliaPackaging/Preferences.jl) is now used to set and store preferences. The following preferences are available:
  - `"define_default_ops"`: if this is set to `false` (default is `true`), the
    "default" operators `a, adag, f, fdag, σx, σy, σz, σp, σm` are not defined
    upon import. Note that changing this value requires restarting the Julia
    session to take effect. The setting can be changed with
    `QuantumAlgebra.set_define_default_ops(true/false)` (which will inform you whether a
    restart is required) or with
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
    ``σ^{xyz}``. This can be changed in a single session by calling
    `QuantumAlgebra.use_σpm(true/false)`, and can be made permanent with
    `QuantumAlgebra.use_σpm(true/false; set_preference=true)` or with
    `Preferences.set_preferences!(QuantumAlgebra,"use_σpm"=>true/false)`.
- Use PrecompileTools.jl instead of SnoopPrecompile.jl for precompilation.


## v1.2.0 (2023-04-22)
This is a minor revision with some bug fixes, and some new features:
- Add the `heisenberg_eom` function to calculate the Heisenberg equations of motions of operators.
- Switch the LaTeX output to use [Latexify.jl](https://github.com/korsbo/Latexify.jl), which also enables support for display of vectors, tuples, etc of `QuExpr`.
- Require at least Julia v1.6.

## v1.1.0 (2021-09-02)
This is a minor revision with some and bug fixes, and one new feature: It is now
possible to use `@anticommuting_fermion_group` to define several fermionic
operators that mutually anticommute, allowing to refer to different kinds of
states for the same "species" (e.g., localized and itinerant electrons):
```julia
julia> @anticommuting_fermion_group c d

julia> normal_form(c()*d() + d()*c())
0
```

## v1.0.0 (2021-07-23)
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
  correlator / cumulant expansion, e.g., ``⟨AB⟩ = ⟨AB⟩_{c}`` - ``⟨A⟩_{c}
  ⟨B⟩_{c}``, with corresponding extensions for products of more operators. Note
  that for a single operator, ``⟨A⟩_{c} = ⟨A⟩``, but we distinguish the two
  formally for clarity.
- There is a new function `julia_expression(A)` that converts a QuantumAlgebra
  object to a julia expression, which helps in using QuantumAlgebra to
  programatically derive codes for numerical implementation. The object `A`
  cannot contain any "bare" operators, but only expectation values or
  correlators. See the documentation for more details.
- QuantumAlgebra expressions are now printed in pretty format in the terminal
  etc.

