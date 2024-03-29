# Release notes

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

