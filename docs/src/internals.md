# Internals

These functions and types are not exported (they are accessible as `QuantumAlgebra.name`). They
mostly form the internal interface of the package and should not be relevant to most users, apart
from the settings (e.g., [`use_σpm`](@ref QuantumAlgebra.use_σpm), [`auto_normal_form`](@ref QuantumAlgebra.auto_normal_form) and
[`use_bracket_display`](@ref QuantumAlgebra.use_bracket_display)) and some functions to inspect or modify expressions (e.g.,
[`replace_inds`](@ref QuantumAlgebra.replace_inds) and [`is_normal_form`](@ref QuantumAlgebra.is_normal_form)).

```@index
Pages = ["internals.md"]
```

```@autodocs
Modules = [QuantumAlgebra]
Public = false
```
