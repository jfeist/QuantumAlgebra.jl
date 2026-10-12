# Migrating code to the new QuantumAlgebra interface

`migrate.jl` converts code that uses the legacy interface of QuantumAlgebra (constructor functions,
e.g., `@boson_ops a` with `a(:i)`, the default operators `a`, `f`, `σx`, ..., and parameters
`Pr"g_i"`) to the interface with operator variables introduced in v1.7 (e.g., `@boson a[:]` with
`a[:i]`, and `@qparam g[:]` with `g[:i]`). It works on Julia files and Jupyter notebooks:

```sh
julia --project=tools/migration -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/migration tools/migration/migrate.jl file.jl notebook.ipynb
```

By default, the result is written to `file_migrated.jl` and `notebook_migrated.ipynb` (use
`--inplace` to overwrite the files, or `--dry-run` to only see what would be done). A file is
migrated together with the files that it includes (with `include("other.jl")`), since operators
are shared between them, and the includes in the migrated files point to the migrated versions.
The outputs of notebook cells are not changed, so rerun migrated notebooks.

Operators declared inside a function, `let` block, `@testset`, ... are only valid there and are
migrated separately. The tool reports what it did, and the code that it could not migrate
automatically, which is left unchanged:
- operators that are used with different numbers of indices (which needs separate names in the new
  interface, since the number of indices is fixed when an operator is defined),
- parameters that are used with different numbers of indices, or both as real and complex
  parameters,
- operators that are used as values (e.g., passed as functions to `map`),
- operators that are declared conditionally (e.g., inside `if`),
- parameters whose names are not valid variable names (e.g., `Pr"f(t)"`), and
- new variables (e.g., `σ` for `σp`/`σm`, or `g` for `Pr"g"`) whose names are already used for
  something else in the code.

Since both interfaces can be used together in QuantumAlgebra 1.x, the migrated code still works in
these cases. Operators that are used with different numbers of indices can be migrated by hand by
giving them different names, e.g., `@boson a b[:]` instead of using `a()` and `a(:i)`.

Note that operators without indices are now created when they are defined. For two-level systems,
whose representation depends on `QuantumAlgebra.use_σpm`, this setting therefore has to be chosen
before the operators are defined (the tool warns about this when necessary).
