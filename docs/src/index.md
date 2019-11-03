# [QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl) - quantum operator algebra in Julia

This package does quantum operator algebra (i.e., algebra with non-commuting
operators) in Julia. It defines an opinionated canonical form (normal ordering
plus some additional rules) that all expressions are automatically transformed
to, which fulfills some invariants that then allows easy use of the resulting
expressions. It is recommended to use an interface that can display LaTeX
formulas (e.g., Jupyter notebooks) for convenient output formatting.

While there is some documentation, it is not always kept fully up to date, and
it is recommended to look at the latest commit messages to get an idea about new
features etc. You can also check out the example notebooks in the `examples`
folder in the [github repository](https://github.com/jfeist/QuantumAlgebra.jl).
You can view them online directly with
[nbviewer](https://nbviewer.jupyter.org/github/jfeist/QuantumAlgebra.jl/blob/master/examples/)
and even try them out interactively with
[Binder](https://mybinder.org/v2/gh/jfeist/QuantumAlgebra.jl/master?filepath=examples).

We define an abstract type that represents an operator, and some concrete subtypes to describe various operators. We have:
- `scal(x)` representing a scalar ``x``
- `param(g,(i,j),state='n')` representing a named scalar parameter ``g_{i,j}``. `state` can be `'r'` for purely real parameters (invariant under complex conjugation), `'n'` for not-conjugated values, and `'c'` for a conjugated parameter ``g_{i,j}^{*}``.
- `a(i)` and `adag(i)` representing ``a_{i}`` and ``a_{i}^{†}``, the annihilation and creation operators for bosonic mode ``i``
- `f(i)` and `fdag(i)` representing ``f_{i}`` and ``f_{i}^{†}``, the annihilation and creation operators for fermionic mode ``i``
- `σ(a,i)` representing the Pauli matrix ``σ_{a,i}`` for two-level system (TLS) ``i``, where ``a ∈ {x=1,y=2,z=3}`` is the type of Pauli matrix.
- `OpProd(A,B)` representing ``A B``, i.e., the product of two operators
- `OpSum(A,B)` representing ``A + B``, i.e., the sum of two operators
- `OpSumAnalytical(i,A)` or `∑(i,A)` representing ``∑_{i} A``, i.e., an analytical sum over an index (assumed to run over all possible values of ``i``).
- `ExpVal(A)` representing the expectation value ``⟨A⟩``
- `Corr(AB)` representing the correlation ``⟨AB⟩_{c} = ⟨AB⟩ - ⟨A⟩⟨B⟩``, with corresponding extensions for products of more operators.

All operations are defined in such a way that the finally created object is automatically transformed to "canonical" form, which is defined by the following:
- Operator sums are expanded fully, such that the final expression is always a sum of operator products. I.e., if we write ``(A + B)(C + D)``, we get ``AC + AD + BC + BD``.
- Operator products are expressed in a well-defined "canonical" order. To achieve this canonical form, relevant commutators etc are used, so that an expression written as a single product can turn into a sum of expressions.
    1. at most one scalar prefactor (i.e., all prefactors collapsed into one)
    1. parameters ordered alphabetically (by string comparison)
    1. expectation values ``⟨A⟩``
    1. many-body correlations ``⟨AB⟩_{c}``
    1. bosonic operators in normal ordering (i.e., first creation, then annihilation operators), ordered by mode index
    1. fermionic operators in normal ordering (i.e., first creation, then annihilation operators), ordered by mode index
    1. Two-level Pauli matrices, ordered by TLS mode index. At most one Pauli matrix per TLS
- Operator sums are ordered first by number of operators (both bare and within expectation values and correlations), and then with the same priority rules.

Some other useful functions that are implemented:
- `comm(A,B)`: calculates the commutator of arbitrary operators ``[A,B] = AB - BA``. This allows, e.g., to calculate Heisenberg equations of motion for the operators.
- `ascorr(x)` takes an expression `x=A B C + D E...` and writes its expectation value in terms of single-body expectation values ``⟨A⟩``, ``⟨B⟩``, ..., and many-body correlations ``⟨AB⟩_{c}``, ``⟨ABC⟩_{c}``, etc. Currently, up to fourth-order correlations (i.e., products of four operators) are supported.
- `Avac(A)` and `vacA(A)` simplify operators by assuming they are applied to the vacuum from the left or right, respectively. To be precise, `Avac(A)` returns ``A'`` such that ``A|0⟩`` = ``A'|0⟩``, while `vacA(A)` does the same for ``⟨0|A``.
- `vacExpVal(A,S=1)` calculates the vacuum expectation value ``⟨0|S^{†}AS|0⟩``, i.e., the expectation value ``⟨ψ|A|ψ⟩`` for the state defined by ``|ψ⟩=S|0⟩``.
