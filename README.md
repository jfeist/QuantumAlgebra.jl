# QuantumAlgebra.jl - quantum operator algebra in Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jfeist.github.io/QuantumAlgebra.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jfeist.github.io/QuantumAlgebra.jl/dev)
[![Build Status](https://travis-ci.com/jfeist/QuantumAlgebra.jl.svg?branch=master)](https://travis-ci.com/jfeist/QuantumAlgebra.jl)
[![Codecov](https://codecov.io/gh/jfeist/QuantumAlgebra.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jfeist/QuantumAlgebra.jl)

This package implements a package to do quantum operator algebra (i.e., algebra
with non-commuting operators) in Julia. It defines an opinionated canonical form
(normal ordering plus some additional rules) that all expressions are
automatically transformed to. The resulting expressions are displayed as LaTeX
formulas in interfaces that support them (e.g., Jupyter notebooks).

We define an abstract type that represents an operator, and some concrete subtypes to describe various operators. We have:
- `scal(x)` representing a scalar _x_
- `param(g)` and `paramconj(g)` representing a named scalar parameter _g_ and its conjugate _g<sup>*</sup>_
- `a(i)` and `adag(i)` representing _a<sub>i</sub>_ and _a<sub>i</sub><sup>†</sup>_, the annihilation and creation operators for bosonic mode _i_
- `σ(a,i)` representing the Pauli matrix _σ<sub>a,i</sub>_ for two-level system (TLS) _i_, where _a ∈ {x=1,y=2,z=3}_ is the type of Pauli matrix.
- `OpProd(A,B)` representing _A B_, i.e., the product of two operators
- `OpSum(A,B)` representing _A + B_, i.e., the sum of two operators
- `OpSumAnalytical(i,A)` representing _∑<sub>i</sub> A_, i.e., an analytical sum over an index (assumed to run over all possible values of _i_).
- `ExpVal(A)` representing the expectation value ⟨_A_⟩
- `Corr(AB)` representing the correlation ⟨_AB_⟩<sub>c</sub> = ⟨_AB_⟩ - ⟨_A_⟩⟨_B_⟩, with corresponding extensions for products of more operators.

All operations are defined in such a way that the finally created object is automatically transformed to "canonical" form, which is defined by the following:
- Operator sums are expanded fully, such that the final expression is always a sum of operator products. I.e., if we write _(A + B)(C + D)_, we get _AC + AD + BC + BD_.
- Operator products are expressed in a well-defined "canonical" order. To achieve this canonical form, relevant commutators etc are used, so that an expression written as a single product can turn into a sum of expressions.
    1. at most one scalar prefactor (i.e., all prefactors collapsed into one)
    1. parameters ordered alphabetically (by string comparison)
    1. expectation values ⟨_A_⟩
    1. many-body correlations ⟨_AB_⟩<sub>c</sub>
    1. bosonic operators in normal ordering (i.e., first creation, then annihilation operators), ordered by mode index
    1. Two-level Pauli matrices, ordered by TLS mode index. At most one Pauli matrix per TLS
- Operator sums are ordered first by number of operators (both bare and within expectation values and correlations), and then with the same priority rules.

Some other useful functions that are implemented:
- `comm(A,B)`: calculates the commutator of arbitrary operators [_A,B_] = _AB - BA_. This allows, e.g., to calculate Heisenberg equations of motion for the operators.
- `ascorr(x)` takes an expression `x=A B C + D E...` and writes its expectation value in terms of single-body expectation values ⟨_A_⟩, ⟨_B_⟩, ..., and many-body correlations ⟨_AB_⟩<sub>c</sub>, ⟨_ABC_⟩<sub>c</sub>, etc. Currently, up to fourth-order correlations (i.e., products of four operators) are supported.
