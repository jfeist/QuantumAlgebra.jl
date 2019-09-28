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
- `scal(x)` representing a scalar x
- `param(g)` and `paramconj(g)` representing a named scalar parameter $g$ and its conjugate $g^*$
- `a(i)` and `adag(i)` representing $a_i$ and $a_i^\dagger$, the annihilation and creation operators for bosonic mode $i$
- `σ(a,i)` representing the Pauli matrix $\sigma_{a,i}$ for TLS $i$, where $a\in\{1=x,2=y,3=z\}$ is the type of Pauli matrix.
- `OpProd(A,B)` representing $A B$, i.e., the product of two operators
- `OpSum(A,B)` representing $A + B$, i.e., the sum of two operators
- `OpSumAnalytical(i,A)` representing $\sum_i A$, i.e., an analytical sum over an index (assumed to run over all possible values of $i$).

The whole idea of the code below is to define all operations in such a way that the finally created object is in "canonical" form, which we want to fulfill the following:
- Operator sums are expanded fully, i.e., if the final expression contains a sum, it is always a sum of operator products. I.e., if we write $(a_1 + a_2)(a_3+a_4)$, we get $a_1 a_3 + a_2 a_3 + a_1 a_4 + a_2 a_4$.
- Operator products of more than two factors are represented in nested form from left to right, i.e., writing $(a_1 a_2) ((a_3 a_4) a_5)$ is rewritten to $a_1 (a_2 (a_3 (a_4 a_5)))$
- Operator products are ordered in the following order (in order to achieve this canonical form, relevant commutators etc are used, so that an expression written as a single product can turn into a sum of expressions):
    1. at most one scalar prefactor (i.e., all prefactors collapsed into one)
    1. conjugated parameters ordered alphabetically (by string comparison)
    1. parameters ordered alphabetically (by string comparison)
    1. expectation values $\langle A \rangle$
    1. many-body correlators $\langle A B\rangle_c$
    1. bosonic creation operators, in ascending order of mode number
    1. bosonic annihilation operators, in ascending order of mode number
    1. Pauli matrices, in ascending order of TLS number. At most one Pauli matrix per TLS
- Operator sums are also ordered with the same priority rules.

Finally, we also implement the `comm(A,B)` function to calculate the commutator of arbitrary operators $A$ and $B$ (which can be arbitrarily nested products and sums of simpler operators). This allows us to calculate the Heisenberg equations of motion in the end.
