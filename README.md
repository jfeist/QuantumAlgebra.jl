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

# Updates in v0.4.0
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
  correlator / cumulant expansion, e.g., ⟨_AB_⟩ = ⟨_AB_⟩<sub>c</sub> - ⟨_A_⟩<sub>c</sub>
  ⟨_B_⟩<sub>c</sub>, with corresponding extensions for products of more
  operators. Note that for a single operator, ⟨_A_⟩<sub>c</sub> = ⟨_A_⟩, but we
  distinguish the two formally for clarity.
- There is a new function `julia_expression(A)` that converts a QuantumAlgebra
  object to a julia expression, which helps in using QuantumAlgebra to
  programatically derive codes for numerical implementation. The object `A`
  cannot contain any "bare" operators, but only expectation values or
  correlators. See the documentation for more details.
- Quantum expressions are now printed in pretty format in the terminal etc.

## Overview

We define an abstract type that represents an operator, and some concrete subtypes to describe various operators. We have:
- `scal(x)` representing a scalar _x_
- `param(g,(i,j),state='n')` representing a named scalar parameter _g<sub>i,j</sub>_. `state` can be `'r'` for purely real parameters (invariant under complex conjugation), `'n'` for not-conjugated values, and `'c'` for a conjugated parameter _g<sub>i,j</sub><sup>*</sup>_.
- `a(i)` and `adag(i)` representing _a<sub>i</sub>_ and _a<sub>i</sub><sup>†</sup>_, the annihilation and creation operators for bosonic mode _i_
- `f(i)` and `fdag(i)` representing _f<sub>i</sub>_ and _f<sub>i</sub><sup>†</sup>_, the annihilation and creation operators for fermionic mode _i_
- `σ(a,i)` representing the Pauli matrix _σ<sub>a,i</sub>_ for two-level system (TLS) _i_, where _a ∈ {x=1,y=2,z=3}_ is the type of Pauli matrix.
- `OpProd(A,B)` representing _A B_, i.e., the product of two operators
- `OpSum(A,B)` representing _A + B_, i.e., the sum of two operators
- `OpSumAnalytical(i,A)` or `∑(i,A)` representing _∑<sub>i</sub> A_, i.e., an analytical sum over an index (assumed to run over all possible values of _i_).
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
    1. fermionic operators in normal ordering (i.e., first creation, then annihilation operators), ordered by mode index
    1. Two-level Pauli matrices, ordered by TLS mode index. At most one Pauli matrix per TLS
- Operator sums are ordered first by number of operators (both bare and within expectation values and correlations), and then with the same priority rules.

Some other useful functions that are implemented:
- `comm(A,B)`: calculates the commutator of arbitrary operators [_A,B_] = _AB - BA_. This allows, e.g., to calculate Heisenberg equations of motion for the operators.
- `ascorr(x)` takes an expression `x=A B C + D E...` and writes its expectation value in terms of single-body expectation values ⟨_A_⟩, ⟨_B_⟩, ..., and many-body correlations ⟨_AB_⟩<sub>c</sub>, ⟨_ABC_⟩<sub>c</sub>, etc. Currently, up to fourth-order correlations (i.e., products of four operators) are supported.
- `Avac(A)` and `vacA(A)` simplify operators by assuming they are applied to the vacuum from the left or right, respectively. To be precise, `Avac(A)` returns _A'_ such that _A_|0⟩ = _A'_|0⟩, while `vacA(A)` does the same for ⟨0|_A_.
- `vacExpVal(A,S=1)` calculates the vacuum expectation value ⟨0|_S<sup>†</sup>AS_|0⟩, i.e., the expectation value ⟨ψ|_A_|ψ⟩ for the state defined by |ψ⟩=_S_|0⟩.

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant references.
