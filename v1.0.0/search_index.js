var documenterSearchIndex = {"docs":
[{"location":"interface/#Interface","page":"Interface","title":"Interface","text":"","category":"section"},{"location":"interface/","page":"Interface","title":"Interface","text":"These functions (mostly) form the internal interface of the package, and should not be relevant to most users.","category":"page"},{"location":"interface/","page":"Interface","title":"Interface","text":"","category":"page"},{"location":"interface/","page":"Interface","title":"Interface","text":"Modules = [QuantumAlgebra, QuantumAlgebra.OpConstructors]","category":"page"},{"location":"interface/#QuantumAlgebra.Avac","page":"Interface","title":"QuantumAlgebra.Avac","text":"Avac(A::QuExpr), vacA(A::QuExpr)\n\nSimplify operator by assuming it is applied to the vacuum from the left or right, respectively. To be precise, Avac(A) returns A such that A0 = A0, while vacA(A) does the same for 0A.\n\n\n\n\n\n","category":"function"},{"location":"interface/#QuantumAlgebra.BosonCreate-Tuple{Any, Vararg{Any, N} where N}","page":"Interface","title":"QuantumAlgebra.BosonCreate","text":"BosonCreate(name,inds): represent bosonic creation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.BosonDestroy-Tuple{Any, Vararg{Any, N} where N}","page":"Interface","title":"QuantumAlgebra.BosonDestroy","text":"BosonDestroy(name,inds): represent bosonic annihilation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.CorrPerm_isless-Tuple{Tuple, Tuple}","page":"Interface","title":"QuantumAlgebra.CorrPerm_isless","text":"CorrPerm_isless(a,b)\n\nisless for Tuples of Tuples representing products of Corr (see above) of a permutation of operators. E.g., ((1,3),(2,)) represents <A1 A3>C <A2>C. It is assumed that the total number of operators in a and b is equal, i.e., that sum(length.(a)) == sum(length.(b)).\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.CorrTup_isless-Union{Tuple{M}, Tuple{N}, Tuple{Tuple{Vararg{Int64, N}}, Tuple{Vararg{Int64, M}}}} where {N, M}","page":"Interface","title":"QuantumAlgebra.CorrTup_isless","text":"CorrTup_isless(a,b)\n\nisless for Tuples of integers that represent Corr of sorted operators (with n representing An such that n<m == An<Am). (n,m,...) ≡ Corr(AnAm...). Defined in such a way that the same order is obtained as with BaseOperator objects\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.FermionCreate-Tuple{Any, Vararg{Any, N} where N}","page":"Interface","title":"QuantumAlgebra.FermionCreate","text":"FermionCreate(name,inds): represent fermionic creation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.FermionDestroy-Tuple{Any, Vararg{Any, N} where N}","page":"Interface","title":"QuantumAlgebra.FermionDestroy","text":"FermionDestroy(name,inds): represent fermionic annihilation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.TLSCreate-Tuple{Any, Vararg{Any, N} where N}","page":"Interface","title":"QuantumAlgebra.TLSCreate","text":"TLSCreate(name,inds): represent TLS creation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.TLSDestroy-Tuple{Any, Vararg{Any, N} where N}","page":"Interface","title":"QuantumAlgebra.TLSDestroy","text":"TLSDestroy(name,inds): represent TLS annihilation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.TLSx-Tuple{Any, Vararg{Any, N} where N}","page":"Interface","title":"QuantumAlgebra.TLSx","text":"TLSx(name,inds): represent TLS x operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.TLSy-Tuple{Any, Vararg{Any, N} where N}","page":"Interface","title":"QuantumAlgebra.TLSy","text":"TLSy(name,inds): represent TLS y operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.TLSz-Tuple{Any, Vararg{Any, N} where N}","page":"Interface","title":"QuantumAlgebra.TLSz","text":"TLSz(name,inds): represent TLS z operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.corr-Tuple{QuExpr}","page":"Interface","title":"QuantumAlgebra.corr","text":"expval(A::QuExpr): replace expression A by its (formal) correlator ⟨A⟩c.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.expval-Tuple{QuExpr}","page":"Interface","title":"QuantumAlgebra.expval","text":"expval(A::QuExpr): replace expression A by its (formal) expectation value ⟨A⟩.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.expval2corrs_inds-Tuple{Int64}","page":"Interface","title":"QuantumAlgebra.expval2corrs_inds","text":"expval2corrs_inds(N::Int)\n\nfor N operators, create an array of tuples of tuples that represents the terms in a sum of products of correlators. Each tuple corresponds to a sum term, see CorrTup_isless and CorrPerm_isless for details of the format. The returned array and terms are sorted such that if the N operators are sorted, the represented expression is also sorted with the conventions of the QuantumAlgebra package. This allows to directly return a normal-ordered form.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.expval_as_corrs-Tuple{QuExpr}","page":"Interface","title":"QuantumAlgebra.expval_as_corrs","text":"expval_as_corrs(expr::QuExpr)\n\nTake an expression expr=A B C + D E... and write its expectation value in terms of correlations A_c B_c AB_c ABC_c ldots. Note that A_c = A.\n\nE.g., expval_as_corrs(adag(:n)*a(:n)) returns a^dagger_n a_n_c + a^dagger_n_c a_n_c (which is equal to a^dagger_n a_n), while expval_as_corrs(adag(:n)*a(:m)*a(:n)) returns langle a_n^dagger a_m a_n rangle_c + langle a_n^dagger rangle_c langle a_m rangle_c langle a_n rangle_c + langle a_n^dagger rangle_c langle a_m a_n rangle_c + langle a_m rangle_c langle a_n^dagger a_n rangle_c + langle a_n rangle_c langle a_n^dagger a_m rangle_c.\n\nSee also: expval, corr\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.extindices-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.extindices","text":"extindices(A) return externally visible indices of an expression\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.symmetric_index_nums-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.symmetric_index_nums","text":"symmetric_index_nums(A) return sequence of numbers of exchange-symmetric indices\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.vacA","page":"Interface","title":"QuantumAlgebra.vacA","text":"Avac(A::QuExpr), vacA(A::QuExpr)\n\nSimplify operator by assuming it is applied to the vacuum from the left or right, respectively. To be precise, Avac(A) returns A such that A0 = A0, while vacA(A) does the same for 0A.\n\n\n\n\n\n","category":"function"},{"location":"interface/#QuantumAlgebra.vacExpVal","page":"Interface","title":"QuantumAlgebra.vacExpVal","text":"vacExpVal(A::QuExpr,S::QuExpr=1)\n\nCalculate the vacuum expectation value 0S^dagger A S0, i.e., the expectation value ψAψ for the state defined by ψ= S0`.\n\n\n\n\n\n","category":"function"},{"location":"interface/#QuantumAlgebra.∑-Tuple{QuantumAlgebra.QuIndex, QuantumAlgebra.QuTerm}","page":"Interface","title":"QuantumAlgebra.∑","text":"∑(ind,A::QuExpr): return (formal) sum of expression A over index ind.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.@boson_ops-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.@boson_ops","text":"@boson_ops name: define functions $name and $(name)dag for creating bosonic annihilation and creation operators with name name\n\n\n\n\n\n","category":"macro"},{"location":"interface/#QuantumAlgebra.@fermion_ops-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.@fermion_ops","text":"@fermion_ops name: define functions $name and $(name)dag for creating fermionic annihilation and creation operators with name name\n\n\n\n\n\n","category":"macro"},{"location":"interface/#QuantumAlgebra.@tlspm_ops-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.@tlspm_ops","text":"@tlspm_ops name: define functions $(name)m and $(name)p creating jump operators for a two-level system with name name.\n\n\n\n\n\n","category":"macro"},{"location":"interface/#QuantumAlgebra.@tlsxyz_ops-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.@tlsxyz_ops","text":"@tlsxyz_ops name: define functions $(name)x, $(name)y, and $(name)z creating Pauli operators for a two-level system with name name.\n\n\n\n\n\n","category":"macro"},{"location":"interface/#QuantumAlgebra.OpConstructors.boson_ops-Tuple{Symbol}","page":"Interface","title":"QuantumAlgebra.OpConstructors.boson_ops","text":"boson_ops(name::Symbol): return functions for creating bosonic annihilation and creation operators with name name (i.e., wrappers of BosonDestroy and BosonCreate)\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.OpConstructors.fermion_ops-Tuple{Symbol}","page":"Interface","title":"QuantumAlgebra.OpConstructors.fermion_ops","text":"fermion_ops(name::Symbol): return functions for creating fermionic annihilation and creation operators with name name (i.e., wrappers of FermionDestroy and FermionCreate)\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.OpConstructors.tlspm_ops-Tuple{Symbol}","page":"Interface","title":"QuantumAlgebra.OpConstructors.tlspm_ops","text":"tlspm_ops(name::Symbol): return functions for creating jump operators for a two-level system with name name. The output of these functions depends on setting of use_σpm.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.OpConstructors.tlsxyz_ops-Tuple{Symbol}","page":"Interface","title":"QuantumAlgebra.OpConstructors.tlsxyz_ops","text":"tlsxyz_ops(name::Symbol): return functions for creating Pauli operators for a two-level system with name name. The output of these functions depends on setting of use_σpm.\n\n\n\n\n\n","category":"method"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = QuantumAlgebra","category":"page"},{"location":"","page":"Home","title":"Home","text":"DocTestSetup = :( using QuantumAlgebra )","category":"page"},{"location":"#[QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl)-quantum-operator-algebra-in-Julia","page":"Home","title":"QuantumAlgebra.jl - quantum operator algebra in Julia","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package does quantum operator algebra (i.e., algebra with non-commuting operators) in Julia, supporting bosonic, fermionic, and two-level system operators, with arbitrary names and indices, as well as sums over any of the indices. It defines an opinionated canonical form (normal ordering plus some additional rules) to automatically simplify expressions. It is recommended to use an interface that can display LaTeX formulas (e.g., Jupyter notebooks) for convenient output formatting. While there is some documentation, it is not always kept fully up to date, and it is recommended to look at the latest commit messages to get an idea about new features etc. You can also check out the notebooks in the examples folder, which can be viewed online with nbviewer and tried out interactively with Binder.","category":"page"},{"location":"#Updates-in-v1.0.0","page":"Home","title":"Updates in v1.0.0","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This is a major revision with some breaking changes. The backend has been almost completely rewritten to make the code more efficient when dealing with large expressions, and the interface has been cleaned up in several places.","category":"page"},{"location":"#Important-changes:","page":"Home","title":"Important changes:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Canonical normal form is not automatically enforced by default. In order to transform expressions to normal form, use normal_form(x). Since automatic conversion to normal form can be convenient for interactive work, it can be enabled with QuantumAlgebra.auto_normal_form(true), or alternatively by setting the environment variable QUANTUMALGEBRA_AUTO_NORMAL_FORM to \"true\" (or any value that parse(Bool,value) parses as true) before using QuantumAlgebra.\nThe function to obtain expectation values is now expval(A) (instead of ExpVal), and expval_as_corrs(A) to express an expectation value through a correlator / cumulant expansion, e.g., AB = AB_c - A_c B_c, with corresponding extensions for products of more operators. Note that for a single operator, A_c = A, but we distinguish the two formally for clarity.\nThere is a new function julia_expression(A) that converts a QuantumAlgebra object to a julia expression, which helps in using QuantumAlgebra to programatically derive codes for numerical implementation. The object A cannot contain any \"bare\" operators, but only expectation values or correlators. See the documentation for more details.\nQuantumAlgebra expressions are now printed in pretty format in the terminal etc.","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The basic functions to create QuantumAlgebra expressions (which are of type QuExpr) are","category":"page"},{"location":"","page":"Home","title":"Home","text":"a(inds...) and adag(inds...) for a and a^, the annihilation and creation operators for a bosonic mode.\nf(inds...) and fdag(inds...) for f and f^, the annihilation and creation operators for a fermionic mode.\nσx(inds...), σy(inds...), σz(inds...) for the Pauli matrices σ^xyz for a two-level system (TLS).\nσp(inds...), σm(inds...) for excitation and deexcitation operators σ^ for a two-level system (TLS).\nIndices: All of these functions take an arbitrary number of indices as arguments, which can be either integers (1,2,...) or symbolic, where symbolic indices must be a single unicode character, with possibly an integer subindex:\njulia> using QuantumAlgebra\n\njulia> a()\na()\n\njulia> adag(:i)\na†(i)\n\njulia> fdag(1,2,:i_9)\nf†(12i₉)\n\njulia> σx(:i_1, 1, :j, :k_2, :μ_2, :◔_1, :😄_121)\nσˣ(i₁1jk₂μ₂◔₁😄₁₂₁)\nYou can define your own bosonic/fermionic/two-level system operators with a set of macros:\n@boson_ops name defines new functions $name() and $(name)dag() for bosonic species name.\n@fermion_ops name defines new functions $name() and $(name)dag() for fermionic species name.\n@tlsxyz_ops name defines new functions $(name)x(), $(name)y() and $(name)z() for the Pauli matrices for two-level system species name.\n@tlspm_ops name defines new functions $(name)p() and $(name)m() for the two-level system excitation and deexcitation operators for species name.\njulia> @boson_ops b\n(QuantumAlgebra.OpConstructors.b, QuantumAlgebra.OpConstructors.bdag)\n\njulia> bdag(:k)*b(:i)\nb†(k) b(i)\nOperators with different names are assumed to belong to different \"species\" and always commute.\nparam(name::Symbol,state='n',inds...) to create a named parameter. state must be one of 'r', 'n', or 'c' for purely real, non-conjugated complex, and conjugated complex parameters. More conveniently, parameters can be entered with string macros Pr\"name_inds...\" and Pc\"name_inds...\" for real and complex parameters:\njulia> Pr\"g_i,j_2,k\"\ng(ij₂k)\n\njulia> Pr\"g_i,j_2,k\" == param(:g,'r',:i,:j_2,:k)\ntrue\n\njulia> Pc\"α_3\" == param(:α,3)\ntrue\nArithmetic operations (*, +, -, ^, adjoint=') are supported (exponents must be nonnegative integers), with any Number types integrating automatically. Division by numbers is also supported.\njulia> 5*adag(:k)*f(3)*σx(3)\n5 a†(k) f(3) σˣ(3)\n\njulia> (5//3+4im) * adag(:k)*f(3)*σx(3) + 9.4\n9.4 + (5//3+4i) a†(k) f(3) σˣ(3)\n\njulia> (a(:i)*f(:k))'\nf†(k) a†(i)\nIf you need a bare number as a QuantumAlgebra expression, you can use x*one(QuExpr) (or one(A), where A is any QuExpr).\n∑(ind,A::QuExpr) to represent an analytic sum over index ind. Since summed indices have no semantic meaning, the index within the expression gets replaced by a special numbered sum index #ᵢ, with i=1,2,....\njulia> ∑(:i,a(:i))\n∑₁ a(#₁)\nnormal_form(A::QuExpr) converts an expression to a well-defined \"canonical\" order. To achieve this canonical form, relevant commutators etc are used, so an expression written as a single product can turn into a sum of expressions. The order is essentially normal ordering (creation before annihilation operators, with σˣʸᶻ in the middle), with some additional conventions to make the normal form (hopefully) unique. In some contexts (e.g., interactive work), it can be convenient to automatically transform all expressions to normal form. This can be enabled by calling QuantumAlgebra.auto_normal_form(true), or alternatively by setting the environment variable QUANTUMALGEBRA_AUTO_NORMAL_FORM to \"true\" (or any value that parse(Bool,value) parses as true) before using QuantumAlgebra.\njulia> normal_form(a(:i)*adag(:j))\nδ(ij)  + a†(j) a(i)\nexpval(A::QuExpr) to represent an expectation value.\njulia> expval(adag(:j)*a(:i))\n⟨a†(j) a(i)⟩\nexpval_as_corrs(A::QuExpr) to represent an expectation value through its correlators, i.e., a cumulant expansion.\njulia> expval_as_corrs(adag(:j)*a(:i))\n⟨a†(j)⟩c ⟨a(i)⟩c  + ⟨a†(j) a(i)⟩c\ncomm(A::QuExpr,B::QuExpr) to calculate the commutator AB = AB - BA.\njulia> comm(a(),adag())\n-a†() a() + a() a†()\n\njulia> normal_form(comm(a(),adag()))\n1\nAvac(A) and vacA(A) simplify operators by assuming they are applied to the vacuum from the left or right, respectively. To be precise, Avac(A) returns A such that A0 = A0, while vacA(A) does the same for 0A. These functions automatically apply normal_form to assure that the operators are simplified as much as possible. Note that \"vacuum\" for two-level systems is interpreted as the lower state, σ^z0 = -0.\njulia> Avac(a())\n0\n\njulia> Avac(a(:i)*adag(:j))\nδ(ij)\n\njulia> Avac(a()*adag()*adag())\n2 a†()\n\njulia> vacA(a()*adag()*adag())\n0\n\njulia> Avac(σx())\nσˣ()\n\njulia> Avac(σz())\n-1\nvacExpVal(A,S=1) calculates the vacuum expectation value 0S^AS0, i.e., the expectation value ψAψ for the state defined by ψ=S0. The result is guaranteed to not contain any operators.\njulia> vacExpVal(adag()*a())\n0\n\njulia> vacExpVal(adag()*a(), adag()^4/sqrt(factorial(4)))\n4.000000000000001\n\njulia> vacExpVal(adag()*a(), adag()^4/sqrt(factorial(big(4))))\n4\n\njulia> vacExpVal(σx())\n0\njulia_expression(A) to obtain a julia expression that can be used to automatically build codes implementing equations derived with QuantumAlgebra. Every expectation value or correlator is treated as a separate array. Daggers are represented as ᴴ, which are valid identifiers that can appear in the array names. Note that expectation values and correlators are not distinguished, so it is best to have all expressions use the same kind.\njulia> julia_expression(expval_as_corrs(adag(:j)*a(:i)))\n:(aᴴ[j] * a[i] + aᴴa[j, i])\nAlso note that expressions are always treated as arrays, even if they have no indices (which gives zero-dimensional arrays). If you are working with scalar quantities exclusively, it might be useful to clean up the resulting expression (e.g., use MacroTools to remove the []).\njulia> julia_expression(expval(adag()*a()*σx()))\n:(aᴴaσˣ[])\nBy default, two-level system operators are represented by the Pauli matrices σˣʸᶻ, and calling σp() and σm() will give results expressed through them:\njulia> σp()\n1//2 σˣ() + 1//2i σʸ()\n\njulia> σm()\n1//2 σˣ() - 1//2i σʸ()\nThis can be changed by calling QuantumAlgebra.use_σpm(true). In this mode, σ⁺ and σ⁻ are the \"fundamental\" operators, and all expressions are written in terms of them. Note that mixing conventions within the same expression is not supported, so it is suggested to set this flag once at the beginning of any calculation.\njulia> QuantumAlgebra.use_σpm(true)\n\njulia> σp()\nσ⁺()\n\njulia> σx()\nσ⁺() + σ⁻()\n\njulia> σz()\n-1 + 2 σ⁺() σ⁻()","category":"page"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you use QuantumAlgebra in academic work, we would appreciate a citation. See CITATION.bib for the relevant references.","category":"page"}]
}
