var documenterSearchIndex = {"docs":
[{"location":"release_notes/#Release-notes","page":"Release notes","title":"Release notes","text":"","category":"section"},{"location":"release_notes/#v1.5.0-(2024-07-19)","page":"Release notes","title":"v1.5.0 (2024-07-19)","text":"","category":"section"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"This is a minor release with some bug fixes and new features:","category":"page"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"Add heisenberg_eom_system function to generate a system of Heisenberg equations of motion for either the expectation values or the cumulants / correlators of operator products, starting from a Hamiltonian and Lindblad operators describing decoherence. Typically, these equation systems are not closed without approximations as equations for products of n operators involve products of m>n operators, so the system has to be truncated. This is achieved with a filter function that removes higher-order terms or rewrites them (approximately) in terms of lower-order expressions. Simple example:\njulia> using QuantumAlgebra\njulia> @boson_ops a;\njulia> H = Pr\"ω\"*a'()*a() + Pr\"χ\"*a'()*(a'()+a())*a();\njulia> Ls = ((Pr\"γ\",a()),);\njulia> EQ = heisenberg_eom_system(H,2,Ls,a())\ndₜ⟨a()⟩ = -1//2 γ ⟨a()⟩  - 1i ω ⟨a()⟩  - 2i χ ⟨a†() a()⟩  - 1i χ ⟨a()²⟩ \ndₜ⟨a†() a()⟩ = -γ ⟨a†() a()⟩ \ndₜ⟨a()²⟩ = -2i χ ⟨a()⟩  - γ ⟨a()²⟩  - 2i ω ⟨a()²⟩  \nAdd additional argument modes_in_vacuum to vacA, Avac, and vacExpVal, to restrict the modes that are considered to be in the vacuum state, with all other operators unaffected.\nFix a LaTeX output error where Latexify misinterpreted a superscript \"^+\" as a sum and replaced \"x^+ - y\" with \"x^- y\".","category":"page"},{"location":"release_notes/#v1.4.0-(2024-06-29)","page":"Release notes","title":"v1.4.0 (2024-06-29)","text":"","category":"section"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"This is a minor release with some bug fixes and new features:","category":"page"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"Allow Number type inputs for most exported functions so that they do not require QuExpr wrapping, such that expval(3) or Avac(1) no longer raise errors (closes #14).\nSupport interoperability with symbolic computer algebra systems (CAS) such as Symbolics.jl or SymPy.jl / SymPyPythonCall.jl. Expressions provided by these systems can be used as the \"scalar\" factors of each term in a QuExpr, and thus provide an alternative to the QuantumAlgebra Param objects. In contrast to Params, they do not support symbolic indices, but provide much more flexibility in terms of the mathematical operations and powerful manipulation functions. For example, one can use\njulia> using QuantumAlgebra, Symbolics\njulia> @boson_ops a;\njulia> @variables x y;\njulia> H = cos(x)^2 * a() * a'() + sin(x)^2 * a'() * a();\njulia> Symbolics.simplify(normal_form(H))\ncos(x)^2 + a†() a()\n\njulia> Symbolics.substitute(H, x => y)\nsin(y)^2 a†() a() + cos(y)^2 a() a†()\n\njulia> Symbolics.substitute(H, x => 2)\n0.826821810431806 a†() a() + 0.17317818956819406 a() a†()\n\njulia> Symbolics.substitute(H, x => 2; fold=false)\nsin(2)^2 a†() a() + cos(2)^2 a() a†()\nWhile overloads for some functions (like simplify and substitute from Symbolics) are available, all other functions/manipulations from the CAS packages can be applied with the new function map_scalar_function which applies a function to each scalar factor in a QuExpr. This can be used as, e.g.,\njulia> map_scalar_function(Symbolics.simplify, normal_form(H))\ncos(x)^2 + a†() a()\nIn order not to add dependencies on several CAS packages to QuantumAlgebra, the interoperability helpers are defined in extension modules, using PackageExtensionCompat.jl to maintain compatibility with Julia versions before v1.9.","category":"page"},{"location":"release_notes/#v1.3.1-(2023-12-22)","page":"Release notes","title":"v1.3.1 (2023-12-22)","text":"","category":"section"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"This is a patch release with some bug fixes and performance improvements.","category":"page"},{"location":"release_notes/#v1.3.0-(2023-09-19)","page":"Release notes","title":"v1.3.0 (2023-09-19)","text":"","category":"section"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"This is a minor release with some bug fixes in LaTeX output, added unit tests for many functions, and some new features and deprecations:","category":"page"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"The functions that generate operators now can be conjugated directly, e.g., a'() is equivalent to adag() (and to a()', which was already the case beforehand). The old _dag functions for bosons and fermions are deprecated and will be removed in a future version.\nPreferences.jl is now used to set and store preferences. The following preferences are available:\n\"define_default_ops\": if this is set to false (default is true), the \"default\" operators a, adag, f, fdag, σx, σy, σz, σp, σm are not defined upon import. Note that changing this value requires restarting the Julia session to take effect. The setting can be changed with QuantumAlgebra.set_define_default_ops(true/false) (which will inform you whether a restart is required) or with Preferences.set_preferences!(QuantumAlgebra,\"define_default_ops\"=>true/false).\n\"auto_normal_form\": Choose whether all expressions are automatically converted to normal form upon creation. The default is false. It can be changed for a single session with QuantumAlgebra.auto_normal_form(true/false), and can be made permanent with QuantumAlgebra.auto_normal_form(true/false; set_preference=true) or with Preferences.set_preferences!(QuantumAlgebra,\"auto_normal_form\"=>true/false). Note that this could previously be set by defining an environment variable \"QUANTUMALGEBRA_AUTO_NORMAL_FORM\", but this usage has been deprecated and will be removed in a future version.\n\"use_σpm\": Choose whether for two-level systems, the \"basic\" operators are excitation/deexcitation operators σ^ or the Pauli matrices σ^xyz. This can be changed in a single session by calling QuantumAlgebra.use_σpm(true/false), and can be made permanent with QuantumAlgebra.use_σpm(true/false; set_preference=true) or with Preferences.set_preferences!(QuantumAlgebra,\"use_σpm\"=>true/false).\nUse PrecompileTools.jl instead of SnoopPrecompile.jl for precompilation.","category":"page"},{"location":"release_notes/#v1.2.0-(2023-04-22)","page":"Release notes","title":"v1.2.0 (2023-04-22)","text":"","category":"section"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"This is a minor revision with some bug fixes, and some new features:","category":"page"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"Add the heisenberg_eom function to calculate the Heisenberg equations of motions of operators.\nSwitch the LaTeX output to use Latexify.jl, which also enables support for display of vectors, tuples, etc of QuExpr.\nRequire at least Julia v1.6.","category":"page"},{"location":"release_notes/#v1.1.0-(2021-09-02)","page":"Release notes","title":"v1.1.0 (2021-09-02)","text":"","category":"section"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"This is a minor revision with some and bug fixes, and one new feature: It is now possible to use @anticommuting_fermion_group to define several fermionic operators that mutually anticommute, allowing to refer to different kinds of states for the same \"species\" (e.g., localized and itinerant electrons):","category":"page"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"julia> @anticommuting_fermion_group c d\n\njulia> normal_form(c()*d() + d()*c())\n0","category":"page"},{"location":"release_notes/#v1.0.0-(2021-07-23)","page":"Release notes","title":"v1.0.0 (2021-07-23)","text":"","category":"section"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"This is a major revision with some breaking changes. The backend has been almost completely rewritten to make the code more efficient when dealing with large expressions, and the interface has been cleaned up in several places.","category":"page"},{"location":"release_notes/#Important-changes:","page":"Release notes","title":"Important changes:","text":"","category":"section"},{"location":"release_notes/","page":"Release notes","title":"Release notes","text":"Canonical normal form is not automatically enforced by default. In order to transform expressions to normal form, use normal_form(x). Since automatic conversion to normal form can be convenient for interactive work, it can be enabled with QuantumAlgebra.auto_normal_form(true), or alternatively by setting the environment variable QUANTUMALGEBRA_AUTO_NORMAL_FORM to \"true\" (or any value that parse(Bool,value) parses as true) before using QuantumAlgebra.\nThe function to obtain expectation values is now expval(A) (instead of ExpVal), and expval_as_corrs(A) to express an expectation value through a correlator / cumulant expansion, e.g., AB = AB_c - A_c B_c, with corresponding extensions for products of more operators. Note that for a single operator, A_c = A, but we distinguish the two formally for clarity.\nThere is a new function julia_expression(A) that converts a QuantumAlgebra object to a julia expression, which helps in using QuantumAlgebra to programatically derive codes for numerical implementation. The object A cannot contain any \"bare\" operators, but only expectation values or correlators. See the documentation for more details.\nQuantumAlgebra expressions are now printed in pretty format in the terminal etc.","category":"page"},{"location":"interface/#Interface","page":"Interface","title":"Interface","text":"","category":"section"},{"location":"interface/","page":"Interface","title":"Interface","text":"These functions (mostly) form the internal interface of the package, and should not be relevant to most users.","category":"page"},{"location":"interface/","page":"Interface","title":"Interface","text":"","category":"page"},{"location":"interface/","page":"Interface","title":"Interface","text":"Modules = [QuantumAlgebra]","category":"page"},{"location":"interface/#QuantumAlgebra.Avac","page":"Interface","title":"QuantumAlgebra.Avac","text":"Avac(A::QuExpr), vacA(A::QuExpr)\n\nSimplify operator by assuming it is applied to the vacuum from the left or right, respectively. To be precise, Avac(A) returns A such that A0 = A0, while vacA(A) does the same for 0A.\n\n\n\n\n\n","category":"function"},{"location":"interface/#QuantumAlgebra.BosonCreate-Tuple{Any, Vararg{Any}}","page":"Interface","title":"QuantumAlgebra.BosonCreate","text":"BosonCreate(name,inds): represent bosonic creation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.BosonDestroy-Tuple{Any, Vararg{Any}}","page":"Interface","title":"QuantumAlgebra.BosonDestroy","text":"BosonDestroy(name,inds): represent bosonic annihilation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.CorrPerm_isless-Tuple{Tuple, Tuple}","page":"Interface","title":"QuantumAlgebra.CorrPerm_isless","text":"CorrPerm_isless(a,b)\n\nisless for Tuples of Tuples representing products of Corr (see above) of a permutation of operators. E.g., ((1,3),(2,)) represents <A1 A3>C <A2>C. It is assumed that the total number of operators in a and b is equal, i.e., that sum(length.(a)) == sum(length.(b)).\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.CorrTup_isless-Union{Tuple{M}, Tuple{N}, Tuple{NTuple{N, Int64}, NTuple{M, Int64}}} where {N, M}","page":"Interface","title":"QuantumAlgebra.CorrTup_isless","text":"CorrTup_isless(a,b)\n\nisless for Tuples of integers that represent Corr of sorted operators (with n representing An such that n<m == An<Am). (n,m,...) ≡ Corr(AnAm...). Defined in such a way that the same order is obtained as with BaseOperator objects\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.FermionCreate-Tuple{Any, Vararg{Any}}","page":"Interface","title":"QuantumAlgebra.FermionCreate","text":"FermionCreate(name,inds): represent fermionic creation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.FermionDestroy-Tuple{Any, Vararg{Any}}","page":"Interface","title":"QuantumAlgebra.FermionDestroy","text":"FermionDestroy(name,inds): represent fermionic annihilation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.TLSCreate-Tuple{Any, Vararg{Any}}","page":"Interface","title":"QuantumAlgebra.TLSCreate","text":"TLSCreate(name,inds): represent TLS creation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.TLSDestroy-Tuple{Any, Vararg{Any}}","page":"Interface","title":"QuantumAlgebra.TLSDestroy","text":"TLSDestroy(name,inds): represent TLS annihilation operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.TLSx-Tuple{Any, Vararg{Any}}","page":"Interface","title":"QuantumAlgebra.TLSx","text":"TLSx(name,inds): represent TLS x operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.TLSy-Tuple{Any, Vararg{Any}}","page":"Interface","title":"QuantumAlgebra.TLSy","text":"TLSy(name,inds): represent TLS y operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.TLSz-Tuple{Any, Vararg{Any}}","page":"Interface","title":"QuantumAlgebra.TLSz","text":"TLSz(name,inds): represent TLS z operator name_inds\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.boson_ops-Tuple{Symbol}","page":"Interface","title":"QuantumAlgebra.boson_ops","text":"boson_ops(name::Symbol): return function for creating bosonic annihilation and creation operators with name name (i.e., wrappers of BosonDestroy and BosonCreate)\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.corr-Tuple{QuExpr}","page":"Interface","title":"QuantumAlgebra.corr","text":"expval(A::QuExpr): replace expression A by its (formal) correlator ⟨A⟩c.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.expval-Tuple{QuExpr}","page":"Interface","title":"QuantumAlgebra.expval","text":"expval(A::QuExpr): replace expression A by its (formal) expectation value ⟨A⟩.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.expval2corrs_inds-Tuple{Int64}","page":"Interface","title":"QuantumAlgebra.expval2corrs_inds","text":"expval2corrs_inds(N::Int)\n\nfor N operators, create an array of tuples of tuples that represents the terms in a sum of products of correlators. Each tuple corresponds to a sum term, see CorrTup_isless and CorrPerm_isless for details of the format. The returned array and terms are sorted such that if the N operators are sorted, the represented expression is also sorted with the conventions of the QuantumAlgebra package. This allows to directly return a normal-ordered form.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.expval_as_corrs-Tuple{QuExpr}","page":"Interface","title":"QuantumAlgebra.expval_as_corrs","text":"expval_as_corrs(expr::QuExpr)\n\nTake an expression expr=A B C + D E... and write its expectation value in terms of correlations A_c B_c AB_c ABC_c ldots. Note that A_c = A.\n\nE.g., expval_as_corrs(a'(:n)*a(:n)) returns a^dagger_n a_n_c + a^dagger_n_c a_n_c (which is equal to a^dagger_n a_n), while expval_as_corrs(a'(:n)*a(:m)*a(:n)) returns langle a_n^dagger a_m a_n rangle_c + langle a_n^dagger rangle_c langle a_m rangle_c langle a_n rangle_c + langle a_n^dagger rangle_c langle a_m a_n rangle_c + langle a_m rangle_c langle a_n^dagger a_n rangle_c + langle a_n rangle_c langle a_n^dagger a_m rangle_c.\n\nSee also: expval, corr\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.extindices-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.extindices","text":"extindices(A) return externally visible indices of an expression\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.fermion_ops-Tuple{Symbol}","page":"Interface","title":"QuantumAlgebra.fermion_ops","text":"fermion_ops(name::Symbol): return function for creating fermionic annihilation and creation operators with name name (i.e., wrappers of FermionDestroy and FermionCreate)\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.heisenberg_eom","page":"Interface","title":"QuantumAlgebra.heisenberg_eom","text":"heisenberg_eom(A,H,Ls=()) calculates dotA = i HA + sum_i (L_i^ A L_i - ½ L_i^ L_i A), where Ls = (L_1,L_2,...) is an iterable of Lindblad operators.\n\n\n\n\n\n","category":"function"},{"location":"interface/#QuantumAlgebra.symmetric_index_nums-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.symmetric_index_nums","text":"symmetric_index_nums(A) return sequence of numbers of exchange-symmetric indices\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.tlspm_ops-Tuple{Symbol}","page":"Interface","title":"QuantumAlgebra.tlspm_ops","text":"tlspm_ops(name::Symbol): return functions for creating jump operators for a two-level system with name name. The output of these functions depends on setting of use_σpm.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.tlsxyz_ops-Tuple{Symbol}","page":"Interface","title":"QuantumAlgebra.tlsxyz_ops","text":"tlsxyz_ops(name::Symbol): return functions for creating Pauli operators for a two-level system with name name. The output of these functions depends on setting of use_σpm.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.vacA","page":"Interface","title":"QuantumAlgebra.vacA","text":"Avac(A::QuExpr), vacA(A::QuExpr)\n\nSimplify operator by assuming it is applied to the vacuum from the left or right, respectively. To be precise, Avac(A) returns A such that A0 = A0, while vacA(A) does the same for 0A.\n\n\n\n\n\n","category":"function"},{"location":"interface/#QuantumAlgebra.vacExpVal","page":"Interface","title":"QuantumAlgebra.vacExpVal","text":"vacExpVal(A::QuExpr,S::QuExpr=1)\n\nCalculate the vacuum expectation value 0S^dagger A S0, i.e., the expectation value ψAψ for the state defined by ψ= S0`.\n\n\n\n\n\n","category":"function"},{"location":"interface/#QuantumAlgebra.∑-Tuple{QuantumAlgebra.QuIndex, QuantumAlgebra.QuTerm}","page":"Interface","title":"QuantumAlgebra.∑","text":"∑(ind,A::QuExpr): return (formal) sum of expression A over index ind.\n\n\n\n\n\n","category":"method"},{"location":"interface/#QuantumAlgebra.@anticommuting_fermion_group-Tuple","page":"Interface","title":"QuantumAlgebra.@anticommuting_fermion_group","text":"@anticommuting_fermion_group name1 name2 ...: define a group of mutually anticommuting fermionic operators\n\n\n\n\n\n","category":"macro"},{"location":"interface/#QuantumAlgebra.@boson_ops-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.@boson_ops","text":"@boson_ops name: define function $name for creating bosonic annihilation operators with name name (also defines deprecated $(name)dag, use $(name)' instead) \n\n\n\n\n\n","category":"macro"},{"location":"interface/#QuantumAlgebra.@fermion_ops-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.@fermion_ops","text":"@fermion_ops name: define function $name for creating fermionic annihilation operators with name name (also defines deprecated $(name)dag, use $(name)' instead)\n\n\n\n\n\n","category":"macro"},{"location":"interface/#QuantumAlgebra.@tlspm_ops-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.@tlspm_ops","text":"@tlspm_ops name: define functions $(name)m and $(name)p creating jump operators for a two-level system with name name.\n\n\n\n\n\n","category":"macro"},{"location":"interface/#QuantumAlgebra.@tlsxyz_ops-Tuple{Any}","page":"Interface","title":"QuantumAlgebra.@tlsxyz_ops","text":"@tlsxyz_ops name: define functions $(name)x, $(name)y, and $(name)z creating Pauli operators for a two-level system with name name.\n\n\n\n\n\n","category":"macro"},{"location":"","page":"Home","title":"Home","text":"CurrentModule = QuantumAlgebra","category":"page"},{"location":"","page":"Home","title":"Home","text":"DocTestSetup = quote\n    using QuantumAlgebra\n    QuantumAlgebra.auto_normal_form(false)\n    QuantumAlgebra.use_σxyz()\n    @static if !QuantumAlgebra._DEFINE_DEFAULT_OPS\n        @boson_ops a\n        @fermion_ops f\n        @tlsxyz_ops σ\n        @tlspm_ops σ\n    end\nend","category":"page"},{"location":"#[QuantumAlgebra.jl](https://github.com/jfeist/QuantumAlgebra.jl)-quantum-operator-algebra-in-Julia","page":"Home","title":"QuantumAlgebra.jl - quantum operator algebra in Julia","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This package does quantum operator algebra (i.e., algebra with non-commuting operators) in Julia, supporting bosonic, fermionic, and two-level system operators, with arbitrary names and indices, as well as sums over any of the indices. It defines an opinionated canonical form (normal ordering plus some additional rules) to automatically simplify expressions. It is recommended to use an interface that can display LaTeX formulas (e.g., Jupyter notebooks) for convenient output formatting. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Starting from v1.4, QuantumAlgebra also interoperates with computer algebra systems (CAS) such as Symbolics.jl or SymPy.jl / SymPyPythonCall.jl, as the \"scalar\" prefactors of each quantum term can be arbitrary expressions provided by these systems. While such expressions do not support symbolic indices in the same way as QuantumAlgebra, they provide much more flexibility in terms of the mathematical operations and powerful manipulation functions possible on the parameters.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Example jupyter notebooks are available in the examples folder and can be viewed online with nbviewer and tried out interactively with Binder.","category":"page"},{"location":"#Release-notes-/-changelog","page":"Home","title":"Release notes / changelog","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Please see the release notes for a summary of changes in each version.","category":"page"},{"location":"#Overview","page":"Home","title":"Overview","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The basic functions to create QuantumAlgebra expressions (which are of type QuExpr) are","category":"page"},{"location":"","page":"Home","title":"Home","text":"a(inds...) and a'(inds...) for a and a^, the annihilation and creation operators for a bosonic mode.\nf(inds...) and f'(inds...) for f and f^, the annihilation and creation operators for a fermionic mode.\nσx(inds...), σy(inds...), σz(inds...) for the Pauli matrices σ^xyz for a two-level system (TLS).\nσp(inds...), σm(inds...) for excitation and deexcitation operators σ^ for a two-level system (TLS).\nIndices: All of these functions take an arbitrary number of indices as arguments, which can be either integers (1,2,...) or symbolic, where symbolic indices must be a single unicode character, with possibly an integer subindex:\njulia> using QuantumAlgebra\n\njulia> a()\na()\n\njulia> a'(:i)\na†(i)\n\njulia> f'(1,2,:i_9)\nf†(12i₉)\n\njulia> σx(:i_1, 1, :j, :k_2, :μ_2, :◔_1, :😄_121)\nσˣ(i₁1jk₂μ₂◔₁😄₁₂₁)\nYou can define your own bosonic/fermionic/two-level system operators with a set of macros:\n@boson_ops name defines new functions $name() and $(name)dag() for bosonic species name.\n@fermion_ops name defines new functions $name() and $(name)dag() for fermionic species name.\n@tlsxyz_ops name defines new functions $(name)x(), $(name)y() and $(name)z() for the Pauli matrices for two-level system species name.\n@tlspm_ops name defines new functions $(name)p() and $(name)m() for the two-level system excitation and deexcitation operators for species name.\nNote that for @boson_ops and @fermion_ops, deprecated $(name)dag() functions are defined for backward compatibility. These will be removed in a future version, as $(name)'() is now the preferred syntax for creating an adjoint.\njulia> @boson_ops b\n(b (QuExpr constructor), b† (QuExpr constructor))\n\njulia> b'(:k)*b(:i)\nb†(k) b(i)\nOperators with different names are assumed to belong to different \"species\" and always commute. For fermions, this is not always desired, since you might want to use different named operators to refer to different kinds of states for the same species (e.g., localized and itinerant electrons). This can be achieved with the macro @anticommuting_fermion_group, which creates several fermionic operators that mutually anticommute:\njulia> @anticommuting_fermion_group c d\n\njulia> normal_form(c()*d() + d()*c())\n0\nparam(name::Symbol,state='n',inds...) to create a named parameter. state must be one of 'r', 'n', or 'c' for purely real, non-conjugated complex, and conjugated complex parameters. More conveniently, parameters can be entered with string macros Pr\"name_inds...\" and Pc\"name_inds...\" for real and complex parameters:\njulia> Pr\"g_i,j_2,k\"\ng(ij₂k)\n\njulia> Pr\"g_i,j_2,k\" == param(:g,'r',:i,:j_2,:k)\ntrue\n\njulia> Pc\"α_3\" == param(:α,3)\ntrue\nArithmetic operations (*, +, -, ^, adjoint=') are supported (exponents must be nonnegative integers), with any Number types integrating automatically. Division by numbers is also supported.\njulia> 5*a'(:k)*f(3)*σx(3)\n5 a†(k) f(3) σˣ(3)\n\njulia> (5//3+4im) * a'(:k)*f(3)*σx(3) + 9.4\n9.4 + (5//3+4i) a†(k) f(3) σˣ(3)\n\njulia> (a(:i)*f(:k))'\nf†(k) a†(i)\nIf you need a bare number as a QuantumAlgebra expression, you can use x*one(QuExpr) (or one(A), where A is any QuExpr).\n∑(ind,A::QuExpr) to represent an analytic sum over index ind. Since summed indices have no semantic meaning, the index within the expression gets replaced by a special numbered sum index #ᵢ, with i=1,2,....\njulia> ∑(:i,a(:i))\n∑₁ a(#₁)\nnormal_form(A::QuExpr) converts an expression to a well-defined \"canonical\" order. To achieve this canonical form, relevant commutators etc are used, so an expression written as a single product can turn into a sum of expressions. The order is essentially normal ordering (creation before annihilation operators, with σˣʸᶻ in the middle), with some additional conventions to make the normal form (hopefully) unique. In some contexts (e.g., interactive work), it can be convenient to automatically transform all expressions to normal form. This can be enabled by calling QuantumAlgebra.auto_normal_form(true). To make the setting permanent, call QuantumAlgebra.auto_normal_form(true; set_preference=true) or alternatively use Preferences.jl directly, i.e., call Preferences.set_preferences!(QuantumAlgebra,\"auto_normal_form\"=>true/false).\njulia> normal_form(a(:i)*a'(:j))\nδ(ij)  + a†(j) a(i)\nexpval(A::QuExpr) to represent an expectation value.\njulia> expval(a'(:j)*a(:i))\n⟨a†(j) a(i)⟩\nexpval_as_corrs(A::QuExpr) to represent an expectation value through its correlators, i.e., a cumulant expansion.\njulia> expval_as_corrs(a'(:j)*a(:i))\n⟨a†(j)⟩c ⟨a(i)⟩c  + ⟨a†(j) a(i)⟩c\ncomm(A::QuExpr,B::QuExpr) to calculate the commutator AB = AB - BA.\njulia> comm(a(),a'())\n-a†() a() + a() a†()\n\njulia> normal_form(comm(a(),a'()))\n1\nAvac(A) and vacA(A) simplify operators by assuming they are applied to the vacuum from the left or right, respectively. To be precise, Avac(A) returns A such that A0 = A0, while vacA(A) does the same for 0A. These functions automatically apply normal_form to assure that the operators are simplified as much as possible. Note that \"vacuum\" for two-level systems is interpreted as the lower state, σ^z0 = -0.\njulia> Avac(a())\n0\n\njulia> Avac(a(:i)*a'(:j))\nδ(ij)\n\njulia> Avac(a()*a'()*a'())\n2 a†()\n\njulia> vacA(a()*a'()*a'())\n0\n\njulia> Avac(σx())\nσˣ()\n\njulia> Avac(σz())\n-1\nBoth functions can also be called with an optional second argument, Avac(A,modes_in_vacuum) or vacA(A,modes_in_vacuum), which is an iterable over operators (or a single operator) that will be assumed to be in the vacuum state, while all others are not. Note that the operators in modes_in_vacuum do not distinguish by index, i.e., if the modes have indices, all modes with the same name are assumed to be in the vacuum state. To avoid confusion, the modes_in_vacuum argument thus does not accept operators with indices.\njulia> Avac(a(),a())\n0\n\njulia> Avac(a(),f())\na()\n\njulia> Avac(a(:i)*a'(:j),f())\nδ(ij)  + a†(j) a(i)\n\njulia> Avac(a'()*a()*f()*f'(),f())\na†() a()\n\njulia> @boson_ops b;\n\njulia> Avac(a'()*a()*b()*b'()^2*f()*f'(),(f(),b()))\n2 a†() b†() a()\nvacExpVal(A,S=1) calculates the vacuum expectation value 0S^AS0, i.e., the expectation value ψAψ for the state defined by ψ=S0. The result is guaranteed to not contain any operators.\njulia> vacExpVal(a'()*a())\n0\n\njulia> vacExpVal(a'()*a(), a'()^4/sqrt(factorial(4)))\n4\n\njulia> vacExpVal(a'()*a(), a'()^4/sqrt(factorial(big(4))))\n4\n\njulia> vacExpVal(σx())\n0\nLike vacA and Avac, vacExpVal also takes an optional modes_in_vacuum argument, vacExpVal(A,S,modes_in_vacuum) (since all arguments are positional, S has to be given explicitly in this case even if it is just the identity operator, i.e., vacExpVal(A,1,a())):\njulia> @boson_ops b;\n\njulia> vacExpVal(a'()*a()*b()^2*b'()^2*f()*f'(), 1, (f(),b()))\n2 a†() a()\nheisenberg_eom(A,H,Ls=()) calculates the Heisenberg equation of motion for operator A under the action of Hamiltonian H and potential Lindblad decay terms Ls, given by fracdAdt = iHA + _i γ_i (L_i^ A L_i - frac12 L_i^ L_iA). The Lindblad decay operators are passed as a tuple (not an array) of tuples, where each inner tuple describes one decay operator. The possible forms are (L,) for decay operator L, (γ,L) for decay operator L with rate γ, and (inds,γ,L) for decay operators summed over the given indices (note that this is different from the operator itself being a sum, seen in the example below). Finally, L can (in all three cases above) be just a single operator or a tuple of two operators L=(X,Y) to represent off-diagonal Lindblad terms L_XYρ = X ρ Y^ - frac12 Y^ Xρ.\njulia> H = Pr\"ω\"*a'()a();\n\njulia> Ls = ((Pr\"γ\",a()),);\n\njulia> normal_form(heisenberg_eom(a(),H,Ls))\n-1//2 γ a() - 1i ω a()\n\njulia> H = QuExpr();\n\njulia> Ls = ((:i,a(:i)),);\n\njulia> normal_form(heisenberg_eom(a(:i),H,Ls))\n-1//2 a(i)\n\njulia> Ls = ((∑(:i,a(:i)),),);\n\njulia> normal_form(heisenberg_eom(a(:i),H,Ls))\n-1//2 ∑₁ a(#₁)\n\njulia> Ls = (((:i,:j),(a(:i),a(:j))),);\n\njulia> normal_form(heisenberg_eom(a(:i),H,Ls))\n-1//2 ∑₁ a(#₁)\nheisenberg_eom_system(H,rhsfilt,Ls=(),ops=nothing) calculates the system of equations of motion for the expectation values of operators appearing in H and Ls (same conventions as for heisenberg_eom above). Typically, these equation systems are not closed without approximations as equations for products of n operators involve products of mn operators, so the system has to be truncated. This is achieved with a filter function that removes higher-order terms or rewrites them (approximately) in terms of lower-order expressions. The function rhsfilt is applied to the right-hand side of the equations to filter them as desired. If rhsfilt(A::QuExpr)::QuExpr is a function, it will be applied to the calculated right-hand side of the equations. QuantumAlgebra comes with two predefined constructors for filter functions, droplen(maxorder::Int), which leads to all terms of order higher than maxorder being neglected, and dropcorr(maxorder::Int), where all terms of order higher than maxorder are rewritten in terms of lower-order expressions up to order maxorder and higher-order correlators, with those correlations being neglected (i.e., dropcorr(1) will replace a^ a = a^ a_c + a^ a  a^ a). If rhsfilt is a number, it will be interpreted as droplen(rhsfilt). Finally, the ops argument can be used to specify the operators that should be used to \"seed\" the system of equations, otherwise all operators appearing in H are used.\njulia> H = Pr\"ω\"*a'()*a() + Pr\"χ\"*a'()*(a'()+a())*a();\n\njulia> Ls = ((Pr\"γ\",a()),);\n\njulia> heisenberg_eom_system(H,2,Ls,a())\ndₜ⟨a()⟩ = -1//2 γ ⟨a()⟩  - 1i ω ⟨a()⟩  - 2i χ ⟨a†() a()⟩  - 1i χ ⟨a()²⟩ \ndₜ⟨a†() a()⟩ = -γ ⟨a†() a()⟩ \ndₜ⟨a()²⟩ = -2i χ ⟨a()⟩  - γ ⟨a()²⟩  - 2i ω ⟨a()²⟩  \nThe heisenberg_eom_system function can also be passed either ExpVal or Corr as a first argument, which will give the equations of motion of the expectation values (the default) or correlators (corresponding to a cumulant expansion) of the operators.\njulia> H = Pr\"ω\"*a'()*a() + Pr\"χ\"*a'()*(a'()+a())*a();\n\njulia> Ls = ((Pr\"γ\",a()),);\n\njulia> heisenberg_eom_system(Corr,H,1,Ls,a())\ndₜ⟨a()⟩c = -1//2 γ ⟨a()⟩c  - 1i ω ⟨a()⟩c  - 2i χ ⟨a†()⟩c ⟨a()⟩c  - 1i χ ⟨a()⟩c² \njulia_expression(A) to obtain a julia expression that can be used to automatically build codes implementing equations derived with QuantumAlgebra. Every expectation value or correlator is treated as a separate array. Daggers are represented as ᴴ, which are valid identifiers that can appear in the array names. Note that expectation values and correlators are not distinguished, so it is best to have all expressions use the same kind.\njulia> julia_expression(expval_as_corrs(a'(:j)*a(:i)))\n:(aᴴ[j] * a[i] + aᴴa[j, i])\nAlso note that expressions are always treated as arrays, even if they have no indices (which gives zero-dimensional arrays). If you are working with scalar quantities exclusively, it might be useful to clean up the resulting expression (e.g., use MacroTools to remove the []).\njulia> julia_expression(expval(a'()*a()*σx()))\n:(aᴴaσˣ[])\nBy default, two-level system operators are represented by the Pauli matrices σ^xyz, and calling σp() and σm() will give results expressed through them:\njulia> σp()\n1//2 σˣ() + 1//2i σʸ()\n\njulia> σm()\n1//2 σˣ() - 1//2i σʸ()\nThis can be changed by calling QuantumAlgebra.use_σpm(true; set_preference=true/false) (where the value of set_preference determines whether this is stored permanently using Preferences.jl). In this mode, σ^+ and σ^- are the \"fundamental\" operators, and all expressions are written in terms of them. Note that mixing conventions within the same expression is not supported, so it is suggested to set this flag once at the beginning of any calculation.\njulia> QuantumAlgebra.use_σpm(true)\n\njulia> σp()\nσ⁺()\n\njulia> σx()\nσ⁺() + σ⁻()\n\njulia> σz()\n-1 + 2 σ⁺() σ⁻()","category":"page"},{"location":"#Preferences","page":"Home","title":"Preferences","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Several preferences changing the behavior of QuantumAlgebra can be set permanently (this uses Preferences.jl):","category":"page"},{"location":"","page":"Home","title":"Home","text":"\"define_default_ops\": if this is set to false (default is true), the \"default\" operators a, adag, f, fdag, σx, σy, σz, σp, σm are not defined upon import. Note that changing this value requires restarting the Julia session to take effect. The setting can be changed with QuantumAlgebra.set_define_default_ops(true/false) (which will inform you whether a restart is required) or with Preferences.set_preferences!(QuantumAlgebra,\"define_default_ops\"=>true/false).\n\"auto_normal_form\": Choose whether all expressions are automatically converted to normal form upon creation. The default is false. It can be changed for a single session with QuantumAlgebra.auto_normal_form(true/false), and can be made permanent with QuantumAlgebra.auto_normal_form(true/false; set_preference=true) or with Preferences.set_preferences!(QuantumAlgebra,\"auto_normal_form\"=>true/false). Note that this could previously be set by defining an environment variable \"QUANTUMALGEBRA_AUTO_NORMAL_FORM\", but this usage has been deprecated and will be removed in a future version.\n\"use_σpm\": Choose whether for two-level systems, the \"basic\" operators are excitation/deexcitation operators σ^ or the Pauli matrices `σ^xyz. This can be changed in a single session by calling QuantumAlgebra.use_σpm(true/false), and can be made permanent with QuantumAlgebra.use_σpm(true/false; set_preference=true) or with Preferences.set_preferences!(QuantumAlgebra,\"use_σpm\"=>true/false).","category":"page"},{"location":"#Citing","page":"Home","title":"Citing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you use QuantumAlgebra in academic work, we would appreciate a citation. See CITATION.bib for the relevant references.","category":"page"}]
}