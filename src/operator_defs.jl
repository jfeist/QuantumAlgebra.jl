using Preferences
using OrderedCollections

export QuExpr
export boson_ops, fermion_ops
export @boson_ops, @fermion_ops, @anticommuting_fermion_group
export tlspm_ops, tlsxyz_ops
export @tlspm_ops, @tlsxyz_ops
export @Pr_str, @Pc_str, ∑
export param, expval, corr
export @qindex, @boson, @fermion, @anticommuting_fermions, @tls, @tlsxyz, @qudit, @spin, @so, @qparam
export qindex, boson, fermion, anticommuting_fermions, tls, tlsxyz, qudit, gellmann, spin, so, qparam
export map_scalar_function

# compile-time options
const _DEFINE_DEFAULT_OPS = @load_preference("define_default_ops", true)

"""
    set_define_default_ops(t::Bool)

Set the preference `"define_default_ops"`, which determines whether the default operators of the
legacy interface (`a`, `adag`, `f`, `fdag`, `σx`, `σy`, `σz`, `σp`, `σm`) are defined when
QuantumAlgebra is loaded (the default is `true`). The default operators are deprecated and will be
removed in version 2.0; define operators with the macros instead (e.g., `@boson a[:]`, see
[`@boson`](@ref)).

Changing this preference requires restarting the Julia session to take effect.
"""
function set_define_default_ops(t::Bool)
    @set_preferences!("define_default_ops" => t)
    if t != _DEFINE_DEFAULT_OPS
        @info("define_default_ops setting changed to $t; restart your Julia session for this change to take effect!")
    end
end

"""
    set_quindices_type(t::String)

Deprecated: indices are now always stored in a compact form that supports an arbitrary
number of indices, so the `"quindices_type"` preference no longer has any effect.
"""
function set_quindices_type(t::String)
    Base.depwarn("`set_quindices_type` is deprecated and has no effect, indices are always stored in a compact form with an arbitrary number of indices.", :set_quindices_type; force=true)
    nothing
end


# dynamically changeable options
const _using_σpm = Ref(false)
"""
    use_σpm(t::Bool=true; set_preference=false)

Choose how the operators of two-level systems are represented: with the raising and lowering
operators ``σ^+`` and ``σ^-`` (for `true`), or with the Pauli matrices ``σ^x``, ``σ^y`` and ``σ^z``
(for `false`, the default). This applies to the operators that are created afterwards, e.g., `σ[i]`
after `@tls σ[:]` (or `σx` after `@tlsxyz σ`, which is created when it is defined), which are then
expressed in the chosen basis. If `set_preference` is `true`, the choice is stored with
`Preferences.jl` (preference `"use_σpm"`) for future sessions.
"""
function use_σpm(t::Bool=true; set_preference=false)
    _using_σpm[] = t
    set_preference && @set_preferences!("use_σpm" => t)
    nothing
end
"""
    use_σxyz(; set_preference=false)

Represent the operators of two-level systems with the Pauli matrices ``σ^x``, ``σ^y`` and ``σ^z``
(the default), equivalent to `use_σpm(false; set_preference)` (see [`use_σpm`](@ref)).
"""
use_σxyz(;set_preference=false) = use_σpm(false;set_preference)
"""
    using_σpm()

Return whether the operators of two-level systems are represented with ``σ^+`` and ``σ^-`` (see
[`use_σpm`](@ref)).
"""
using_σpm() = _using_σpm[]

const _auto_normal_form = Ref(false)
"""
    auto_normal_form(t::Bool=true; set_preference=false)

Enable (or disable, with `t = false`) automatic normal ordering, which brings every expression into
normal form when it is created (e.g., `a*a'` then directly gives `1 + a†() a()`). It is disabled by
default, so that [`normal_form`](@ref) has to be called explicitly. If `set_preference` is `true`,
the choice is stored with `Preferences.jl` (preference `"auto_normal_form"`) for future sessions.
"""
function auto_normal_form(t::Bool=true; set_preference=false)
    _auto_normal_form[] = t
    set_preference && @set_preferences!("auto_normal_form" => t)
    nothing
end
"""
    using_auto_normal_form()

Return whether automatic normal ordering is enabled (see [`auto_normal_form`](@ref)).
"""
using_auto_normal_form() = _auto_normal_form[]

const _bracket_display = Ref(false)
"""
    use_bracket_display(t::Bool=true; set_preference=false)

Select whether operators and parameters are printed in the syntax of the interface with operator
variables (see [`@boson`](@ref)), e.g., `a†[i,j]`, `g[i]` and `a` (for `true`), or with parentheses,
e.g., `a†(ij)`, `g(i)` and `a()` (for `false`, the default in v1.x). If `set_preference` is `true`, the
choice is stored with `Preferences.jl` for future sessions.
"""
function use_bracket_display(t::Bool=true; set_preference=false)
    _bracket_display[] = t
    set_preference && @set_preferences!("bracket_display" => t)
    nothing
end
"""
    using_bracket_display()

Return whether operators are printed in the syntax with brackets (see [`use_bracket_display`](@ref)).
"""
using_bracket_display() = _bracket_display[]

const IndexInt = Int32

@concrete struct QuIndex
    # sym is the index character, apart from two special values:
    # '\0' for integer indices, which ensures they are ordered before any other indices
    # '#' for sum indices, which ensures that they are ordered before letters etc, and
    sym::Char
    # num is integer subindex, special value typemin(IndexInt) means no subindex
    num::IndexInt
    QuIndex(sym::Char,num::Integer=typemin(IndexInt)) = new(sym,num)
end

"""
    QuIndex(i)

An index of an operator or parameter: an integer (`QuIndex(2)`), a symbolic index, given by a single
unicode character with an optional integer subindex (`QuIndex(:i)`, or `QuIndex(:j_2)` and
`QuIndex("j_2")`, printed as `j₂`), or a summation index (see [`sumindex`](@ref)).

Indices are usually not created directly: operators and parameters accept integers and symbols as
indices (e.g., `a[:i,2]`), or index variables defined with [`@qindex`](@ref) (e.g., `a[i,j[2]]`).

See also [`isintindex`](@ref), [`issumindex`](@ref).
"""
QuIndex
QuIndex(ii::QuIndex) = ii
QuIndex(ii::Symbol) = QuIndex(string(ii))
function QuIndex(ii::String)
    s = split(ii,"_")
    length(s) <= 2 || throw(ArgumentError("index can have at most one subindex, got $ii with subindices $(s[2:end])"))
    length(s[1]) == 1 || throw(ArgumentError("index names must be single character, got $ii."))
    sym = s[1][1]
    # we do not allow the ascii control characters as symbolic indices to ensure ordering
    sym > '#' || throw(ArgumentError("Symbolic index name cannot be character less than '#' (codepoint 35), passed Int(sym) = $(Int(sym))."))
    if length(s)==2
        num = try
            parse(IndexInt,s[2])
        catch
            throw(ArgumentError("Only integer literals are supported as index subscripts. Got \"$(s[2])\"."))
        end
        QuIndex(sym,num)
    else
        QuIndex(sym)
    end
end
QuIndex(ii::Integer) = QuIndex('\0',ii)

"""
    sumindex(i)

Return the summation index number `i` (printed as `#ᵢ`). Sums (see [`∑`](@ref)) are stored with
summation indices for the indices that are summed over, so that they cannot clash with other indices.
"""
sumindex(ii) = QuIndex('#',ii)
# by convention, these should never be present in the input or output of any function,
# but just as temporary variables that cannot conflict with any other within a computation
tmpindex(ii) = QuIndex('!',ii)

"""
    isintindex(i::QuIndex)

Check if the index `i` is an integer index (created with `QuIndex(i::Integer)`).
"""
isintindex(ii::QuIndex) = ii.sym=='\0'

"""
    issumindex(i::QuIndex)

Check if the index `i` is a summation index (created with `sumindex(i)`).
"""
issumindex(ii::QuIndex) = ii.sym=='#'
const NoIndex = QuIndex(typemin(IndexInt))
isnoindex(ii::QuIndex) = ii == NoIndex

@inline Base.isless(i1::QuIndex,i2::QuIndex) = isless((i1.sym,i1.num),(i2.sym,i2.num))

"`make_indices(inds...)`: convert the arguments (or a single tuple or vector of them) to a vector of `QuIndex`"
make_indices(inds...) = QuIndex[QuIndex(i) for i in inds]
make_indices(inds::Union{AbstractVector,Tuple}) = make_indices(inds...)
make_indices(inds::AbstractVector{QuIndex}) = inds
# all indices are "assigned" (i.e., not NoIndex) since they are stored with variable length
assignedinds(inds) = inds

# Operator names are stored as integer ids into a table of names. The table is never modified, but
# replaced (under a lock) when a name is added, so that it can be read from several threads at once.
# (This is a mutable struct with constant fields so that it is stored as a reference, which can be
# loaded atomically by the hardware, while larger immutable structs would need a lock.)
mutable struct _NameTableData
    const ids::Dict{Symbol,IndexInt}
    const syms::Vector{Symbol}
    # hashes of the names that only depend on their content (FNV-1a), used to number sum indices canonically
    const stablehash::Vector{UInt64}
    # names are ordered alphabetically, rank[i] is the position of the name with id i in that order
    const rank::Vector{IndexInt}
end
const _NameTable = _AtomicRef(_NameTableData(Dict{Symbol,IndexInt}(), Symbol[], UInt64[], IndexInt[]))
const _NameTableLock = ReentrantLock()
_fnv1a(s::AbstractString) = foldl((h,b) -> (h ⊻ b) * 0x00000100000001b3, codeunits(s); init=0xcbf29ce484222325)

function _name_id(name::Symbol)
    i = get(_NameTable[].ids, name, zero(IndexInt))
    iszero(i) || return i
    @lock _NameTableLock begin
        old = _NameTable[]
        i = get(old.ids, name, zero(IndexInt))
        iszero(i) || return i
        n = length(old.syms) + 1
        n < 2^24 || error("too many operator names")
        syms = push!(copy(old.syms), name)
        ranks = similar(old.rank, n)
        ranks[sortperm(syms)] .= 1:n
        ids = copy(old.ids)
        ids[name] = n
        _NameTable[] = _NameTableData(ids, syms, push!(copy(old.stablehash), _fnv1a(string(name))), ranks)
        IndexInt(n)
    end
end

@concrete struct QuOpName
    i::IndexInt
    QuOpName(name::Symbol) = new(_name_id(name))
    QuOpName(name::QuOpName) = name
    # create from the id of an existing name
    QuOpName(i::IndexInt, ::Val{:id}) = new(i)
end

"""
    QuOpName(name::Symbol)

The name of an operator or parameter. Names are stored as integer ids into a global table of names,
so that they can be compared and hashed quickly, and they are ordered alphabetically.
"""
QuOpName
sym(ind::QuOpName) = _NameTable[].syms[ind.i]
Base.print(io::IO, ind::QuOpName) = print(io, sym(ind))
_namestablehash(ind::QuOpName) = @inbounds _NameTable[].stablehash[ind.i]
# (read the table only once, since another thread can replace it in between. Adding names does not
# change the relative order of existing names, so it does not matter which version is used.)
Base.cmp(i1::QuOpName,i2::QuOpName) = (rank = _NameTable[].rank; @inbounds cmp(rank[i1.i],rank[i2.i]))
Base.isless(i1::QuOpName,i2::QuOpName) = cmp(i1,i2) < 0
const NoName = QuOpName(Symbol())

# the enum also directly defines a natural ordering,so choose this directly how we later want it
# start the counting at 1 so we can index into the tuples defined below with Int(BaseOpType)
# (for qudits and so(n) generators, the symbols are determined from the levels or the indices a,b of L_ab)
"""
    BaseOpType

The type of a [`BaseOperator`](@ref), e.g., `BosonCreate_` or `TLSz_`. Operators are ordered by their
type (in the order of the values of `BaseOpType`), then by name, additional data and indices. Normal
ordering brings products of operators into this order as far as their algebra allows, so that, e.g.,
creation operators end up to the left of annihilation operators.
"""
@enum BaseOpType        BosonCreate_=1 FermionCreate_ TLSCreate_ QuditCreate_ SpinPlus_ TLSx_  TLSy_  TLSz_  SpinMinus_ SpinZ_ SOGen_ QuditDestroy_ TLSDestroy_ FermionDestroy_ BosonDestroy_
const BaseOpType_adj   = (BosonDestroy_, FermionDestroy_, TLSDestroy_, QuditDestroy_, SpinMinus_, TLSx_, TLSy_, TLSz_, SpinPlus_, SpinZ_, SOGen_, QuditCreate_, TLSCreate_, FermionCreate_, BosonCreate_)
const BaseOpType_sym   = ("†", "†", "⁺", "", "⁺", "ˣ", "ʸ", "ᶻ", "⁻", "ᶻ", "", "", "⁻", "", "")
const BaseOpType_expr  = ("ᴴ", "ᴴ", "⁺", "", "⁺", "ˣ", "ʸ", "ᶻ", "⁻", "ᶻ", "", "", "⁻", "", "")
const BaseOpType_latex = ("^{\\dagger}", "^{\\dagger}", "^{+}", "", "^{+}", "^{x}", "^{y}", "^{z}", "^{-}", "^{z}", "", "", "^{-}", "", "")

# a struct representing a delta function with unequal indices
@concrete struct δ
    iA::QuIndex
    iB::QuIndex
    δ(iA,iB) = (@assert !isnoindex(iA) && !isnoindex(iB); new(iA,iB))
end

"""
    δ(iA, iB)

The Kronecker delta ``δ_{iA,iB}`` of two indices in a term. It is only stored for two different
indices that can be equal (i.e., not for two different integers). Kronecker deltas appear, e.g., when
normal ordering, such as in `normal_form(b[:i]*b'[:j])`, which gives `δ(ij) + b†(j) b(i)`.

`δ(Ainds, Binds)`, with two vectors of indices, returns the Kronecker deltas for the pairwise
equality of `Ainds` and `Binds` (only for pairs of different indices), or `nothing` if not all of
them can be equal.
"""
δ
function δ(Ainds::AbstractVector{QuIndex},Binds::AbstractVector{QuIndex})
    length(Ainds) == length(Binds) || return nothing

    res = Vector{δ}(undef,length(Ainds))
    jj = 0
    for (iA,iB) in zip(Ainds,Binds)
        if iA != iB
            (isnoindex(iA) || isnoindex(iB)) && return nothing
            # if they are integer indices and different, the result is zero
            isintindex(iA) && isintindex(iB) && return nothing
            res[jj+=1] = iB<iA ? δ(iB,iA) : δ(iA,iB)
        end
    end
    resize!(res,jj)
    res
end

include("term_storage.jl")

"""
    BaseOperator

A single operator in a term, given by its type (see [`BaseOpType`](@ref)), its name, its indices,
and additional data for some kinds of operators (e.g., the levels of a qudit transition, or the
value of a spin). Operators are usually not created directly, but with the macros that define
operators (e.g., [`@boson`](@ref)), which give `QuExpr` objects. For bosons, fermions and two-level
systems, `BaseOperator`s can also be created with `BosonDestroy(name, inds...)`,
`BosonCreate(name, inds...)`, `FermionDestroy`, `FermionCreate`, `TLSDestroy`, `TLSCreate`, `TLSx`,
`TLSy` and `TLSz` (e.g., `BosonCreate(:a, :i)` is ``a_i^†``).
"""
BaseOperator

for (op,desc) in (
    (:BosonDestroy,"bosonic annihilation"),
    (:BosonCreate,"bosonic creation"),
    (:FermionDestroy,"fermionic annihilation"),
    (:FermionCreate,"fermionic creation"),
    (:TLSDestroy,"TLS annihilation"),
    (:TLSCreate,"TLS creation"),
    (:TLSx,"TLS x"),
    (:TLSy,"TLS y"),
    (:TLSz,"TLS z"))
    @eval begin
        "`$($op)(name,inds)`: represent $($desc) operator ``name_{inds}``"
        $op(name,inds...) = BaseOperator($(Symbol(op,:_)),name,make_indices(inds...))
    end
end

"""
    ExpVal(ops::BaseOpProduct)

The expectation value ``⟨ops⟩`` of a product of operators in a term. Expectation values are created
with [`expval`](@ref), e.g., `expval(a'*a)` gives `⟨a†() a()⟩`. As the first argument of
[`heisenberg_eom_system`](@ref), `ExpVal` selects equations for expectation values.

See also [`Corr`](@ref).
"""
ExpVal

"""
    Corr(ops::BaseOpProduct)

The correlator (or cumulant) ``⟨ops⟩_c`` of a product of operators in a term, i.e., the connected
part of the expectation value ``⟨ops⟩`` (see [`expval_as_corrs`](@ref)). Correlators are created with
[`corr`](@ref), e.g., `corr(a'*a)` gives `⟨a†() a()⟩c`. As the first argument of
[`heisenberg_eom_system`](@ref), `Corr` selects equations for correlators.

See also [`ExpVal`](@ref).
"""
Corr

"""
    Param

A parameter in a term, given by its name, its state (`'r'` for a real parameter, `'n'` for a complex
parameter, and `'c'` for the complex conjugate of a complex parameter) and its indices. Parameters
are defined with [`@qparam`](@ref) or [`qparam`](@ref) (or, in the legacy interface, with
[`param`](@ref) and the string macros [`@Pr_str`](@ref) and [`@Pc_str`](@ref)).
"""
Param


"""
    QuExpr

A quantum expression: a sum of terms ([`QuTerm`](@ref)) with numerical coefficients, which are
stored in the dictionary `A.terms` (from terms to coefficients). The coefficients can be numbers of
any type, including symbolic numbers (e.g., from Symbolics.jl or SymPy.jl), and numbers of different
types can appear in the same expression.

Expressions are built from operators and parameters (see [`@boson`](@ref) and [`@qparam`](@ref)),
which are `QuExpr`s themselves, with `+`, `-`, `*`, `^`, the adjoint `'`, sums ([`∑`](@ref)), and
functions such as [`normal_form`](@ref), [`comm`](@ref) or [`expval`](@ref). `QuExpr()` (or
`zero(QuExpr)`) is the empty expression, i.e., zero, and `QuExpr(x)` is the constant `x` (e.g.,
`QuExpr(1) == one(QuExpr)`).
"""
struct QuExpr
    # A QuantumAlgebra Expression is saved as a Dictionary of QuTerms with scalar prefactors
    terms::Dict{QuTerm,Number}
    QuExpr() = new(Dict{QuTerm,Number}())
    QuExpr(terms::Dict{QuTerm,Number}) = new(terms)
end

function QuExpr(itr)
    A = QuExpr()
    for (t,s) in itr
        _withnum(s -> _add_with_auto_order!(A,t,simplify_number(s)), s)
    end
    A
end
QuExpr(A::Union{BaseOperator,Param,Corr,ExpVal}) = QuExpr(QuTerm(A))
QuExpr(A::QuTerm) = QuExpr(((A,1),))
QuExpr(s::Number) = QuExpr(((QuTerm(),s),))
Base.isempty(A::QuExpr) = isempty(A.terms)
Base.copy(A::QuExpr) = QuExpr(copy(A.terms))

_add_with_auto_order!(A::QuExpr,B::QuTerm,sB) = using_auto_normal_form() ? _add_with_normal_order!(A,B,sB) : _add_sum_term!(A,B,sB)

# (inlined into the callers, and calls the implementation with the type of sB known, see `_withnum`)
@inline _add_sum_term!(A::QuExpr,oB::QuTerm,sB) = _withnum(sB -> _add_sum_term_impl!(A,oB,sB), sB)
function _add_sum_term_impl!(A::QuExpr,oB::QuTerm,sB)
    iszero(sB) && return A
    sold = get(A.terms,oB,nothing)
    if sold === nothing
        A.terms[oB] = simplify_number(sB)
    else
        _withnum(sold -> _add_to_term!(A,oB,sB,sold), sold)
    end
    A
end
function _add_to_term!(A::QuExpr,oB::QuTerm,sB,sold)
    snew = sB + sold
    if iszero(snew)
        delete!(A.terms,oB)
    else
        A.terms[oB] = simplify_number(snew)
    end
    nothing
end
_map_quexpr_ops(f,A::QuExpr) = QuExpr((f(t),s) for (t,s) in A.terms)

"""
    map_scalar_function(f, A::QuExpr)

Return the expression `A` with the function `f` applied to the coefficient of each term, e.g., to
round numerical coefficients or to simplify symbolic ones.

# Examples
```julia
julia> @boson a;

julia> map_scalar_function(x -> round(x, digits=3), a/3 + a'/7)
0.143 a†() + 0.333 a()
```
"""
map_scalar_function(f,A::QuExpr) = QuExpr((t,f(s)) for (t,s) in A.terms)

"""
    QuExprConstructor

A callable object that creates operators (or parameters) as `QuExpr`s, e.g., `a` in `a(:i)` for the
operators of the legacy interface (see [`boson_ops`](@ref)). Its adjoint `a'` creates the adjoint
operators. The variables of operators with indices (e.g., `a` after `@boson a[:]`, see
[`QuExprCstrctArrAlias`](@ref)) call it with the indices.
"""
struct QuExprConstructor{F,Fdag}
    name::String
    f::F
    namedag::String
    fdag::Fdag
    QuExprConstructor(name,f,namedag,fdag) = new{Core.Typeof(f),Core.Typeof(fdag)}(name,f,namedag,fdag)
    QuExprConstructor(name,f) = QuExprConstructor(name,f,name,f)
end
(q::QuExprConstructor)(args...)::QuExpr = q.f(args...)
Base.adjoint(q::QuExprConstructor) = QuExprConstructor(q.namedag,q.fdag,q.name,q.f)
Base.print(io::IO, q::QuExprConstructor) = print(io, q.name, " (QuExpr constructor)")
Base.show(io::IO, q::QuExprConstructor) = print(io, q)
Base.show(io::IO, ::MIME"text/plain", q::QuExprConstructor) = show(io, q)

wrapdagdeprecated(f) = (args...) -> (Base.depwarn("`xdag()` constructors are deprecated, use `x'()` instead", :xdag; force=true); f(args...))

#################################################################
## "external" functions that always construct QuExpr           ##
#################################################################

"""
    boson_ops(name::Symbol)

Return the constructors `(a, adag)` for the bosonic operators with the given name in the legacy
interface: `a(inds...)` is the annihilation operator with these indices, and `a'(inds...)` the
creation operator (`adag` is deprecated, use `a'` instead). See also [`@boson_ops`](@ref).

!!! note "Legacy interface"
    This function is part of the legacy interface, which is deprecated and will be removed in
    version 2.0. Use [`boson`](@ref) or [`@boson`](@ref) instead.
"""
function boson_ops(name::Symbol)
    op_name = QuOpName(name)
    c = QuExprConstructor(string(name),     (args...) -> QuExpr(BosonDestroy(op_name,args...)),
                          string(name,"†"), (args...) -> QuExpr(BosonCreate( op_name,args...)))
    cdag = QuExprConstructor(c.namedag, wrapdagdeprecated(c.fdag), c.name, wrapdagdeprecated(c.f))
    c, cdag
end

"""
    fermion_ops(name::Symbol)

Return the constructors `(f, fdag)` for the fermionic operators with the given name in the legacy
interface: `f(inds...)` is the annihilation operator with these indices, and `f'(inds...)` the
creation operator (`fdag` is deprecated, use `f'` instead). See also [`@fermion_ops`](@ref).

!!! note "Legacy interface"
    This function is part of the legacy interface, which is deprecated and will be removed in
    version 2.0. Use [`fermion`](@ref) or [`@fermion`](@ref) instead.
"""
function fermion_ops(name::Symbol)
    op_name = QuOpName(name)
    c = QuExprConstructor(string(name),     (args...) -> QuExpr(FermionDestroy(op_name,args...)),
                          string(name,"†"), (args...) -> QuExpr(FermionCreate( op_name,args...)))
    cdag = QuExprConstructor(c.namedag, wrapdagdeprecated(c.fdag), c.name, wrapdagdeprecated(c.f))
    c, cdag
end

"""
    tlspm_ops(name::Symbol)

Return the constructors `(σm, σp)` for the lowering and raising operators ``σ^-`` and ``σ^+`` of the
two-level system `name` in the legacy interface. Depending on [`use_σpm`](@ref) (when the operators
are created), they are represented as jump operators or through the Pauli matrices.

!!! note "Legacy interface"
    This function is part of the legacy interface, which is deprecated and will be removed in
    version 2.0. Use [`tls`](@ref) or [`@tls`](@ref) instead.
"""
function tlspm_ops(name::Symbol)
    op_name = QuOpName(name)
    namem = string(name,"m")
    namep = string(name,"p")
    fm = (args...) -> using_σpm() ? QuExpr(TLSDestroy(op_name,args...)) : QuExpr((QuTerm(TLSx(op_name,args...))=>1//2, QuTerm(TLSy(op_name,args...))=>-1im//2))
    fp = (args...) -> using_σpm() ? QuExpr(TLSCreate( op_name,args...)) : QuExpr((QuTerm(TLSx(op_name,args...))=>1//2, QuTerm(TLSy(op_name,args...))=>1im//2))
    (QuExprConstructor(namem, fm, namep, fp), QuExprConstructor(namep, fp, namem, fm))
end

"""
    tlsxyz_ops(name::Symbol)

Return the constructors `(σx, σy, σz)` for the Pauli matrices of the two-level system `name` in the
legacy interface. Depending on [`use_σpm`](@ref) (when the operators are created), they are
represented as Pauli matrices or through the jump operators ``σ^±``.

!!! note "Legacy interface"
    This function is part of the legacy interface, which is deprecated and will be removed in
    version 2.0. Use [`tlsxyz`](@ref) or [`@tlsxyz`](@ref) instead.
"""
function tlsxyz_ops(name::Symbol)
    op_name = QuOpName(name)
    namex = string(name,"x")
    namey = string(name,"y")
    namez = string(name,"z")
    tlsx = QuExprConstructor(namex, (args...) -> using_σpm() ? QuExpr((QuTerm(TLSDestroy(op_name,args...))=>1, QuTerm(TLSCreate(op_name,args...))=>1)) : QuExpr(TLSx(op_name,args...)))
    tlsy = QuExprConstructor(namey, (args...) -> using_σpm() ? QuExpr((QuTerm(TLSDestroy(op_name,args...))=>1im, QuTerm(TLSCreate(op_name,args...))=>-1im)) : QuExpr(TLSy(op_name,args...)))
    tlsz = QuExprConstructor(namez, (args...) -> using_σpm() ? QuExpr((QuTerm(BaseOpProduct([TLSCreate(op_name,args...),TLSDestroy(op_name,args...)]))=>2, QuTerm()=>-1)) : QuExpr(TLSz(op_name,args...)))
    tlsx, tlsy, tlsz
end

# constructor for the transition operators of a qudit, q(i, j, inds...) = |i⟩⟨j| (see `qudit`)
struct QuditOps
    name::QuOpName
    d::Int
end
function _qudit_ops(name::Symbol, d::Integer)
    2 ≤ d ≤ 255 || throw(ArgumentError("the number of levels of a qudit has to be between 2 and 255, got $d"))
    QuditOps(QuOpName(name), d)
end
function (q::QuditOps)(i::Integer, j::Integer, inds...)
    (1 ≤ i ≤ q.d && 1 ≤ j ≤ q.d) || throw(ArgumentError("levels have to be between 1 and $(q.d), got ($i,$j)"))
    ii = make_indices(inds...)
    # |i⟩⟨j| = |i⟩⟨1|1⟩⟨j|
    ij(i,j) = i == 1 ? QuTerm(qudit_op(q.name,q.d,1,j,ii)) : j == 1 ? QuTerm(qudit_op(q.name,q.d,i,1,ii)) :
                       QuTerm(BaseOpProduct([qudit_op(q.name,q.d,i,1,ii), qudit_op(q.name,q.d,1,j,ii)]))
    if (i,j) == (1,1)
        # |1⟩⟨1| = 1 - ∑_{k≥2} |k⟩⟨k|
        QuExpr(Iterators.flatten((((QuTerm(),1),), ((ij(k,k),-1) for k=2:q.d))))
    else
        QuExpr(((ij(i,j),1),))
    end
end
# q(M, inds...) = ∑_ij M_ij |i⟩⟨j|, where the diagonal is ∑_i M_ii |i⟩⟨i| = M_11 + ∑_{k≥2} (M_kk - M_11) |k⟩⟨k|
# (so that every operator only appears once). The entries can be numbers or expressions (e.g., parameters).
function (q::QuditOps)(M::AbstractMatrix, inds...)
    Base.require_one_based_indexing(M)
    size(M) == (q.d, q.d) || throw(DimensionMismatch("expected a $(q.d)×$(q.d) matrix for a $(q.d)-level qudit, got size $(size(M))"))
    M11 = M[1,1]
    terms = QuExpr[]
    iszero(M11) || push!(terms, M11*QuExpr(1))
    for j = 1:q.d, i = 1:q.d
        (i,j) == (1,1) && continue
        c = i == j ? M[i,i] - M11 : M[i,j]
        iszero(c) || push!(terms, c*q(i,j,inds...))
    end
    # (sum adds pairwise, which avoids copying large intermediate sums)
    isempty(terms) ? QuExpr() : sum(terms)
end

# constructors for the generalized Gell-Mann matrices of a qudit (see `gellmann`)
function _gellmann_ops(q::QuditOps, sqrt)
    cstrs = QuExprConstructor[]
    λname() = "λ$(subscript(length(cstrs)+1))"
    for k = 2:q.d
        for j = 1:k-1
            push!(cstrs, QuExprConstructor(λname(), (inds...) -> q(j,k,inds...) + q(k,j,inds...)))
            push!(cstrs, QuExprConstructor(λname(), (inds...) -> -1im*q(j,k,inds...) + 1im*q(k,j,inds...)))
        end
        l = k-1
        c = _sqrt_rational(2//(l*(l+1)), sqrt)
        push!(cstrs, QuExprConstructor(λname(), (inds...) -> c*(sum(q(m,m,inds...) for m=1:l) - l*q(l+1,l+1,inds...))))
    end
    cstrs
end
# the square root of the rational number x: exact if x is the square of a rational number, and sqrt(x) otherwise
function _sqrt_rational(x::Rational, sqrt)
    n, d = isqrt(numerator(x)), isqrt(denominator(x))
    (n^2 == numerator(x) && d^2 == denominator(x)) || return sqrt(x)
    d == 1 ? n : n//d
end

"`@boson_ops name`: define the function `\$name` that creates bosonic annihilation operators with name `name`, e.g., `\$name(:i)` (and the deprecated `\$(name)dag`, use `\$(name)'` instead). Legacy interface (deprecated, will be removed in version 2.0), use [`@boson`](@ref) instead."
macro boson_ops(name)
    :( ($(esc(name)), $(esc(Symbol(name,:dag)))) = boson_ops($(Meta.quot(name))) )
end

"`@fermion_ops name`: define the function `\$name` that creates fermionic annihilation operators with name `name`, e.g., `\$name(:i)` (and the deprecated `\$(name)dag`, use `\$(name)'` instead). Legacy interface (deprecated, will be removed in version 2.0), use [`@fermion`](@ref) instead."
macro fermion_ops(name)
    :( ($(esc(name)), $(esc(Symbol(name,:dag)))) = fermion_ops($(Meta.quot(name))) )
end

"`@anticommuting_fermion_group name1 name2 ...`: define the functions `name1`, `name2`, … that create fermionic annihilation operators that mutually anticommute. Legacy interface (deprecated, will be removed in version 2.0), use [`@anticommuting_fermions`](@ref) instead."
macro anticommuting_fermion_group(names...)
    # start the groupname (which is the "internal" species name) by concatenating all names,
    # to ensure reasonable sorting relative to other fermionic species
    code = quote
        groupname = Symbol($names...,gensym())
        add_groupaliases(groupname,$names)
        ann, cre = fermion_ops(groupname)
    end
    for (ii,name) in enumerate(names)
        push!(code.args,:( $(esc(name)) = QuExprConstructor($(string(name)), (args...) -> ann($ii,args...),
                                                            $(string(name,"†")), (args...) -> ann'($ii,args... )) ))
        push!(code.args,:( $(esc(Symbol(name,:dag))) = ($(esc(name)))' ))
    end
    push!(code.args, :( nothing ))
    code
end

# functions for constructing Pauli operators depending on whether we use (σ+,σ-) or (σx,σy,σz) as the "basic" operators

"`@tlspm_ops name`: define the functions `\$(name)m` and `\$(name)p` that create the lowering and raising operators ``σ^-`` and ``σ^+`` of the two-level system `name` (see [`tlspm_ops`](@ref)). Legacy interface (deprecated, will be removed in version 2.0), use [`@tls`](@ref) instead."
macro tlspm_ops(name)
    :( ($(esc(Symbol(name,:m))), $(esc(Symbol(name,:p)))) = tlspm_ops($(Meta.quot(name))) )
end

"`@tlsxyz_ops name`: define the functions `\$(name)x`, `\$(name)y`, and `\$(name)z` that create the Pauli matrices of the two-level system `name` (see [`tlsxyz_ops`](@ref)). Legacy interface (deprecated, will be removed in version 2.0), use [`@tlsxyz`](@ref) instead."
macro tlsxyz_ops(name)
    :( ($(esc(Symbol(name,:x))), $(esc(Symbol(name,:y))), $(esc(Symbol(name,:z)))) = tlsxyz_ops($(Meta.quot(name))) )
end

# constructors (J₊, J₋, J_z) for the operators of a spin (see `spin`)
function _spin_ops(name::Symbol, j=nothing)
    op_name = QuOpName(name)
    extra = _spin_extra(j)
    cstr(t) = (args...) -> QuExpr(BaseOperator(t, op_name, make_indices(args...), extra))
    Jp = QuExprConstructor(string(name,"⁺"), cstr(SpinPlus_), string(name,"⁻"), cstr(SpinMinus_))
    (Jp, Jp', QuExprConstructor(string(name,"ᶻ"), cstr(SpinZ_)))
end

# constructor for the generators of so(n), L(a, b, inds...) = L_ab (see `so`)
struct SOOps
    name::QuOpName
    n::Int
end
function _so_ops(name::Symbol, n::Integer)
    2 ≤ n ≤ 255 || throw(ArgumentError("so(n) is only supported for 2 ≤ n ≤ 255, got n = $n"))
    SOOps(QuOpName(name), n)
end
function (L::SOOps)(a::Integer, b::Integer, inds...)
    (1 ≤ a ≤ L.n && 1 ≤ b ≤ L.n) || throw(ArgumentError("the indices of so($(L.n)) generators have to be between 1 and $(L.n), got ($a,$b)"))
    # L_aa = 0, and L_ba = -L_ab
    a == b && return QuExpr()
    A = QuExpr(so_op(L.name, L.n, min(a,b), max(a,b), make_indices(inds...)))
    a < b ? A : -A
end

"""
    QuExprCstrctArrAlias

The variable of an operator or parameter with indices, e.g., `a` after `@boson a[:,:]` (or
`a = boson(:a; ndims=2)`): indexing it with the declared number of indices gives the operator with
these indices, `a[i,j]`, and its adjoint `a'` is the variable of the adjoint operator, `a'[i,j]`.
It is also used for indices with a subindex, e.g., `j` after `@qindex j[:]` (with `j[2]` giving `j₂`).
"""
struct QuExprCstrctArrAlias{T}
    constructor::T
    numdims::Int
end
function Base.getindex(q::QuExprCstrctArrAlias,is...)
    length(is) == q.numdims || throw(ArgumentError("Expected $(q.numdims) indices, got $(length(is))."))
    q.constructor(is...)
end
Base.adjoint(q::QuExprCstrctArrAlias) = QuExprCstrctArrAlias(q.constructor',q.numdims)
_aliasname(c::QuExprConstructor) = c.name
_aliasname(c::Base.Fix1{Type{QuIndex}}) = string(c.x)
_aliasname(c) = string(c)
Base.print(io::IO, q::QuExprCstrctArrAlias) = print(io, _aliasname(q.constructor), "[", join(fill(":",q.numdims),","), "]")
Base.show(io::IO, q::QuExprCstrctArrAlias) = print(io, q)

function parse_name_array(expr)
    if expr isa Symbol
        expr, 0
    elseif expr isa Expr && expr.head == :ref
        name, args = expr.args[1], expr.args[2:end]
        name isa Symbol || throw(ArgumentError("Expected symbol, got $(name)."))
        all(args .== :(:)) || throw(ArgumentError("Expected all indices to be `:`, got $(args)."))
        name, length(args)
    else
        throw(ArgumentError("Expected symbol or array statement, got $(expr)."))
    end
end

"""
    qindex(name::Symbol; ndims=0)

Return the index `name` (the function version of [`@qindex`](@ref)): with `ndims=0`, the index
itself, and with `ndims=1`, an object `j` such that `j[n]` is the index with name `j` and numeric
subindex `n` (printed as `jₙ`).
"""
function qindex(name::Symbol; ndims::Integer=0)
    0 ≤ ndims ≤ 1 || throw(ArgumentError("indices can have at most one subscript, got ndims = $ndims"))
    ind = QuIndex(name)
    iszero(ndims) ? ind : QuExprCstrctArrAlias(Base.Fix1(QuIndex,ind.sym), 1)
end

"""
    @qindex i j[:] ...

Define the indices `i`, `j`, … (see [`qindex`](@ref)), where an index declared as `j[:]` has a
numeric subindex, `j[2]` (printed as `j₂`). Indices can then be used directly, e.g., `a[i,j[2]]`
(instead of `a[:i,:j_2]`).
"""
macro qindex(exprs...)
    names = [first(parse_name_array(expr)) for expr in exprs]
    calls = [:( qindex($(QuoteNode(name)); ndims=$(last(parse_name_array(expr)))) ) for (name,expr) in zip(names,exprs)]
    :( ($(esc.(names)...),) = ($(calls...),) )
end

# the variable name, the name of the operators (as code) and the number of indices in a declaration `x`,
# `x[:,:]`, or `x = name` (or `x[:,:] = name`), where `name` (a string or symbol) is the name that is printed
function parse_decl(expr)
    if expr isa Expr && expr.head == :(=) && length(expr.args) == 2
        var, ndims = parse_name_array(expr.args[1])
        var, esc(expr.args[2]), ndims
    else
        var, ndims = parse_name_array(expr)
        var, QuoteNode(var), ndims
    end
end
# whether an expression is a string literal (also raw"…" or with interpolation)
_isstringliteral(x) = x isa String || (x isa Expr && (x.head == :string || (x.head == :macrocall && x.args[1] == Symbol("@raw_str"))))

# the variable for the operators created by the constructor c: the operator itself (for ndims = 0), or an
# object that gives the operators when indexed with ndims indices, e.g., a[i,j]
function _opvariable(c, ndims::Integer)
    ndims ≥ 0 || throw(ArgumentError("the number of indices cannot be negative, got ndims = $ndims"))
    iszero(ndims) ? c() : QuExprCstrctArrAlias(c, Int(ndims))
end

# code for a macro `@name x y[:,:] z = "name" ...` that calls `f(:x, args...; ndims)` for each declaration (with
# the name that is given, if any) and assigns the result to `x` (or to `x\$p` for each p in postfixes, if f
# returns several variables). The macro returns the variables (as a tuple for several declarations).
function _opvariable_macro(f, exprs, args...; postfixes=nothing)
    isempty(exprs) && throw(ArgumentError("no names given to define with $f"))
    code = Expr(:block)
    results = Any[]
    for expr in exprs
        var, name, ndims = parse_decl(expr)
        call = :( $f($name, $(map(esc,args)...); ndims=$ndims) )
        if postfixes === nothing
            push!(code.args, :( $(esc(var)) = $call ))
            push!(results, esc(var))
        else
            vars = [esc(Symbol(var,p)) for p in postfixes]
            push!(code.args, :( ($(vars...),) = $call ))
            push!(results, :( ($(vars...),) ))
        end
    end
    push!(code.args, length(results) == 1 ? only(results) : :( ($(results...),) ))
    code
end

# names of operators and parameters can be given as symbols or strings
const _Name = Union{Symbol,AbstractString}

"""
    boson(name; ndims=0)

Return the bosonic annihilation operator `name` (the function version of [`@boson`](@ref)), where
`name` is a `Symbol` or a string (which can be any name, e.g., `"c↑"`). With `ndims=0`, this is the
operator itself, and otherwise an object `a` such that `a[i,j,…]` (with `ndims` indices) is the
operator with these indices. The creation operator is obtained as `a'` (or `a'[i,j,…]`).
"""
boson(name::_Name; ndims::Integer=0) = _opvariable(first(boson_ops(Symbol(name))), ndims)

"""
    fermion(name; ndims=0)

Return the fermionic annihilation operator `name` (the function version of [`@fermion`](@ref)), see
[`boson`](@ref). Fermions with different names commute, see [`anticommuting_fermions`](@ref) for
fermions that anticommute.
"""
fermion(name::_Name; ndims::Integer=0) = _opvariable(first(fermion_ops(Symbol(name))), ndims)

"""
    tls(name; ndims=0)

Return the lowering operator ``σ^-`` of the two-level system `name` (the function version of
[`@tls`](@ref)), see [`boson`](@ref). The raising operator ``σ^+`` is obtained as `σ'`. The operators are
represented as jump operators or through Pauli matrices depending on [`use_σpm`](@ref) (at the time
the operator is created).
"""
tls(name::_Name; ndims::Integer=0) = _opvariable(first(tlspm_ops(Symbol(name))), ndims)

"""
    tlsxyz(name; ndims=0)

Return the Pauli matrices ``(σ^x, σ^y, σ^z)`` of the two-level system `name` (the function version of
[`@tlsxyz`](@ref)), see [`boson`](@ref) and [`tls`](@ref). With an empty name, the Pauli matrices are
printed as `X`, `Y` and `Z`, as in quantum computing (e.g., `X, Y, Z = tlsxyz(""; ndims=1)`).
"""
tlsxyz(name::_Name; ndims::Integer=0) = map(c -> _opvariable(c, ndims), tlsxyz_ops(Symbol(name)))

"""
    anticommuting_fermions(names...; ndims=0)

Return fermionic annihilation operators with the given names that mutually anticommute (the function
version of [`@anticommuting_fermions`](@ref)), see [`fermion`](@ref).
"""
function anticommuting_fermions(names::_Name...; ndims::Integer=0)
    names = Symbol.(names)
    # start the groupname (which is the "internal" species name) by concatenating all names,
    # to ensure reasonable sorting relative to other fermionic species
    groupname = Symbol(names...,gensym())
    add_groupaliases(groupname,names)
    ann = first(fermion_ops(groupname))
    map(Tuple(enumerate(names))) do (ii,name)
        c = QuExprConstructor(string(name), (args...) -> ann(ii,args...), string(name,"†"), (args...) -> ann'(ii,args...))
        _opvariable(c, ndims)
    end
end

"""
    @boson a
    @boson a[:,:]
    @boson a b[:] ...
    @boson au="a↑" ad[:]="a↓" ...

Define the variable `a` for a bosonic annihilation operator (see [`boson`](@ref)). With `a[:,:]` (with
any number of `:`), the operator has as many indices, given as `a[i,j]` (with indices declared by
[`@qindex`](@ref), or with symbols and integers, e.g., `a[:i,2]`). The creation operator is `a'`
(or `a'[i,j]`). Several operators can be defined at once, and with `au = "a↑"`, the variable `au` is
an operator with the name `a↑` (also for the other macros that define operators or parameters). This
can be any name, e.g., one that is not a valid variable name. Note that unicode characters in names
are converted to LaTeX for LaTeX output (e.g., `@boson b̄` gives ``\\bar{b}``), so a LaTeX name
(which is also used for text output) is usually not needed.
"""
macro boson(exprs...) _opvariable_macro(boson, exprs) end

"""
    @fermion f
    @fermion f[:,:] g cu="c↑" ...

Define the variable `f` for a fermionic annihilation operator (see [`fermion`](@ref) and [`@boson`](@ref)).
"""
macro fermion(exprs...) _opvariable_macro(fermion, exprs) end

"""
    @tls σ
    @tls σ[:] τ ...

Define the variable `σ` for the lowering operator ``σ^-`` of a two-level system, with ``σ^+`` = `σ'`
(see [`tls`](@ref) and [`@boson`](@ref)).
"""
macro tls(exprs...) _opvariable_macro(tls, exprs) end

"""
    @tlsxyz σ
    @tlsxyz σ[:] τ ...

Define the variables `σx`, `σy` and `σz` for the Pauli matrices of a two-level system (see
[`tlsxyz`](@ref) and [`@boson`](@ref)), and similarly for each further name. With `@tlsxyz σ = ""`,
the Pauli matrices are printed as `X`, `Y` and `Z`.
"""
macro tlsxyz(exprs...) _opvariable_macro(tlsxyz, exprs; postfixes=(:x,:y,:z)) end

"""
    @anticommuting_fermions c d[:] ...
    @anticommuting_fermions cu="c↑" cd="c↓" ...

Define fermionic annihilation operators that mutually anticommute (see [`anticommuting_fermions`](@ref)
and [`@boson`](@ref)). All of them must have the same number of indices.
"""
macro anticommuting_fermions(exprs...)
    decls = parse_decl.(exprs)
    ndims = last(first(decls))
    all(d -> last(d) == ndims, decls) || throw(ArgumentError("all fermionic species in a group must have the same number of indices, got $(join(exprs,", "))."))
    :( ($(esc.(first.(decls))...),) = anticommuting_fermions($(getindex.(decls, 2)...); ndims=$ndims) )
end

"""
    QuditVariable

A qudit created by [`qudit`](@ref) or [`@qudit`](@ref): `q(i,j)` is the transition operator ``|i⟩⟨j|``
(for a qudit with indices, an object that is indexed with the indices, `q(i,j)[k]`), and `q(M)` the
operator with matrix `M`, ``\\sum_{ij} M_{ij} |i⟩⟨j|`` (or `q(M)[k]`).
"""
struct QuditVariable
    ops::QuditOps
    ndims::Int
end
"""
    SOVariable

The generators of so(n) created by [`so`](@ref) or [`@so`](@ref): `L(a,b)` is the generator ``L_{ab}``
(for generators with indices, an object that is indexed with the indices, `L(a,b)[k]`).
"""
struct SOVariable
    ops::SOOps
    ndims::Int
end
const _PairVariable = Union{QuditVariable,SOVariable}
(v::_PairVariable)(a::Integer, b::Integer) = _pairvariable_op(v, "$(v.ops.name)($a,$b)", a, b)
(q::QuditVariable)(M::AbstractMatrix) = _pairvariable_op(q, "$(q.ops.name)(M)", M)
# the operator v.ops(args...), or for variables with indices, an object that is indexed with the indices
function _pairvariable_op(v::_PairVariable, name, args...)
    iszero(v.ndims) && return v.ops(args...)
    f = (inds...) -> v.ops(args...,inds...)
    QuExprCstrctArrAlias(QuExprConstructor(name, f, name*"†", (inds...) -> f(inds...)'), v.ndims)
end
_dimsstring(v::_PairVariable) = iszero(v.ndims) ? "" : "[$(join(fill(":",v.ndims),","))]"
Base.print(io::IO, q::QuditVariable) = print(io, q.ops.name, _dimsstring(q), " (", q.ops.d, "-level qudit)")
Base.print(io::IO, L::SOVariable) = print(io, L.ops.name, _dimsstring(L), " (so(", L.ops.n, ") generators)")
Base.show(io::IO, v::_PairVariable) = print(io, v)
_checkndims(ndims) = ndims ≥ 0 ? Int(ndims) : throw(ArgumentError("the number of indices cannot be negative, got ndims = $ndims"))

"""
    qudit(name, d::Integer; ndims=0)

Return the qudit (`d`-level system) `name` (the function version of [`@qudit`](@ref)), such that
`q(i,j)` is the transition operator ``|i⟩⟨j|`` (with levels `1 ≤ i,j ≤ d`). With `ndims > 0`, the
qudit has `ndims` indices, given as `q(i,j)[k,…]`. Level 1 is the ground state, i.e., the vacuum.

Internally, all operators are expressed through the "creation" operators ``|i⟩⟨1|`` and
"annihilation" operators ``|1⟩⟨i|`` (with ``i ≥ 2``), using ``|i⟩⟨j| = |i⟩⟨1|1⟩⟨j|`` and
``|1⟩⟨1| = 1 - \\sum_{k≥2} |k⟩⟨1|1⟩⟨k|``, which gives a unique normal form. For `d=2`, `q(2,1)`
and `q(1,2)` are equivalent to ``σ^{+}`` and ``σ^{-}``. When printing, products ``|i⟩⟨1|1⟩⟨j|``
acting on the same qudit are shown as ``|i⟩⟨j|`` (written `qⁱʲ`).

With a `d×d` matrix `M`, `q(M)` (or `q(M)[k,…]`) is the operator ``\\sum_{ij} M_{ij} |i⟩⟨j|`` with this
matrix. The entries of `M` can be numbers or expressions, e.g., parameters. This gives, e.g., the
Hamiltonian of a `d`-level system, or the generators of a Lie algebra in a `d`-dimensional
representation (which then satisfy all the relations of this representation).

# Examples
```julia
q = qudit(:q, 3)
normal_form(q(1,2)*q(2,3)) # |1⟩⟨3|
q(1,1)                     # 1 - |2⟩⟨2| - |3⟩⟨3|
r = qudit(:r, 3; ndims=1)
r(2,3)[:i]                 # |2⟩⟨3| = |2⟩⟨1|1⟩⟨3| acting on qudit i
# so(3) in the vector representation, (L_k)_ab = -i ϵ_kab, with L² = l(l+1) = 2
Lx, Ly, Lz = q([0 0 0; 0 0 -im; 0 im 0]), q([0 0 im; 0 0 0; -im 0 0]), q([0 -im 0; im 0 0; 0 0 0])
normal_form(Lx^2 + Ly^2 + Lz^2)  # 2
```
"""
qudit(name::_Name, d::Integer; ndims::Integer=0) = QuditVariable(_qudit_ops(Symbol(name), d), _checkndims(ndims))

"""
    gellmann(q::QuditVariable; sqrt=Base.sqrt)

Return the `d²-1` generalized Gell-Mann matrices ``λ_1, …, λ_{d²-1}`` of the `d`-level qudit `q` (see
[`qudit`](@ref)), with the standard normalization ``\\mathrm{Tr}(λ_a λ_b) = 2δ_{ab}`` (for a qudit
with indices, as objects that are indexed with the indices, `λ[1][k]`). For `k = 2, …, d`, these
are (in this order) the pairs ``|j⟩⟨k| + |k⟩⟨j|`` and ``-i|j⟩⟨k| + i|k⟩⟨j|`` for ``j = 1, …, k-1``,
followed by the diagonal matrix ``\\sqrt{2/(l(l+1))} (\\sum_{m≤l} |m⟩⟨m| - l |l+1⟩⟨l+1|)`` with
``l = k-1``. For `d = 3`, these are the Gell-Mann matrices ``λ_1, …, λ_8`` in their usual order.
Note that since level 1 is the ground state, for `d = 2` they correspond to ``(σ^x, -σ^y, -σ^z)``.

The normalization factors ``\\sqrt{2/(l(l+1))}`` of the diagonal matrices are computed by calling
`sqrt` with the rational number `2//(l*(l+1))` (unless the result is rational, e.g., for `l = 1`).
By default, this gives `Float64` numbers. To keep them exact, pass a function that computes the
square root symbolically, e.g., `sqrt = x -> sqrt(SymPyPythonCall.Sym(x))`.

# Examples
```julia
q = qudit(:q, 3)
λ = gellmann(q)
normal_form(comm(λ[1], λ[2]))  # 2i λ₃
λ[8]                           # 0.57735 - 1.73205 q³³()
```
"""
function gellmann(q::QuditVariable; sqrt=Base.sqrt)
    cstrs = _gellmann_ops(q.ops, sqrt)
    iszero(q.ndims) ? [c() for c in cstrs] : [QuExprCstrctArrAlias(c, q.ndims) for c in cstrs]
end

"""
    spin(name, j=nothing; ndims=0)

Return the operators `(Jp, Jm, Jz)`, i.e., ``J_+``, ``J_-`` and ``J_z``, of the spin `name` (the function
version of [`@spin`](@ref)), see [`boson`](@ref) for `ndims`. They generate the Lie algebra
su(2) ≅ so(3), with ``[J_z,J_±] = ±J_±`` and ``[J_+,J_-] = 2J_z``. If `j` is not given, only these
commutation relations are used, which hold in any representation (e.g., for collective spins). If
the spin `j` is given (any number, e.g., `1//2`, `3`, or a symbolic variable), products ``J_+ J_-``
acting on the same spin are replaced by ``j(j+1) - J_z^2 + J_z``, and the vacuum is the
lowest-weight state with ``J_-|0⟩ = 0`` and ``J_z|0⟩ = -j|0⟩``. If `j` is a number (and not
symbolic), all relations of the ``(2j+1)``-dimensional representation are used, i.e., also
``J_+^a ∏_{m=-j}^{j-a} (J_z - m) = 0`` and ``J_-^b ∏_{m=b-j}^{j} (J_z - m) = 0`` (e.g., ``J_+^2 = 0`` and
``J_z^2 = 1/4`` for `j = 1//2`), so that the products on one site are ``J_+^a J_z^n`` with
``n ≤ 2j-a`` and ``J_-^b J_z^n`` with ``n ≤ 2j-b``.
"""
spin(name::_Name, j=nothing; ndims::Integer=0) = map(c -> _opvariable(c, ndims), _spin_ops(Symbol(name), j))

"""
    so(name, n::Integer; ndims=0)

Return the generators `L` of the Lie algebra so(n) of rotations in `n` dimensions (the function
version of [`@so`](@ref)), such that `L(a,b)` is the generator ``L_{ab} = -L_{ba}`` (``a,b = 1,…,n``).
With `ndims > 0`, the generators have `ndims` indices, given as `L(a,b)[k,…]`. The generators are
Hermitian, with
```math
[L_{ab}, L_{cd}] = i (δ_{ac} L_{bd} - δ_{ad} L_{bc} - δ_{bc} L_{ad} + δ_{bd} L_{ac}),
```
as, e.g., for the angular momenta ``L_{ab} = x_a p_b - x_b p_a`` of a particle in `n` dimensions (for
`n = 3`, ``(L_x, L_y, L_z) = (L_{23}, L_{31}, L_{12})``). Only these commutation relations are used,
so that results hold in any representation. Operators are stored and printed as ``L_{ab}`` with
``a < b``.

# Examples
```julia
L = so(:L, 4)
L(2,1)                             # -L¹²()
normal_form(comm(L(1,2), L(2,3)))  # -i L¹³()
```
"""
so(name::_Name, n::Integer; ndims::Integer=0) = SOVariable(_so_ops(Symbol(name), n), _checkndims(ndims))

"""
    @qudit d q
    @qudit d q[:,:]
    @qudit d q r[:] ...

Define the variable `q` for a qudit with `d` levels (see [`qudit`](@ref)), such that `q(i,j)`
is the transition operator ``|i⟩⟨j|``. With `q[:,:]` (with any number of `:`), the qudit has indices,
which are given as `q(i,j)[k,l]`. With several names, all qudits have `d` levels.

# Examples
```julia
@qudit 3 q
q(2,1) * q(1,3)        # |2⟩⟨3|
@qudit 3 m[:]
∑(:k, m(2,1)[:k] * m(1,2)[:k])
```
"""
macro qudit(d, exprs...) _opvariable_macro(qudit, exprs, d) end

"""
    @spin J
    @spin J K[:] ...
    @spin 1//2 J K[:] ...
    @spin j=jval J K[:] ...

Define the variables `Jp`, `Jm` and `Jz` for the operators ``J_+``, ``J_-`` and ``J_z`` of a spin
`J` (see [`spin`](@ref)), and similarly for each further name. With `J[:,:]` (with any number of
`:`), the spin has indices, given as, e.g., `Jp[i,j]`. A fixed spin `j` (the same for all names)
is given as the first argument if it is a number (e.g., `1//2` or `3`), and as `j=jval` otherwise
(e.g., for a variable `jval` or a symbolic value), since a name as the first argument is
interpreted as the name of a spin. (`j = "name"`, with a string, is the spin `j` with that name.)
"""
macro spin(args...)
    isempty(args) && throw(ArgumentError("usage: @spin J K[:] ..., @spin 1//2 J ..., or @spin j=jval J ..."))
    a1 = first(args)
    if a1 isa Expr && a1.head == :(=) && length(a1.args) == 2 && a1.args[1] == :j && !_isstringliteral(a1.args[2])
        _opvariable_macro(spin, args[2:end], a1.args[2]; postfixes=(:p,:m,:z))
    elseif a1 isa Symbol || (a1 isa Expr && a1.head in (:ref, :(=)))
        _opvariable_macro(spin, args; postfixes=(:p,:m,:z))
    else
        _opvariable_macro(spin, args[2:end], a1; postfixes=(:p,:m,:z))
    end
end

"""
    @so n L
    @so n L[:,:]
    @so n L M[:] ...

Define the variable `L` for the generators of so(n) (see [`so`](@ref)), such that `L(a,b)` is the
generator ``L_{ab}``. With `L[:,:]` (with any number of `:`), the generators have indices, which are
given as `L(a,b)[k,l]`. With several names, all of them are generators of so(n) for the same `n`.
"""
macro so(n, exprs...) _opvariable_macro(so, exprs, n) end

"""
    qparam(name, T=Real; ndims=0)

Return the parameter `name` (the function version of [`@qparam`](@ref)), which is real for `T = Real`
and complex for `T = Complex`. The name is a `Symbol` or a string, which can be any name, e.g.,
`"f(t)"` (which is also how the parameter appears in [`julia_expression`](@ref)). With `ndims=0`, this
is the parameter itself, and otherwise an object `g` such that `g[i,j,…]` (with `ndims` indices) is
the parameter with these indices. For a complex parameter `α`, `α'` (or `α'[i]`) is its complex
conjugate.
"""
function qparam(name::_Name, T::Type=Real; ndims::Integer=0)
    name = Symbol(name)
    T in (Real, Complex) || throw(ArgumentError("parameters are either Real or Complex, got $T"))
    c = if T == Complex
        QuExprConstructor(string(name), (inds...) -> param(name,'n',inds...), string(name,"*"), (inds...) -> param(name,'c',inds...))
    else
        QuExprConstructor(string(name), (inds...) -> param(name,'r',inds...))
    end
    _opvariable(c, ndims)
end

"""
    @qparam g
    @qparam g[:,:] α::Complex ft="f(t)" ...

Define the variable `g` for a real parameter, or several parameters at once (see [`qparam`](@ref)).
With `g[:,:]` (with any number of `:`), the parameter has indices, given as `g[i,j]`. Parameters
declared as `α::Complex` (or `α[:]::Complex`) are complex, and `α'` (or `α'[i]`) is their complex
conjugate. With `ft = "f(t)"`, the variable `ft` is the parameter with the name `f(t)` (see
[`@boson`](@ref)).
"""
macro qparam(exprs...)
    isempty(exprs) && throw(ArgumentError("usage: @qparam g or @qparam g[:,:] α::Complex ..."))
    names = Symbol[]
    calls = Expr[]
    for expr in exprs
        pname = nothing
        if expr isa Expr && expr.head == :(=) && length(expr.args) == 2
            expr, pname = expr.args
        end
        T = :Real
        if expr isa Expr && expr.head == :(::)
            length(expr.args) == 2 || throw(ArgumentError("Expected `name::Complex` or `name::Real`, got `$expr`."))
            expr, T = expr.args
            T ∈ (:Real, :Complex) || throw(ArgumentError("Expected `::Complex` or `::Real`, got `::$T`."))
        end
        name, ndims = parse_name_array(expr)
        push!(names, name)
        push!(calls, :( qparam($(pname === nothing ? QuoteNode(name) : esc(pname)), $(T == :Complex ? Complex : Real); ndims=$ndims) ))
    end
    length(names) == 1 ? :( $(esc(only(names))) = $(only(calls)) ) : :( ($(esc.(names)...),) = ($(calls...),) )
end

## functions for constructing `param`s with string macros,
# Pc"ω_i,j" = param(:ω,'n',:i,:j) (complex parameter)
# Pr"ω_i,j" = param(:ω,'r',:i,:j) (real parameter)

"""
    param(name, [state], inds...)

Return the parameter `name` with the given indices as a `QuExpr`, where `state` is `'r'` for a real
parameter, `'n'` (the default) for a complex parameter, and `'c'` for its complex conjugate. E.g.,
`param(:g,'r',:i)` is the real parameter ``g_i``.

!!! note "Legacy interface"
    This function is part of the legacy interface, which is deprecated and will be removed in
    version 2.0. Use [`@qparam`](@ref) or [`qparam`](@ref) instead, e.g., `qparam(:g; ndims=1)[:i]`
    for `param(:g,'r',:i)`, or `qparam(:α, Complex)'` for `param(:α,'c')`.
"""
param(name::Symbol,args...) = param(QuOpName(name),args...)
param(name::QuOpName,args...) = param(name,'n',args...)
function param(name::QuOpName,state::Char,inds...)
    state ∈ ('r','n','c') || throw(ArgumentError("state has to be one of n,r,c"))
    QuExpr(Param(name,state,make_indices(inds...)))
end

function parse_paramstr(s)
    s = split(s, "_"; limit=2)
    par = Symbol(s[1])
    if length(s) == 1
        inds = ()
    else
        indstr = s[2]
        # allow {} around the index expression
        if startswith(indstr,"{") && endswith(indstr,"}")
            indstr = chop(indstr,head=1,tail=1)
        end
        inds = Meta.parse.(split(indstr,","))
    end
    par, inds
end

"""
    Pc"name_inds"

Return the complex parameter `name` with the indices `inds` (separated by commas, optionally in
braces), e.g., `Pc"α"` or `Pc"g_i,j_2"` (for ``g_{i j_2}``). See also [`@Pr_str`](@ref).

!!! note "Legacy interface"
    This string macro is part of the legacy interface, which is deprecated and will be removed in
    version 2.0. Use parameters defined with [`@qparam`](@ref) instead, e.g., `g[:i,:j_2]` after
    `@qparam g[:,:]::Complex`.
"""
macro Pc_str(s)
    par, inds = parse_paramstr(s)
    param(par,'n',inds)
end
"""
    Pr"name_inds"

Return the real parameter `name` with the indices `inds` (separated by commas, optionally in braces),
e.g., `Pr"ω"` or `Pr"g_i,j_2"` (for ``g_{i j_2}``). See also [`@Pc_str`](@ref).

!!! note "Legacy interface"
    This string macro is part of the legacy interface, which is deprecated and will be removed in
    version 2.0. Use parameters defined with [`@qparam`](@ref) instead, e.g., `g[:i,:j_2]` after
    `@qparam g[:,:]`.
"""
macro Pr_str(s)
    par, inds = parse_paramstr(s)
    param(par,'r',inds)
end

"""
    expval(A)

Return the expectation value of `A`, i.e., the operators in each term are wrapped in an expectation
value (the coefficients and parameters stay outside), e.g., `expval(2g[:i]*a)` gives `2 g(i) ⟨a()⟩`.
`A` can be a `QuExpr`, a `QuTerm` or a number (which is returned as a constant expression).

See also [`corr`](@ref), [`expval_as_corrs`](@ref).
"""
function expval(A::QuTerm)
    if isempty(bares(A))
        A
    else
        QuTerm(nsuminds(A),deltas(A),params(A),[collect(expvals(A)); ExpVal(bares(A))],corrs(A),BaseOpProduct())
    end
end
"""
    corr(A)

Return the correlator (cumulant) of `A`, i.e., the operators in each term are wrapped in a
correlator ``⟨…⟩_c`` (the coefficients and parameters stay outside), e.g., `corr(a'*a)` gives
`⟨a†() a()⟩c`. `A` can be a `QuExpr`, a `QuTerm` or a number (which is returned as a constant
expression).

See also [`expval`](@ref), [`corr_as_expvals`](@ref).
"""
function corr(A::QuTerm)
    if isempty(bares(A))
        A
    else
        QuTerm(nsuminds(A),deltas(A),params(A),expvals(A),[collect(corrs(A)); Corr(bares(A))],BaseOpProduct())
    end
end

expval(A::QuExpr) = _map_quexpr_ops(expval,A)
corr(A::QuExpr) = _map_quexpr_ops(corr,A)

expval(A::Number) = QuExpr(A)
corr(A::Number) = QuExpr(A)

"""
    ∑(ind, A)
    ∑((ind1, ind2, …), A)

Return the sum of the expression `A` over the index `ind` (a symbol such as `:i`, or an index
variable, see [`@qindex`](@ref)), or over several indices. (`∑` is entered as `\\sum<tab>`.) The
sum is formal, i.e., the range of the index is not specified. The summed index is replaced by a
summation index (see [`sumindex`](@ref)), so that `∑(:i, g[:i]*b[:i])` is shown as
`∑₁ g(#₁) b(#₁)`. When normal ordering, terms with a Kronecker delta of a summation index are summed,
e.g., `normal_form(∑(:i, b[:i])*b'[:j])` gives `1 + ∑₁ b†(j) b(#₁)`.
"""
function ∑(ind::QuIndex,A::QuTerm)
    (issumindex(ind) || isintindex(ind)) && error("Index $ind to be summed over needs to be symbolic!")
    sumind = sumindex(nsuminds(A)+1)
    f = replace_inds(ind=>sumind)
    g = reorder_suminds()
    # use the form that also sets nsuminds
    g(f(A,sumind.num))
end
∑(ind::QuIndex,A::QuExpr) = _map_quexpr_ops(t->∑(ind,t),A)
∑(ind::Symbol,A::QuExpr) = ∑(QuIndex(ind),A)
# foldr calculates ∑(i1,∑(i2,∑(i3,A))), i.e., the sum over all indices in inds
∑(inds::T,A::QuExpr) where T<:NTuple{N,Union{Symbol,QuIndex}} where N = foldr(∑,inds,init=A)
∑(ind,A::Number) = ∑(ind,QuExpr(A))

const EqDict{T<:Union{ExpVal,Corr}} = OrderedDict{T,QuExpr}

"""
    QuEqSys

A system of equations of motion for expectation values (`QuEqSys{ExpVal}`) or correlators
(`QuEqSys{Corr}`), as returned by [`heisenberg_eom_system`](@ref). `EQ.eqs` is an ordered dictionary
from the expectation values (or correlators) to the right-hand sides of their equations, `EQ[k]` is
the right-hand side of the `k`-th equation, and `EQ[A]` the one for the operator `A` (e.g., `EQ[a]`
for ``⟨a⟩``) or the expectation value (or correlator) `A`. The system is printed as the list of
equations, and can be converted to LaTeX with [`latex`](@ref).
"""
struct QuEqSys{T<:Union{ExpVal,Corr}}
    eqs::EqDict{T}
end

const QuantumObject = Union{QuIndex,QuOpName,BaseOperator,Param,BaseOpProduct,OpSeq,ExpVal,Corr,QuTerm,QuExpr,QuEqSys}

@static if _DEFINE_DEFAULT_OPS
    # do not use the macros here so that we can define the constructors as const
    const a, adag = boson_ops(:a)
    const f, fdag = fermion_ops(:f)
    const σx, σy, σz = tlsxyz_ops(:σ)
    const σm, σp = tlspm_ops(:σ)

    export σx, σy, σz, σp, σm
    export a, adag, f, fdag

    @doc """
        a(inds...)
        a'(inds...)

    The default bosonic operator of the legacy interface: `a(inds...)` is the annihilation operator
    with the given indices (e.g., `a(:i)` for ``a_i`` and `a()` without indices), and `a'(inds...)` the
    creation operator. The default operators are deprecated and will be removed in version 2.0;
    define operators with the macros instead, e.g., `@boson a[:]` (see [`@boson`](@ref)).
    """ a
    @doc """
        f(inds...)
        f'(inds...)

    The default fermionic operator of the legacy interface: `f(inds...)` is the annihilation operator
    with the given indices, and `f'(inds...)` the creation operator. The default operators are
    deprecated and will be removed in version 2.0; define operators with the macros instead, e.g.,
    `@fermion f[:]` (see [`@fermion`](@ref)).
    """ f
    @doc """
        adag(inds...)
        fdag(inds...)

    Deprecated, use `a'(inds...)` and `f'(inds...)` instead (see [`a`](@ref) and [`f`](@ref)).
    """ adag
    @doc (@doc adag) fdag
    @doc """
        σx(inds...)
        σy(inds...)
        σz(inds...)

    The Pauli matrices of the default two-level system `σ` of the legacy interface, with the given
    indices. Depending on [`use_σpm`](@ref) (when they are created), they are represented as Pauli
    matrices or through ``σ^±``. The default operators are deprecated and will be removed in version
    2.0; define operators with the macros instead, e.g., `@tlsxyz σ[:]` (see [`@tlsxyz`](@ref)).
    """ σx
    @doc (@doc σx) σy
    @doc (@doc σx) σz
    @doc """
        σp(inds...)
        σm(inds...)

    The raising and lowering operators ``σ^+`` and ``σ^-`` of the default two-level system `σ` of the
    legacy interface, with the given indices. Depending on [`use_σpm`](@ref) (when they are created),
    they are represented as jump operators or through Pauli matrices. The default operators are
    deprecated and will be removed in version 2.0; define operators with the macros instead, e.g.,
    `@tls σ[:]`, with `σ'[i]` for ``σ^+_i`` and `σ[i]` for ``σ^-_i`` (see [`@tls`](@ref)).
    """ σp
    @doc (@doc σp) σm
end
