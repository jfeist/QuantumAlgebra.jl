using Preferences
using OrderedCollections

export QuExpr
export boson_ops, fermion_ops
export @boson_ops, @fermion_ops, @anticommuting_fermion_group
export tlspm_ops, tlsxyz_ops, qudit_ops, gellmann_ops, spin_ops, so_ops
export @tlspm_ops, @tlsxyz_ops
export @Pr_str, @Pc_str, ∑
export param, expval, corr
export @boson, @fermion, @anticommuting_fermions, @tls, @tlsxyz, @param
export map_scalar_function

# compile-time options
const _DEFINE_DEFAULT_OPS = @load_preference("define_default_ops", false)

"""
    set_define_default_ops(t::Bool)

Set the preference for whether to define default operators (`a`, `adag`, `f`, `fdag`, `σx`, `σy`, `σz`, `σp`, `σm`) upon import.
Default is `true`.

Note that changing this preference requires restarting the Julia session to take effect.
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

Select whether TLS operators are expressed in the `(σ₊, σ₋)` basis (`true`) or
in the `(σx, σy, σz)` basis (`false`). If `set_preference` is `true`, the choice
is persisted via `Preferences.jl` for future sessions.
"""
function use_σpm(t::Bool=true; set_preference=false)
    _using_σpm[] = t
    set_preference && @set_preferences!("use_σpm" => t)
    nothing
end
"""
    use_σxyz(; set_preference=false)

Convenience helper to switch to the `(σx, σy, σz)` basis. Equivalent to
`use_σpm(false; set_preference)`.
"""
use_σxyz(;set_preference=false) = use_σpm(false;set_preference)
"""
    using_σpm()

Return whether TLS operators are currently configured to use the `(σ₊, σ₋)`
basis.
"""
using_σpm() = _using_σpm[]

const _auto_normal_form = Ref(false)
"""
    auto_normal_form(t::Bool=true; set_preference=false)

Enable or disable automatic normal-ordering of expressions on construction.
When `set_preference` is `true`, the setting is stored with `Preferences.jl`.
"""
function auto_normal_form(t::Bool=true; set_preference=false)
    _auto_normal_form[] = t
    set_preference && @set_preferences!("auto_normal_form" => t)
    nothing
end
"""
    using_auto_normal_form()

Return whether automatic normal-ordering is currently enabled.
"""
using_auto_normal_form() = _auto_normal_form[]

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
    QuIndex

Represents an index in a quantum operator expression.

A `QuIndex` can be:
- An integer index (created with `QuIndex(i)` for integer `i`)
- A symbolic index like `QuIndex(:i)` or `QuIndex("i_2")`
- A sum index created with `sumindex(n)`

# Examples
```julia
i = QuIndex(1)              # Integer index 1
i = QuIndex(:i)             # Symbolic index 'i'
i = QuIndex("i_2")          # Symbolic index 'i' with subscript 2
si = sumindex(1)            # Sum index #₁
```

# See also
[`sumindex`](@ref), [`issumindex`](@ref), [`isintindex`](@ref)
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

Create a summation index with index `i`. These indices are used internally for
implicit summation (Einstein summation convention).
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

Represents a Kronecker delta function δ_{iA,iB} in a quantum expression.

The delta function is zero if indices are different, and one if they are the same.

# Examples
```julia
δ(QuIndex(1), QuIndex(2))    # δ₁₂
δ(QuIndex(:i), QuIndex(:j))  # δᵢⱼ
```
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

Represents a single quantum operator in an expression.

A `BaseOperator` is one of:
- Bosonic operators: `BosonCreate`, `BosonDestroy`
- Fermionic operators: `FermionCreate`, `FermionDestroy`
- Two-level system (TLS) operators: `TLSCreate`, `TLSDestroy`, `TLSx`, `TLSy`, `TLSz`

Each operator is associated with a name and indices.

# Examples
```julia
a = BosonCreate(:a)        # Bosonic creation operator a†
f = FermionDestroy(:f, 1)  # Fermionic annihilation operator f₁
σx = TLSx(:σ, :i)          # TLS x operator σ_i
```
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

Represents an expectation value ⟨ops⟩ in a quantum expression.

Expectation values contain products of quantum operators that should be evaluated as a single unit.

# Examples
```julia
ev = ExpVal(BaseOpProduct([a, a_dag]))  # ⟨a† a⟩
```

# See also
[`Corr`](@ref), [`expval`](@ref)
"""
ExpVal

"""
    Corr(ops::BaseOpProduct)

Represents a correlator ⟨ops⟩c in a quantum expression.

Correlators (or cumulants) represent connected correlations in the expectation value.

# Examples
```julia
corr = Corr(BaseOpProduct([a, a_dag]))  # ⟨a† a⟩c
```

# See also
[`ExpVal`](@ref), [`corr`](@ref)
"""
Corr

"""
    Param(name, state, inds)

Represents a parameter in a quantum expression.

Parameters are scalar values that can appear in quantum expressions, with optional indices
and conjugation state.

# Arguments
- `name`: Symbol identifying the parameter
- `state`: 'r' for real, 'n' for normal, 'c' for conjugate
- `inds`: Quantum indices (optional)

# Examples
```julia
p = Param(:λ, 'r')        # Real parameter λ
p = Param(:ω, 'r', 1)     # Real parameter ω₁
```
"""
Param


"""
    QuExpr

A quantum algebra expression represented as a sum of terms.

A `QuExpr` is a dictionary mapping `QuTerm`s to numeric coefficients stored as `Number`.
Using `Number` allows heterogeneous numeric types (e.g., `Int`, `Float64`, `Rational`,
`Complex`) to coexist without forced type conversions.

# Constructors

- The quantum operator constructors generated with `@boson_ops`, `@fermion_ops`, etc., return `QuExpr` objects.
- `QuExpr()`: Create an empty expression (zero) - `zero(QuExpr)` works as well.
- `QuExpr(x)`: Create a constant expression from a number `x`. `QuExpr(1)` is equivalent to `one(QuExpr)`.

# See also
[`QuTerm`](@ref), [`normal_form`](@ref)
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
        _add_with_auto_order!(A,t,simplify_number(s))
    end
    A
end
QuExpr(A::Union{BaseOperator,Param,Corr,ExpVal}) = QuExpr(QuTerm(A))
QuExpr(A::QuTerm) = QuExpr(((A,1),))
QuExpr(s::Number) = QuExpr(((QuTerm(),s),))
Base.isempty(A::QuExpr) = isempty(A.terms)
Base.copy(A::QuExpr) = QuExpr(copy(A.terms))

_add_with_auto_order!(A::QuExpr,B::QuTerm,sB) = using_auto_normal_form() ? _add_with_normal_order!(A,B,sB) : _add_sum_term!(A,B,sB)

function _add_sum_term!(A::QuExpr,oB::QuTerm,sB)
    iszero(sB) && return A
    sold = get(A.terms,oB,zero(sB))
    # function barrier to have concrete types
    _add_sum_term!(A,oB,sB,sold)
end
function _add_sum_term!(A::QuExpr,oB::QuTerm,sB,sold)
    snew = sB + sold
    if iszero(snew)
        delete!(A.terms,oB)
    else
        A.terms[oB] = simplify_number(snew)
    end
    A
end
_map_quexpr_ops(f,A::QuExpr) = QuExpr((f(t),s) for (t,s) in A.terms)

"""
    map_scalar_function(f, A::QuExpr)

Apply a function to all scalar coefficients in a quantum expression.

# Arguments
- `f`: Function to apply to each coefficient
- `A`: Quantum expression

# Returns
A new `QuExpr` with function applied to all coefficients

# Examples
```julia
A = QuExpr(a) + 2 * QuExpr(a_dag)
B = map_scalar_function(x -> 2x, A)  # Double all coefficients
```
"""
map_scalar_function(f,A::QuExpr) = QuExpr((t,f(s)) for (t,s) in A.terms)

"""
    QuExprConstructor

Helper struct for creating operator expressions with automatic conjugate support.

Used internally by macros like `@boson_ops` and `@fermion_ops` to provide
both normal and daggered operator constructors.
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

Create constructors for bosonic operators with the given name.

Returns a tuple of `(annihilation, creation)` operators that can be called
to create quantum expressions containing these operators.

# Arguments
- `name`: Symbol for the operator name (e.g., `:a`, `:b`)

# Returns
- Tuple `(a, a†)` of `QuExprConstructor` objects

# Examples
```julia
a, a_dag = boson_ops(:a)
expr = a(1) + a_dag(2)  # a₁ + a₂†
```

# See also
[`fermion_ops`](@ref), [`@boson_ops`](@ref)
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

Create constructors for fermionic operators with the given name.

Returns a tuple of `(annihilation, creation)` operators that can be called
to create quantum expressions containing these operators.

# Arguments
- `name`: Symbol for the operator name (e.g., `:f`, `:c`)

# Returns
- Tuple `(f, f†)` of `QuExprConstructor` objects

# Examples
```julia
f, f_dag = fermion_ops(:f)
expr = f(1) + f_dag(2)  # f₁ + f₂†
```

# See also
[`boson_ops`](@ref), [`@fermion_ops`](@ref)
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

Return constructors for TLS ladder operators `(σ₋, σ₊)` for the given `name`.

The concrete output depends on `use_σpm`: when true, returns annihilation/creation
expressions; when false, returns linear combinations of `σx` and `σy`.
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

Return constructors for Pauli operators `(σₓ, σ_y, σ_z)` for the given `name`.

The concrete output depends on `use_σpm`: when true, returns combinations of ladder
operators; when false, returns the direct `σx/σy/σz` operators.
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

"""
    QuditOps

Constructor for the transition operators of a qudit (see [`qudit_ops`](@ref)).
"""
struct QuditOps
    name::QuOpName
    d::Int
end

"""
    qudit_ops(name::Symbol, d::Integer)

Return a constructor `q` for the operators of a `d`-level system (qudit) with the given `name`.
`q(i, j, inds...)` gives the transition operator ``|i⟩⟨j|`` (with levels `1 ≤ i,j ≤ d`)
acting on the qudit with indices `inds`. Level 1 is the ground state, i.e., the vacuum.

Internally, all operators are expressed through the "creation" operators ``|i⟩⟨1|`` and
"annihilation" operators ``|1⟩⟨i|`` (with ``i ≥ 2``), using ``|i⟩⟨j| = |i⟩⟨1|1⟩⟨j|`` and
``|1⟩⟨1| = 1 - \\sum_{k≥2} |k⟩⟨1|1⟩⟨k|``, which gives a unique normal form. For `d=2`, `q(2,1)`
and `q(1,2)` are equivalent to ``σ^{+}`` and ``σ^{-}``. When printing, products ``|i⟩⟨1|1⟩⟨j|``
acting on the same qudit are shown as ``|i⟩⟨j|`` (written `qⁱʲ`).

# Examples
```julia
q = qudit_ops(:q, 3)
q(2,1,:i)                  # |2⟩⟨1| acting on qudit i
q(2,3,:i)                  # |2⟩⟨3| = |2⟩⟨1|1⟩⟨3|
normal_form(q(1,2)*q(2,3)) # |1⟩⟨3|
q(1,1)                     # 1 - |2⟩⟨2| - |3⟩⟨3|
```
"""
function qudit_ops(name::Symbol, d::Integer)
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
Base.print(io::IO, q::QuditOps) = print(io, q.name, " (constructor for ", q.d, "-level qudit)")
Base.show(io::IO, q::QuditOps) = print(io, q)

"""
    gellmann_ops(q::QuditOps; sqrt=Base.sqrt)

Return constructors for the `d²-1` generalized Gell-Mann matrices ``λ_1, …, λ_{d²-1}`` of the
`d`-level qudit `q` (see [`qudit_ops`](@ref)), with the standard normalization
``\\mathrm{Tr}(λ_a λ_b) = 2δ_{ab}``. For `k = 2, …, d`, these are (in this order) the pairs
``|j⟩⟨k| + |k⟩⟨j|`` and ``-i|j⟩⟨k| + i|k⟩⟨j|`` for ``j = 1, …, k-1``, followed by the diagonal matrix
``\\sqrt{2/(l(l+1))} (\\sum_{m≤l} |m⟩⟨m| - l |l+1⟩⟨l+1|)`` with ``l = k-1``. For `d = 3`, these are
the Gell-Mann matrices ``λ_1, …, λ_8`` in their usual order. Note that since level 1 is the ground
state, for `d = 2` they correspond to ``(σ^x, -σ^y, -σ^z)``.

The normalization factors ``\\sqrt{2/(l(l+1))}`` of the diagonal matrices are computed by calling
`sqrt` with the rational number `2//(l*(l+1))` (unless the result is rational, e.g., for `l = 1`).
By default, this gives `Float64` numbers. To keep them exact, pass a function that computes the
square root symbolically, e.g., `sqrt = x -> sqrt(SymPyPythonCall.Sym(x))`.

# Examples
```julia
q = qudit_ops(:q, 3)
λ = gellmann_ops(q)
normal_form(comm(λ[1](:i), λ[2](:i)))  # 2i λ₃(i)
λ[8](:i)                               # 0.57735 - 1.73205 q³³(i)
```
"""
function gellmann_ops(q::QuditOps; sqrt=Base.sqrt)
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

"`@boson_ops name`: define function `\$name` for creating bosonic annihilation operators with name `name` (also defines deprecated `\$(name)dag`, use `\$(name)'` instead) "
macro boson_ops(name)
    :( ($(esc(name)), $(esc(Symbol(name,:dag)))) = boson_ops($(Meta.quot(name))) )
end

"`@fermion_ops name`: define function `\$name` for creating fermionic annihilation operators with name `name` (also defines deprecated `\$(name)dag`, use `\$(name)'` instead)"
macro fermion_ops(name)
    :( ($(esc(name)), $(esc(Symbol(name,:dag)))) = fermion_ops($(Meta.quot(name))) )
end

"`@anticommuting_fermion_group name1 name2 ...`: define a group of mutually anticommuting fermionic operators"
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

"`@tlspm_ops name`: define functions `\$(name)m` and `\$(name)p` creating jump operators for a two-level system with name `name`."
macro tlspm_ops(name)
    :( ($(esc(Symbol(name,:m))), $(esc(Symbol(name,:p)))) = tlspm_ops($(Meta.quot(name))) )
end

"`@tlsxyz_ops name`: define functions `\$(name)x`, `\$(name)y`, and `\$(name)z` creating Pauli operators for a two-level system with name `name`."
macro tlsxyz_ops(name)
    :( ($(esc(Symbol(name,:x))), $(esc(Symbol(name,:y))), $(esc(Symbol(name,:z)))) = tlsxyz_ops($(Meta.quot(name))) )
end

"""
    spin_ops(name::Symbol, j=nothing)

Return constructors `(Jp, Jm, Jz)` for the operators ``J_+``, ``J_-`` and ``J_z`` of a spin, i.e.,
of the Lie algebra su(2) ≅ so(3), with ``[J_z,J_±] = ±J_±`` and ``[J_+,J_-] = 2J_z``. If `j` is not
given, only these commutation relations are used, which hold in any representation (e.g., for
collective spins). If the spin `j` is given (any number, e.g., `1//2`, `3`, or a symbolic variable),
products ``J_+ J_-`` acting on the same spin are replaced by ``j(j+1) - J_z^2 + J_z``, and the vacuum
is the lowest-weight state with ``J_-|0⟩ = 0`` and ``J_z|0⟩ = -j|0⟩``. Note that other relations that
hold for a specific `j` (e.g., ``J_+^2 = 0`` for `j = 1//2`) are not used.
"""
function spin_ops(name::Symbol, j=nothing)
    op_name = QuOpName(name)
    extra = _spin_extra(j)
    cstr(t) = (args...) -> QuExpr(BaseOperator(t, op_name, make_indices(args...), extra))
    Jp = QuExprConstructor(string(name,"⁺"), cstr(SpinPlus_), string(name,"⁻"), cstr(SpinMinus_))
    (Jp, Jp', QuExprConstructor(string(name,"ᶻ"), cstr(SpinZ_)))
end

"""
    SOOps

Constructor for the generators of so(n) (see [`so_ops`](@ref)).
"""
struct SOOps
    name::QuOpName
    n::Int
end

"""
    so_ops(name::Symbol, n::Integer)

Return a constructor `L` for the generators ``L_{ab} = -L_{ba}`` (``a,b = 1,…,n``) of the Lie algebra
so(n) of rotations in `n` dimensions, called as `L(a, b, inds...)`. The generators are Hermitian, with
```math
[L_{ab}, L_{cd}] = i (δ_{ac} L_{bd} - δ_{ad} L_{bc} - δ_{bc} L_{ad} + δ_{bd} L_{ac}),
```
as, e.g., for the angular momenta ``L_{ab} = x_a p_b - x_b p_a`` of a particle in `n` dimensions (for
`n = 3`, ``(L_x, L_y, L_z) = (L_{23}, L_{31}, L_{12})``). Only these commutation relations are used,
so that results hold in any representation. Operators are stored and printed as ``L_{ab}`` with
``a < b``.

# Examples
```julia
L = so_ops(:L, 4)
L(2,1,:i)                          # -L¹²(i)
normal_form(comm(L(1,2), L(2,3)))  # -i L¹³()
```
"""
function so_ops(name::Symbol, n::Integer)
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
Base.print(io::IO, L::SOOps) = print(io, L.name, " (constructor for so(", L.n, ") generators)")
Base.show(io::IO, L::SOOps) = print(io, L)

struct QuExprCstrctArrAlias{T}
    constructor::T
    numdims::Int
end
function Base.getindex(q::QuExprCstrctArrAlias,is...)
    length(is) == q.numdims || throw(ArgumentError("Expected $(q.numdims) indices, got $(length(is))."))
    q.constructor(is...)
end
Base.adjoint(q::QuExprCstrctArrAlias) = QuExprCstrctArrAlias(q.constructor',q.numdims)

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

function quexpr_variable_code(expr, opfunc, postfixes=(Symbol(),))
    name, ndims = parse_name_array(expr)
    funs = ntuple(i->gensym(), length(postfixes))
    vars = esc.(Symbol.(name,postfixes))
    code = Expr(:block, :( ($(funs...),) = $(opfunc)($(Meta.quot(name))) ))
    for (var,fun) in zip(vars,funs)
        rhs = iszero(ndims) ? :( $fun() ) : :( QuExprCstrctArrAlias($fun,$ndims) )
        push!(code.args, :( $var = $rhs ))
    end
    push!(code.args, :( ($(vars...),) ))
    code
end

macro boson(name) quexpr_variable_code(name, boson_ops) end
macro fermion(name) quexpr_variable_code(name, fermion_ops) end
macro tls(name) quexpr_variable_code(name, tlspm_ops) end
macro tlsxyz(name) quexpr_variable_code(name, tlsxyz_ops, (:x,:y,:z)) end

macro anticommuting_fermions(exprs...)
    namesndims = parse_name_array.(exprs)
    names = first.(namesndims)
    # start the groupname (which is the "internal" species name) by concatenating all names,
    # to ensure reasonable sorting relative to other fermionic species
    code = quote
        groupname = Symbol($names...,gensym())
        add_groupaliases(groupname,$names)
        ann, cre = fermion_ops(groupname)
    end
    ndims = -1
    for (ii,(name,ndims_i)) in enumerate(namesndims)
        ann_i = Symbol(:ann,ii)
        ndims == -1 && (ndims = ndims_i)
        ndims_i == ndims || throw(ArgumentError("All fermionic species in a group must have the same number of indices, got $(ndims_i) for $(name), expected $(ndims)."))
        push!(code.args,:( $ann_i = QuExprConstructor($(string(name)), (args...) -> ann($ii,args...), $(string(name,"†")), (args...) -> ann'($ii,args... )) ))
        if ndims==0
            push!(code.args, :( $(esc(Symbol(name))) = $(ann_i)() ))
        else
            push!(code.args, :( $(esc(Symbol(name))) = QuExprCstrctArrAlias($(ann_i),$ndims) ))
        end
    end
    push!(code.args, :( nothing ))
    code
end

macro param(name)
    partype = :Real
    if name isa Expr && name.head == :(::)
        @assert length(name.args) == 2
        name, partype = name.args
        partype ∈ (:Real, :Complex) || throw(ArgumentError("Expected `::Complex` or `::Real`, got `::$(partype)`."))
    end
    name, ndims = parse_name_array(name)
    state = partype == :Complex ? 'n' : 'r'
    if ndims == 0
        :( $(esc(Symbol(name))) = param($(Meta.quot(name)),$state) )
    else
        cstrctfun(state) = :( (args...) -> param($(Meta.quot(name)), $state,  args...) )
        cstrctargs = (string(name), cstrctfun(state))
        if state == 'n'
            cstrctargs = (cstrctargs..., string(name,"*"), cstrctfun('c'))
        end
        :( $(esc(Symbol(name))) = QuExprCstrctArrAlias(QuExprConstructor($(cstrctargs...)),$ndims) )
    end
end

## functions for constructing `param`s with string macros,
# Pc"ω_i,j" = param(:ω,'n',:i,:j) (complex parameter)
# Pr"ω_i,j" = param(:ω,'r',:i,:j) (real parameter)

"""
    param(name, [state], inds...)

Construct a parameter as a `QuExpr`. `state` is one of `n` (normal/complex),
`r` (real), or `c` (conjugate). Indices are optional.
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

macro Pc_str(s)
    par, inds = parse_paramstr(s)
    param(par,'n',inds)
end
macro Pr_str(s)
    par, inds = parse_paramstr(s)
    param(par,'r',inds)
end

"""
    expval(A)

Wrap bare operators inside `A` into a formal expectation value ⟨A⟩.
Accepts `QuTerm` or `QuExpr`; numeric inputs are returned as constant expressions.
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

Wrap bare operators inside `A` into a formal connected correlator ⟨A⟩₍c₎.
Accepts `QuTerm` or `QuExpr`; numeric inputs are returned as constant expressions.
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

"`∑(ind,A::QuExpr)`: return (formal) sum of expression A over index ind."
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
end
