# export scal,param,a,adag,f,fdag,OpSumAnalytic,ExpVal,Corr
# export boson_ops, fermion_ops
# export @boson_ops,@fermion_ops
# export σx,σy,σz,σp,σm
# export @Pr_str, @Pc_str, ∑

# define for σx, σy, σz
@enum SpatialIndex x=1 y=2 z=3
SpatialIndex(a::SpatialIndex) = a

const IndexInt = Int32

@concrete struct OpIndex
    # sym is the index character, apart from two special values:
    # '#' for sum indices, which ensures that they are ordered before letters etc, and
    # '\U10ffff' (unicode noncharacter, highest allowed value) for integer indices, which ensures they are ordered after any other indices
    sym::Char
    # num is integer subindex, special value typemin(IndexInt) means no subindex
    num::IndexInt
    OpIndex(sym::Char,num::Integer=typemin(IndexInt)) = new(sym,num)
end
OpIndex(ii::OpIndex) = ii
OpIndex(ii::Symbol) = OpIndex(string(ii))
function OpIndex(ii::String)
    s = split(ii,"_")
    length(s) <= 2 || error("index can have at most one subindex")
    length(s[1]) == 1 || error("index names must be single character.")
    sym = s[1][1]
    # we do not allow the ascii control characters as symbolic indices to ensure ordering
    sym > '#' || error("Symbolic index name cannot be character below '#' (codepoint 35), passed Int(sym) = $(Int(sym)).")
    if length(s)==2
        OpIndex(sym,parse(IndexInt,s[2]))
    else
        OpIndex(sym)
    end
end
OpIndex(ii::Integer) = OpIndex('\U10ffff',ii)
sumindex(ii) = OpIndex('#',ii)
isintindex(ii::OpIndex) = ii.sym=='\U10ffff'
issumindex(ii::OpIndex) = ii.sym=='#'

@inline Base.isless(i1::OpIndex,i2::OpIndex) = isless((i1.sym,i1.num),(i2.sym,i2.num))

const OpIndices = Vector{OpIndex}

make_indices(inds::OpIndices) = inds
make_indices(inds::Vector)::OpIndices = OpIndex.(inds)
make_indices(inds::Tuple)::OpIndices = make_indices(inds...)
make_indices(inds...)::OpIndices = OpIndex[OpIndex(i) for i in inds]

# the enum also directly defines a natural ordering,so choose this directly how we later want it
@enum OpType BosonCreate_ BosonDestroy_ FermionCreate_ FermionDestroy_ σplus_ σminus_ σ_

const _NameTable = Dict{Symbol,IndexInt}()
const _NameTableInv = Symbol[]

@concrete struct NameIndex
    i::IndexInt
    function NameIndex(name::Symbol)
        i = get!(_NameTable,name) do
            push!(_NameTableInv,name)
            length(_NameTableInv)
        end
        new(i)
    end
end
sym(ind::NameIndex) = _NameTableInv[ind.i]
Base.print(io::IO, ind::NameIndex) = print(io, sym(ind))
Base.isless(i1::NameIndex,i2::NameIndex) = isless(sym(i1),sym(i2))
const NoName = NameIndex(Symbol())

@concrete struct BaseOperator
    t::OpType
    a::SpatialIndex
    name::NameIndex
    inds::OpIndices    
end
# provide outer constructors useful for the specific operators
BaseOperator(t,inds::OpIndices) = BaseOperator(t,x,NoName,inds)
BaseOperator(t,a::SpatialIndex,inds::OpIndices) = BaseOperator(t,a,NoName,inds)
BaseOperator(t,name::NameIndex,inds::OpIndices) = BaseOperator(t,x,name,inds)
BaseOperator(t,name::Symbol,inds::OpIndices) = BaseOperator(t,x,NameIndex(name),inds)

for (op,desc) in (
    (:BosonDestroy,"bosonic annihilation"),
    (:BosonCreate,"bosonic creation"),
    (:FermionDestroy,"fermionic annihilation"),
    (:FermionCreate,"fermionic creation"))
    @eval begin
        "`$($op)(name,inds)`: represent $($desc) operator ``name_{inds}``"
        $op(name::Union{Symbol,NameIndex},inds...) = BaseOperator($(Symbol(op,:_)),name,make_indices(inds...))
    end
end

"`boson_ops(name::Symbol)`: return functions for creating bosonic annihilation and creation operators with name `name` (i.e., wrappers of [`BosonDestroy`](@ref) and [`BosonCreate`](@ref))"
function boson_ops(name::Symbol)
    op_ind = NameIndex(name)
    ann(args...) = BosonDestroy(op_ind,args...)
    cre(args...) = BosonCreate(op_ind,args...)
    ann,cre
end

"`@boson_ops name`: define functions `\$name` and `\$(name)dag` for creating bosonic annihilation and creation operators with name `name` (i.e., wrappers of [`BosonDestroy`](@ref) and [`BosonCreate`](@ref))"
macro boson_ops(name)
    :( ($(esc(name)), $(esc(Symbol(name,:dag)))) = boson_ops($(Meta.quot(name))) )
end

"`fermion_ops(name::Symbol)`: return functions for creating fermionic annihilation and creation operators with name `name` (i.e., wrappers of [`FermionDestroy`](@ref) and [`FermionCreate`](@ref))"
function fermion_ops(name::Symbol)
    op_ind = NameIndex(name)
    ann(args...) = FermionDestroy(op_ind,args...)
    cre(args...) = FermionCreate(op_ind,args...)
    ann,cre
end

"`@fermion_ops name`: define functions `\$name` and `\$(name)dag` for creating fermionic annihilation and creation operators with name `name` (i.e., wrappers of [`FermionDestroy`](@ref) and [`FermionCreate`](@ref))"
macro fermion_ops(name)
    :( ($(esc(name)), $(esc(Symbol(name,:dag)))) = fermion_ops($(Meta.quot(name))) )
end

@boson_ops a
@fermion_ops f

for (op,desc,sym) in (
    (:σminus,"TLS annihilation","σ^-"),
    (:σplus,"TLS creation","σ^+"))
    @eval begin
        "`$($op)(inds)`: represent $($desc) operator ``$($sym)_{inds}``"
        $op(inds...) = BaseOperator($(Symbol(op,:_)),make_indices(inds...))
    end
end

"`σ(a,inds)`: represent Pauli matrix ``σ_{a,inds}`` for two-level system (TLS), where ``a ∈ \\{x,y,z\\}``."
σ(a,inds...) = BaseOperator(σ_,SpatialIndex(a),make_indices(inds...))

@concrete struct BaseOpProduct
    v::Vector{BaseOperator}
end
Base.isempty(A::BaseOpProduct) = isempty(A.v)
BaseOpProduct() = BaseOpProduct(BaseOperator[])

## functions for constructing Pauli operators depending on whether we use (σ+,σ-) or (σx,σy,σz) as the "basic" operators

const _using_σpm = Ref(false)
use_σpm(t::Bool=true) = (_using_σpm[] = t; nothing)
use_σxyz() = use_σpm(false)
using_σpm() = _using_σpm[]

"`σp(n)`: construct ``σ^+_n = \\frac12 σ_{x,n} + \\frac{i}{2} σ_{y,n}``"
σp(n...) = using_σpm() ? σplus(n...) : OpSum((OpTerm(σ(x,n...))=>1//2, OpTerm(σ(y,n...))=>1im//2))
"`σm(n)`: construct ``σ^-_n = \\frac12 σ_{x,n} - \\frac{i}{2} σ_{y,n}``"
σm(n...) = using_σpm() ? σminus(n...) : OpSum((OpTerm(σ(x,n...))=>1//2, OpTerm(σ(y,n...))=>-1im//2))
"`σx(n)`: construct ``σ_{x,n}``"
σx(n...) = using_σpm() ? OpSum((OpTerm(σminus(n...))=>1, OpTerm(σplus(n...))=>1)) : σ(x,n...)
"`σy(n)`: construct ``σ_{y,n}``"
σy(n...) = using_σpm() ? OpSum((OpTerm(σminus(n...))=>1im, OpTerm(σplus(n...))=>-1im)) : σ(y,n...)
"`σz(n)`: construct ``σ_{z,n}``"
σz(n...) = using_σpm() ? OpSum((OpTerm(BaseOpProduct([σplus(n...),σminus(n...)]))=>2, OpTerm()=>-1)) : σ(z,n...)

# a struct representing a delta function with unequal indices
@concrete struct δ
    iA::OpIndex
    iB::OpIndex
end
function δ(Ainds::OpIndices,Binds::OpIndices)
    length(Ainds) == length(Binds) || return nothing

    res = Vector{δ}(undef,length(Ainds))
    jj = 0
    for (iA,iB) in zip(Ainds,Binds)
        if iA != iB
            # if they are integer indices and different, the result is zero
            isintindex(iA) && isintindex(iB) && return nothing
            res[jj+=1] = iB<iA ? δ(iB,iA) : δ(iA,iB)
        end
    end
    resize!(res,jj)
    res
end

@concrete struct ExpVal
    ops::BaseOpProduct
end
@concrete struct Corr
    ops::BaseOpProduct
end

@concrete struct Param
    name::NameIndex
    state::Char
    inds::OpIndices
end
Param(name::Symbol,args...) = Param(NameIndex(name),args...)

## functions for constructing `param`s with string macros,
# Pc"ω_i,j" = param(:ω,'n',:i,:j) (complex parameter)
# Pr"ω_i,j" = param(:ω,'r',:i,:j) (real parameter)

function parse_paramstr(s)
    s = split(s,"_")
    par = Symbol(s[1])
    @assert length(s) <= 2
    inds = length(s)==1 ? () : Meta.parse.(split(s[2],","))
    par,inds
end

macro Pc_str(s)
    par, inds = parse_paramstr(s)
    Param(par,'n',make_indices(inds))
end
macro Pr_str(s)
    par, inds = parse_paramstr(s)
    Param(par,'r',make_indices(inds))
end

struct OpTerm
    nsuminds::IndexInt # have a sum over n indices, represented by OpIndex with issumindex(ind)==true
    δs::Vector{δ}
    params::Vector{Param}
    expvals::Vector{ExpVal}
    corrs::Vector{Corr}
    bares::BaseOpProduct
end
OpTerm() = OpTerm(BaseOpProduct())
OpTerm(op::BaseOperator) = OpTerm(BaseOpProduct([op]))
OpTerm(ops::BaseOpProduct) = OpTerm(0,δ[],Param[],ExpVal[],Corr[],ops)
OpTerm(δs::Vector{δ},ops::BaseOpProduct) = OpTerm(0,δs,Param[],ExpVal[],Corr[],ops)

_switch_bares(A::OpTerm,newbares::BaseOpProduct) = OpTerm(A.nsuminds,A.δs,A.params,A.expvals,A.corrs,newbares)

Base.isempty(A::OpTerm) = A.nsuminds == 0 && isempty(A.δs) && isempty(A.params) && isempty(A.expvals) && isempty(A.corrs) && isempty(A.bares)

function expval(A::OpTerm)
    if isempty(A.bares)
        A
    else
        OpTerm(A.nsuminds,A.δs,A.params,[A.expvals; ExpVal(A.bares)],A.corrs,BaseOpProduct())
    end
end
function corr(A::OpTerm)
    if isempty(A.bares)
        A
    else
        OpTerm(A.nsuminds,A.δs,A.params,A.expvals,[A.corrs; Corr(A.bares)],BaseOpProduct())
    end
end

PREFAC_TYPES = Union{Int,Float64,Rational{Int},ComplexF64,Complex{Int},Complex{Rational{Int}}}
struct OpSum
    # sum of Operators is saved as Dictionary of operators with scalar prefactors
    terms::Dict{OpTerm,PREFAC_TYPES}
    OpSum() = new(Dict{OpTerm,PREFAC_TYPES}())
end
function OpSum(itr)
    A = OpSum()
    for (t,s) in itr
        _add_sum_term!(A,t,simplify_number(s))
    end
    A
end

function _add_sum_term!(A::OpSum,oB::OpTerm,sB)
    iszero(sB) && return A
    sold = get(A.terms,oB,zero(sB))
    # function barrier to have concrete types
    _add_sum_term!(A,oB,sB,sold)
end
function _add_sum_term!(A::OpSum,oB::OpTerm,sB,sold)
    snew = sB + sold
    if iszero(snew)
        delete!(A.terms,oB)
    else
        A.terms[oB] = simplify_number(snew)
    end
    A
end
_map_opsum_ops(f,A::OpSum) = OpSum((f(t),s) for (t,s) in A.terms)
expval(A::OpSum) = _map_opsum_ops(expval,A)
corr(A::OpSum) = _map_opsum_ops(corr,A)
