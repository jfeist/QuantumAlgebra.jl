export boson_ops, fermion_ops
export @boson_ops, @fermion_ops
export tlspm_ops, tlsxyz_ops
export @tlspm_ops, @tlsxyz_ops
export @Pr_str, @Pc_str, ∑
export param, expval, corr

export σx, σy, σz, σp, σm
export a, adag, f, fdag

const IndexInt = Int32

@concrete struct OpIndex
    # sym is the index character, apart from two special values:
    # '\0' for integer indices, which ensures they are ordered before any other indices
    # '#' for sum indices, which ensures that they are ordered before letters etc, and
    sym::Char
    # num is integer subindex, special value typemin(IndexInt) means no subindex
    num::IndexInt
    OpIndex(sym::Char,num::Integer=typemin(IndexInt)) = new(sym,num)
end
OpIndex(ii::OpIndex) = ii
OpIndex(ii::Symbol) = OpIndex(string(ii))
function OpIndex(ii::String)
    s = split(ii,"_")
    length(s) <= 2 || throw(ArgumentError("index can have at most one subindex, got $ii with subindices $(s[2:end])"))
    length(s[1]) == 1 || throw(ArgumentError("index names must be single character, got $ii."))
    sym = s[1][1]
    # we do not allow the ascii control characters as symbolic indices to ensure ordering
    sym > '#' || throw(ArgumentError("Symbolic index name cannot be character less than '#' (codepoint 35), passed Int(sym) = $(Int(sym))."))
    if length(s)==2
        OpIndex(sym,parse(IndexInt,s[2]))
    else
        OpIndex(sym)
    end
end
OpIndex(ii::Integer) = OpIndex('\0',ii)
sumindex(ii) = OpIndex('#',ii)
isintindex(ii::OpIndex) = ii.sym=='\0'
issumindex(ii::OpIndex) = ii.sym=='#'
const NoIndex = OpIndex(typemin(IndexInt))
isnoindex(ii::OpIndex) = ii == NoIndex

@inline Base.isless(i1::OpIndex,i2::OpIndex) = isless((i1.sym,i1.num),(i2.sym,i2.num))

const OpIndices = Vector{OpIndex}
assignedinds(inds::OpIndices) = inds
# const OpIndices = NTuple{5,OpIndex}
# assignedinds(inds::OpIndices) = filter(!isnoindex,inds)

make_indices(inds::OpIndices) = inds
make_indices(inds::Union{Vector,Tuple}) = make_indices(inds...)
make_indices(inds...)::OpIndices = [OpIndex(i) for i in inds]
#make_indices(i1=NoIndex,i2=NoIndex,i3=NoIndex,i4=NoIndex,i5=NoIndex)::OpIndices = (OpIndex(i1),OpIndex(i2),OpIndex(i3),OpIndex(i4),OpIndex(i5))

const _NameTable = Dict{Symbol,IndexInt}()
const _NameTableInv = Symbol[]

@concrete struct OpName
    i::IndexInt
    function OpName(name::Symbol)
        i = get!(_NameTable,name) do
            push!(_NameTableInv,name)
            length(_NameTableInv)
        end
        new(i)
    end
    OpName(name::OpName) = name
end
sym(ind::OpName) = _NameTableInv[ind.i]
Base.print(io::IO, ind::OpName) = print(io, sym(ind))
Base.isless(i1::OpName,i2::OpName) = isless(sym(i1),sym(i2))
const NoName = OpName(Symbol())

# the enum also directly defines a natural ordering,so choose this directly how we later want it
# start the counting at 1 so we can index into the tuples defined below with Int(OpType)
@enum OpType        BosonCreate_=1 FermionCreate_   TLSCreate_   TLSx_  TLSy_  TLSz_  TLSDestroy_ FermionDestroy_ BosonDestroy_
const OpType_adj = (BosonDestroy_, FermionDestroy_, TLSDestroy_, TLSx_, TLSy_, TLSz_, TLSCreate_, FermionCreate_, BosonCreate_)
const OpType_sym = ("†", "†", "⁺", "ˣ", "ʸ", "ᶻ", "⁻", "", "")
const OpType_latex = ("^\\dagger", "^\\dagger", "^+", "^x", "^y", "^z", "^-", "", "")

@concrete struct BaseOperator
    t::OpType
    name::OpName
    inds::OpIndices
    BaseOperator(t,name,inds...) = new(t,OpName(name),make_indices(inds...))
end

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

@concrete struct BaseOpProduct
    v::Vector{BaseOperator}
end
Base.isempty(A::BaseOpProduct) = isempty(A.v)
BaseOpProduct() = BaseOpProduct(BaseOperator[])

# a struct representing a delta function with unequal indices
@concrete struct δ
    iA::OpIndex
    iB::OpIndex
    δ(iA,iB) = (@assert !isnoindex(iA) && !isnoindex(iB); new(iA,iB))
end
function δ(Ainds::OpIndices,Binds::OpIndices)
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

@concrete struct ExpVal
    ops::BaseOpProduct
end
@concrete struct Corr
    ops::BaseOpProduct
end

@concrete struct Param
    name::OpName
    state::Char
    inds::OpIndices
end

@concrete struct OpTerm
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
OpTerm(p::Param)  = OpTerm(0,δ[],Param[p],ExpVal[],Corr[],BaseOpProduct())
OpTerm(E::ExpVal) = OpTerm(0,δ[],Param[],ExpVal[E],Corr[],BaseOpProduct())
OpTerm(C::Corr)   = OpTerm(0,δ[],Param[],ExpVal[],Corr[C],BaseOpProduct())

Base.isempty(A::OpTerm) = A.nsuminds == 0 && isempty(A.δs) && isempty(A.params) && isempty(A.expvals) && isempty(A.corrs) && isempty(A.bares)

const PREFAC_TYPES = Union{Int,Float64,Rational{Int},ComplexF64,Complex{Int},Complex{Rational{Int}}}
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
OpSum(A::Union{BaseOperator,Param,Corr,ExpVal}) = OpSum(OpTerm(A))
OpSum(A::OpTerm) = OpSum(((A,1),))
Base.isempty(A::OpSum) = isempty(A.terms)

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


#################################################################
## Here the "external" functions that always construct OpSum   ##
#################################################################

"`boson_ops(name::Symbol)`: return functions for creating bosonic annihilation and creation operators with name `name` (i.e., wrappers of [`BosonDestroy`](@ref) and [`BosonCreate`](@ref))"
function boson_ops(name::Symbol)
    op_name = OpName(name)
    ann(args...) = OpSum(BosonDestroy(op_name,args...))
    cre(args...) = OpSum(BosonCreate(op_name,args...))
    ann,cre
end

"`@boson_ops name`: define functions `\$name` and `\$(name)dag` for creating bosonic annihilation and creation operators with name `name` (i.e., wrappers of [`BosonDestroy`](@ref) and [`BosonCreate`](@ref))"
macro boson_ops(name)
    :( ($(esc(name)), $(esc(Symbol(name,:dag)))) = boson_ops($(Meta.quot(name))) )
end

"`fermion_ops(name::Symbol)`: return functions for creating fermionic annihilation and creation operators with name `name` (i.e., wrappers of [`FermionDestroy`](@ref) and [`FermionCreate`](@ref))"
function fermion_ops(name::Symbol)
    op_name = OpName(name)
    ann(args...) = OpSum(FermionDestroy(op_name,args...))
    cre(args...) = OpSum(FermionCreate(op_name,args...))
    ann,cre
end

"`@fermion_ops name`: define functions `\$name` and `\$(name)dag` for creating fermionic annihilation and creation operators with name `name` (i.e., wrappers of [`FermionDestroy`](@ref) and [`FermionCreate`](@ref))"
macro fermion_ops(name)
    :( ($(esc(name)), $(esc(Symbol(name,:dag)))) = fermion_ops($(Meta.quot(name))) )
end

# functions for constructing Pauli operators depending on whether we use (σ+,σ-) or (σx,σy,σz) as the "basic" operators

const _using_σpm = Ref(false)
use_σpm(t::Bool=true) = (_using_σpm[] = t; nothing)
use_σxyz() = use_σpm(false)
using_σpm() = _using_σpm[]

"`tlspm_ops(name::Symbol)`: return functions for creating jump operators for a two-level system with name `name`. The output of these functions depends on setting of use_σpm."
function tlspm_ops(name::Symbol)
    op_name = OpName(name)
    tlsm(args...) = using_σpm() ? OpSum(TLSDestroy(op_name,args...)) : OpSum((OpTerm(TLSx(op_name,args...))=>1//2, OpTerm(TLSy(op_name,args...))=>-1im//2))
    tlsp(args...) = using_σpm() ? OpSum(TLSCreate( op_name,args...)) : OpSum((OpTerm(TLSx(op_name,args...))=>1//2, OpTerm(TLSy(op_name,args...))=>1im//2))
    tlsm,tlsp
end

"`@tlspm_ops name`: define functions `\$(name)m` and `\$(name)p` creating jump operators for a two-level system with name `name`."
macro tlspm_ops(name)
    :( ($(esc(Symbol(name,:m))), $(esc(Symbol(name,:p)))) = tlspm_ops($(Meta.quot(name))) )
end

"`tlsxyz_ops(name::Symbol)`: return functions for creating Pauli operators for a two-level system with name `name`. The output of these functions depends on setting of use_σpm."
function tlsxyz_ops(name::Symbol)
    op_name = OpName(name)
    tlsx(args...) = using_σpm() ? OpSum((OpTerm(TLSDestroy(op_name,args...))=>1, OpTerm(TLSCreate(op_name,args...))=>1)) : OpSum(TLSx(op_name,args...))
    tlsy(args...) = using_σpm() ? OpSum((OpTerm(TLSDestroy(op_name,args...))=>1im, OpTerm(TLSCreate(op_name,args...))=>-1im)) : OpSum(TLSy(op_name,args...))
    tlsz(args...) = using_σpm() ? OpSum((OpTerm(BaseOpProduct([TLSCreate(op_name,args...),TLSDestroy(op_name,args...)]))=>2, OpTerm()=>-1)) : OpSum(TLSz(op_name,args...))
    tlsx,tlsy,tlsz
end

"`@tlsxyz_ops name`: define functions `\$(name)x`, `\$(name)y`, and `\$(name)z` creating Pauli operators for a two-level system with name `name`."
macro tlsxyz_ops(name)
    :( ($(esc(Symbol(name,:x))), $(esc(Symbol(name,:y))), $(esc(Symbol(name,:z)))) = tlsxyz_ops($(Meta.quot(name))) )
end

module DefaultOps
    using ..QuantumAlgebra
    export a, adag, f, fdag, σm, σp, σx, σy, σz
    @boson_ops a
    @fermion_ops f
    @tlspm_ops σ
    @tlsxyz_ops σ
end
using .DefaultOps

## functions for constructing `param`s with string macros,
# Pc"ω_i,j" = param(:ω,'n',:i,:j) (complex parameter)
# Pr"ω_i,j" = param(:ω,'r',:i,:j) (real parameter)

param(name::Symbol,args...) = param(OpName(name),args...)
function param(name::OpName,state::Char,inds...)
    state ∈ ('r','n','c') || throw(ArgumentError("state has to be one of n,r,c"))
    OpSum(Param(name,state,make_indices(inds...)))
end

function parse_paramstr(s)
    s = split(s,"_")
    par = Symbol(s[1])
    @assert length(s) <= 2
    inds = length(s)==1 ? () : Meta.parse.(split(s[2],","))
    par,inds
end

macro Pc_str(s)
    par, inds = parse_paramstr(s)
    param(par,'n',inds)
end
macro Pr_str(s)
    par, inds = parse_paramstr(s)
    param(par,'r',inds)
end

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

expval(A::OpSum) = _map_opsum_ops(expval,A)
corr(A::OpSum) = _map_opsum_ops(corr,A)

function ∑(ind::OpIndex,A::OpTerm)
    (issumindex(ind) || isintindex(ind)) && error("Index $ind to be summed over needs to be symbolic!")
    sumind = sumindex(A.nsuminds+one(A.nsuminds))
    f = replace_inds(ind=>sumind)
    g = reorder_suminds()
    # use the form that also sets nsuminds
    g(f(A,sumind.num))
end
∑(ind::OpIndex,A::OpSum) = _map_opsum_ops(t->∑(ind,t),A)
∑(ind::Symbol,A) = ∑(OpIndex(ind),A)