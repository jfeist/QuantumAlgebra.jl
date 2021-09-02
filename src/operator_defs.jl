export QuExpr
export boson_ops, fermion_ops
export @boson_ops, @fermion_ops, @anticommuting_fermion_group
export tlspm_ops, tlsxyz_ops
export @tlspm_ops, @tlsxyz_ops
export @Pr_str, @Pc_str, ∑
export param, expval, corr
export σx, σy, σz, σp, σm
export a, adag, f, fdag

# externally changeable options
const _using_σpm = Ref(false)
use_σpm(t::Bool=true) = (_using_σpm[] = t; nothing)
use_σxyz() = use_σpm(false)
using_σpm() = _using_σpm[]

const _auto_normal_form = Ref(false)
auto_normal_form(t::Bool=true) = (_auto_normal_form[] = t; nothing)
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
sumindex(ii) = QuIndex('#',ii)
isintindex(ii::QuIndex) = ii.sym=='\0'
issumindex(ii::QuIndex) = ii.sym=='#'
const NoIndex = QuIndex(typemin(IndexInt))
isnoindex(ii::QuIndex) = ii == NoIndex

@inline Base.isless(i1::QuIndex,i2::QuIndex) = isless((i1.sym,i1.num),(i2.sym,i2.num))

const QuIndices = Vector{QuIndex}
Base.tail(inds::QuIndices) = inds[2:end]
assignedinds(inds::QuIndices) = inds
# const QuIndices = NTuple{5,QuIndex}
# # no need to define Base.tail for NTuple
# assignedinds(inds::QuIndices) = filter(!isnoindex,inds)

make_indices(inds::QuIndices) = inds
make_indices(inds::Union{Vector,Tuple}) = make_indices(inds...)
make_indices(inds...)::QuIndices = [QuIndex(i) for i in inds]
#make_indices(i1=NoIndex,i2=NoIndex,i3=NoIndex,i4=NoIndex,i5=NoIndex)::QuIndices = (QuIndex(i1),QuIndex(i2),QuIndex(i3),QuIndex(i4),QuIndex(i5))

const _NameTable = Dict{Symbol,IndexInt}()
const _NameTableInv = Symbol[]

@concrete struct QuOpName
    i::IndexInt
    function QuOpName(name::Symbol)
        i = get!(_NameTable,name) do
            push!(_NameTableInv,name)
            length(_NameTableInv)
        end
        new(i)
    end
    QuOpName(name::QuOpName) = name
end
sym(ind::QuOpName) = _NameTableInv[ind.i]
Base.print(io::IO, ind::QuOpName) = print(io, sym(ind))
Base.isless(i1::QuOpName,i2::QuOpName) = isless(sym(i1),sym(i2))
const NoName = QuOpName(Symbol())

# the enum also directly defines a natural ordering,so choose this directly how we later want it
# start the counting at 1 so we can index into the tuples defined below with Int(BaseOpType)
@enum BaseOpType        BosonCreate_=1 FermionCreate_   TLSCreate_   TLSx_  TLSy_  TLSz_  TLSDestroy_ FermionDestroy_ BosonDestroy_
const BaseOpType_adj = (BosonDestroy_, FermionDestroy_, TLSDestroy_, TLSx_, TLSy_, TLSz_, TLSCreate_, FermionCreate_, BosonCreate_)
const BaseOpType_sym  = ("†", "†", "⁺", "ˣ", "ʸ", "ᶻ", "⁻", "", "")
const BaseOpType_expr = ("ᴴ", "ᴴ", "⁺", "ˣ", "ʸ", "ᶻ", "⁻", "", "")
const BaseOpType_latex = ("^\\dagger", "^\\dagger", "^+", "^x", "^y", "^z", "^-", "", "")

@concrete struct BaseOperator
    t::BaseOpType
    name::QuOpName
    inds::QuIndices
    BaseOperator(t,name,inds...) = new(t,QuOpName(name),make_indices(inds...))
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
    iA::QuIndex
    iB::QuIndex
    δ(iA,iB) = (@assert !isnoindex(iA) && !isnoindex(iB); new(iA,iB))
end
function δ(Ainds::QuIndices,Binds::QuIndices)
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
    name::QuOpName
    state::Char
    inds::QuIndices
end

@concrete struct QuTerm
    nsuminds::IndexInt # have a sum over n indices, represented by QuIndex with issumindex(ind)==true
    δs::Vector{δ}
    params::Vector{Param}
    expvals::Vector{ExpVal}
    corrs::Vector{Corr}
    bares::BaseOpProduct
end
QuTerm() = QuTerm(BaseOpProduct())
QuTerm(op::BaseOperator) = QuTerm(BaseOpProduct([op]))
QuTerm(ops::BaseOpProduct) = QuTerm(0,δ[],Param[],ExpVal[],Corr[],ops)
QuTerm(δs::Vector{δ},ops::BaseOpProduct) = QuTerm(0,δs,Param[],ExpVal[],Corr[],ops)
QuTerm(p::Param)  = QuTerm(0,δ[],Param[p],ExpVal[],Corr[],BaseOpProduct())
QuTerm(E::ExpVal) = QuTerm(0,δ[],Param[],ExpVal[E],Corr[],BaseOpProduct())
QuTerm(C::Corr)   = QuTerm(0,δ[],Param[],ExpVal[],Corr[C],BaseOpProduct())
QuTerm(Es::Vector{ExpVal}) = QuTerm(0,δ[],Param[],Es,Corr[],BaseOpProduct())
QuTerm(Cs::Vector{Corr}) = QuTerm(0,δ[],Param[],ExpVal[],Cs,BaseOpProduct())
QuTerm(δs::Vector{δ},Es::Vector{ExpVal}) = QuTerm(0,δs,Param[],Es,Corr[],BaseOpProduct())
QuTerm(δs::Vector{δ},Cs::Vector{Corr})   = QuTerm(0,δs,Param[],ExpVal[],Cs,BaseOpProduct())

Base.isempty(A::QuTerm) = A.nsuminds == 0 && isempty(A.δs) && isempty(A.params) && isempty(A.expvals) && isempty(A.corrs) && isempty(A.bares)

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


#################################################################
## For the "external" functions that always construct QuExpr    ##
#################################################################
module OpConstructors
    using ..QuantumAlgebra: using_σpm, QuOpName, QuTerm, QuExpr, BaseOpProduct,
                            BosonDestroy, BosonCreate,
                            FermionDestroy, FermionCreate,
                            TLSDestroy, TLSCreate, TLSx, TLSy, TLSz

    "`boson_ops(name::Symbol)`: return functions for creating bosonic annihilation and creation operators with name `name` (i.e., wrappers of [`BosonDestroy`](@ref) and [`BosonCreate`](@ref))"
    function boson_ops(name::Symbol)
        op_name = QuOpName(name)
        namedag = Symbol(name,:dag)
        ann = @eval $name(   args...) = QuExpr(BosonDestroy($op_name,args...))
        cre = @eval $namedag(args...) = QuExpr(BosonCreate( $op_name,args...))
        ann, cre
    end

    "`fermion_ops(name::Symbol)`: return functions for creating fermionic annihilation and creation operators with name `name` (i.e., wrappers of [`FermionDestroy`](@ref) and [`FermionCreate`](@ref))"
    function fermion_ops(name::Symbol)
        op_name = QuOpName(name)
        namedag = Symbol(name,:dag)
        ann = @eval $name(   args...) = QuExpr(FermionDestroy($op_name,args...))
        cre = @eval $namedag(args...) = QuExpr(FermionCreate( $op_name,args...))
        ann, cre
    end

    "`tlspm_ops(name::Symbol)`: return functions for creating jump operators for a two-level system with name `name`. The output of these functions depends on setting of use_σpm."
    function tlspm_ops(name::Symbol)
        op_name = QuOpName(name)
        namep = Symbol(name,:p)
        namem = Symbol(name,:m)
        tlsm = @eval $namem(args...) = using_σpm() ? QuExpr(TLSDestroy($op_name,args...)) : QuExpr((QuTerm(TLSx($op_name,args...))=>1//2, QuTerm(TLSy($op_name,args...))=>-1im//2))
        tlsp = @eval $namep(args...) = using_σpm() ? QuExpr(TLSCreate( $op_name,args...)) : QuExpr((QuTerm(TLSx($op_name,args...))=>1//2, QuTerm(TLSy($op_name,args...))=>1im//2))
        tlsm, tlsp
    end

    "`tlsxyz_ops(name::Symbol)`: return functions for creating Pauli operators for a two-level system with name `name`. The output of these functions depends on setting of use_σpm."
    function tlsxyz_ops(name::Symbol)
        op_name = QuOpName(name)
        namex = Symbol(name,:x)
        namey = Symbol(name,:y)
        namez = Symbol(name,:z)
        tlsx = @eval $namex(args...) = using_σpm() ? QuExpr((QuTerm(TLSDestroy($op_name,args...))=>1, QuTerm(TLSCreate($op_name,args...))=>1)) : QuExpr(TLSx($op_name,args...))
        tlsy = @eval $namey(args...) = using_σpm() ? QuExpr((QuTerm(TLSDestroy($op_name,args...))=>1im, QuTerm(TLSCreate($op_name,args...))=>-1im)) : QuExpr(TLSy($op_name,args...))
        tlsz = @eval $namez(args...) = using_σpm() ? QuExpr((QuTerm(BaseOpProduct([TLSCreate($op_name,args...),TLSDestroy($op_name,args...)]))=>2, QuTerm()=>-1)) : QuExpr(TLSz($op_name,args...))
        tlsx, tlsy, tlsz
    end

    # default operators
    boson_ops(:a)
    fermion_ops(:f)
    tlspm_ops(:σ)
    tlsxyz_ops(:σ)
end
using .OpConstructors: boson_ops, fermion_ops, tlspm_ops, tlsxyz_ops
using .OpConstructors: a, adag, f, fdag, σx, σy, σz, σp, σm

"`@boson_ops name`: define functions `\$name` and `\$(name)dag` for creating bosonic annihilation and creation operators with name `name`"
macro boson_ops(name)
    :( ($(esc(name)), $(esc(Symbol(name,:dag)))) = boson_ops($(Meta.quot(name))) )
end

"`@fermion_ops name`: define functions `\$name` and `\$(name)dag` for creating fermionic annihilation and creation operators with name `name`"
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
        push!(code.args,:( $(esc(name))(args...) = ann($ii,args...) ))
        push!(code.args,:( $(esc(Symbol(name,:dag)))(args...) = cre($ii,args...) ))
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

## functions for constructing `param`s with string macros,
# Pc"ω_i,j" = param(:ω,'n',:i,:j) (complex parameter)
# Pr"ω_i,j" = param(:ω,'r',:i,:j) (real parameter)

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

function expval(A::QuTerm)
    if isempty(A.bares)
        A
    else
        QuTerm(A.nsuminds,A.δs,A.params,[A.expvals; ExpVal(A.bares)],A.corrs,BaseOpProduct())
    end
end
function corr(A::QuTerm)
    if isempty(A.bares)
        A
    else
        QuTerm(A.nsuminds,A.δs,A.params,A.expvals,[A.corrs; Corr(A.bares)],BaseOpProduct())
    end
end

"`expval(A::QuExpr)`: replace expression A by its (formal) expectation value ⟨A⟩."
expval(A::QuExpr) = _map_quexpr_ops(expval,A)
"`expval(A::QuExpr)`: replace expression A by its (formal) correlator ⟨A⟩c."
corr(A::QuExpr) = _map_quexpr_ops(corr,A)

"`∑(ind,A::QuExpr)`: return (formal) sum of expression A over index ind."
function ∑(ind::QuIndex,A::QuTerm)
    (issumindex(ind) || isintindex(ind)) && error("Index $ind to be summed over needs to be symbolic!")
    sumind = sumindex(A.nsuminds+one(A.nsuminds))
    f = replace_inds(ind=>sumind)
    g = reorder_suminds()
    # use the form that also sets nsuminds
    g(f(A,sumind.num))
end
∑(ind::QuIndex,A::QuExpr) = _map_quexpr_ops(t->∑(ind,t),A)
∑(ind::Symbol,A) = ∑(QuIndex(ind),A)

const QuantumObject = Union{QuIndex,QuOpName,BaseOperator,Param,BaseOpProduct,ExpVal,Corr,QuTerm,QuExpr}
