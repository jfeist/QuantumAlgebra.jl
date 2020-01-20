using DataStructures

export scal,param,a,adag,f,fdag,OpSumAnalytic,ExpVal,Corr
export @boson_ops,@fermion_ops
export σx,σy,σz,σp,σm
export @Pr_str, @Pc_str, ∑

# define for σx, σy, σz
@enum SpatialIndex x=1 y=2 z=3
SpatialIndex(a::SpatialIndex) = a

const OpIndex = Union{Int,Symbol}
const OpIndices = NTuple{N,OpIndex} where N

abstract type Operator end
abstract type Scalar <: Operator; end
"Represent a scalar value (i.e., a number)"
struct scal{T<:Number} <: Scalar; v::T; end
"`param(g,(:i,:j),'n')`: represent a scalar named parameter ``g_{i,j}``. state can be purely real (`'r'`), not conjugated (`'n'`), or conjugated (`'c'`)"
struct param{T<:OpIndices} <: Scalar
    name::Symbol
    state::Char
    inds::T
    function param(name,state::Char,inds::OpIndex...)
        state in ('n','r','c') || throw(ArgumentError("state has to be one of n,r,c"))
        new{typeof(inds)}(name,state,inds)
    end
    param(name,inds::OpIndex...) = param(name,'n',inds...)
    param(name,state::Char,inds::Tuple) = param(name,state,inds...)
    param(name,inds::Tuple) = param(name,'n',inds...)
end
struct δ <: Scalar
    inds::Tuple{OpIndex,OpIndex}
    δ(inds) = δ(inds...)
    function δ(iA::OpIndex,iB::OpIndex)
        # sort indices
        if iA == iB
            scal(1)
        elseif iA isa Integer && iB isa Integer
            # e.g., δ_1,3 = 0 (integers are not symbolic!)
            scal(0)
        elseif sortsentinel(iB) < sortsentinel(iA)
            new((iB,iA))
        else
            new((iA,iB))
        end
    end
end
function δ(Ainds::OpIndices,Binds::OpIndices)::Operator
    if length(Ainds) == length(Binds)
        prod(collect(Operator,δ.(Ainds,Binds)))
    else
        scal(0)
    end
end

abstract type BaseOperator <: Operator; end
for (op,desc) in (
    (:BosonDestroy,"bosonic annihilation"),
    (:BosonCreate,"bosonic creation"),
    (:FermionDestroy,"fermionic annihilation"),
    (:FermionCreate,"fermionic creation"))
    @eval begin
        "`$($op)(name,inds)`: represent $($desc) operator ``name_{inds}``"
        struct $op{T<:OpIndices} <: BaseOperator
            name::Symbol
            inds::T
            $op(name::Symbol,inds::OpIndex...) = new{typeof(inds)}(name,inds)
            $op(name::Symbol,inds::Tuple) = $op(name,inds...)
        end
    end
end

"`@boson_ops name`: return functions for creating bosonic annihilation and creation operators with name `name` (i.e., wrappers of [`BosonDestroy`](@ref) and [`BosonCreate`](@ref))"
macro boson_ops(name)
    quote
        ann(args...) = BosonDestroy($(Meta.quot(name)),args...)
        cre(args...) = BosonCreate($(Meta.quot(name)),args...)
        ann,cre
    end
end

"`@fermion_ops name`: return functions for creating fermionic annihilation and creation operators with name `name` (i.e., wrappers of [`FermionDestroy`](@ref) and [`FermionCreate`](@ref))"
macro fermion_ops(name)
    quote
        ann(args...) = FermionDestroy($(Meta.quot(name)),args...)
        cre(args...) = FermionCreate($(Meta.quot(name)),args...)
        ann,cre
    end
end

a, adag = @boson_ops a
f, fdag = @fermion_ops f

for (op,desc,sym) in (
    (:σminus,"TLS annihilation","σ^-"),
    (:σplus,"TLS creation","σ^+"))
    @eval begin
        "`$($op)(inds)`: represent $($desc) operator ``$($sym)_{inds}``"
        struct $op{T<:OpIndices} <: BaseOperator
            inds::T
            $op(inds::OpIndex...) = new{typeof(inds)}(inds)
            $op(inds::Tuple) = $op(inds...)
        end
    end
end


"`σ(a,inds)`: represent Pauli matrix ``σ_{a,inds}`` for two-level system (TLS), where ``a ∈ \\{x,y,z\\}``."
struct σ{T<:OpIndices} <: BaseOperator
    a::SpatialIndex
    inds::T
    σ(a,inds::OpIndex...) = new{typeof(inds)}(SpatialIndex(a),inds)
    σ(a,inds::Tuple) = σ(a,inds...)
end

struct OpProd <: Operator; A::Operator; B::Operator; end

struct OpSum  <: Operator
    # sum of Operators is saved as Dictionary of operators with scalar prefactors
    terms::SortedDict{Operator,scal}
    OpSum(A::Operator) = new(SortedDict(A => scal(1)))
    OpSum(A::scal) = new(SortedDict(scal(1) => A))
    OpSum(A::OpProd) = new(SortedDict(A.A isa scal ? A.B => A.A : A => scal(1)))
    OpSum(A::OpSum) = A
    function OpSum(itr)
        A = new(SortedDict{Operator,scal}())
        for (t,s) in itr
            _add_sum_term!(A,t,s)
        end
        _clean_sum(A)
    end
end
function _add_sum_term!(A::OpSum,oB::Operator,sB::scal)
    sB == scal(0) && return A
    sB += get(A.terms,oB,scal(0))
    if sB == scal(0)
        delete!(A.terms,oB)
    else
        A.terms[oB] = sB
    end
    A
end
# all scalar terms go into prefactor of operator scal(1)
_add_sum_term!(A::OpSum,oB::scal,sB::scal) = invoke(_add_sum_term!,Tuple{OpSum,Operator,scal},A,scal(1),oB*sB)

function _clean_sum(A)
    length(A.terms) == 0 && return scal(0)
    length(A.terms) == 1 && ((t,s) = first(A.terms); return s*t)
    A
end

# map function where scalar prefactors are unchanged
_map_opsum_ops(f,A::OpSum) = OpSum(f(t)=>s for (t,s) in A.terms)

"`OpSumAnalytic(i::Symbol,A::Operator)` or `∑(i,A)`: represent ``\\sum_{i} A``, with all possible values of ``i`` assumed to be included"
struct OpSumAnalytic <: Operator
    ind::Symbol
    A::Operator
    OpSumAnalytic(ind::Symbol,A::scal) = A.v == 0 ? A : new(ind,A)
    OpSumAnalytic(ind::Symbol,A::Operator) = begin
        if A isa OpProd && A.A isa scal
            A.A*OpSumAnalytic(ind,A.B)
        elseif A isa OpSumAnalytic && A.ind < ind
            # reorder nested sums by indices
            OpSumAnalytic(A.ind,OpSumAnalytic(ind,A.A))
        else
            for t in preftuple(A)
                # if there is a δ in the expression and it has the sum index, the sum disappears
                if t isa δ && ind in t.inds
                    irepl = ind==t.inds[1] ? t.inds[2] : t.inds[1]
                    return replace_index(A,ind,irepl)
                end
            end
            return new(ind,A)
        end
    end
    OpSumAnalytic(ind::Symbol,A::OpSum) = _map_opsum_ops(t->OpSumAnalytic(ind,t), A)
end

const ∑ = OpSumAnalytic

"`ExpVal(A::Operator)`: represent expectation value ``⟨A⟩``"
struct ExpVal <: Scalar
    A::Operator
    ExpVal(A::Scalar) = A
    ExpVal(A::OpSum) = _map_opsum_ops(ExpVal, A)
    ExpVal(A::OpProd) = A.A isa Scalar ? A.A*ExpVal(A.B) : new(A)
    ExpVal(A::OpSumAnalytic) = OpSumAnalytic(A.ind,ExpVal(A.A))
    ExpVal(A::Operator) = new(A)
end

"""
    Corr(A::Operator)

Represents correlations of operator ``A``. ``A`` should be a product for this to
make sense, in which case ``⟨AB⟩_c = ⟨AB⟩ - ⟨A⟩⟨B⟩``, with corresponding
extensions for products of more operators.

See also: [`ascorr`](@ref)
"""
struct Corr <: Scalar
    A::Operator
    Corr(A::Scalar) = A
    Corr(A::OpSum) = _map_opsum_ops(Corr, A)
    Corr(A::OpProd) = A.A isa Scalar ? A.A*Corr(A.B) : new(A)
    Corr(A::Operator) = new(A)
end

## functions for constructing Pauli operators depending on whether we use (σ+,σ-) or (σx,σy,σz) as the "basic" operators

const _using_σpm = Ref(false)
use_σpm(t::Bool=true) = (_using_σpm[] = t; nothing)
use_σxyz() = use_σpm(false)
using_σpm() = _using_σpm[]

"`σp(n)`: construct ``σ^+_n = \\frac12 σ_{x,n} + \\frac{i}{2} σ_{y,n}``"
σp(n...) = using_σpm() ? σplus(n...) : scal(1//2)*σ(x,n...) + scal(1im//2)*σ(y,n...)
"`σm(n)`: construct ``σ^-_n = \\frac12 σ_{x,n} - \\frac{i}{2} σ_{y,n}``"
σm(n...) = using_σpm() ? σminus(n...) : scal(1//2)*σ(x,n...) - scal(1im//2)*σ(y,n...)
"`σx(n)`: construct ``σ_{x,n}``"
σx(n...) = using_σpm() ? σminus(n...) + σplus(n...) : σ(x,n...)
"`σy(n)`: construct ``σ_{y,n}``"
σy(n...) = using_σpm() ? scal(1im) * (σminus(n...) - σplus(n...)) : σ(y,n...)
"`σz(n)`: construct ``σ_{z,n}``"
σz(n...) = using_σpm() ? scal(2)*σplus(n...)*σminus(n...) - scal(1) : σ(z,n...)

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
    param(par,'n',inds...)
end
macro Pr_str(s)
    par, inds = parse_paramstr(s)
    param(par,'r',inds...)
end
