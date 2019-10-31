module QuantumAlgebra

export scal,param,a,adag,f,fdag,OpSumAnalytic,ExpVal,Corr
export σx,σy,σz,σp,σm,comm,latex
export Avac,vacA,vacExpVal
export CorrOrExp,ascorr
export @Pr_str, @Pc_str, ∑
#export σ,preftuple,exptuple,optuple,prodtuples,sumtuples

using Printf

# we will want to overload these operators and functions for our custom types
import Base: ==, ≈, *, +, -, isless, length, adjoint, print, zero, one

# define for σx, σy, σz
@enum SpatialIndex x=1 y=2 z=3
SpatialIndex(a::SpatialIndex) = a

# we sometimes want to sort integers and symbols together, with integers coming first
sortsentinel(x::Integer) = (1,x)
sortsentinel(x::Symbol) = (2,x)

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
        state in ('n','r','c') || error("state has to be one of n,r,c")
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
for (op,desc,sym) in (
    (:a,   "bosonic annihilation","a"),
    (:adag,"bosonic creation","a^†"),
    (:f,   "fermionic annihilation","f"),
    (:fdag,"fermionic creation","f^†"),
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

"`σ(a,n)`: represent Pauli matrix ``σ_{a,n}`` for two-level system (TLS) ``n``, where ``a ∈ \\{x,y,z\\}`` or ``\\{1,2,3\\}`` is the type of Pauli matrix."
struct σ{T<:OpIndices} <: BaseOperator
    a::SpatialIndex
    inds::T
    σ(a,inds::OpIndex...) = new{typeof(inds)}(SpatialIndex(a),inds)
    σ(a,inds::Tuple) = σ(a,inds...)
end

zero(::Type{<:Operator}) = scal(0)
zero(::Operator) = scal(0)
one(::Type{<:Operator}) = scal(1)
one(::Operator) = scal(1)

using_σpm = false
use_σpm(t::Bool=true) = eval(:( using_σpm = $t; nothing ))
use_σxyz() = use_σpm(false)

"`σp(n)`: construct ``σ^+_n = \\frac12 σ_{x,n} + \\frac{i}{2} σ_{y,n}``"
σp(n...) = using_σpm ? σplus(n...) : scal(1//2)*σ(x,n...) + scal(1im//2)*σ(y,n...)
"`σm(n)`: construct ``σ^-_n = \\frac12 σ_{x,n} - \\frac{i}{2} σ_{y,n}``"
σm(n...) = using_σpm ? σminus(n...) : scal(1//2)*σ(x,n...) - scal(1im//2)*σ(y,n...)
"`σx(n)`: construct ``σ_{x,n}``"
σx(n...) = using_σpm ? σminus(n...) + σplus(n...) : σ(x,n...)
"`σy(n)`: construct ``σ_{y,n}``"
σy(n...) = using_σpm ? scal(1im) * (σminus(n...) - σplus(n...)) : σ(y,n...)
"`σz(n)`: construct ``σ_{z,n}``"
σz(n...) = using_σpm ? scal(2)*σplus(n...)*σminus(n...) - scal(1) : σ(z,n...)

struct OpProd <: Operator; A::Operator; B::Operator; end
struct OpSum  <: Operator; A::Operator; B::Operator; end

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
    OpSumAnalytic(ind::Symbol,A::OpSum) = OpSumAnalytic(ind,A.A) + OpSumAnalytic(ind,A.B)
end

const ∑ = OpSumAnalytic

"`ExpVal(A::Operator)`: represent expectation value ``⟨A⟩``"
struct ExpVal <: Scalar
    A::Operator
    ExpVal(A::Scalar) = A
    ExpVal(A::OpSum) = ExpVal(A.A) + ExpVal(A.B)
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
    Corr(A::OpSum) = Corr(A.A) + Corr(A.B)
    Corr(A::OpProd) = A.A isa Scalar ? A.A*Corr(A.B) : new(A)
    Corr(A::Operator) = new(A)
end

==(A::scal,B::scal) = A.v == B.v
==(A::OpProd,B::OpProd) = A.A == B.A && A.B == B.B
==(A::OpSum, B::OpSum)  = A.A == B.A && A.B == B.B
==(A::OpSumAnalytic, B::OpSumAnalytic) = A.A == B.A && A.ind == B.ind
==(A::ExpVal, B::ExpVal) = A.A == B.A
==(A::Corr, B::Corr) = A.A == B.A
==(A::param, B::param) = (A.name,A.inds,A.state) == (B.name,B.inds,B.state)

≈(A::Operator,B::Operator) = A == B
≈(A::scal,B::scal) = A.v ≈ B.v
≈(A::OpProd,B::OpProd) = A.A ≈ B.A && A.B ≈ B.B
≈(A::OpSum, B::OpSum)  = A.A ≈ B.A && A.B ≈ B.B
≈(A::OpSumAnalytic, B::OpSumAnalytic) = A.A ≈ B.A && A.ind == B.ind
≈(A::ExpVal, B::ExpVal) = A.A ≈ B.A
≈(A::Corr, B::Corr) = A.A ≈ B.A

prodtuple(A::Operator) = (A,)
prodtuple(A::OpProd) = (prodtuple(A.A)...,prodtuple(A.B)...)

# make a tuple from a product, containing only either prefactor types, expectation value types, or operators
for (name,types) in [(:pref,(scal,param,δ)),(:exp,(ExpVal,Corr)),(:op,(adag,a,f,fdag,σ,σplus,σminus))]
    name = Symbol(name,:tuple)
    types = Union{types...}
    @eval $name(A::$types) = (A,)
    @eval $name(A::Operator) = ()
    @eval $name(A::OpSum) = error("cannot get $($name) for OpSum!")
    @eval $name(A::OpSumAnalytic) = ($name(A.A)...,)
    @eval $name(A::OpProd) = ($name(A.A)...,$name(A.B)...)
end
prodtuples(A::Operator) = (preftuple(A), exptuple(A), optuple(A))

sumtuple(A::Operator) = (A,)
sumtuple(A::OpSum) = (sumtuple(A.A)...,sumtuple(A.B)...)

OpOrder = (scal,δ,param,ExpVal,Corr,adag,a,fdag,f,σplus,σminus,σ,OpProd,OpSumAnalytic,OpSum)
for (ii,op1) in enumerate(OpOrder)
    for op2 in OpOrder[ii+1:end]
        @eval isless(::$op1,::$op2) = true
        @eval isless(::$op2,::$op1) = false
    end
end
isless(A::scal,B::scal) = (real(A.v),imag(A.v)) < (real(B.v),imag(B.v))

# do not order parameters by whether they are purely real, or conjugated or not
isless(A::param,B::param) = (A.name,sortsentinel.(A.inds)) < (B.name,sortsentinel.(B.inds))
isless(A::σ,B::σ) = (sortsentinel.(A.inds),A.a) < (sortsentinel.(B.inds),B.a)
for op in (δ,a,adag,f,fdag,σminus,σplus)
    @eval isless(A::$op,B::$op) = sortsentinel.(A.inds) < sortsentinel.(B.inds)
end
for op in (ExpVal,Corr,OpSumAnalytic)
    @eval isless(A::$op,B::$op) = A.A < B.A
end
function isless(A::OpProd,B::OpProd)
    # order operator products first by number of operators (also within expectation values),
    # then by operators, expectation values and correlations, then by reversed prefactors (to order params, not scalars). Again use tuples to write this easily
    Aprefs,Aexps,Aops = prodtuples(A)
    Bprefs,Bexps,Bops = prodtuples(B)
    (length(A),Aops,Aexps,reverse(Aprefs)) < (length(B),Bops,Bexps,reverse(Bprefs))
end

# prefactors do not count for length calculation
for op in [scal,δ,param]
    @eval length(::$op) = 0
end
for op in [adag,a,fdag,f,σ,σminus,σplus]
    @eval length(::$op) = 1
end
for op in [ExpVal,Corr,OpSumAnalytic]
    @eval length(A::$op) = length(A.A)
end
length(A::OpProd) = length(A.A) + length(A.B)

*(A::Operator) = A # needed to make prod((A,)) work
*(A::OpProd,B::Operator) = A.A*(A.B*B)
*(A::OpSum,B::Operator)      = A.A*B + A.B*B
*(A::OpSum,B::OpSum)         = A.A*B + A.B*B # resolve ambiguity
*(A::OpSum,B::OpSumAnalytic) = A.A*B + A.B*B # resolve ambiguity
*(A::Operator,     B::OpSum) = A*B.A + A*B.B
*(A::OpProd,       B::OpSum) = A*B.A + A*B.B # resolve ambiguity
*(A::OpSumAnalytic,B::OpSum) = A*B.A + A*B.B # resolve ambiguity

# allow addition, substraction, and multiplication with a number x by promoting it to scal(x) operator
*(x::Number,A::Operator) = scal(x)*A
*(A::Operator,x::Number) = scal(x)*A
+(x::Number,A::Operator) = scal(x)+A
+(A::Operator,x::Number) = scal(x)+A
-(x::Number,A::Operator) = scal(x)-A
-(A::Operator,x::Number) = A-scal(x)

function *(A::Operator,B::Operator)::Operator
    if A isa scal && A.v==0
        A
    elseif A isa scal && A.v==1
        B
    elseif A isa scal && B isa scal
        scal(A.v*B.v)
    elseif A isa scal && B isa OpProd && B.A isa scal
        scal(A.v*B.A.v)*B.B
    elseif B isa OpProd && (A*B.A) != OpProd(A,B.A)
        # if A*B.A is not ordered as we want, evaluate A*B.A first
        (A*B.A)*B.B
    elseif A isa δ && A.inds[2] in indextuple(B)
        # if second index of δ_iA,iB shows up on RHS, replace by iA
        # (relies on indices in δ being sorted)
        A * replace_index(B,A.inds[2],A.inds[1])
    elseif A isa σ && B isa σ && A.inds == B.inds
        # σa σb = δab + i ϵabc σc = δab + 1/2 [σa,σb]
        A.a==B.a ? scal(1) : scal(1//2)*comm(A,B)
    elseif A isa Union{σplus,σminus,f,fdag} && typeof(B)==typeof(A) && A.inds == B.inds
        scal(0)
    elseif A>B
        B*A + comm(A,B)
    else
        OpProd(A,B)
    end
end

# these are the only overrides for minus, everything else falls back to +
-(A::Operator) = scal(-1)*A
-(A::Operator,B::Operator) = A + scal(-1)*B

separate_prefac(A::Operator) = (1,A)
separate_prefac(A::OpProd) = A.A isa scal ? (A.A.v,A.B) : (1,A)
function combinable(A::Operator,B::Operator)
    sA, oA = separate_prefac(A)
    sB, oB = separate_prefac(B)
    if oA==oB
        true, scal(sA+sB)*oA
    else
        false, scal(0)
    end
end

+(A::OpSum,B::Operator) = A.A + (A.B + B)
function +(A::Operator,B::Operator)::Operator
    A isa scal && A.v == 0 && return B
    A isa scal && B isa scal && return scal(A.v+B.v)
    A isa scal && B isa OpSum && B.A isa scal && return scal(A.v+B.A.v) + B.B

    # ignore scalars in ordering
    sA, oA = separate_prefac(A)
    if B isa OpSum
        sB, oB = separate_prefac(B.A)
        oA>oB && return B.A + (A + B.B)
    else
        sB, oB = separate_prefac(B)
        oA>oB && return B + A
    end

    flag, S = combinable(A,B)
    flag && return S

    if B isa OpSum
        flag, S = combinable(A,B.A)
        flag && return S + B.B
    end
    return OpSum(A,B)
end

adjoint(A::scal) = scal(conj(A.v))
adjoint(A::param) = A.state=='r' ? A : param(A.name,A.state=='n' ? 'c' : 'n',A.inds)
adjoint(A::a) = adag(A.inds)
adjoint(A::adag) = a(A.inds)
adjoint(A::f) = fdag(A.inds)
adjoint(A::fdag) = f(A.inds)
adjoint(A::σminus) = σplus(A.inds)
adjoint(A::σplus) = σminus(A.inds)
adjoint(A::σ) = A
adjoint(A::OpSum) = A.A' + A.B'
adjoint(A::OpProd) = A.B' * A.A'
adjoint(A::ExpVal) = ExpVal(A.A')
adjoint(A::Corr) = Corr(A.A')
adjoint(A::OpSumAnalytic) = OpSumAnalytic(A.ind,A.A')

comm(A::Scalar,B::Operator) = scal(0)
comm(A::Scalar,B::OpProd)   = scal(0)
comm(A::Scalar,B::OpSum)    = scal(0)
comm(A::Scalar,B::OpSumAnalytic) = scal(0)
comm(A::Operator,B::Scalar) = scal(0)
comm(A::OpSum,   B::Scalar) = scal(0)
comm(A::OpProd,  B::Scalar) = scal(0)
comm(A::OpSumAnalytic,B::Scalar) = scal(0)
comm(A::Scalar,B::Scalar)   = scal(0)
comm(A::OpProd,B::Operator) = comm(A.A,B)*A.B + A.A*comm(A.B,B)
comm(A::OpProd,B::OpSumAnalytic) = comm(A.A,B)*A.B + A.A*comm(A.B,B)
comm(A::OpProd,B::OpProd)   = comm(A.A,B)*A.B + A.A*comm(A.B,B)
comm(A::Operator,B::OpProd) = comm(A,B.A)*B.B + B.A*comm(A,B.B)
comm(A::OpSumAnalytic, B::OpProd) = comm(A,B.A)*B.B + B.A*comm(A,B.B)
comm(A::OpSum,B::Operator)  = comm(A.A,B) + comm(A.B,B)
comm(A::OpSum,B::OpSum)     = comm(A.A,B) + comm(A.B,B)
comm(A::OpSum,B::OpProd)    = comm(A.A,B) + comm(A.B,B)
comm(A::Operator,B::OpSum)  = comm(A,B.A) + comm(A,B.B)
comm(A::OpProd,  B::OpSum)  = comm(A,B.A) + comm(A,B.B)
comm(A::OpSumAnalytic, B::OpSum) = comm(A,B.A) + comm(A,B.B)
# different types of operators commute
commgroups = (Union{a,adag},Union{f,fdag},Union{σ,σminus,σplus})
for (ii,op) in enumerate(commgroups)
    for op2 in commgroups[ii+1:end]
        @eval comm(A::$op,B::$op2) = scal(0)
        @eval comm(A::$op2,B::$op) = scal(0)
    end
end
# these operators commute with themselves
for op in (a,adag,σplus,σminus)
    @eval comm(A::$op,B::$op) = scal(0)
end
comm(A::a,B::adag) = δ(A.inds,B.inds)
comm(A::adag,B::a) = -comm(B,A)
comm(A::σplus,B::σminus) = δ(A.inds,B.inds)*σz(A.inds)
comm(A::σminus,B::σplus) = -comm(B,A)

# {f_i, fdag_j} = f_i fdag_j + fdag_j f_i = δ_{i,j}
# f_i fdag_j = δ_{i,j} - fdag_j f_i
# [f_i, fdag_j] = f_i fdag_j - fdag_j f_i = δ_{i,j} - fdag_j f_i - fdag_j f_i = δ_{i,j} - 2 fdag_j f_i
comm(A::f,B::fdag) = δ(A.inds,B.inds) - 2*B*A
comm(A::fdag,B::f) = -comm(B,A)
# {f_i,f_j} = f_i f_j + f_j f_i = 0
# we assume that if we need the commutator, it's because of exchanging order
# [f_i,f_j] = f_i f_j - f_j f_i = -f_j f_i - f_j f_i
comm(A::f,B::f) = -2*B*A
comm(A::fdag,B::fdag) = -2*B*A

# levicivita_lut[a,b] contains the Levi-Cevita symbol ϵ_abc
# for c=6-a-b, i.e, when a,b,c is a permutation of 1,2,3
const levicivita_lut = [0 1 -1; -1 0 1; 1 -1 0]

function comm(A::σ,B::σ)
    if A.a == B.a
        scal(0)
    else
        a = Int(A.a)
        b = Int(B.a)
        # a+b+c == 6 (since a,b,c is a permutation of 1,2,3)
        c = 6 - a - b
        s = levicivita_lut[a,b]
        δ(A.inds,B.inds)*scal(2im*s)*σ(c,A.inds)
    end
end

replace_index(A::scal,iold,inew) = A
replace_index(A::param,iold,inew) = param(A.name,A.state,(n-> n==iold ? inew : n).(A.inds))
replace_index(A::ExpVal,iold,inew) = ExpVal(replace_index(A.A,iold,inew))
replace_index(A::Corr,iold,inew) = Corr(replace_index(A.A,iold,inew))
for op in (δ,a,adag,f,fdag,σminus,σplus)
    @eval replace_index(A::$op,iold,inew) = $op((n-> n==iold ? inew : n).(A.inds))
end
replace_index(A::σ,iold,inew) = σ(A.a,(n-> n==iold ? inew : n).(A.inds))
replace_index(A::OpProd,iold,inew) = replace_index(A.A,iold,inew)*replace_index(A.B,iold,inew)
replace_index(A::OpSum,iold,inew) = replace_index(A.A,iold,inew) + replace_index(A.B,iold,inew)
replace_index(A::OpSumAnalytic,iold,inew) = begin
    (A.ind==iold || A.ind==inew) && error("replace_index in OpSumAnalytic cannot have iold ($iold) or inew ($inew) be the same as the sum index ($(A.ind))!")
    tmp = replace_index(A.A,iold,inew)
    # we have to be careful here - when replacing an index inside the expression, there might be reorderings that
    # have be done with the same index explicitly to make sure that the commutator shows up
    # so take out the term with the same index explicitly, and redo that one term
    OpSumAnalytic(A.ind,tmp) - replace_index(tmp,A.ind,inew) + replace_index(replace_index(A.A,A.ind,iold),iold,inew)
end
replace_index(A::Operator,reps) = begin
    Anew = A
    for (iold,inew) in reps
        Anew = replace_index(Anew,iold,inew)
    end
    Anew
end

function exchange_inds(A::Operator,i1,i2)
    ss = gensym()
    A = replace_index(A,i1,ss)
    A = replace_index(A,i2,i1)
    A = replace_index(A,ss,i2)
    A
end

distribute_indices!(inds,A::scal) = A
distribute_indices!(inds,A::param) = param(A.name,A.state,(popfirst!(inds) for _ in A.inds)...)
distribute_indices!(inds,A::ExpVal) = ExpVal(distribute_indices!(inds,A.A))
distribute_indices!(inds,A::Corr) = Corr(distribute_indices!(inds,A.A))
for op in (a,adag,f,fdag,σminus,σplus)
    @eval distribute_indices!(inds,A::$op) = $op((popfirst!(inds) for _ in A.inds)...)
end
distribute_indices!(inds,A::σ) = σ(A.a,(popfirst!(inds) for _ in A.inds)...)
distribute_indices!(inds,A::OpProd) = distribute_indices!(inds,A.A)*distribute_indices!(inds,A.B)
distribute_indices!(inds,A::OpSum) = distribute_indices!(inds,A.A) + distribute_indices!(inds,A.B)
# on purpose do not define this for OpSumAnalytic or δ

"sum indices have no semantic meaning, so rename them in case they happen to occur in the other expression"
function ensure_compatible_sumind(S::OpSumAnalytic,A::Operator)
    Ainds = indextuple(A)
    if S.ind in Ainds
        oldinds = Set{OpIndex}((Ainds...,indextuple(S)...))
        m = match(r"(.*)_([0-9]+)",string(S.ind))
        indstem, ii = (m === nothing) ? (string(S.ind), 1) : (m.captures[1], 1+parse(Int,m.captures[2]))
        while (newind = Symbol(indstem,:_,ii)) in oldinds
            ii += 1
        end
        OpSumAnalytic(newind,replace_index(S.A,S.ind,newind))
    else
        S
    end
end

# when multiplying (or commuting) with an operator with an index, take into account that the term in the sum with equal index has to be treated specially
*(A::Union{param,ExpVal,Corr,BaseOperator},B::OpSumAnalytic) = (B = ensure_compatible_sumind(B,A); OpSumAnalytic(B.ind,A*B.A))
# no need to check indices here since we just dispatch to another routine
*(A::OpSumAnalytic,B::OpProd) = (A*B.A)*B.B
*(A::OpSumAnalytic,B::Operator) = (A = ensure_compatible_sumind(A,B); OpSumAnalytic(A.ind,A.A*B))

comm(A::Union{param,ExpVal,Corr,BaseOperator},B::OpSumAnalytic) = (B = ensure_compatible_sumind(B,A); OpSumAnalytic(B.ind,comm(A,B.A)))
comm(A::OpSumAnalytic,B::OpSumAnalytic) = (B = ensure_compatible_sumind(B,A); OpSumAnalytic(B.ind,comm(A,B.A)))
comm(A::OpSumAnalytic,B::Operator) = -comm(B,A)

print(io::IO,A::a) = print(io,"a($(A.inds...))")
print(io::IO,A::adag) = print(io,"a†($(A.inds...))")
print(io::IO,A::f) = print(io,"f($(A.inds...))")
print(io::IO,A::fdag) = print(io,"f†($(A.inds...))")
print(io::IO,A::σ) = print(io,"σ$(A.a)($(A.inds...))")
print(io::IO,A::σminus) = print(io,"σ-($(A.inds...))")
print(io::IO,A::σplus) = print(io,"σ+($(A.inds...))")
print(io::IO,A::scal) = print(io,"$(A.v)")
print(io::IO,A::param) = print(io,string(A.name, A.state=='c' ? "'" : "",length(A.inds)==0 ? "" : "($(A.inds...))"))
print(io::IO,A::OpProd) = print(io,"$(A.A) $(A.B)")
print(io::IO,A::OpSum) = print(io,"$(A.A) + $(A.B)")
print(io::IO,A::ExpVal) = print(io,"⟨$(A.A)⟩")
print(io::IO,A::Corr) = print(io,"⟨$(A.A)⟩c")
print(io::IO,A::OpSumAnalytic) = print(io,"∑_$(A.ind) $(A.A)")

mystring(x::Number) = @sprintf "%g" x
mystring(x::Rational) = denominator(x)==1 ? "$(numerator(x))" : "\\frac{$(numerator(x))}{$(denominator(x))}"
mystring(x::Complex) = @sprintf "(%g%+gi)" real(x) imag(x)
mystring(x::Complex{Rational{T}}) where T = @sprintf "\\left(%s%s%si\\right)" mystring(real(x)) (imag(x)>=0 ? "+" : "-") mystring(abs(imag(x)))

Base.show(io::IO, ::MIME"text/latex", A::Operator) = print(io,"\$",latex(A),"\$")
latex(A::σ) = string("\\sigma_{$(A.a)",length(A.inds)>0 ? ",$(A.inds...)}" : "}")
latexindstr(inds::OpIndices) = length(inds)==0 ? "" : "_{$(inds...)}"
latex(A::δ) = "δ" * latexindstr(A.inds)
latex(A::a) = "a" * latexindstr(A.inds)
latex(A::adag) = "a$(latexindstr(A.inds))^\\dagger"
latex(A::f) = "f" * latexindstr(A.inds)
latex(A::fdag) = "f$(latexindstr(A.inds))^\\dagger"
latex(A::σminus) = "\\sigma^-" * latexindstr(A.inds)
latex(A::σplus) = "\\sigma^+" * latexindstr(A.inds)
latex(A::scal) = imag(A.v)==0 ? mystring(real(A.v)) : (real(A.v)==0 ? mystring(imag(A.v))*"i" : mystring(A.v))
latex(A::param)  = string(A.name, latexindstr(A.inds), A.state=='c' ? "^*" : "")
latex(A::OpProd) = string(latex(A.A)," ",latex(A.B))
latex(A::OpSum) = string(latex(A.A)," + ",latex(A.B))
latex(A::ExpVal) = "\\langle $(latex(A.A)) \\rangle"
latex(A::Corr) = "\\langle $(latex(A.A)) \\rangle_{c}"
latex(A::OpSumAnalytic) = string("\\sum_{$(A.ind)}",latex(A.A))

indextuple(A::scal)::OpIndices = ()
indextuple(A::Union{param,δ,a,adag,f,fdag,σ,σminus,σplus})::OpIndices = A.inds
indextuple(A::Union{OpProd,OpSum})::OpIndices = (indextuple(A.A)...,indextuple(A.B)...)
indextuple(A::Union{ExpVal,Corr})::OpIndices = indextuple(A.A)
indextuple(A::OpSumAnalytic)::OpIndices = (A.ind,indextuple(A.A)...)
indexset(A) = Set{OpIndex}(indextuple(A))
sumindextuple(A::Operator)::OpIndices = ()
sumindextuple(A::Union{OpProd,OpSum})::OpIndices = (sumindextuple(A.A)...,sumindextuple(A.B)...)
sumindextuple(A::Union{ExpVal,Corr})::OpIndices = sumindextuple(A.A)
sumindextuple(A::OpSumAnalytic)::OpIndices = (A.ind,sumindextuple(A.A)...)
sumindexset(A) = Set{OpIndex}(sumindextuple(A))

"`extindices(A::Operator)` return externally visible indices of an expression"
extindices(A::Operator) = [ind for ind in indextuple(A) if !in(ind,sumindextuple(A))]

"`symmetric_index_nums(A::Operator)` return sequence of numbers of exchange-symmetric indices"
function symmetric_index_nums(A::Operator)
    inds = extindices(A)
    Nsyms = [1]
    for ii=2:length(inds)
        if A == exchange_inds(A,inds[ii-1],inds[ii])
            Nsyms[end] += 1
        else
            push!(Nsyms,1)
        end
    end
    Nsyms
end

"""
    ascorr(expr::Operator)

Take an expression `expr=A B C + D E...` and write its expectation value in
terms of single-body expectation values ``⟨A⟩, ⟨B⟩, \\ldots``, and many-body
correlations ``⟨AB⟩_c, ⟨ABC⟩_c``, etc. Currently, up to fourth-order
correlations (i.e., products of four operators) are supported.

E.g., `ascorr(adag(:n)*a(:n))` returns ``⟨a^\\dagger_n a_n⟩_c + ⟨a^\\dagger_n⟩
⟨a_n⟩`` (which is equal to ``⟨a^\\dagger_n a_n⟩``), while
`ascorr(adag(:n)*a(:m)*a(:n))` returns ``\\langle a_{n}^\\dagger a_{m} a_{n}
\\rangle_{c} + \\langle a_{n}^\\dagger \\rangle \\langle a_{m} \\rangle \\langle
a_{n} \\rangle + \\langle a_{n}^\\dagger \\rangle \\langle a_{m} a_{n}
\\rangle_{c} + \\langle a_{m} \\rangle \\langle a_{n}^\\dagger a_{n}
\\rangle_{c} + \\langle a_{n} \\rangle \\langle a_{n}^\\dagger a_{m}
\\rangle_{c}``.

See also: [`ExpVal`](@ref), [`Corr`](@ref)"""
function ascorr end

ascorr(A::Scalar) = A
for op in (adag,a,f,fdag,σ,σminus,σplus)
    @eval ascorr(A::$op) = ExpVal(A)
end
ascorr(A::OpSum) = ascorr(A.A) + ascorr(A.B)
ascorr(A::OpSumAnalytic) = begin
    # first calculate the correlation for the term in the sum with the "bare" indices, which means that the sum index
    # is assumed to be distinct from the indices of the expressions
    # then, calculate the correction expression for the term where the sum index is the same as any of the non-summed indices,
    # and add the correction between that and the "incorrect" term where index identity is not taken into account
    tmp = ascorr(A.A)
    res = OpSumAnalytic(A.ind,tmp)
    for ii in setdiff(indexset(A),sumindexset(A))
        res += ascorr(replace_index(A.A,A.ind,ii)) - replace_index(tmp,A.ind,ii)
    end
    res
end
function ascorr(A::OpProd)::Operator
    A.A isa OpSumAnalytic && error("should not occurr!")
    A.A isa scal && A.B isa OpSumAnalytic && return A.A*ascorr(A.B)
    preftup, exptup, optup = prodtuples(A)
    pref = prod((scal(1),preftup...,exptup...))
    if length(optup)==0
        pref
    elseif length(optup)==1
        pref*ExpVal(optup[1])
    elseif length(optup)==2
        A,B = optup
        pref*(Corr(A*B) + ExpVal(A)*ExpVal(B))
    elseif length(optup)==3
        A,B,C = optup
        pref*(Corr(A*B*C) + ExpVal(A)*Corr(B*C) + ExpVal(B)*Corr(A*C) + ExpVal(C)*Corr(A*B) + ExpVal(A)*ExpVal(B)*ExpVal(C))
    elseif length(optup)==4
        A,B,C,D = optup
        pref*(Corr(A*B*C*D) + ExpVal(A)*ExpVal(B)*ExpVal(C)*ExpVal(D)
            + ExpVal(A)*Corr(B*C*D) + ExpVal(B)*Corr(A*C*D) + ExpVal(C)*Corr(A*B*D) + ExpVal(D)*Corr(A*B*C)
            + ExpVal(A)*ExpVal(B)*Corr(C*D) + Corr(A*B)*ExpVal(C)*ExpVal(D) + Corr(A*B)*Corr(C*D)
            + ExpVal(A)*ExpVal(C)*Corr(B*D) + Corr(A*C)*ExpVal(B)*ExpVal(D) + Corr(A*C)*Corr(B*D)
            + ExpVal(A)*ExpVal(D)*Corr(B*C) + Corr(A*D)*ExpVal(B)*ExpVal(C) + Corr(A*D)*Corr(B*C))
    else
        error("ERROR: Only correlations up to fourth order are implemented for now")
    end
end
CorrOrExp(A::Operator) = length(A)==1 ? ExpVal(A) : Corr(A)


Avac(A::a) = scal(0)
Avac(A::adag) = A
Avac(A::f) = scal(0)
Avac(A::fdag) = A
Avac(A::σminus) = scal(0)
Avac(A::σplus) = A
# vacuum is an eigenstate of σz
Avac(A::σ) = (A.a == z) ? scal(-1) : A
Avac(A::OpSum) = Avac(A.A) + Avac(A.B)
Avac(A::OpSumAnalytic) = OpSumAnalytic(A.ind,Avac(A.A))
Avac(A::OpProd) = begin
    Bv = Avac(A.B)
    Bv != A.B && return Avac(A.A*Bv)
    # if A.A and A.B commute, try what happens with Avac(A.B*A.A)
    if comm(A.A,A.B) == scal(0)
        Av = Avac(A.A)
        Av != A.A && return Avac(A.B*Av)
    end
    # did not find a way to simplify
    A
end
Avac(A::Scalar) = A

vacA(A::a) = A
vacA(A::adag) = scal(0)
vacA(A::f) = A
vacA(A::fdag) = scal(0)
vacA(A::σminus) = A
vacA(A::σplus) = scal(0)
# vacuum is an eigenstate of σz
vacA(A::σ) = (A.a == z) ? scal(-1) : A
vacA(A::OpSum) = vacA(A.A) + vacA(A.B)
vacA(A::OpSumAnalytic) = OpSumAnalytic(A.ind,vacA(A.A))
vacA(A::OpProd) = begin
    if A.A isa Scalar
        A.A*vacA(A.B)
    else
        # for applying ⟨0|AB, see if ⟨0|A changes A, and if so,
        # keep going with the new operator
        vA = vacA(A.A)
        vA==A.A ?  A : vacA(vA * A.B)
    end
end
vacA(A::Scalar) = A

"""
    Avac(A::Operator), vacA(A::Operator)

Simplify operator by assuming it is applied to the vacuum from the left or
right, respectively. To be precise, `Avac(A)` returns ``A'`` such that ``A'|0⟩ =
A|0⟩``, while `vacA(A)` does the same for ``⟨0|A``."""
Avac, vacA

"""
    vacExpVal(A::Operator,S::Operator=scal(1))

Calculate the vacuum expectation value ``⟨0|S^\\dagger A S|0⟩``, i.e., the
expectation value ``⟨ψ|A|ψ⟩`` for the state defined by ``|ψ⟩= S|0⟩```.
"""
function vacExpVal(A::Operator,stateop::Operator=scal(1))
    # simplify down as much as possible by applying vacuum from left and right
    vsAsv = vacA(Avac(stateop' * A * stateop))
    # only operators that should survive here as operators are σs
    _vacExpVal(vsAsv)
end
_vacExpVal(A::Scalar) = A
# the terms that survive until here have at most a single σ of any particle
# so it does not matter if we are inside a product
_vacExpVal(A::σ) = A.a == z ? scal(-1) : scal(0)
_vacExpVal(A::OpSum) = _vacExpVal(A.A) + _vacExpVal(A.B)
_vacExpVal(A::OpSumAnalytic) = OpSumAnalytic(A.ind,_vacExpVal(A.A))
# we know that the operators here commute (all a and a† have disappeared, and at most a single σ remaining for each particle)
_vacExpVal(A::OpProd) = _vacExpVal(A.A) * _vacExpVal(A.B)

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

include("precompile.jl")

end # module
