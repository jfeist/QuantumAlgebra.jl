module QuantumAlgebra

export scal,param,a,adag,σ,σp,σm,OpSumAnalytic,ExpVal,Corr
export x,y,z,comm,latex
export Avac,vacA,vacExpVal
export CorrOrExp,ascorr
#export preftuple,exptuple,optuple,prodtuples,sumtuples

using Combinatorics
using Printf

# we will want to overload these operators and functions for our custom types
import Base: ==, *, +, -, isless, length, adjoint, repr

# define for σx, σy, σz
@enum SpatialIndex x=1 y=2 z=3
SpatialIndex(a::SpatialIndex) = a

abstract type Operator end
abstract type Scalar <: Operator; end
struct scal{T<:Number} <: Scalar; v::T; end
struct param{T<:Tuple} <: Scalar
    name::Symbol
    inds::T
    state::Char
    function param(name,inds::T=(),state='n') where T<:Tuple
        state in ('n','r','c') || error("state has to be one of n,r,c")
        new{T}(name,inds,state)
    end
    param(name,n,state='n') = param(name,(n,),state)
end
struct a{T}    <: Operator; n::T; end
struct adag{T} <: Operator; n::T; end
struct σ{T}    <: Operator
    a::SpatialIndex
    n::T
    σ(a,n::T) where T = new{T}(SpatialIndex(a),n)
end
struct OpProd <: Operator; A::Operator; B::Operator; end
struct OpSum  <: Operator; A::Operator; B::Operator; end

# represents a sum over the symbol "ind", with all values assumed to be included
struct OpSumAnalytic <: Operator
    ind::Symbol
    A::Operator
    OpSumAnalytic(ind::Symbol,A::Operator) = A == scal(0) ? A : new(ind,A)
    OpSumAnalytic(ind::Symbol,A::OpProd) = A.A isa scal ? A.A*OpSumAnalytic(ind,A.B) : new(ind,A)
    OpSumAnalytic(ind::Symbol,A::OpSum) = OpSumAnalytic(ind,A.A) + OpSumAnalytic(ind,A.B)
end

struct ExpVal <: Scalar
    A::Operator
    ExpVal(A::Scalar) = A
    ExpVal(A::OpSum) = ExpVal(A.A) + ExpVal(A.B)
    ExpVal(A::OpProd) = A.A isa Scalar ? A.A*ExpVal(A.B) : new(A)
    ExpVal(A::OpSumAnalytic) = OpSumAnalytic(A.ind,ExpVal(A.A))
    ExpVal(A::Operator) = new(A)
end

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
==(A::param, B::param) = (A.name,A.inds) == (B.name,B.inds) && (A.state==B.state=='c' || A.state !='c' && B.state != 'c')

prodtuple(A::Operator) = (A,)
prodtuple(A::OpProd) = (prodtuple(A.A)...,prodtuple(A.B)...)

# make a tuple from a product, containing only either prefactor types, expectation value types, or operators
for (name,types) in [(:pref,(scal,param)),(:exp,(ExpVal,Corr)),(:op,(adag,a,σ))]
    name = Symbol(name,:tuple)
    types = Union{types...}
    @eval $name(A::$types) = (A,)
    @eval $name(A::Operator) = ()
    @eval $name(A::OpSumAnalytic) = ($name(A.A)...,)
    @eval $name(A::OpProd) = ($name(A.A)...,$name(A.B)...)
end
prodtuples(A::Operator) = (preftuple(A), exptuple(A), optuple(A))

sumtuple(A::Operator) = (A,)
sumtuple(A::OpSum) = (sumtuple(A.A)...,sumtuple(A.B)...)

OpOrder = (scal,param,ExpVal,Corr,adag,a,σ,OpProd,OpSumAnalytic,OpSum)
for (ii,op1) in enumerate(OpOrder)
    for op2 in OpOrder[ii+1:end]
        @eval isless(::$op1,::$op2) = true
        @eval isless(::$op2,::$op1) = false
    end
end
isless(A::scal,B::scal) = (real(A.v),imag(A.v)) < (real(B.v),imag(B.v))
# do not order parameters by whether they are purely real, or conjugated or not
isless(A::param,B::param) = (A.name,A.inds) < (B.name,B.inds)
isless(A::σ{T},B::σ{T}) where T = (A.n,A.a) < (B.n,B.a)
isless(A::a{T},B::a{T}) where T = A.n < B.n
isless(A::adag{T},B::adag{T}) where T = A.n < B.n
for op in (a,adag,σ)
    @eval isless(A::$op{Symbol},B::$op{Int}) = false
    @eval isless(A::$op{Int},B::$op{Symbol}) = true
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

# define σ+ and σ- for some checks below
σp(n) = scal(1/2)*σ(x,n) + scal(1im/2)*σ(y,n)
σm(n) = scal(1/2)*σ(x,n) - scal(1im/2)*σ(y,n);

# prefactors do not count for length calculation
for op in [scal,param]
    @eval length(::$op) = 0
end
for op in [adag,a,σ]
    @eval length(::$op) = 1
end
for op in [ExpVal,Corr,OpSumAnalytic]
    @eval length(A::$op) = length(A.A)
end
length(A::OpProd) = length(A.A) + length(A.B)

*(A::Operator) = A # needed to make prod((A,)) work
*(A::OpProd,B::Operator) = A.A*(A.B*B)
*(A::OpSum,B::Operator) = A.A*B + A.B*B
*(A::OpSum,B::OpSum)    = A.A*B + A.B*B # resolve ambiguity
*(A::Operator,B::OpSum) = A*B.A + A*B.B
*(A::OpProd,  B::OpSum) = A*B.A + A*B.B

function *(A::Operator,B::Operator)::Operator
    if A isa scal && A.v==0
        A
    elseif A isa scal && A.v==1
        B
    elseif A isa scal && B isa scal
        scal(A.v*B.v)
    elseif A isa scal && B isa OpProd && B.A isa scal
        scal(A.v*B.A.v)*B.B
    elseif B isa OpProd && (A>B.A)
         # the (A*B.A) will get replaced by B.A*A+comm(A,B.A)
        (A*B.A)*B.B
    elseif A isa σ && B isa σ && A.n == B.n
        # σa σb = δab + i ϵabc σc = δab + 1/2 [σa,σb]
        A.a==B.a ? scal(1) : scal(1//2)*comm(A,B)
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
adjoint(A::param) = A.state=='r' ? A : param(A.name,A.inds,A.state=='n' ? 'c' : 'n')
adjoint(A::a) = adag(A.n)
adjoint(A::adag) = a(A.n)
adjoint(A::σ) = A
adjoint(A::OpSum) = A.A' + A.B'
adjoint(A::OpProd) = A.B' * A.A'
adjoint(A::ExpVal) = ExpVal(A.A')
adjoint(A::Corr) = Corr(A.A')
# since the term inside is fully simplified, we have normal ordering or single two-level operators
# the order reversion of the adjoint and subsequent simplification thus never needs commutations of non-commuting operators
# (which would be problematic because we would have to take into account that the sum index can be equal to any other index)
adjoint(A::OpSumAnalytic) = OpSumAnalytic(A.ind,A.A')

comm(A::Scalar,B::Operator) = scal(0)
comm(A::Scalar,B::OpProd)   = scal(0)
comm(A::Scalar,B::OpSum)    = scal(0)
comm(A::Operator,B::Scalar) = scal(0)
comm(A::OpSum,   B::Scalar) = scal(0)
comm(A::OpProd,  B::Scalar) = scal(0)
comm(A::Scalar,B::Scalar)   = scal(0)
comm(A::OpProd,B::Operator) = comm(A.A,B)*A.B + A.A*comm(A.B,B)
comm(A::OpProd,B::OpProd)   = comm(A.A,B)*A.B + A.A*comm(A.B,B)
comm(A::Operator,B::OpProd) = comm(A,B.A)*B.B + B.A*comm(A,B.B)
comm(A::OpSum,B::Operator)  = comm(A.A,B) + comm(A.B,B)
comm(A::OpSum,B::OpSum)     = comm(A.A,B) + comm(A.B,B)
comm(A::OpSum,B::OpProd)    = comm(A.A,B) + comm(A.B,B)
comm(A::Operator,B::OpSum)  = comm(A,B.A) + comm(A,B.B)
comm(A::OpProd,  B::OpSum)  = comm(A,B.A) + comm(A,B.B)
for op in (a,adag)
    @eval comm(A::$op,B::$op) = scal(0)
    @eval comm(A::$op,B::σ)   = scal(0)
    @eval comm(A::σ,B::$op)   = scal(0)
end
comm(A::a,B::adag) = scal(A.n==B.n ? 1 : 0)
comm(A::adag,B::a) = scal(A.n==B.n ? -1 : 0)

function comm(A::σ,B::σ)
    if A.n != B.n
        scal(0)
    elseif A.a == B.a
        scal(0)
    else
        a = Int(A.a)
        b = Int(B.a)
        # a+b+c == 6 (since a,b,c is a permutation of 1,2,3)
        c = 6 - a - b
        s = levicivita([a,b,c])
        scal(2im*s)*σ(c,A.n)
    end
end

replace_index(A::scal,iold,inew) = A
replace_index(A::param,iold,inew) = param(A.name,(n-> n==iold ? inew : n).(A.inds),A.state)
replace_index(A::ExpVal,iold,inew) = ExpVal(replace_index(A.A,iold,inew))
replace_index(A::Corr,iold,inew) = Corr(replace_index(A.A,iold,inew))
replace_index(A::a,iold,inew) = A.n==iold ? a(inew) : A
replace_index(A::adag,iold,inew) = A.n==iold ? adag(inew) : A
replace_index(A::σ,iold,inew) = A.n==iold ? σ(A.a,inew) : A
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

distribute_indices!(inds,A::scal) = A
distribute_indices!(inds,A::param) = param(A.name,((popfirst!(inds) for _ in A.inds)...,),A.state)
distribute_indices!(inds,A::ExpVal) = ExpVal(distribute_indices!(inds,A.A))
distribute_indices!(inds,A::Corr) = Corr(distribute_indices!(inds,A.A))
distribute_indices!(inds,A::a) = a(popfirst!(inds))
distribute_indices!(inds,A::adag) = adag(popfirst!(inds))
distribute_indices!(inds,A::σ) = σ(A.a,popfirst!(inds))
distribute_indices!(inds,A::OpProd) = distribute_indices!(inds,A.A)*distribute_indices!(inds,A.B)
distribute_indices!(inds,A::OpSum) = distribute_indices!(inds,A.A) + distribute_indices!(inds,A.B)
# on purpose do not define this for OpSumAnalytic

distribute_indices!([:a,:b,:c,:d,:e,:f,:g,:h],param(:ω,:y)*a(1)*adag(1)*a(3)*adag(:a))

# when multiplying (or commuting) with an operator with an index, take into account that the term in the sum with equal index has to be treated specially
*(A::Union{a,adag,σ},B::OpSumAnalytic) = (tmp = A*B.A; OpSumAnalytic(B.ind,tmp) - replace_index(tmp,B.ind,A.n) + A*replace_index(B.A,B.ind,A.n))
*(A::Union{param,ExpVal,Corr},B::OpSumAnalytic) = OpSumAnalytic(B.ind,A*B.A)
*(A::OpSumAnalytic,B::Union{a,adag,σ}) = (tmp = A.A*B; OpSumAnalytic(A.ind,tmp) - replace_index(tmp,A.ind,B.n) + replace_index(A.A,A.ind,B.n)*B)
*(A::OpSumAnalytic,B::OpProd) = (A*B.A)*B.B
*(A::OpSumAnalytic,B::Operator) = OpSumAnalytic(A.ind,A.A*B)

comm(A::Union{a,adag,σ},B::OpSumAnalytic) = (tmp = comm(A,B.A); OpSumAnalytic(B.ind,tmp) - replace_index(tmp,B.ind,A.n) + comm(A,replace_index(B.A,B.ind,A.n)))
comm(A::OpSumAnalytic,B::Union{a,adag,σ}) = (tmp = comm(A.A,B); OpSumAnalytic(A.ind,tmp) - replace_index(tmp,A.ind,B.n) + comm(replace_index(A.A,A.ind,B.n),B))

repr(A::a) = "a($(A.n))"
repr(A::adag) = "adag($(A.n))"
repr(A::σ) = "σ($(A.a),$(A.n))"
repr(A::scal) = "scal($(A.v))"
repr(A::param) = "param($(A.name),$(A.inds),$(A.state))"
repr(A::OpProd) = "OpProd($(A.A),$(A.B))"
repr(A::OpSum) = "OpSum($(A.A),$(A.B))"
repr(A::ExpVal) = "ExpVal($(A.A))"
repr(A::Corr) = "Corr($(A.A))"

mystring(x::Number) = @sprintf "%g" x
mystring(x::Rational) = denominator(x)==1 ? "$(numerator(x))" : "\\frac{$(numerator(x))}{$(denominator(x))}"
mystring(x::Complex) = @sprintf "(%g%+gi)" real(x) imag(x)
mystring(x::Complex{Rational{T}}) where T = @sprintf "\\left(%s%s%si\\right)" mystring(real(x)) (imag(x)>=0 ? "+" : "-") mystring(abs(imag(x)))

Base.show(io::IO, ::MIME"text/latex", A::Operator) where {T} = print(io,"\$",latex(A),"\$")
latex(A::a) = "a_{$(A.n)}"
latex(A::adag) = "a_{$(A.n)}^\\dagger"
latex(A::σ) = "\\sigma_{$(A.a),$(A.n)}"
latex(A::scal) = imag(A.v)==0 ? mystring(real(A.v)) : (real(A.v)==0 ? mystring(imag(A.v))*"i" : mystring(A.v))
latex(A::param)  = string(A.name, length(A.inds)==0 ? "" : "_{$(A.inds...)}", A.state=='c' ? "^*" : "")
latex(A::OpProd) = string(latex(A.A)," ",latex(A.B))
latex(A::OpSum) = string(latex(A.A)," + ",latex(A.B))
latex(A::ExpVal) = "\\langle $(latex(A.A)) \\rangle"
latex(A::Corr) = "\\langle $(latex(A.A)) \\rangle_{c}"
latex(A::OpSumAnalytic) = string("\\sum_{$(A.ind)}",latex(A.A))

indextuple(A::scal) = ()
indextuple(A::param) = A.inds
indextuple(A::Union{a,adag,σ}) = (A.n,)
indextuple(A::Union{OpProd,OpSum}) = (indextuple(A.A)...,indextuple(A.B)...)
indextuple(A::Union{ExpVal,Corr}) = indextuple(A.A)
indextuple(A::OpSumAnalytic) = (A.ind,indextuple(A.A)...)
indexset(A) = Set(indextuple(A))
sumindextuple(A::Operator) = ()
sumindextuple(A::Union{OpProd,OpSum}) = (sumindextuple(A.A)...,sumindextuple(A.B)...)
sumindextuple(A::Union{ExpVal,Corr}) = sumindextuple(A.A)
sumindextuple(A::OpSumAnalytic) = (A.ind,sumindextuple(A.A)...)
sumindexset(A) = Set(sumindextuple(A))

ascorr(A::Scalar) = A
for op in (adag,a,σ)
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
vacA(A::σ) = (A.a == z) ? scal(-1) : A
vacA(A::OpSum) = vacA(A.A) + vacA(A.B)
vacA(A::OpSumAnalytic) = OpSumAnalytic(A.ind,vacA(A.A))
vacA(A::OpProd) = begin
    if A.A isa Scalar
        A.A*vacA(A.B)
    else
        vA = vacA(A.A)
        vA==A.A ?  A : vacA(vA * A.B)
    end
end
vacA(A::Scalar) = A

# corresponds to |state> = stateop|0>, <state|A|state> = <0|state† A state|0>
function vacExpVal(stateop::Operator,A::Operator=scal(1))
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

end # module
