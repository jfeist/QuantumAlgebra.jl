# we will want to overload these operators and functions for our custom types
import Base: ==, ≈, *, +, -, isless, length, adjoint, print, zero, one
export comm

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

# we sometimes want to sort integers and symbols together, with integers coming first
sortsentinel(x::Integer) = (1,x)
sortsentinel(x::Symbol) = (2,x)

OpOrder = (scal,δ,param,ExpVal,Corr,BosonCreate,BosonDestroy,FermionCreate,FermionDestroy,σplus,σminus,σ,OpProd,OpSumAnalytic,OpSum)
for (ii,op1) in enumerate(OpOrder)
    for op2 in OpOrder[ii+1:end]
        @eval isless(::$op1,::$op2) = true
        @eval isless(::$op2,::$op1) = false
    end
end
isless(A::scal,B::scal) = (real(A.v),imag(A.v)) < (real(B.v),imag(B.v))

for op in (BosonCreate,BosonDestroy,FermionCreate,FermionDestroy)
    @eval isless(A::$op,B::$op) = (A.name,sortsentinel.(A.inds)) < (B.name,sortsentinel.(B.inds))
end
# use conjugation state as last order parameter
isless(A::param,B::param) = (A.name,sortsentinel.(A.inds),A.state) < (B.name,sortsentinel.(B.inds),B.state)
isless(A::σ,B::σ) = (sortsentinel.(A.inds),A.a) < (sortsentinel.(B.inds),B.a)
for op in (δ,σminus,σplus)
    @eval isless(A::$op,B::$op) = sortsentinel.(A.inds) < sortsentinel.(B.inds)
end
for op in (ExpVal,Corr)
    @eval isless(A::$op,B::$op) = A.A < B.A
end
isless(A::OpSumAnalytic,B::OpSumAnalytic) = (A.A,A.ind) < (B.A,B.ind)

# return -1 if A<B, 0 if A==B, 1 if B<A
function iterlesseq(A,B)::Int
    resa = iterate(A)
    resb = iterate(B)
    while resa !== nothing && resb !== nothing
        vala, statea = resa
        valb, stateb = resb
        vala == valb || return isless(vala,valb) ? -1 : 1
        resa = iterate(A,statea)
        resb = iterate(B,stateb)
    end
    resa === nothing && resb === nothing && return 0
    # the longer iterator is larger
    resa === nothing ? 1 : -1
end

function isless(A::OpProd,B::OpProd)
    # only evaluate each part that we need for each step of the comparison to avoid unnecessary work
    # order operator products first by number of operators (also within expectation values)
    lA,lB = length(A), length(B)
    lA == lB || return lA < lB
    # then by operators
    c1 = iterlesseq(opiter(A,true),opiter(B,true))
    c1 == 0 || return c1 < 0
    # then by expectation values
    c2 = iterlesseq(expiter(A,true),expiter(B,true))
    c2 == 0 || return c2 < 0
    # then by reversed prefactors (to order params, not scalars)
    rpA = reverse(preftuple(A))
    rpB = reverse(preftuple(B))
    rpA == rpB || return rpA < rpB
    # finally without swallowing opsumanalytic to make sure we always have definite ordering
    iterlesseq(proditer(A,false),proditer(B,false)) < 0
end

# prefactors do not count for length calculation
length(::Union{scal,δ,param}) = 0
length(::BaseOperator) = 1
length(A::Union{ExpVal,Corr,OpSumAnalytic}) = length(A.A)
length(A::OpProd) = length(A.A) + length(A.B)

adjoint(A::scal) = scal(conj(A.v))
adjoint(A::param) = A.state=='r' ? A : param(A.name,A.state=='n' ? 'c' : 'n',A.inds)
adjoint(A::δ) = A
adjoint(A::BosonDestroy) = BosonCreate(A.name,A.inds)
adjoint(A::BosonCreate) = BosonDestroy(A.name,A.inds)
adjoint(A::FermionDestroy) = FermionCreate(A.name,A.inds)
adjoint(A::FermionCreate) = FermionDestroy(A.name,A.inds)
adjoint(A::σminus) = σplus(A.inds)
adjoint(A::σplus) = σminus(A.inds)
adjoint(A::σ) = A
adjoint(A::OpSum) = A.A' + A.B'
adjoint(A::OpProd) = A.B' * A.A'
adjoint(A::ExpVal) = ExpVal(A.A')
adjoint(A::Corr) = Corr(A.A')
adjoint(A::OpSumAnalytic) = OpSumAnalytic(A.ind,A.A')

zero(::Type{<:Operator}) = scal(0)
zero(::Operator) = scal(0)
one(::Type{<:Operator}) = scal(1)
one(::Operator) = scal(1)

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
    elseif A isa Union{σplus,σminus,FermionDestroy,FermionCreate} && A == B
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

+(A::OpSum,B::Operator) = A.A + (A.B + B)
function +(A::Operator,B::Operator)::Operator
    A isa scal && A.v == 0 && return B
    A isa scal && B isa scal && return scal(A.v+B.v)
    A isa scal && B isa OpSum && B.A isa scal && return scal(A.v+B.A.v) + B.B

    # ignore scalars in ordering
    sA, oA = separate_prefac(A)
    if B isa OpSum
        sB, oB = separate_prefac(B.A)
        if oA > oB
            return B.A + (A + B.B)
        elseif oA == oB
            return scal(sA+sB)*oA + B.B
        end
    else
        sB, oB = separate_prefac(B)
        if oA > oB
            return B + A
        elseif oA == oB
            return scal(sA+sB)*oA 
        end
    end
    oA < oB || error("invariant oA < oB violated! oA: $oA, oB: $oB, oA<oB: $(oA<oB), oA==oB: $(oA==oB), oA>oB: $(oA>oB)")

    return OpSum(A,B)
end

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
commgroups = (Union{BosonDestroy,BosonCreate},Union{FermionDestroy,FermionCreate},Union{σ,σminus,σplus})
for (ii,op) in enumerate(commgroups)
    for op2 in commgroups[ii+1:end]
        @eval comm(A::$op,B::$op2) = scal(0)
        @eval comm(A::$op2,B::$op) = scal(0)
    end
end
# these operators commute with themselves
for op in (BosonDestroy,BosonCreate,σplus,σminus)
    @eval comm(A::$op,B::$op) = scal(0)
end
comm(A::BosonDestroy,B::BosonCreate) = A.name == B.name ? δ(A.inds,B.inds) : scal(0)
comm(A::BosonCreate,B::BosonDestroy) = -comm(B,A)
comm(A::σplus,B::σminus) = δ(A.inds,B.inds)*σz(A.inds)
comm(A::σminus,B::σplus) = -comm(B,A)

# {f_i, fdag_j} = f_i fdag_j + fdag_j f_i = δ_{i,j}
# f_i fdag_j = δ_{i,j} - fdag_j f_i
# [f_i, fdag_j] = f_i fdag_j - fdag_j f_i = δ_{i,j} - fdag_j f_i - fdag_j f_i = δ_{i,j} - 2 fdag_j f_i
# if they correspond to operators on different systems (i.e., with different names), also Fermion operators commute
comm(A::FermionDestroy,B::FermionCreate) = A.name == B.name ? δ(A.inds,B.inds) - 2*B*A : scal(0)
comm(A::FermionCreate,B::FermionDestroy) = -comm(B,A)
# {f_i,f_j} = f_i f_j + f_j f_i = 0
# we assume that if we need the commutator, it's because of exchanging order
# [f_i,f_j] = f_i f_j - f_j f_i = -f_j f_i - f_j f_i
comm(A::FermionDestroy,B::FermionDestroy) = A.name == B.name ? -2*B*A : scal(0)
comm(A::FermionCreate,B::FermionCreate) = A.name == B.name ? -2*B*A : scal(0)

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

# when multiplying (or commuting) with an operator with an index, take into account that the term in the sum with equal index has to be treated specially
*(A::Union{param,ExpVal,Corr,BaseOperator},B::OpSumAnalytic) = (B = ensure_compatible_sumind(B,A); OpSumAnalytic(B.ind,A*B.A))
# no need to check indices here since we just dispatch to another routine
*(A::OpSumAnalytic,B::OpProd) = (A*B.A)*B.B
*(A::OpSumAnalytic,B::Operator) = (A = ensure_compatible_sumind(A,B); OpSumAnalytic(A.ind,A.A*B))

comm(A::Union{param,ExpVal,Corr,BaseOperator},B::OpSumAnalytic) = (B = ensure_compatible_sumind(B,A); OpSumAnalytic(B.ind,comm(A,B.A)))
comm(A::OpSumAnalytic,B::OpSumAnalytic) = (B = ensure_compatible_sumind(B,A); OpSumAnalytic(B.ind,comm(A,B.A)))
comm(A::OpSumAnalytic,B::Operator) = -comm(B,A)
