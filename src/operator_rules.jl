# Algebraic rules for the different kinds of operators.
#
# Every operator kind (BaseOpType) belongs to a family (bosons, fermions, two-level systems).
# Operators from different families always commute. Within a family, the rules are given by
# the methods (for a family F, with operators A ≤ B in the normal order)
#
#   exchange!(::Val{F}, emit, B, A)::Int
#       rewrite B A = x A B + R, where the "remainder" R is a sum of terms each passed as
#       emit(c, δs, ops...) (meaning c δs ops...), and x is returned
#   contract!(::Val{F}, emit, A, B)
#       simplify the product A B of two operators that have been brought next to each other by
#       normal ordering. Returns `nothing` if the product does not simplify, and otherwise
#       (fac, op) such that A B = fac op + R (op === nothing stands for the identity), where R is
#       again passed term by term to emit(c, δs, ops...)
#   commutes(::Val{F}, A, B)::Bool
#       whether A and B commute (i.e., exchange! would return 1 and not emit any terms)
#
# In addition, families can define (with the defaults given below)
#   lookahead_contraction(::Val{F}, t::BaseOpType)
#       whether operators of kind t can contract with operators that are not next to them in
#       the normal order (e.g., σˣ₁ σˣ₂ σʸ₁), and operators in between need to be checked
#   ladder(::Val{F}, A)
#       +1 for creation, -1 for annihilation and 0 for other operators
#   vacuum_eigenvalue_right(::Val{F}, A) / vacuum_eigenvalue_left(::Val{F}, A)
#       λ if A|0⟩ = λ|0⟩ (⟨0|A = λ⟨0|), including λ=0, and `nothing` otherwise

@enum OpFamily::UInt8 BosonFamily FermionFamily TLSFamily

const ComplexInt = Complex{Int}

#                           BosonCreate_ FermionCreate_ TLSCreate_ TLSx_     TLSy_     TLSz_     TLSDestroy_ FermionDestroy_ BosonDestroy_
const BaseOpType_family = (BosonFamily, FermionFamily, TLSFamily, TLSFamily, TLSFamily, TLSFamily, TLSFamily, FermionFamily, BosonFamily)

opfamily(t::BaseOpType) = @inbounds BaseOpType_family[Int(t)]
opfamily(A::BaseOperator) = opfamily(optype(A))

# call f(Val(fam), args...), with the family known at compile time in each branch
@inline function _with_family(f::F, fam::OpFamily, args...) where F
    if fam == BosonFamily
        f(Val(BosonFamily), args...)
    elseif fam == FermionFamily
        f(Val(FermionFamily), args...)
    elseif fam == TLSFamily
        f(Val(TLSFamily), args...)
    else
        error("unknown operator family $fam")
    end
end

#################################################################
## generic entry points                                        ##
#################################################################

"`exchange!(emit, B, A)`: rewrite `B A = x A B + R`, with `A < B`. Passes the terms of `R` to `emit(c, δs, ops...)` and returns `x`."
@inline function exchange!(emit, B::BaseOperator, A::BaseOperator)
    fam = opfamily(A)
    opfamily(B) == fam || return 1
    _with_family(exchange!, fam, emit, B, A)
end

"`contract!(emit, A, B)`: simplify the product `A B`. Returns `nothing` if it does not simplify, or `(fac, op)` with `A B = fac op + R` (see `exchange!` for `R`)."
@inline function contract!(emit, A::BaseOperator, B::BaseOperator)
    fam = opfamily(A)
    opfamily(B) == fam || return nothing
    _with_family(contract!, fam, emit, A, B)
end

"`commutes(A, B)`: whether operators `A` and `B` commute."
@inline function commutes(A::BaseOperator, B::BaseOperator)
    fam = opfamily(A)
    opfamily(B) == fam || return true
    _with_family(commutes, fam, A, B)
end

lookahead_contraction(A::BaseOperator) = _with_family(lookahead_contraction, opfamily(A), optype(A))
ladder(A::BaseOperator) = _with_family(ladder, opfamily(A), A)
vacuum_eigenvalue_right(A::BaseOperator) = _with_family(vacuum_eigenvalue_right, opfamily(A), A)
vacuum_eigenvalue_left(A::BaseOperator) = _with_family(vacuum_eigenvalue_left, opfamily(A), A)

"`annihilates_vacuum_left(A)`: whether ⟨0|A = 0"
annihilates_vacuum_left(A::BaseOperator) = (λ = vacuum_eigenvalue_left(A); λ !== nothing && iszero(λ))

# defaults
lookahead_contraction(::Val, ::BaseOpType) = false
function ladder(::Val, A::BaseOperator)
    t = optype(A)
    t in (BosonCreate_, FermionCreate_, TLSCreate_) ? 1 : t in (BosonDestroy_, FermionDestroy_, TLSDestroy_) ? -1 : 0
end
# for ladder operators, the vacuum is annihilated by the annihilation operators (and their adjoints from the left)
vacuum_eigenvalue_right(F::Val, A::BaseOperator) = ladder(F, A) < 0 ? 0 : nothing
vacuum_eigenvalue_left(F::Val, A::BaseOperator) = ladder(F, A) > 0 ? 0 : nothing

# whether two operators can act on the same mode (i.e., all their indices can be equal)
function maybe_same_site(A::BaseOperator, B::BaseOperator)
    opname(A) == opname(B) || return false
    Ainds, Binds = inds(A), inds(B)
    length(Ainds) == length(Binds) || return false
    for (iA,iB) in zip(Ainds,Binds)
        iA == iB && continue
        (isnoindex(iA) || isnoindex(iB)) && return false
        # if they are integer indices and different, they cannot be equal
        isintindex(iA) && isintindex(iB) && return false
    end
    true
end
# the δs that make two operators act on the same mode (only call if maybe_same_site(A,B))
site_deltas(A::BaseOperator, B::BaseOperator) = δ(inds(A), inds(B))::Vector{δ}

#################################################################
## bosons: [a_i, a†_j] = δ_ij                                  ##
#################################################################

function exchange!(::Val{BosonFamily}, emit, B, A)
    # a_i a†_j = a†_j a_i + δ_ij
    if optype(A) == BosonCreate_ && optype(B) == BosonDestroy_ && maybe_same_site(A,B)
        emit(1, site_deltas(A,B))
    end
    return 1
end
contract!(::Val{BosonFamily}, emit, A, B) = nothing
commutes(::Val{BosonFamily}, A, B) = optype(A) == optype(B) || !maybe_same_site(A,B)

#################################################################
## fermions: {f_i, f†_j} = δ_ij                                ##
#################################################################
# operators with different names commute, while operators with the same name anticommute
# (use @anticommuting_fermion_group to create several species that anticommute)

function exchange!(::Val{FermionFamily}, emit, B, A)
    opname(A) == opname(B) || return 1
    # f_i f†_j = -f†_j f_i + δ_ij
    if optype(A) != optype(B) && maybe_same_site(A,B)
        emit(1, site_deltas(A,B))
    end
    return -1
end
# f_i f_i = f†_i f†_i = 0
contract!(::Val{FermionFamily}, emit, A, B) = A == B ? (zero(ComplexInt), nothing) : nothing
commutes(::Val{FermionFamily}, A, B) = opname(A) != opname(B)

#################################################################
## two-level systems                                           ##
#################################################################
# either described by σ⁺ and σ⁻ (TLSCreate_ and TLSDestroy_), or by the
# Pauli matrices σˣ, σʸ, σᶻ (TLSx_, TLSy_, TLSz_). These cannot be mixed.

ispauli(t::BaseOpType) = t in (TLSx_,TLSy_,TLSz_)

# levicivita_lut[a][b] contains the Levi-Civita symbol ϵ_abc
# for c=6-a-b, i.e, when a,b,c is a permutation of 1,2,3
const levicivita_lut = ((0,1,-1), (-1,0,1), (1,-1,0))
"`ϵ_ab(A,B)`: for Pauli matrices σᵃ, σᵇ (with a≠b), return the third one, σᶜ, and the Levi-Civita symbol ϵ_abc"
function ϵ_ab(A::BaseOperator,B::BaseOperator)
    # a+b+c == 6 (since a,b,c is a permutation of 1,2,3)
    a = Int(optype(A)) - Int(TLSx_) + 1
    b = Int(optype(B)) - Int(TLSx_) + 1
    c = BaseOpType(Int(TLSx_) - 1 + (6 - a - b))
    s = @inbounds levicivita_lut[a][b]
    c, s
end

function exchange!(::Val{TLSFamily}, emit, B, A)
    ta, tb = optype(A), optype(B)
    ta == tb && return 1
    ispauli(ta) == ispauli(tb) || throw(ArgumentError("QuantumAlgebra currently does not support mixing 'normal' (x,y,z) Pauli matrices with jump operators (+,-)."))
    maybe_same_site(A,B) || return 1
    dd = site_deltas(A,B)
    if ta == TLSCreate_ && tb == TLSDestroy_
        if isempty(dd)
            # σ⁻ σ⁺ = 1 - σ⁺ σ⁻
            emit(1, dd)
            return -1
        else
            # σ⁻_i σ⁺_j = σ⁺_j σ⁻_i - δ_ij σᶻ_i = σ⁺_j σ⁻_i + δ_ij (1 - 2 σ⁺_i σ⁻_i)
            emit(-2, dd, TLSCreate(opname(A),inds(A)), TLSDestroy(opname(A),inds(A)))
            emit(1, dd)
            return 1
        end
    else
        # both are Pauli matrices, we need ϵ_bac
        c, s = ϵ_ab(B,A)
        σc = BaseOperator(c,opname(A),inds(A))
        if isempty(dd)
            # indices are all the same, σᵇ σᵃ = δ_ab + i ϵ_bac σᶜ
            # since A < B, we know that a != b
            emit(im*s, dd, σc)
            return 0
        else
            # [σᵇ_j,σᵃ_i] = 2i ϵ_bac δ_ij σᶜ_i => σᵇ_j σᵃ_i = σᵃ_i σᵇ_j + 2i ϵ_bac δ_ij σᶜ_i
            emit(2im*s, dd, σc)
            return 1
        end
    end
end

function contract!(::Val{TLSFamily}, emit, A, B)
    ta, tb = optype(A), optype(B)
    if ta in (TLSCreate_,TLSDestroy_)
        # σ⁺ σ⁺ = σ⁻ σ⁻ = 0
        A == B ? (zero(ComplexInt), nothing) : nothing
    elseif ispauli(ta) && ispauli(tb) && opname(A) == opname(B) && inds(A) == inds(B)
        # σᵃ σᵇ = δ_ab + i ϵ_abc σᶜ
        if ta == tb
            (one(ComplexInt), nothing)
        else
            c, s = ϵ_ab(A,B)
            (ComplexInt(0,s), BaseOperator(c,opname(A),inds(A)))
        end
    else
        nothing
    end
end

commutes(::Val{TLSFamily}, A, B) = optype(A) == optype(B) || !maybe_same_site(A,B)

lookahead_contraction(::Val{TLSFamily}, t::BaseOpType) = ispauli(t)

# σᶻ|0⟩ = -|0⟩, where |0⟩ is the ground state
vacuum_eigenvalue_right(F::Val{TLSFamily}, A::BaseOperator) = optype(A) == TLSz_ ? -1 : ladder(F, A) < 0 ? 0 : nothing
vacuum_eigenvalue_left(F::Val{TLSFamily}, A::BaseOperator) = optype(A) == TLSz_ ? -1 : ladder(F, A) > 0 ? 0 : nothing
