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
#   op_adjoint(::Val{F}, A)
#       the adjoint of the operator A
#   opsymbol(::Val{F}, A, style)
#       the symbol appended to the operator name when printing (style is :text, :expr or :latex)

@enum OpFamily::UInt8 BosonFamily FermionFamily TLSFamily QuditFamily

const ComplexInt = Complex{Int}

#                         BosonCreate_ FermionCreate_ TLSCreate_ QuditCreate_ TLSx_      TLSy_      TLSz_      QuditDestroy_ TLSDestroy_ FermionDestroy_ BosonDestroy_
const BaseOpType_family = (BosonFamily, FermionFamily, TLSFamily, QuditFamily, TLSFamily, TLSFamily, TLSFamily, QuditFamily,   TLSFamily,  FermionFamily,  BosonFamily)

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
    elseif fam == QuditFamily
        f(Val(QuditFamily), args...)
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

op_adjoint(A::BaseOperator) = _with_family(op_adjoint, opfamily(A), A)
opsymbol(A::BaseOperator, style::Symbol) = _with_family(opsymbol, opfamily(A), A, style)

"`annihilates_vacuum_left(A)`: whether ⟨0|A = 0"
annihilates_vacuum_left(A::BaseOperator) = (λ = vacuum_eigenvalue_left(A); λ !== nothing && iszero(λ))

# defaults
lookahead_contraction(::Val, ::BaseOpType) = false
function ladder(::Val, A::BaseOperator)
    t = optype(A)
    t in (BosonCreate_, FermionCreate_, TLSCreate_, QuditCreate_) ? 1 : t in (BosonDestroy_, FermionDestroy_, TLSDestroy_, QuditDestroy_) ? -1 : 0
end
op_adjoint(::Val, A::BaseOperator) = BaseOperator(BaseOpType_adj[Int(optype(A))],opname(A),inds(A),opextra(A))
function opsymbol(::Val, A::BaseOperator, style::Symbol)
    t = Int(optype(A))
    style === :text ? BaseOpType_sym[t] : style === :expr ? BaseOpType_expr[t] : BaseOpType_latex[t]
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
# Pauli matrices σˣ, σʸ, σᶻ (TLSx_, TLSy_, TLSz_). These cannot be mixed for
# operators that can act on the same two-level system.

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
    maybe_same_site(A,B) || return 1
    ispauli(ta) == ispauli(tb) || throw(ArgumentError("QuantumAlgebra currently does not support mixing 'normal' (x,y,z) Pauli matrices with jump operators (+,-) acting on the same two-level system."))
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

#################################################################
## qudits (d-level systems)                                    ##
#################################################################
# The operators of a d-level system are written in terms of the "creation" operators
# b_i† = |i⟩⟨1| and the "annihilation" operators b_i = |1⟩⟨i| (i = 2…d), where level 1 is
# the ground state (the vacuum). All other transition operators are products of these,
# |i⟩⟨j| = b_i† b_j for i,j ≥ 2, and |1⟩⟨1| = 1 - ∑_{k≥2} b_k† b_k. On the same qudit,
# b_i† b_j† = b_i b_j = 0 and b_j b_i† = δ_ij |1⟩⟨1|, so the normal-ordered products of
# operators on one qudit are 1, b_i†, b_j and b_i† b_j, i.e., d² operators that form a basis.
# For d = 2, b_2† and b_2 are equivalent to σ⁺ and σ⁻.
# The extra information stored with each operator is d:8 | i:8 | j:8 for |i⟩⟨j|.

_qudit_extra(d::Integer, i::Integer, j::Integer) = (d << 16) | (i << 8) | j
"`qudit_levels(A)`: (d,i,j) for a qudit operator A = |i⟩⟨j| (with i == 1 or j == 1) of a d-level system"
function qudit_levels(A::BaseOperator)
    x = opextra(A)
    (Int(x >> 16), Int((x >> 8) & 0xff), Int(x & 0xff))
end
# the level i ≥ 2 of b_i† = |i⟩⟨1| or b_i = |1⟩⟨i|
_qudit_level(A::BaseOperator) = ((_, i, j) = qudit_levels(A); max(i,j))

"`qudit_op(name, d, i, j, inds)`: the operator |i⟩⟨j| of a d-level system, with either i == 1 or j == 1 (but not both)"
function qudit_op(name::QuOpName, d::Integer, i::Integer, j::Integer, inds::AbstractVector{QuIndex})
    2 ≤ d ≤ 255 || throw(ArgumentError("the number of levels of a qudit has to be between 2 and 255, got $d"))
    (1 ≤ i ≤ d && 1 ≤ j ≤ d) || throw(ArgumentError("levels have to be between 1 and $d, got ($i,$j)"))
    (i == 1) ⊻ (j == 1) || throw(ArgumentError("only |i⟩⟨1| and |1⟩⟨i| (i ≥ 2) are basis operators, got |$i⟩⟨$j|"))
    BaseOperator(j == 1 ? QuditCreate_ : QuditDestroy_, name, inds, _qudit_extra(d,i,j))
end
# b_k† and b_k on the same qudit (and with the same indices) as A
_qudit_create_like(A::BaseOperator, k) = qudit_op(opname(A), qudit_levels(A)[1], k, 1, inds(A))
_qudit_destroy_like(A::BaseOperator, k) = qudit_op(opname(A), qudit_levels(A)[1], 1, k, inds(A))

function exchange!(::Val{QuditFamily}, emit, B, A)
    # operators of the same kind commute (on the same qudit, both orders give zero)
    optype(A) == optype(B) && return 1
    maybe_same_site(A,B) || return 1
    qudit_levels(A)[1] == qudit_levels(B)[1] || throw(ArgumentError("qudit $(opname(A)) is used with different numbers of levels"))
    # A = b_i† and B = b_j, and we need B A = b_j b_i†
    i, j = _qudit_level(A), _qudit_level(B)
    d = qudit_levels(A)[1]
    dd = site_deltas(A,B)
    if isempty(dd)
        # on the same qudit, b_j b_i† = δ_ij |1⟩⟨1| = δ_ij (1 - ∑_k b_k† b_k)
        i == j || return 0
        emit(1, dd)
        for k = 2:d
            k == i || emit(-1, dd, _qudit_create_like(A,k), _qudit_destroy_like(A,k))
        end
        return -1
    else
        # b_j b_i† = b_i† b_j + δ_site (δ_ij |1⟩⟨1| - |i⟩⟨j|) = b_i† b_j + δ_site (δ_ij (1 - ∑_k b_k† b_k) - b_i† b_j)
        if i == j
            emit(1, dd)
            for k = 2:d
                emit(k == i ? -2 : -1, dd, _qudit_create_like(A,k), _qudit_destroy_like(A,k))
            end
        else
            emit(-1, dd, _qudit_create_like(A,i), _qudit_destroy_like(A,j))
        end
        return 1
    end
end

# on the same qudit, b_i† b_j† = b_i b_j = 0
function contract!(::Val{QuditFamily}, emit, A, B)
    (optype(A) == optype(B) && opname(A) == opname(B) && inds(A) == inds(B)) || return nothing
    (zero(ComplexInt), nothing)
end

commutes(::Val{QuditFamily}, A, B) = optype(A) == optype(B) || !maybe_same_site(A,B)

# operators on the same qudit are not necessarily next to each other in the normal order (e.g., b₂†(1) b₂†(2) b₃†(1))
lookahead_contraction(::Val{QuditFamily}, ::BaseOpType) = true

op_adjoint(::Val{QuditFamily}, A::BaseOperator) = ((d,i,j) = qudit_levels(A); qudit_op(opname(A), d, j, i, inds(A)))

opsymbol(::Val{QuditFamily}, A::BaseOperator, style::Symbol) = ((_, i, j) = qudit_levels(A); _qudit_symbol(i, j, style))
function _qudit_symbol(i, j, style::Symbol)
    if style === :latex
        "^{$i$(max(i,j) ≥ 10 ? "," : "")$j}"
    elseif style === :expr
        "_$(i)_$(j)"
    else
        max(i,j) ≥ 10 ? "⁽$(superscript(i))'$(superscript(j))⁾" : "$(superscript(i))$(superscript(j))"
    end
end
