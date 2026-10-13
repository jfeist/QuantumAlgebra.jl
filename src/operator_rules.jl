# Algebraic rules for the different kinds of operators.
#
# Every operator kind (BaseOpType) belongs to a family (bosons, fermions, two-level systems, ...).
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
#   lookahead_contraction(::Val{F}, A)
#       whether A can contract with operators that are not next to it in the normal order (e.g.,
#       σˣ₁ σˣ₂ σʸ₁). By default, the operators in between have to commute with the later one.
#   lookahead_meet(::Val{F}, A)
#       if true, the operators in between can also commute with A instead (so that A and the
#       other operator meet somewhere in between)
#   lookahead_exchange(::Val{F}, A)
#       if true, the later operator is also moved past operators in between that it does not
#       commute with (with exchange!, which has to return x = 1 for them, and the remainders are
#       added as further terms), so that operators that can contract are always contracted
#   ladder(::Val{F}, A)
#       +1 for creation, -1 for annihilation and 0 for other operators
#   vacuum_eigenvalue_right(::Val{F}, A) / vacuum_eigenvalue_left(::Val{F}, A)
#       λ if A|0⟩ = λ|0⟩ (⟨0|A = λ⟨0|), including λ=0, and `nothing` otherwise. λ is either an
#       Int, or an object x (of a concrete type) for which _vacuum_value(x) gives the number
#   op_adjoint(::Val{F}, A)
#       the adjoint of the operator A
#   opsymbol(::Val{F}, A, style)
#       the symbol appended to the operator name when printing (style is :text, :expr or :latex)

"""
    OpFamily

The family of an operator (bosons, fermions, two-level systems, qudits, spins or so(n) generators, see
[`opfamily`](@ref)). Operators of different families always commute, and the algebra within a family
`F` is defined by methods of [`exchange!`](@ref), [`contract!`](@ref) and [`commutes`](@ref) with the
first argument `::Val{F}` (see the comments at the beginning of `src/operator_rules.jl`), together
with [`ladder`](@ref), [`vacuum_eigenvalue_right`](@ref), [`op_adjoint`](@ref), [`opsymbol`](@ref),
[`lookahead_contraction`](@ref), [`lookahead_meet`](@ref) and [`lookahead_exchange`](@ref), which have
default methods.
"""
@enum OpFamily::UInt8 BosonFamily FermionFamily TLSFamily QuditFamily SpinFamily SOFamily

const ComplexInt = Complex{Int}

#                         BosonCreate_ FermionCreate_ TLSCreate_ QuditCreate_ SpinPlus_   TLSx_      TLSy_      TLSz_      SpinMinus_  SpinZ_      SOGen_    QuditDestroy_ TLSDestroy_ FermionDestroy_ BosonDestroy_
const BaseOpType_family = (BosonFamily, FermionFamily, TLSFamily, QuditFamily, SpinFamily, TLSFamily, TLSFamily, TLSFamily, SpinFamily, SpinFamily, SOFamily, QuditFamily,   TLSFamily,  FermionFamily,  BosonFamily)

"`opfamily(A)`: the [`OpFamily`](@ref) of the operator `A` (or of the operator type `A::BaseOpType`)."
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
    elseif fam == SpinFamily
        f(Val(SpinFamily), args...)
    elseif fam == SOFamily
        f(Val(SOFamily), args...)
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

"`lookahead_contraction(A)`: whether `A` can contract with an operator that is not next to it in normal order (e.g., σˣ₁ in σˣ₁ σˣ₂ σʸ₁), if the operators in between commute with the later one."
lookahead_contraction(A::BaseOperator) = _with_family(lookahead_contraction, opfamily(A), A)
"`lookahead_meet(A)`: whether (for [`lookahead_contraction`](@ref)) the operators in between can also commute with `A` instead, so that the two operators meet in between."
lookahead_meet(A::BaseOperator) = _with_family(lookahead_meet, opfamily(A), A)
"`lookahead_exchange(A)`: whether (for [`lookahead_contraction`](@ref)) the later operator can also be moved past operators in between that it does not commute with (which gives additional terms with their commutators), so that operators that can contract are always contracted."
lookahead_exchange(A::BaseOperator) = _with_family(lookahead_exchange, opfamily(A), A)
"`ladder(A)`: `+1` for a creation operator, `-1` for an annihilation operator, and `0` for other operators."
ladder(A::BaseOperator) = _with_family(ladder, opfamily(A), A)
"`vacuum_eigenvalue_right(A)`: `λ` if ``A|0⟩ = λ|0⟩`` (including `λ = 0`), and `nothing` otherwise."
vacuum_eigenvalue_right(A::BaseOperator) = _with_family(vacuum_eigenvalue_right, opfamily(A), A)
"`vacuum_eigenvalue_left(A)`: `λ` if ``⟨0|A = λ⟨0|`` (including `λ = 0`), and `nothing` otherwise."
vacuum_eigenvalue_left(A::BaseOperator) = _with_family(vacuum_eigenvalue_left, opfamily(A), A)

"`op_adjoint(A)`: the adjoint of the operator `A`."
op_adjoint(A::BaseOperator) = _with_family(op_adjoint, opfamily(A), A)
"`opsymbol(A, style)`: the symbol that is appended to the name of the operator `A` when printing it, for `style` `:text`, `:expr` (in [`julia_expression`](@ref)) or `:latex`."
opsymbol(A::BaseOperator, style::Symbol) = _with_family(opsymbol, opfamily(A), A, style)

"`annihilates_vacuum_left(A)`: whether ⟨0|A = 0"
# (eigenvalues are usually integers, but can be any number, e.g., for spins)
annihilates_vacuum_left(A::BaseOperator) = (λ = vacuum_eigenvalue_left(A); λ isa Int && iszero(λ))

# defaults
lookahead_contraction(::Val, ::BaseOperator) = false
lookahead_meet(::Val, ::BaseOperator) = false
lookahead_exchange(::Val, ::BaseOperator) = false
function ladder(::Val, A::BaseOperator)
    t = optype(A)
    t in (BosonCreate_, FermionCreate_, TLSCreate_, QuditCreate_, SpinPlus_) ? 1 : t in (BosonDestroy_, FermionDestroy_, TLSDestroy_, QuditDestroy_, SpinMinus_) ? -1 : 0
end
op_adjoint(::Val, A::BaseOperator) = BaseOperator(BaseOpType_adj[Int(optype(A))],opname(A),inds(A),opextra(A))
function opsymbol(::Val, A::BaseOperator, style::Symbol)
    t = Int(optype(A))
    style === :text ? BaseOpType_sym[t] : style === :expr ? BaseOpType_expr[t] : BaseOpType_latex[t]
end
# for ladder operators, the vacuum is annihilated by the annihilation operators (and their adjoints from the left)
vacuum_eigenvalue_right(F::Val, A::BaseOperator) = ladder(F, A) < 0 ? 0 : nothing
vacuum_eigenvalue_left(F::Val, A::BaseOperator) = ladder(F, A) > 0 ? 0 : nothing

const _NO_DELTAS = δ[]

"`maybe_same_site(A, B)`: whether the operators `A` and `B` can act on the same mode, i.e., have the same name and indices that can all be equal."
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
"`site_deltas(A, B)`: the Kronecker deltas that make the operators `A` and `B` act on the same mode (only for [`maybe_same_site`](@ref)`(A, B)`)."
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

lookahead_contraction(::Val{TLSFamily}, A::BaseOperator) = ispauli(optype(A))
# e.g., σˣ(i) σʸ(j) σᶻ(i) = σˣ(i) σᶻ(i) σʸ(j) + σˣ(i) [σʸ(j),σᶻ(i)] = -i σʸ(i) σʸ(j) - 2 δ(ij), so that the normal form
# of products of Pauli matrices with symbolic indices is unique (Pauli matrices with different indices that do
# not commute have to be on the same site, so their exchange gives x = 1)
lookahead_exchange(::Val{TLSFamily}, A::BaseOperator) = ispauli(optype(A))

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
lookahead_contraction(::Val{QuditFamily}, ::BaseOperator) = true

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

#################################################################
## spins (su(2) ≅ so(3), in the Cartan-Weyl basis)              ##
#################################################################
# Operators J₊, J₋ and J_z with [J_z,J±] = ±J± and [J₊,J₋] = 2J_z. Since these are the commutation
# relations of the Lie algebra, they hold in any representation (e.g., for collective spins). If
# the spin j of the representation is fixed, the Casimir J² = J₊J₋ + J_z² - J_z = j(j+1) is used to
# eliminate products J₊J₋ (on the same site), so that the normal-ordered products on one site are
# J₊ᵃ J_zⁿ and J₋ᵇ J_zⁿ. The vacuum is the lowest-weight state, with J₋|0⟩ = 0 and J_z|0⟩ = -j|0⟩.
# The normal order is J₊ < J₋ < J_z, so that J₊ and J₋ on the same site can always be brought
# next to each other (by moving them past operators of the same kind).
# The extra information stored with each operator is the index of j in _SPIN_J (0 if j is not fixed).
# (The table of values of j is replaced instead of modified, so that it can be read from several threads.)

const _SPIN_J = _AtomicRef(Any[])
const _SPIN_J_LOCK = ReentrantLock()
_spin_extra(j::Nothing) = 0
function _spin_extra(j)
    k = findfirst(isequal(j), _SPIN_J[])
    k === nothing || return k
    @lock _SPIN_J_LOCK begin
        js = _SPIN_J[]
        k = findfirst(isequal(j), js)
        k === nothing || return k
        _SPIN_J[] = push!(copy(js), j)
        length(js) + 1
    end
end
"`spin_j(A)`: the spin j of the spin operator A, or `nothing` if it is not fixed"
spin_j(A::BaseOperator) = (x = opextra(A); x == 0 ? nothing : _SPIN_J[][x])

_spin_like(A::BaseOperator, t::BaseOpType) = BaseOperator(t, opname(A), inds(A), opextra(A))

function exchange!(::Val{SpinFamily}, emit, B, A)
    ta, tb = optype(A), optype(B)
    ta == tb && return 1
    maybe_same_site(A,B) || return 1
    opextra(A) == opextra(B) || throw(ArgumentError("spin $(opname(A)) is used with different values of j"))
    dd = site_deltas(A,B)
    if ta == SpinPlus_ && tb == SpinMinus_
        # J₋ J₊ = J₊ J₋ - 2J_z
        emit(-2, dd, _spin_like(A,SpinZ_))
    elseif ta == SpinPlus_
        # J_z J₊ = J₊ J_z + J₊
        emit(1, dd, A)
    else
        # J_z J₋ = J₋ J_z - J₋
        emit(-1, dd, A)
    end
    return 1
end

# J₊ J₋ = j(j+1) - J_z² + J_z on the same site (if j is fixed)
function contract!(::Val{SpinFamily}, emit, A, B)
    (optype(A) == SpinPlus_ && optype(B) == SpinMinus_ && opname(A) == opname(B) && inds(A) == inds(B)) || return nothing
    opextra(A) == opextra(B) || throw(ArgumentError("spin $(opname(A)) is used with different values of j"))
    j = spin_j(A)
    j === nothing && return nothing
    Jz = _spin_like(A,SpinZ_)
    emit(j*(j+1), _NO_DELTAS)
    emit(-1, _NO_DELTAS, Jz, Jz)
    emit(1, _NO_DELTAS, Jz)
    (zero(ComplexInt), nothing)
end

commutes(::Val{SpinFamily}, A, B) = optype(A) == optype(B) || !maybe_same_site(A,B)

lookahead_contraction(::Val{SpinFamily}, A::BaseOperator) = optype(A) in (SpinPlus_, SpinMinus_) && opextra(A) != 0
lookahead_meet(::Val{SpinFamily}, ::BaseOperator) = true

# the eigenvalue -j of J_z for the vacuum (since j can be any number, it is wrapped so that the
# return type of vacuum_eigenvalue_* stays concrete for all families)
struct _SpinVacuumEigenvalue
    j::Any
end
_vacuum_value(λ::_SpinVacuumEigenvalue) = -λ.j
function _spin_vacuum_eigenvalue(A::BaseOperator)
    j = spin_j(A)
    j === nothing && throw(ArgumentError("the vacuum (lowest-weight state) is only defined for spins with fixed j, got $A"))
    _SpinVacuumEigenvalue(j)
end
vacuum_eigenvalue_right(::Val{SpinFamily}, A::BaseOperator) = (t = optype(A); t == SpinMinus_ ? 0 : t == SpinZ_ ? _spin_vacuum_eigenvalue(A) : nothing)
vacuum_eigenvalue_left(::Val{SpinFamily}, A::BaseOperator) = (t = optype(A); t == SpinPlus_ ? 0 : t == SpinZ_ ? _spin_vacuum_eigenvalue(A) : nothing)

#################################################################
## so(n): generators of rotations in n dimensions               ##
#################################################################
# Hermitian generators L_ab = -L_ba (a,b = 1…n) with
#   [L_ab, L_cd] = i (δ_ac L_bd - δ_ad L_bc - δ_bc L_ad + δ_bd L_ac).
# Only the commutation relations are used (the algebra is not tied to a representation), so the
# normal-ordered products on one site are the products of generators sorted by (a,b), which form a
# basis (Poincaré-Birkhoff-Witt theorem). Each generator is stored as L_ab with a < b, with the
# extra information n:8 | a:8 | b:8.

_so_extra(n::Integer, a::Integer, b::Integer) = (n << 16) | (a << 8) | b
"`so_generator(A)`: (n,a,b) for the generator A = L_ab (with a < b) of so(n)"
function so_generator(A::BaseOperator)
    x = opextra(A)
    (Int(x >> 16), Int((x >> 8) & 0xff), Int(x & 0xff))
end

"`so_op(name, n, a, b, inds)`: the generator L_ab of so(n), with a < b"
function so_op(name::QuOpName, n::Integer, a::Integer, b::Integer, inds::AbstractVector{QuIndex})
    2 ≤ n ≤ 255 || throw(ArgumentError("so(n) is only supported for 2 ≤ n ≤ 255, got n = $n"))
    1 ≤ a < b ≤ n || throw(ArgumentError("so($n) generators are stored as L_ab with 1 ≤ a < b ≤ $n, got ($a,$b)"))
    BaseOperator(SOGen_, name, inds, _so_extra(n,a,b))
end

# emit c L_xy (with L_xy = -L_yx) on the same site as A
function _emit_so(emit, c::ComplexInt, dd, A::BaseOperator, n, x, y)
    x < y ? emit(c, dd, so_op(opname(A), n, x, y, inds(A))) : emit(-c, dd, so_op(opname(A), n, y, x, inds(A)))
end

function exchange!(::Val{SOFamily}, emit, B, A)
    opextra(A) == opextra(B) && return 1
    maybe_same_site(A,B) || return 1
    n, c, d = so_generator(A)
    nB, a, b = so_generator(B)
    n == nB || throw(ArgumentError("so(n) generators $(opname(A)) are used with different values of n"))
    # B A = A B + [L_ab, L_cd], where a<b and c<d, and (a,b) ≠ (c,d) share at most one index
    dd = site_deltas(A,B)
    a == c && _emit_so(emit, ComplexInt(0, 1), dd, A, n, b, d)
    a == d && _emit_so(emit, ComplexInt(0,-1), dd, A, n, b, c)
    b == c && _emit_so(emit, ComplexInt(0,-1), dd, A, n, a, d)
    b == d && _emit_so(emit, ComplexInt(0, 1), dd, A, n, a, c)
    return 1
end

# products of generators do not simplify (without fixing a representation)
contract!(::Val{SOFamily}, emit, A, B) = nothing

function commutes(::Val{SOFamily}, A, B)
    (opextra(A) == opextra(B) || !maybe_same_site(A,B)) && return true
    n, c, d = so_generator(A)
    nB, a, b = so_generator(B)
    # (with different n, exchange! throws an error)
    n == nB && a != c && a != d && b != c && b != d
end

opsymbol(::Val{SOFamily}, A::BaseOperator, style::Symbol) = ((_, a, b) = so_generator(A); _qudit_symbol(a, b, style))

function _so_no_vacuum(A::BaseOperator)
    throw(ArgumentError("vacuum expectation values are not defined for so(n) generators, got $A"))
end
vacuum_eigenvalue_right(::Val{SOFamily}, A::BaseOperator) = _so_no_vacuum(A)
vacuum_eigenvalue_left(::Val{SOFamily}, A::BaseOperator) = _so_no_vacuum(A)
