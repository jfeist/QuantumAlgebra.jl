# Prototype of a compact representation for expressions on a fixed set of concrete modes (operators
# with integer indices, no sums, δs, expectation values or correlators), for fast products and
# commutators. All kinds of operators are supported.
#
# Each term is a monomial stored as a fixed number of 64-bit words (the "key"), with fields for:
# - bosons: 16 bits per mode, the exponents m (low 8 bits) and n (high 8 bits) of a†ᵐ aⁿ. Products of
#   monomials follow Wick's theorem, aⁿ¹ a†ᵐ² = ∑ₖ k! C(n₁,k) C(m₂,k) a†ᵐ²⁻ᵏ aⁿ¹⁻ᵏ, mode by mode.
# - two-level systems with Pauli matrices: 2 bits per site, the Pauli matrix 1, σᶻ, σˣ, σʸ (labels 0, 1,
#   2, 3). The product of two Pauli matrices is a single Pauli matrix (XOR of the labels) times a phase.
# - all other sites (fermions, two-level systems with σ±, qudits, spins, so(n) generators): 16 bits per
#   site, a label for a product of operators on that site in QuantumAlgebra's normal form (0 is the
#   identity). Products of two such local monomials are computed once with QuantumAlgebra and cached,
#   so that the results agree with QuantumAlgebra by construction. Fermions with the same name
#   anticommute, which gives a sign when operators of different sites are reordered.
# - parameters: 8 bits per parameter (with name, indices and state 'r', 'n' or 'c'), the exponent.
# The modes and parameters are collected from the expressions passed to `concrete`, which share one
# `ConcreteBasis`. Operators on modes that do not appear in a monomial are the identity.
module ConcreteProto

using QuantumAlgebra
const QA = QuantumAlgebra
using QuantumAlgebra: QuExpr, QuTerm, QuIndex, QuOpName, BaseOperator, BaseOpProduct, Param, ExpVal, Corr

export ConcreteBasis, ConcreteExpr, concrete

const IndsTuple = Tuple{Vararg{QuIndex}}

# a product of operators on one site (in normal order), given by the types and extra information of the operators
const LocalMono = Vector{Tuple{QA.BaseOpType,UInt32}}

"The monomials on single sites of one family of operators (shared by all its sites), labeled 1, 2, ..."
struct LocalAlgebra
    family::QA.OpFamily
    monos::Vector{LocalMono}
    labels::Dict{LocalMono,Int}
    odd::Vector{Bool}   # fermions: odd number of operators
end
LocalAlgebra(fam) = LocalAlgebra(fam, LocalMono[], Dict{LocalMono,Int}(), Bool[])
function _label!(alg::LocalAlgebra, m::LocalMono)
    get!(alg.labels, m) do
        push!(alg.monos, m)
        push!(alg.odd, alg.family == QA.FermionFamily && isodd(length(m)))
        length(alg.monos)
    end
end

"The modes and parameters of a set of expressions, and where their fields are stored in the keys."
struct ConcreteBasis
    bosons::Vector{Tuple{QuOpName,IndsTuple}}
    paulis::Vector{Tuple{QuOpName,IndsTuple}}
    locals::Vector{Tuple{QuOpName,IndsTuple}}
    localalg::Vector{Int}                 # algebra of each local site
    algebras::Vector{LocalAlgebra}
    fermions::Vector{Tuple{Int,QuOpName}} # local sites with fermions, and their names
    params::Vector{Tuple{QuOpName,IndsTuple,Char}}
    lookup::Dict{Any,Int}
    nwb::Int   # words for bosons (4 modes per word)
    nwt::Int   # words for two-level systems with Pauli matrices (32 sites per word)
    nwl::Int   # words for other sites (4 sites per word)
    nwp::Int   # words for parameters (8 per word)
    caches::Dict{DataType,Any}   # products of local monomials, for each coefficient type
end
nwords(b::ConcreteBasis) = max(1, b.nwb + b.nwt + b.nwl + b.nwp)

_concreteinds(X) = (ii = Tuple(QA.inds(X)); all(QA.isintindex, ii) || throw(ArgumentError("only integer indices are supported, got $X")); ii)
_ispauli(O::BaseOperator) = QA.optype(O) in (QA.TLSx_, QA.TLSy_, QA.TLSz_)
_sitekind(O::BaseOperator) = (fam = QA.opfamily(O); fam == QA.BosonFamily ? :boson : fam == QA.TLSFamily && _ispauli(O) ? :pauli : :local)

function ConcreteBasis(exprs)
    bosons, paulis, locals = (Tuple{QuOpName,IndsTuple}[] for _ in 1:3)
    localalg = Int[]
    algebras = LocalAlgebra[]
    fermions = Tuple{Int,QuOpName}[]
    params = Tuple{QuOpName,IndsTuple,Char}[]
    lookup = Dict{Any,Int}()
    for A in exprs, t in keys(A.terms)
        (QA.nsuminds(t) == 0 && isempty(QA.deltas(t)) && isempty(QA.expvals(t)) && isempty(QA.corrs(t))) ||
            throw(ArgumentError("sums, δs, expectation values and correlators are not supported"))
        for p in QA.params(t)
            k = (:param, QA.opname(p), _concreteinds(p), QA.parstate(p))
            haskey(lookup, k) || (push!(params, k[2:4]); lookup[k] = length(params))
        end
        for O in QA.bares(t)
            fam, kind = QA.opfamily(O), _sitekind(O)
            k = (fam, QA.opname(O), _concreteinds(O))
            if haskey(lookup, k)
                # (a two-level system cannot have both Pauli matrices and σ±, as in QuantumAlgebra)
                haskey(lookup, (kind, k)) || throw(ArgumentError("two-level system $(QA.opname(O)) with both Pauli matrices and σ±"))
                continue
            end
            v = kind == :boson ? bosons : kind == :pauli ? paulis : locals
            push!(v, k[2:3])
            lookup[k] = length(v)
            lookup[(kind, k)] = 1
            if kind == :local
                ia = findfirst(alg -> alg.family == fam, algebras)
                ia === nothing && (push!(algebras, LocalAlgebra(fam)); ia = length(algebras))
                push!(localalg, ia)
                fam == QA.FermionFamily && push!(fermions, (length(locals), QA.opname(O)))
            end
        end
    end
    ConcreteBasis(bosons, paulis, locals, localalg, algebras, fermions, params, lookup,
                  cld(length(bosons), 4), cld(length(paulis), 32), cld(length(locals), 4), cld(length(params), 8), Dict{DataType,Any}())
end

# position (word, shift) of the field of boson mode i (m at shift, n at shift+8), Pauli site i, local site i, and parameter i
@inline _bpos(b::ConcreteBasis, i) = (((i-1) >> 2) + 1, 16*((i-1) & 3))
@inline _tpos(b::ConcreteBasis, i) = (b.nwb + ((i-1) >> 5) + 1, 2*((i-1) & 31))
@inline _lpos(b::ConcreteBasis, i) = (b.nwb + b.nwt + ((i-1) >> 2) + 1, 16*((i-1) & 3))
@inline _ppos(b::ConcreteBasis, i) = (b.nwb + b.nwt + b.nwl + ((i-1) >> 3) + 1, 8*((i-1) & 7))
@inline _localfield(b::ConcreteBasis, key, i) = ((w, sh) = _lpos(b, i); Int((key[w] >> sh) & 0xffff))

"`ConcreteExpr{W,T}`: a sum of monomials (keys of `W` words) with coefficients of type `T`."
struct ConcreteExpr{W,T}
    basis::ConcreteBasis
    terms::Dict{NTuple{W,UInt64},T}
end
Base.length(A::ConcreteExpr) = length(A.terms)

@inline _zerokey(::Val{W}) where W = ntuple(_ -> UInt64(0), Val(W))
@inline _setfield(key, (w, sh), x) = Base.setindex(key, key[w] + (UInt64(x) << sh), w)

@inline function _add!(d::Dict{K,T}, k::K, s) where {K,T}
    iszero(s) && return nothing
    snew = get(d, k, zero(T)) + s
    iszero(snew) ? delete!(d, k) : (d[k] = snew)
    nothing
end

#################################################################
## products of local monomials (computed with QuantumAlgebra) ##
#################################################################

const _CacheKey = Tuple{Int,Int,Int}
_cache(b::ConcreteBasis, ::Type{T}) where T = get!(() -> Dict{_CacheKey,Vector{Tuple{Int,T}}}(), b.caches, T)::Dict{_CacheKey,Vector{Tuple{Int,T}}}

const _REPNAME = Ref{QuOpName}()
_repexpr(alg::LocalAlgebra, l) = QuExpr(QuTerm(BaseOpProduct([BaseOperator(t, _REPNAME[], QuIndex[], x) for (t, x) in alg.monos[l]])))

# the product of local monomials la and lb (on the same site) as a list of (label, coefficient)
function _localproduct(alg::LocalAlgebra, la, lb, ::Type{T}) where T
    isassigned(_REPNAME) || (_REPNAME[] = QuOpName(:ξ_concrete))
    R = normal_form(_repexpr(alg, la) * _repexpr(alg, lb))
    res = Tuple{Int,T}[]
    for (t, s) in R.terms
        @assert QA.nsuminds(t) == 0 && isempty(QA.deltas(t)) && isempty(QA.params(t))
        m = LocalMono([(QA.optype(O), QA.opextra(O)) for O in QA.bares(t)])
        push!(res, (isempty(m) ? 0 : _label!(alg, m), convert(T, s)))
    end
    res
end
@inline _localprod(cache, b::ConcreteBasis, ia, la, lb, ::Type{T}) where T =
    get!(() -> _localproduct(b.algebras[ia], la, lb, T), cache, (ia, la, lb))

#################################################################
## products of monomials                                       ##
#################################################################

const _LOW2 = 0x5555555555555555
const _LOW16 = 0x0001000100010001
# low bit of each 2-bit (16-bit) field set if the field is nonzero
@inline _nz2(t::UInt64) = (t | (t >> 1)) & _LOW2
@inline function _nz16(x::UInt64)
    x |= x >> 8; x |= x >> 4; x |= x >> 2; x |= x >> 1
    x & _LOW16
end
# position of the Pauli matrices (labels 1, 2, 3 = σᶻ, σˣ, σʸ) in the cyclic order σˣ, σʸ, σᶻ
@inline _cyc(l) = l == 2 ? 0 : l == 3 ? 1 : 2
const _IPOW = (Complex(1,0), Complex(0,1), Complex(-1,0), Complex(0,-1))

# whether the operators of the local site i in monomial k are fermionic with an odd number of operators
@inline _fodd(b::ConcreteBasis, key, i) = (l = _localfield(b, key, i); l != 0 && b.algebras[b.localalg[i]].odd[l])

# the sign from bringing ka kb = (∏ᵢ aᵢ)(∏ᵢ bᵢ) into the form ∏ᵢ (aᵢ bᵢ), i.e., moving each bᵢ past the
# aₖ with k > i (only fermions with the same name anticommute)
function _fermionsign(b::ConcreteBasis, ka, kb)
    neg = false
    fs = b.fermions
    for x in eachindex(fs)
        i, name = fs[x]
        _fodd(b, kb, i) || continue
        for y in x+1:length(fs)
            k, namek = fs[y]
            namek == name && _fodd(b, ka, k) && (neg = !neg)
        end
    end
    neg
end

# call f(key, coefficient) for each monomial of s * ka * kb
@inline function _mulmono!(f::F, b::ConcreteBasis, cache, ka::NTuple{W,UInt64}, kb::NTuple{W,UInt64}, s) where {F,W}
    t0, t1, l1 = b.nwb, b.nwb + b.nwt, b.nwb + b.nwt + b.nwl
    # bosons (without contractions) and parameters: add the exponents, Pauli matrices: XOR of the labels,
    # local sites: the label of the monomial that is not the identity (products are done in _local!)
    key = ntuple(w -> (t0 < w <= t1) ? (ka[w] ⊻ kb[w]) : (t1 < w <= l1) ? (ka[w] | kb[w]) : (ka[w] + kb[w]), Val(W))
    # phase from the Pauli matrices: σᵃσᵇ = iϵ_abc σᶜ for a ≠ b
    phase = 0
    for w in t0+1:t1
        ta, tb = ka[w], kb[w]
        both = _nz2(ta) & _nz2(tb) & _nz2(ta ⊻ tb)
        while both != 0
            sh = trailing_zeros(both)
            d = _cyc((tb >> sh) & 3) - _cyc((ta >> sh) & 3)
            phase += (d == 1 || d == -2) ? 1 : 3
            both &= both - 1
        end
    end
    isempty(b.fermions) || (_fermionsign(b, ka, kb) && (phase += 2))
    phase & 3 != 0 && (s *= _IPOW[(phase & 3) + 1])
    _local!(f, b, cache, ka, kb, key, s, t1 + 1, 0)
end

# go through the local sites (starting at word w, slot `slot`) and multiply the monomials of ka and kb on
# the sites where both are not the identity, then continue with the bosons
function _local!(f::F, b::ConcreteBasis, cache::Dict{_CacheKey,Vector{Tuple{Int,T}}}, ka, kb, key, s, w, slot) where {F,T}
    l1 = b.nwb + b.nwt + b.nwl
    while w <= l1
        wa, wb = ka[w], kb[w]
        if wa != 0 && wb != 0
            while slot < 4
                sh = 16slot
                la, lb = Int((wa >> sh) & 0xffff), Int((wb >> sh) & 0xffff)
                if la != 0 && lb != 0
                    i = 4*(w - b.nwb - b.nwt - 1) + slot + 1
                    cleared = key[w] & ~(UInt64(0xffff) << sh)
                    for (l, c) in _localprod(cache, b, b.localalg[i], la, lb, T)
                        _local!(f, b, cache, ka, kb, Base.setindex(key, cleared | (UInt64(l) << sh), w), s*c, w, slot+1)
                    end
                    return nothing
                end
                slot += 1
            end
        end
        w += 1
        slot = 0
    end
    _wick!(f, b, ka, kb, key, s, 1, 0)
end

# go through the boson modes (starting at word w, slot `slot`) and add the contractions of the
# annihilation operators of ka with the creation operators of kb
function _wick!(f::F, b::ConcreteBasis, ka, kb, key, s, w, slot) where F
    while w <= b.nwb
        wa, wb = ka[w], kb[w]
        if wa != 0 && wb != 0
            while slot < 4
                sh = 16slot
                n1 = Int((wa >> (sh+8)) & 0xff)
                m2 = Int((wb >> sh) & 0xff)
                if n1 > 0 && m2 > 0
                    _wick!(f, b, ka, kb, key, s, w, slot+1)
                    c = 1
                    for k in 1:min(n1,m2)
                        # k! C(n1,k) C(m2,k), computed incrementally
                        c = c * (n1-k+1) * (m2-k+1) ÷ k
                        d = (UInt64(k) << sh) | (UInt64(k) << (sh+8))
                        _wick!(f, b, ka, kb, Base.setindex(key, key[w] - d, w), s*c, w, slot+1)
                    end
                    return nothing
                end
                slot += 1
            end
        end
        w += 1
        slot = 0
    end
    f(key, s)
    nothing
end

# the fermionic parity of the operators with the given name in monomial k
function _fparity(b::ConcreteBasis, key, name)
    p = false
    for (i, n) in b.fermions
        n == name && _fodd(b, key, i) && (p = !p)
    end
    p
end

# whether two monomials do not necessarily commute: they act on a common site, or both have an odd
# number of fermion operators with the same name
function _overlap(b::ConcreteBasis, ka, kb)
    for w in 1:b.nwb
        _nz16(ka[w]) & _nz16(kb[w]) != 0 && return true
    end
    for w in b.nwb+1:b.nwb+b.nwt
        _nz2(ka[w]) & _nz2(kb[w]) != 0 && return true
    end
    for w in b.nwb+b.nwt+1:b.nwb+b.nwt+b.nwl
        _nz16(ka[w]) & _nz16(kb[w]) != 0 && return true
    end
    for (_, name) in b.fermions
        _fparity(b, ka, name) && _fparity(b, kb, name) && return true
    end
    false
end

#################################################################
## arithmetic                                                   ##
#################################################################

_check(A::ConcreteExpr, B::ConcreteExpr) = A.basis === B.basis || throw(ArgumentError("expressions have different bases, create them together with `concrete`"))

function Base.:*(A::ConcreteExpr{W,T}, B::ConcreteExpr{W,T}) where {W,T}
    _check(A, B)
    b = A.basis
    cache = _cache(b, T)
    C = Dict{NTuple{W,UInt64},T}()
    add = (k, x) -> _add!(C, k, x)
    for (kb, sb) in B.terms, (ka, sa) in A.terms
        _mulmono!(add, b, cache, ka, kb, sa*sb)
    end
    ConcreteExpr{W,T}(b, C)
end
Base.:*(x::Number, A::ConcreteExpr{W,T}) where {W,T} = ConcreteExpr{W,T}(A.basis, Dict(k => convert(T, x*s) for (k, s) in A.terms if !iszero(x*s)))
function Base.:+(A::ConcreteExpr{W,T}, B::ConcreteExpr{W,T}) where {W,T}
    _check(A, B)
    C = copy(A.terms)
    for (k, s) in B.terms
        _add!(C, k, s)
    end
    ConcreteExpr{W,T}(A.basis, C)
end
Base.:-(A::ConcreteExpr, B::ConcreteExpr) = A + (-1)*B

"`comm(A, B)`: the commutator AB - BA (only pairs of monomials that do not commute trivially contribute)."
function QA.comm(A::ConcreteExpr{W,T}, B::ConcreteExpr{W,T}) where {W,T}
    _check(A, B)
    b = A.basis
    cache = _cache(b, T)
    C = Dict{NTuple{W,UInt64},T}()
    add = (k, x) -> _add!(C, k, x)
    sub = (k, x) -> _add!(C, k, -x)
    for (ka, sa) in A.terms, (kb, sb) in B.terms
        _overlap(b, ka, kb) || continue
        s = sa*sb
        _mulmono!(add, b, cache, ka, kb, s)
        _mulmono!(sub, b, cache, kb, ka, s)
    end
    ConcreteExpr{W,T}(b, C)
end

#################################################################
## conversion from and to QuExpr                                ##
#################################################################

"""
    concrete(exprs...; T=nothing)

Convert the expressions `exprs` (which can only contain operators and parameters with integer indices)
to `ConcreteExpr`s with a common basis. The coefficient type `T` is by default the promotion of
`Complex{Int}`, the types of the coefficients, and the type of j(j+1) for spins with fixed j.
"""
function concrete(exprs::QuExpr...; T=nothing)
    b = ConcreteBasis(exprs)
    W = nwords(b)
    if T === nothing
        T = Complex{Int}
        for A in exprs, (t, s) in A.terms
            T = promote_type(T, typeof(s))
            for O in QA.bares(t)
                QA.opfamily(O) == QA.SpinFamily && (j = QA.spin_j(O); j === nothing || (T = promote_type(T, typeof(j*(j+1)))))
            end
        end
    end
    map(A -> _convert(ConcreteExpr{W,T}, b, A), exprs)
end

# the monomials (key => coefficient) of a single operator
function _opmonos(::Type{ConcreteExpr{W,T}}, b::ConcreteBasis, O::BaseOperator) where {W,T}
    t = QA.optype(O)
    k0 = _zerokey(Val(W))
    i = b.lookup[(QA.opfamily(O), QA.opname(O), Tuple(QA.inds(O)))]
    kind = _sitekind(O)
    if kind == :boson
        w, sh = _bpos(b, i)
        return [(_setfield(k0, (w, t == QA.BosonCreate_ ? sh : sh+8), 1), one(T))]
    elseif kind == :pauli
        l = t == QA.TLSz_ ? 1 : t == QA.TLSx_ ? 2 : 3
        return [(_setfield(k0, _tpos(b, i), l), one(T))]
    else
        l = _label!(b.algebras[b.localalg[i]], LocalMono([(t, QA.opextra(O))]))
        return [(_setfield(k0, _lpos(b, i), l), one(T))]
    end
end

function _convert(::Type{ConcreteExpr{W,T}}, b::ConcreteBasis, A::QuExpr) where {W,T}
    cache = _cache(b, T)
    C = Dict{NTuple{W,UInt64},T}()
    for (t, s) in A.terms
        key = _zerokey(Val(W))
        for p in QA.params(t)
            key = _setfield(key, _ppos(b, b.lookup[(:param, QA.opname(p), Tuple(QA.inds(p)), QA.parstate(p))]), 1)
        end
        cur = Dict(key => convert(T, s))
        for O in QA.bares(t)
            next = Dict{NTuple{W,UInt64},T}()
            add = (k, x) -> _add!(next, k, x)
            for (k1, s1) in cur, (k2, s2) in _opmonos(ConcreteExpr{W,T}, b, O)
                _mulmono!(add, b, cache, k1, k2, s1*s2)
            end
            cur = next
        end
        for (k, x) in cur
            _add!(C, k, x)
        end
    end
    ConcreteExpr{W,T}(b, C)
end

const _PAULI_TYPES = (QA.TLSz_, QA.TLSx_, QA.TLSy_)

# the sign of sorting ops (as a product of operators) into normal order: -1 for each pair of fermions
# with the same name whose order is exchanged
function _sortsign(ops::Vector{BaseOperator})
    neg = false
    for x in eachindex(ops), y in x+1:length(ops)
        A, B = ops[x], ops[y]
        QA.opfamily(A) == QA.FermionFamily && QA.opfamily(B) == QA.FermionFamily && QA.opname(A) == QA.opname(B) &&
            isless(B, A) && (neg = !neg)
    end
    neg
end

# σˣ = σ⁺ + σ⁻, σʸ = i(σ⁻ - σ⁺), σᶻ = 2σ⁺σ⁻ - 1 (for use_σpm(true))
function _pauli_pm(name, ii, l)
    p, m = QuExpr(BaseOperator(QA.TLSCreate_, name, QuIndex[ii...])), QuExpr(BaseOperator(QA.TLSDestroy_, name, QuIndex[ii...]))
    l == 2 ? p + m : l == 3 ? 1im*(m - p) : 2*p*m - 1
end

"Convert back to a `QuExpr` (in normal form)."
function QA.QuExpr(C::ConcreteExpr{W,T}) where {W,T}
    b = C.basis
    R = QuExpr()
    pm = QA.using_σpm()
    for (key, s) in C.terms
        ps = Param[]
        for (i, (name, ii, st)) in enumerate(b.params)
            w, sh = _ppos(b, i)
            for _ in 1:Int((key[w] >> sh) & 0xff)
                push!(ps, Param(name, st, QuIndex[ii...]))
            end
        end
        ops = BaseOperator[]
        for (i, (name, ii)) in enumerate(b.bosons)
            w, sh = _bpos(b, i)
            append!(ops, fill(BaseOperator(QA.BosonCreate_, name, QuIndex[ii...]), Int((key[w] >> sh) & 0xff)))
            append!(ops, fill(BaseOperator(QA.BosonDestroy_, name, QuIndex[ii...]), Int((key[w] >> (sh+8)) & 0xff)))
        end
        for (i, (name, ii)) in enumerate(b.locals)
            l = _localfield(b, key, i)
            l == 0 && continue
            for (t, x) in b.algebras[b.localalg[i]].monos[l]
                push!(ops, BaseOperator(t, name, QuIndex[ii...], x))
            end
        end
        _sortsign(ops) && (s = -s)
        paulis = Tuple{QuOpName,IndsTuple,Int}[]
        for (i, (name, ii)) in enumerate(b.paulis)
            w, sh = _tpos(b, i)
            l = Int((key[w] >> sh) & 3)
            l == 0 && continue
            pm ? push!(paulis, (name, ii, l)) : push!(ops, BaseOperator(_PAULI_TYPES[l], name, QuIndex[ii...]))
        end
        t = QuTerm(0, QA.δ[], sort!(ps), ExpVal[], Corr[], BaseOpProduct(sort!(ops)))
        if isempty(paulis)
            QA._add_sum_term!(R, t, s)
        else
            X = s * QuExpr(t)
            for (name, ii, l) in paulis
                X = X * _pauli_pm(name, ii, l)
            end
            R = R + normal_form(X)
        end
    end
    R
end

end
