# Prototype of a compact representation for expressions on a fixed set of concrete modes (operators
# with integer indices, no sums, δs, expectation values or correlators), for fast products and
# commutators. Supported are bosons, two-level systems and parameters.
#
# Each term is a monomial stored as a fixed number of 64-bit words (the "key"):
# - bosons: 16 bits per mode, the exponents m (low 8 bits) and n (high 8 bits) of a†ᵐ aⁿ. Products of
#   monomials follow Wick's theorem, aⁿ¹ a†ᵐ² = ∑ₖ k! C(n₁,k) C(m₂,k) a†ᵐ²⁻ᵏ aⁿ¹⁻ᵏ, mode by mode.
# - two-level systems: 2 bits per site, the Pauli matrix 1, σᶻ, σˣ, σʸ (labels 0, 1, 2, 3). The product
#   of two Pauli matrices is a single Pauli matrix (the XOR of the labels) times a phase, so the
#   product of two monomials never has more terms than the bosonic contractions give.
# - parameters: 8 bits per parameter (with name, indices and state 'r', 'n' or 'c'), the exponent.
# The modes and parameters are collected from the expressions passed to `concrete`, which share one
# `ConcreteBasis`. Operators on modes that do not appear in a monomial are the identity.
module ConcreteProto

using QuantumAlgebra
const QA = QuantumAlgebra
using QuantumAlgebra: QuExpr, QuTerm, QuIndex, QuOpName, BaseOperator, BaseOpProduct, Param, ExpVal, Corr

export ConcreteBasis, ConcreteExpr, concrete

const IndsTuple = Tuple{Vararg{QuIndex}}

"The modes and parameters of a set of expressions, and where their fields are stored in the keys."
struct ConcreteBasis
    bosons::Vector{Tuple{QuOpName,IndsTuple}}
    tlss::Vector{Tuple{QuOpName,IndsTuple}}
    params::Vector{Tuple{QuOpName,IndsTuple,Char}}
    lookup::Dict{Any,Int}
    nwb::Int   # words for bosons (4 modes per word)
    nwt::Int   # words for two-level systems (32 sites per word)
    nwp::Int   # words for parameters (8 per word)
end
nwords(b::ConcreteBasis) = max(1, b.nwb + b.nwt + b.nwp)

_concreteinds(X) = (ii = Tuple(QA.inds(X)); all(QA.isintindex, ii) || throw(ArgumentError("only integer indices are supported, got $X")); ii)

function ConcreteBasis(exprs)
    bosons = Tuple{QuOpName,IndsTuple}[]
    tlss = Tuple{QuOpName,IndsTuple}[]
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
            fam = QA.opfamily(O)
            v = fam == QA.BosonFamily ? bosons : fam == QA.TLSFamily ? tlss : throw(ArgumentError("unsupported operator $O"))
            k = (fam, QA.opname(O), _concreteinds(O))
            haskey(lookup, k) || (push!(v, k[2:3]); lookup[k] = length(v))
        end
    end
    ConcreteBasis(bosons, tlss, params, lookup, cld(length(bosons), 4), cld(length(tlss), 32), cld(length(params), 8))
end

# position (word, shift) of the field of boson mode i (m at shift, n at shift+8), two-level system i, and parameter i
@inline _bpos(b::ConcreteBasis, i) = (((i-1) >> 2) + 1, 16*((i-1) & 3))
@inline _tpos(b::ConcreteBasis, i) = (b.nwb + ((i-1) >> 5) + 1, 2*((i-1) & 31))
@inline _ppos(b::ConcreteBasis, i) = (b.nwb + b.nwt + ((i-1) >> 3) + 1, 8*((i-1) & 7))

"`ConcreteExpr{W,T}`: a sum of monomials (keys of `W` words) with coefficients of type `T`."
struct ConcreteExpr{W,T}
    basis::ConcreteBasis
    terms::Dict{NTuple{W,UInt64},T}
end
ConcreteExpr{W,T}(b::ConcreteBasis) where {W,T} = ConcreteExpr{W,T}(b, Dict{NTuple{W,UInt64},T}())
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

# call f(key, coefficient) for each monomial of s * ka * kb
@inline function _mulmono!(f::F, b::ConcreteBasis, ka::NTuple{W,UInt64}, kb::NTuple{W,UInt64}, s) where {F,W}
    t0, t1 = b.nwb, b.nwb + b.nwt
    # bosons (without contractions) and parameters: add the exponents, Pauli matrices: XOR of the labels
    key = ntuple(w -> (t0 < w <= t1) ? (ka[w] ⊻ kb[w]) : (ka[w] + kb[w]), Val(W))
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
    phase & 3 != 0 && (s *= _IPOW[(phase & 3) + 1])
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

# whether two monomials act on a common mode (otherwise they commute)
@inline function _overlap(b::ConcreteBasis, ka, kb)
    for w in 1:b.nwb
        _nz16(ka[w]) & _nz16(kb[w]) != 0 && return true
    end
    for w in b.nwb+1:b.nwb+b.nwt
        _nz2(ka[w]) & _nz2(kb[w]) != 0 && return true
    end
    false
end

#################################################################
## arithmetic                                                   ##
#################################################################

_check(A::ConcreteExpr, B::ConcreteExpr) = A.basis === B.basis || throw(ArgumentError("expressions have different bases, create them together with `concrete`"))

function Base.:*(A::ConcreteExpr{W,T}, B::ConcreteExpr{W,T}) where {W,T}
    _check(A, B)
    C = Dict{NTuple{W,UInt64},T}()
    add = (k, x) -> _add!(C, k, x)
    for (kb, sb) in B.terms, (ka, sa) in A.terms
        _mulmono!(add, A.basis, ka, kb, sa*sb)
    end
    ConcreteExpr{W,T}(A.basis, C)
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

"`comm(A, B)`: the commutator AB - BA (only pairs of monomials that act on a common mode contribute)."
function QA.comm(A::ConcreteExpr{W,T}, B::ConcreteExpr{W,T}) where {W,T}
    _check(A, B)
    b = A.basis
    C = Dict{NTuple{W,UInt64},T}()
    add = (k, x) -> _add!(C, k, x)
    sub = (k, x) -> _add!(C, k, -x)
    for (ka, sa) in A.terms, (kb, sb) in B.terms
        _overlap(b, ka, kb) || continue
        s = sa*sb
        _mulmono!(add, b, ka, kb, s)
        _mulmono!(sub, b, kb, ka, s)
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
`Complex{Int}` and the types of the coefficients (and `Rational{Int}` if there are σ⁺ or σ⁻).
"""
function concrete(exprs::QuExpr...; T=nothing)
    b = ConcreteBasis(exprs)
    W = nwords(b)
    if T === nothing
        T = Complex{Int}
        for A in exprs, (t, s) in A.terms
            T = promote_type(T, typeof(s))
            any(O -> QA.optype(O) in (QA.TLSCreate_, QA.TLSDestroy_), QA.bares(t)) && (T = promote_type(T, Rational{Int}))
        end
    end
    map(A -> _convert(ConcreteExpr{W,T}, b, A), exprs)
end

# the monomials (key => coefficient) of a single operator
function _opmonos(::Type{ConcreteExpr{W,T}}, b::ConcreteBasis, O::BaseOperator) where {W,T}
    t = QA.optype(O)
    k0 = _zerokey(Val(W))
    if QA.opfamily(O) == QA.BosonFamily
        w, sh = _bpos(b, b.lookup[(QA.BosonFamily, QA.opname(O), Tuple(QA.inds(O)))])
        return [(_setfield(k0, (w, t == QA.BosonCreate_ ? sh : sh+8), 1), one(T))]
    else
        pos = _tpos(b, b.lookup[(QA.TLSFamily, QA.opname(O), Tuple(QA.inds(O)))])
        X, Y, Z = _setfield(k0, pos, 2), _setfield(k0, pos, 3), _setfield(k0, pos, 1)
        t == QA.TLSx_ && return [(X, one(T))]
        t == QA.TLSy_ && return [(Y, one(T))]
        t == QA.TLSz_ && return [(Z, one(T))]
        # σ± = (σˣ ± iσʸ)/2
        sgn = t == QA.TLSCreate_ ? 1 : -1
        return [(X, convert(T, 1//2)), (Y, convert(T, sgn*im//2))]
    end
end

function _convert(::Type{ConcreteExpr{W,T}}, b::ConcreteBasis, A::QuExpr) where {W,T}
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
                _mulmono!(add, b, k1, k2, s1*s2)
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

# σˣ = σ⁺ + σ⁻, σʸ = i(σ⁻ - σ⁺), σᶻ = 2σ⁺σ⁻ - 1 (for use_σpm(true))
function _pauli_pm(name, ii, l)
    p, m = QuExpr(BaseOperator(QA.TLSCreate_, name, collect(ii))), QuExpr(BaseOperator(QA.TLSDestroy_, name, collect(ii)))
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
                push!(ps, Param(name, st, collect(ii)))
            end
        end
        ops = BaseOperator[]
        for (i, (name, ii)) in enumerate(b.bosons)
            w, sh = _bpos(b, i)
            append!(ops, fill(BaseOperator(QA.BosonCreate_, name, collect(ii)), Int((key[w] >> sh) & 0xff)))
            append!(ops, fill(BaseOperator(QA.BosonDestroy_, name, collect(ii)), Int((key[w] >> (sh+8)) & 0xff)))
        end
        paulis = Tuple{QuOpName,IndsTuple,Int}[]
        for (i, (name, ii)) in enumerate(b.tlss)
            w, sh = _tpos(b, i)
            l = Int((key[w] >> sh) & 3)
            l == 0 && continue
            pm ? push!(paulis, (name, ii, l)) : push!(ops, BaseOperator(_PAULI_TYPES[l], name, collect(ii)))
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
