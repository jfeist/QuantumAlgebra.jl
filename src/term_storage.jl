# Storage of terms and their constituents.
#
# A QuTerm is stored as a flat vector of 64-bit words (plus a cached hash). The parts of a term
# (Kronecker deltas, parameters, expectation values, correlators, and "bare" operators) are views
# into this vector, so that accessing them does not allocate. Operators, parameters etc. that are
# created on their own (e.g., by the user) are stored in their own short word vector.
#
# Layout of a term:
#   w[1]: nsuminds:16 | number of δs:16 | number of params:16 | number of expvals:16
#   w[2]: number of corrs:16 | number of bare operators:16 | total number of operators:16 | 0:16
#   w[3]: positions (in w) of the first param:16 | expval:16 | corr:16 | bare operator:16
#   then the δs (2 index words each), params, expvals, corrs and bare operators
# Each operator is stored as a header word followed by one word per index:
#   header: kind:8 | name id:24 | number of indices:8 | extra:24
# where `extra` is available for additional (family-specific) information.
# Parameters are stored in the same way, with header
#   header: 0:8 | name id:24 | number of indices:8 | state (as UInt8):24
# Expectation values and correlators are stored as a header word followed by the operators:
#   header: number of operators:32 | number of words of the operators:32
# Indices are stored such that comparing the words gives the same order as comparing the indices.

@static if VERSION >= v"1.11"
    const Words = Memory{UInt64}
else
    const Words = Vector{UInt64}
end
_words(n::Integer) = Words(undef, n)

const _TERM_HEADER_WORDS = 3

#################################################################
## words for indices and headers                              ##
#################################################################

@inline idxword(i::QuIndex) = (UInt64(reinterpret(UInt32, i.sym)) << 32) | UInt64(reinterpret(UInt32, i.num) ⊻ 0x80000000)
@inline wordidx(w::UInt64) = QuIndex(reinterpret(Char, UInt32(w >> 32)), reinterpret(IndexInt, UInt32(w & 0xffffffff) ⊻ 0x80000000))

@inline function _itemhead(kind::Integer, name::QuOpName, nind::Integer, extra::Integer)
    nind ≤ 0xff || throw(ArgumentError("at most 255 indices are supported, got $nind"))
    (UInt64(kind) << 56) | (UInt64(name.i) << 32) | (UInt64(nind) << 24) | UInt64(extra)
end
@inline _headkind(h::UInt64) = reinterpret(BaseOpType, Int32(h >> 56))
@inline _headname(h::UInt64) = QuOpName(IndexInt((h >> 32) & 0xffffff), Val(:id))
@inline _headnind(h::UInt64) = Int((h >> 24) & 0xff)
@inline _headextra(h::UInt64) = UInt32(h & 0xffffff)

@inline _grouphead(nops::Integer, nwords::Integer) = (UInt64(nops) << 32) | UInt64(nwords)
@inline _groupnops(h::UInt64) = Int(h >> 32)
@inline _groupnwords(h::UInt64) = Int(h & 0xffffffff)

#################################################################
## views of indices and Kronecker deltas                      ##
#################################################################

"`IndexView`: the indices of an operator or parameter (a read-only view into its storage)"
struct IndexView <: AbstractVector{QuIndex}
    w::Words
    off::Int
    n::Int
end
Base.size(v::IndexView) = (v.n,)
Base.@propagate_inbounds Base.getindex(v::IndexView, i::Int) = (@boundscheck checkbounds(v,i); wordidx(@inbounds v.w[v.off+i-1]))
Base.IndexStyle(::Type{IndexView}) = IndexLinear()
Base.tail(v::IndexView) = (@assert v.n > 0; IndexView(v.w, v.off+1, v.n-1))

"`DeltaView`: the Kronecker deltas of a term (a read-only view into its storage)"
struct DeltaView <: AbstractVector{δ}
    w::Words
    off::Int
    n::Int
end
Base.size(v::DeltaView) = (v.n,)
Base.@propagate_inbounds function Base.getindex(v::DeltaView, i::Int)
    @boundscheck checkbounds(v,i)
    k = v.off + 2(i-1)
    @inbounds δ(wordidx(v.w[k]), wordidx(v.w[k+1]))
end
Base.IndexStyle(::Type{DeltaView}) = IndexLinear()

#################################################################
## operators, parameters, expectation values, correlators     ##
#################################################################

@concrete struct BaseOperator
    w::Words
    off::Int
end
@concrete struct Param
    w::Words
    off::Int
end
@concrete struct ExpVal
    w::Words
    off::Int
end
@concrete struct Corr
    w::Words
    off::Int
end

const _IndexedItem = Union{BaseOperator,Param}
const _GroupItem = Union{ExpVal,Corr}
const _Item = Union{_IndexedItem,_GroupItem}

@inline _head(A::_Item) = @inbounds A.w[A.off]
# number of words used to store an item
@inline _nwords(A::_IndexedItem) = 1 + _headnind(_head(A))
@inline _nwords(A::_GroupItem) = 1 + _groupnwords(_head(A))
# number of words used by the item stored at w[off]
@inline _nwords(::Type{<:_IndexedItem}, w::Words, off::Int) = 1 + _headnind(@inbounds w[off])
@inline _nwords(::Type{<:_GroupItem}, w::Words, off::Int) = 1 + _groupnwords(@inbounds w[off])

"`ItemSeq{T}`: `n` consecutive items of type `T` (stored in `nw` words starting at `w[off]`)"
struct ItemSeq{T<:_Item}
    w::Words
    off::Int
    n::Int
    nw::Int
end
const OpSeq = ItemSeq{BaseOperator}

Base.length(s::ItemSeq) = s.n
Base.isempty(s::ItemSeq) = s.n == 0
Base.eltype(::Type{ItemSeq{T}}) where T = T
# the iteration state is (position of the next item, number of items already visited)
@inline Base.iterate(s::ItemSeq) = iterate(s, (s.off, 0))
@inline function Base.iterate(s::ItemSeq{T}, state::Tuple{Int,Int}) where T
    off, k = state
    k == s.n && return nothing
    (T(s.w,off), (off + _nwords(T,s.w,off), k+1))
end
function Base.getindex(s::ItemSeq{T}, i::Int) where T
    1 ≤ i ≤ s.n || throw(BoundsError(s,i))
    off = s.off
    for _ = 1:i-1
        off += _nwords(T,s.w,off)
    end
    T(s.w,off)
end
Base.firstindex(s::ItemSeq) = 1
Base.lastindex(s::ItemSeq) = s.n
Base.collect(s::ItemSeq{T}) where T = collect(T, s)
Base.broadcastable(s::ItemSeq) = collect(s)

#################################################################
## comparison and hashing of stored items                     ##
#################################################################

# compare words w1[off1:off1+n-1] and w2[off2:off2+n-1]
@inline function _words_equal(w1::Words, off1::Int, w2::Words, off2::Int, n::Int)
    @inbounds for k = 0:n-1
        w1[off1+k] == w2[off2+k] || return false
    end
    true
end

# a fast hash of the words w[off:off+n-1]
@inline function _hash_words(w::Words, off::Int, n::Int, h::UInt)
    h = UInt64(h)
    @inbounds for k = off:off+n-1
        h = _hashmix(h ⊻ w[k])
    end
    UInt(h)
end
@inline function _hashmix(x::UInt64)
    p = widemul(x, 0x9e3779b97f4a7c15)
    (p % UInt64) ⊻ ((p >> 64) % UInt64)
end

Base.:(==)(A::T,B::T) where T<:_Item = _head(A) == _head(B) && _words_equal(A.w, A.off, B.w, B.off, _nwords(A))
Base.hash(A::T,h::UInt) where T<:_Item = _hash_words(A.w, A.off, _nwords(A), hash(T,h))

Base.:(==)(A::ItemSeq{T},B::ItemSeq{T}) where T = A.n == B.n && A.nw == B.nw && _words_equal(A.w, A.off, B.w, B.off, A.nw)
Base.hash(A::ItemSeq{T},h::UInt) where T = _hash_words(A.w, A.off, A.nw, hash(A.n,h))

#################################################################
## accessors                                                   ##
#################################################################

"`optype(A::BaseOperator)`: the `BaseOpType` of an operator"
optype(A::BaseOperator) = _headkind(_head(A))
"`opname(A)`: the name (as a `QuOpName`) of an operator or parameter"
opname(A::_IndexedItem) = _headname(_head(A))
"`inds(A)`: the indices of an operator or parameter"
inds(A::_IndexedItem) = IndexView(A.w, A.off+1, _headnind(_head(A)))
"`opextra(A::BaseOperator)`: additional (family-specific) information stored with an operator"
opextra(A::BaseOperator) = _headextra(_head(A))
"`parstate(A::Param)`: one of 'r' (real), 'n' (normal/complex) or 'c' (complex conjugate)"
parstate(A::Param) = Char(_headextra(_head(A)))
"`ops(A::Union{ExpVal,Corr})`: the product of operators inside the expectation value or correlator"
ops(A::_GroupItem) = (h = _head(A); OpSeq(A.w, A.off+1, _groupnops(h), _groupnwords(h)))

#################################################################
## construction of standalone items                           ##
#################################################################

function BaseOperator(t::BaseOpType, name::QuOpName, inds::AbstractVector{QuIndex}, extra::Integer=0)
    0 ≤ extra < 2^24 || throw(ArgumentError("extra has to fit in 24 bits, got $extra"))
    n = length(inds)
    w = _words(1+n)
    @inbounds w[1] = _itemhead(Int(t), name, n, extra)
    for (k,ind) in enumerate(inds)
        @inbounds w[1+k] = idxword(ind)
    end
    BaseOperator(w, 1)
end
BaseOperator(t::BaseOpType, name, inds...) = BaseOperator(t, QuOpName(name), make_indices(inds...))

function Param(name::QuOpName, state::Char, inds::AbstractVector{QuIndex})
    state ∈ ('r','n','c') || throw(ArgumentError("state has to be one of n,r,c"))
    n = length(inds)
    w = _words(1+n)
    @inbounds w[1] = _itemhead(0, name, n, UInt8(state))
    for (k,ind) in enumerate(inds)
        @inbounds w[1+k] = idxword(ind)
    end
    Param(w, 1)
end
Param(name, state::Char, inds...) = Param(QuOpName(name), state, make_indices(inds...))

#################################################################
## operator products                                           ##
#################################################################

# a mutable product of operators (used when modifying operator products, e.g., in normal ordering)
@concrete struct BaseOpProduct
    v::Vector{BaseOperator}
end
BaseOpProduct() = BaseOpProduct(BaseOperator[])
BaseOpProduct(s::OpSeq) = BaseOpProduct(collect(s))

Base.isempty(A::BaseOpProduct) = isempty(A.v)
Base.length(A::BaseOpProduct) = length(A.v)
Base.iterate(A::BaseOpProduct, args...) = iterate(A.v, args...)
Base.eltype(::Type{BaseOpProduct}) = BaseOperator
Base.getindex(A::BaseOpProduct, i) = A.v[i]
Base.firstindex(A::BaseOpProduct) = firstindex(A.v)
Base.lastindex(A::BaseOpProduct) = lastindex(A.v)
Base.eachindex(A::BaseOpProduct) = eachindex(A.v)
Base.collect(A::BaseOpProduct) = collect(A.v)
Base.copy(A::BaseOpProduct) = BaseOpProduct(copy(A.v))

"any (stored or mutable) product of operators"
const OpProduct = Union{BaseOpProduct,OpSeq}

Base.:(==)(A::OpProduct,B::OpProduct) = length(A) == length(B) && all(splat(==), zip(A,B))
Base.:(==)(A::OpSeq,B::OpSeq) = A.n == B.n && A.nw == B.nw && _words_equal(A.w, A.off, B.w, B.off, A.nw)
Base.hash(A::BaseOpProduct,h::UInt) = foldl((h,O) -> hash(O,h), A.v; init=hash(length(A),h))

# number of words used to store the operators
_opwords(A::BaseOpProduct) = sum(_nwords, A.v; init=0)
_opwords(A::OpSeq) = A.nw

function _group(::Type{T}, ops::OpProduct) where T<:_GroupItem
    nw = _opwords(ops)
    w = _words(1+nw)
    @inbounds w[1] = _grouphead(length(ops), nw)
    pos = 2
    for O in ops
        pos = _write!(w, pos, O, nothing)
    end
    T(w, 1)
end
ExpVal(ops::OpProduct) = _group(ExpVal, ops)
Corr(ops::OpProduct) = _group(Corr, ops)

#################################################################
## writing items into word vectors                            ##
#################################################################

# write item x to w starting at w[pos], applying the index function f
# to all indices (unless f === nothing). returns the position after it.
@inline function _write!(w::Words, pos::Int, x::_IndexedItem, f)
    n = _nwords(x)
    if f === nothing
        unsafe_copyto!(w, pos, x.w, x.off, n)
    else
        @inbounds w[pos] = _head(x)
        @inbounds for k = 1:n-1
            w[pos+k] = idxword(f(wordidx(x.w[x.off+k])))
        end
    end
    pos + n
end
@inline function _write!(w::Words, pos::Int, d::δ, f)
    if f === nothing
        @inbounds w[pos], w[pos+1] = idxword(d.iA), idxword(d.iB)
    else
        dn = f(d)
        @inbounds w[pos], w[pos+1] = idxword(dn.iA), idxword(dn.iB)
    end
    pos + 2
end
# expectation values and correlators can be given either as such or as operator products
_groupops(x::_GroupItem) = ops(x)
_groupops(x::OpProduct) = x
_groupwords(x) = 1 + _opwords(_groupops(x))
function _write_group!(w::Words, pos::Int, x, f)
    xops = _groupops(x)
    @inbounds w[pos] = _grouphead(length(xops), _opwords(xops))
    _write_ops!(w, pos+1, xops, f)
end
@inline function _write_ops!(w::Words, pos::Int, xops::OpSeq, f)
    if f === nothing
        unsafe_copyto!(w, pos, xops.w, xops.off, xops.nw)
        pos + xops.nw
    else
        for O in xops
            pos = _write!(w, pos, O, f)
        end
        pos
    end
end
@inline function _write_ops!(w::Words, pos::Int, xops, f)
    for O in xops
        pos = _write!(w, pos, O, f)
    end
    pos
end

#################################################################
## terms                                                       ##
#################################################################

"""
    QuTerm

Represents a single term in a quantum algebra expression.

A term consists of:
- Sum indices (summed-over indices)
- Delta functions
- Parameters
- Expectation values
- Correlators/cumulants
- Bare operators

# See also
[`QuExpr`](@ref), [`BaseOperator`](@ref)
"""
struct QuTerm
    w::Words
    hash::UInt
    QuTerm(w::Words) = new(w, _hash_words(w, 1, length(w), zero(UInt)))
end

@inline _termhead(t::QuTerm, k) = @inbounds t.w[k]
@inline _field16(h::UInt64, k) = Int((h >> (48 - 16k)) & 0xffff)

"`nsuminds(t::QuTerm)`: number of summation indices in the term"
nsuminds(t::QuTerm) = _field16(_termhead(t,1), 0)
_ndeltas(t::QuTerm) = _field16(_termhead(t,1), 1)
_nparams(t::QuTerm) = _field16(_termhead(t,1), 2)
_nexpvals(t::QuTerm) = _field16(_termhead(t,1), 3)
_ncorrs(t::QuTerm) = _field16(_termhead(t,2), 0)
_nbares(t::QuTerm) = _field16(_termhead(t,2), 1)
_nops(t::QuTerm) = _field16(_termhead(t,2), 2)
_offparams(t::QuTerm) = _field16(_termhead(t,3), 0)
_offexpvals(t::QuTerm) = _field16(_termhead(t,3), 1)
_offcorrs(t::QuTerm) = _field16(_termhead(t,3), 2)
_offbares(t::QuTerm) = _field16(_termhead(t,3), 3)

"`deltas(t::QuTerm)`: the Kronecker deltas `δ` in the term"
deltas(t::QuTerm) = DeltaView(t.w, _TERM_HEADER_WORDS+1, _ndeltas(t))
"`params(t::QuTerm)`: the parameters `Param` in the term"
params(t::QuTerm) = ItemSeq{Param}(t.w, _offparams(t), _nparams(t), _offexpvals(t) - _offparams(t))
"`expvals(t::QuTerm)`: the expectation values `ExpVal` in the term"
expvals(t::QuTerm) = ItemSeq{ExpVal}(t.w, _offexpvals(t), _nexpvals(t), _offcorrs(t) - _offexpvals(t))
"`corrs(t::QuTerm)`: the correlators `Corr` in the term"
corrs(t::QuTerm) = ItemSeq{Corr}(t.w, _offcorrs(t), _ncorrs(t), _offbares(t) - _offcorrs(t))
"`bares(t::QuTerm)`: the product of \"bare\" operators (outside expectation values) in the term"
bares(t::QuTerm) = OpSeq(t.w, _offbares(t), _nbares(t), length(t.w) + 1 - _offbares(t))

Base.:(==)(A::QuTerm,B::QuTerm) = A.hash == B.hash && length(A.w) == length(B.w) && _words_equal(A.w, 1, B.w, 1, length(A.w))
Base.hash(A::QuTerm,h::UInt) = hash(A.hash,h)
Base.isempty(A::QuTerm) = length(A.w) == _TERM_HEADER_WORDS && nsuminds(A) == 0
Base.length(A::QuTerm) = _nops(A)

_check16(n, what) = (0 ≤ n ≤ 0xffff || throw(ArgumentError("too many $what in term ($n)")); UInt64(n))

"""
    _encode_term(nsuminds, δs, params, expvals, corrs, bares, fδ, f)

Create a `QuTerm` from its parts. Each part is a tuple of collections whose concatenation gives
the respective part of the term. The index function `fδ` is applied to the indices of the δs and
`f` to all other indices (`nothing` stands for no change). Index functions are applied in the order
in which the indices are stored, i.e., δs, params, expvals, corrs, bares.
"""
function _encode_term(nsuminds::Integer, δs::Tuple, params::Tuple, expvals::Tuple, corrs::Tuple, bares::Tuple, fδ, f)
    nδ = sum(length, δs)
    np = sum(length, params)
    ne = sum(length, expvals)
    nc = sum(length, corrs)
    nb = sum(length, bares)
    nwp = sum(ps -> sum(_nwords, ps; init=0), params)
    nwe = sum(es -> sum(_groupwords, es; init=0), expvals)
    nwc = sum(cs -> sum(_groupwords, cs; init=0), corrs)
    nwb = sum(_opwords_any, bares)
    nops = nb + sum(es -> sum(e -> length(_groupops(e)), es; init=0), expvals) + sum(cs -> sum(c -> length(_groupops(c)), cs; init=0), corrs)

    offp = _TERM_HEADER_WORDS + 1 + 2nδ
    offe = offp + nwp
    offc = offe + nwe
    offb = offc + nwc
    w = _words(offb + nwb - 1)
    @inbounds w[1] = (_check16(nsuminds,"sum indices") << 48) | (_check16(nδ,"δs") << 32) | (_check16(np,"params") << 16) | _check16(ne,"expvals")
    @inbounds w[2] = (_check16(nc,"corrs") << 48) | (_check16(nb,"operators") << 32) | (_check16(nops,"operators") << 16)
    @inbounds w[3] = (_check16(offp,"words") << 48) | (_check16(offe,"words") << 32) | (_check16(offc,"words") << 16) | _check16(offb,"words")
    pos = _TERM_HEADER_WORDS + 1
    for ds in δs, d in ds
        pos = _write!(w, pos, d, fδ)
    end
    for ps in params, p in ps
        pos = _write!(w, pos, p, f)
    end
    for es in expvals, e in es
        pos = _write_group!(w, pos, e, f)
    end
    for cs in corrs, c in cs
        pos = _write_group!(w, pos, c, f)
    end
    for bs in bares
        pos = _write_ops!(w, pos, bs, f)
    end
    @assert pos == length(w) + 1
    QuTerm(w)
end
_opwords_any(x::OpProduct) = _opwords(x)
_opwords_any(x) = sum(_nwords, x; init=0)

QuTerm(nsuminds::Integer,δs,params,expvals,corrs,bares) = _encode_term(nsuminds,(δs,),(params,),(expvals,),(corrs,),(bares,),nothing,nothing)

QuTerm() = QuTerm(BaseOpProduct())
QuTerm(op::BaseOperator) = QuTerm(BaseOpProduct([op]))
QuTerm(ops::OpProduct) = QuTerm(0,δ[],Param[],ExpVal[],Corr[],ops)
QuTerm(δs::AbstractVector{δ},ops::OpProduct) = QuTerm(0,δs,Param[],ExpVal[],Corr[],ops)
QuTerm(p::Param)  = QuTerm(0,δ[],Param[p],ExpVal[],Corr[],BaseOpProduct())
QuTerm(E::ExpVal) = QuTerm(0,δ[],Param[],ExpVal[E],Corr[],BaseOpProduct())
QuTerm(C::Corr)   = QuTerm(0,δ[],Param[],ExpVal[],Corr[C],BaseOpProduct())
QuTerm(Es::Union{AbstractVector{ExpVal},ItemSeq{ExpVal}}) = QuTerm(0,δ[],Param[],Es,Corr[],BaseOpProduct())
QuTerm(Cs::Union{AbstractVector{Corr},ItemSeq{Corr}}) = QuTerm(0,δ[],Param[],ExpVal[],Cs,BaseOpProduct())
QuTerm(δs::AbstractVector{δ},Es::Union{AbstractVector{ExpVal},ItemSeq{ExpVal}}) = QuTerm(0,δs,Param[],Es,Corr[],BaseOpProduct())
QuTerm(δs::AbstractVector{δ},Cs::Union{AbstractVector{Corr},ItemSeq{Corr}}) = QuTerm(0,δs,Param[],ExpVal[],Cs,BaseOpProduct())

# apply the index function f to all indices of t (in the order in which they are stored), and set the number of sum indices
_map_term_indices(f, t::QuTerm, nsuminds::Integer) =
    _encode_term(nsuminds, (deltas(t),), (params(t),), (expvals(t),), (corrs(t),), (bares(t),), f, f)

# fold op over all indices of t (in the order in which they are stored), i.e., compute
# op(...op(op(init, ind1), ind2)..., indN), optionally leaving out the indices of the δs
@inline function _foldl_indices(op::F, acc, t::QuTerm, withδs::Bool=true) where F
    w = t.w
    if withδs
        @inbounds for k = _TERM_HEADER_WORDS+1 : _offparams(t)-1
            acc = op(acc, wordidx(w[k]))
        end
    end
    # params, expvals, corrs and bares: skip the headers of items and groups
    pos = _offparams(t)
    @inbounds for _ = 1:_nparams(t)
        n = _headnind(w[pos])
        for k = pos+1:pos+n
            acc = op(acc, wordidx(w[k]))
        end
        pos += n+1
    end
    @inbounds for _ = 1:_nexpvals(t)+_ncorrs(t)
        nops = _groupnops(w[pos])
        pos += 1
        for _ = 1:nops
            n = _headnind(w[pos])
            for k = pos+1:pos+n
                acc = op(acc, wordidx(w[k]))
            end
            pos += n+1
        end
    end
    @inbounds while pos ≤ length(w)
        n = _headnind(w[pos])
        for k = pos+1:pos+n
            acc = op(acc, wordidx(w[k]))
        end
        pos += n+1
    end
    acc
end

"""
    WorkTerm

Mutable working form of a `QuTerm` used while normal ordering. All vectors are owned by
the `WorkTerm` and can be modified freely. Convert back with `QuTerm(w::WorkTerm)`.
"""
struct WorkTerm
    nsuminds::IndexInt
    δs::Vector{δ}
    params::Vector{Param}
    expvals::Vector{BaseOpProduct}
    corrs::Vector{BaseOpProduct}
    bares::BaseOpProduct
end
WorkTerm(t::QuTerm) = WorkTerm(nsuminds(t), collect(deltas(t)), collect(params(t)), [BaseOpProduct(ops(E)) for E in expvals(t)],
                               [BaseOpProduct(ops(C)) for C in corrs(t)], BaseOpProduct(bares(t)))
QuTerm(w::WorkTerm) = QuTerm(w.nsuminds, w.δs, w.params, w.expvals, w.corrs, w.bares)
