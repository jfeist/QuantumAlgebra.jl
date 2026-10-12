# Canonical numbering of sum indices.
#
# Terms that only differ in the numbering of their sum indices (which are dummy variables) are equal,
# and should be represented in the same way so that they can be combined. To ensure this, the sum
# indices of normal-ordered terms are numbered in a canonical way: Each sum index gets a signature
# that only depends on how it appears in the term (but not on the numbering of the sum indices), and
# sum indices are numbered in the order of their signatures. If several sum indices have the same
# signature, all their possible numberings are tried, and the one that gives the "smallest" term is
# chosen. After renumbering, the operators in a term have to be sorted again. Since this only
# reorders operators of the same type (with the same name and extra information), it does not
# generate additional terms, but can give a sign (for fermions). If renumbering a term gives the same
# term with the opposite sign, the term vanishes.

# signature contributions for sum indices that appear in the same item as the sum index under consideration
const _SIG_SELF = 0x5bd1e9955bd1e995
const _SIG_OTHER = 0x27d4eb2f165667c5

# a hash of the header word of an item that does not depend on the order in which names were registered
# (or on the numerical values of operator types)
const _KindStableHash = Tuple(_fnv1a(string(t)) for t in instances(BaseOpType))
@inline function _stable_head(h::UInt64)
    k = Int(h >> 56)
    kh = k == 0 ? zero(UInt64) : @inbounds(_KindStableHash[k])
    (h & ~(UInt64(0xffffffff) << 32)) ⊻ kh ⊻ _namestablehash(_headname(h))
end

# add the contributions of the occurrences of sum indices in the item with index words w[off+1:off+nind] and
# label-independent description desc to the signatures in newsig (using the signatures sig of the previous round)
@inline function _add_occurrences!(newsig::Vector{UInt64}, sig::Vector{UInt64}, w::Words, off::Int, nind::Int, desc::UInt64)
    n = length(newsig)
    @inbounds for p = 1:nind
        ind = wordidx(w[off+p])
        issumindex(ind) || continue
        k = Int(ind.num)
        1 ≤ k ≤ n || throw(ArgumentError("invalid sum index $ind in term with $n sum indices"))
        h = _hashmix(desc ⊻ UInt64(p))
        for q = 1:nind
            q == p && continue
            iq = wordidx(w[off+q])
            x = !issumindex(iq) ? w[off+q] : iq.num == k ? _SIG_SELF : _SIG_OTHER ⊻ (1 ≤ iq.num ≤ n ? sig[iq.num] : throw(ArgumentError("invalid sum index $iq")))
            h = _hashmix(h ⊻ x ⊻ (UInt64(q) << 56))
        end
        newsig[k] += h
    end
end

# one round of refinement: hashes newc of the sum indices from their old colors and how they appear in the term
function _refine_colors_round!(newc::Vector{UInt64}, t::QuTerm, colors::Vector{UInt64})
    newc .= _hashmix.(colors .⊻ _SIG_SELF)
    w = t.w
    for k = _TERM_HEADER_WORDS+1 : 2 : _offparams(t)-1
        _add_occurrences!(newc, colors, w, k-1, 2, UInt64(0x01))
    end
    pos = _offparams(t)
    for _ = 1:_nparams(t)
        nind = _headnind(w[pos])
        _add_occurrences!(newc, colors, w, pos, nind, _stable_head(w[pos]) ⊻ 0x02)
        pos += nind+1
    end
    for (ngroups, seg) in ((_nexpvals(t), 0x03), (_ncorrs(t), 0x04))
        for _ = 1:ngroups
            nops = _groupnops(w[pos])
            # describe the group by the (unordered) operators it contains, without their indices
            gh = UInt64(nops)
            p = pos+1
            for _ = 1:nops
                gh += _hashmix(_stable_head(w[p]))
                p += _headnind(w[p])+1
            end
            p = pos+1
            for _ = 1:nops
                nind = _headnind(w[p])
                _add_occurrences!(newc, colors, w, p, nind, _hashmix(_stable_head(w[p]) ⊻ gh) ⊻ seg)
                p += nind+1
            end
            pos = p
        end
    end
    while pos ≤ length(w)
        nind = _headnind(w[pos])
        _add_occurrences!(newc, colors, w, pos, nind, _stable_head(w[pos]) ⊻ 0x05)
        pos += nind+1
    end
    newc
end

# set r to the (dense) ranks of the values in v, i.e., r[k] = 1 + the number of distinct values smaller than v[k],
# using the buffer `sorted` (the number of sum indices is small, so this is quadratic in their number)
function _ranks!(r::Vector{UInt64}, v::AbstractVector, sorted::AbstractVector)
    copyto!(sorted, v)
    sort!(sorted)
    for k in eachindex(v)
        rank = 1
        for j = 2:searchsortedfirst(sorted, v[k])
            sorted[j] != sorted[j-1] && (rank += 1)
        end
        r[k] = rank
    end
    r
end
_ranks(v::Vector{UInt64}) = _ranks!(similar(v), v, similar(v))
# the number of different colors (for dense ranks)
_ncolors(colors::Vector{UInt64}) = isempty(colors) ? 0 : Int(maximum(colors))
_allunique_colors(colors::Vector{UInt64}) = _ncolors(colors) == length(colors)

# Refine the colors (one per sum index, given as dense ranks) until the partition of the sum indices into colors
# does not change anymore. Refinement only splits groups of sum indices with equal colors, i.e., the order of the
# groups is kept.
function _refine_colors(t::QuTerm, colors::Vector{UInt64})
    n = length(colors)
    newc = Vector{UInt64}(undef, n)
    pairs = Vector{Tuple{UInt64,UInt64}}(undef, n)
    sorted = similar(pairs)
    colors = copy(colors)
    ncolors = _ncolors(colors)
    while ncolors < n
        _refine_colors_round!(newc, t, colors)
        pairs .= tuple.(colors, newc)
        _ranks!(colors, pairs, sorted)
        _ncolors(colors) == ncolors && break
        ncolors = _ncolors(colors)
    end
    colors
end

"`_suminds_signatures(t::QuTerm)`: signatures of the sum indices of `t` that do not depend on their numbering"
function _suminds_signatures(t::QuTerm)
    n = nsuminds(t)
    _refine_colors(t, ones(UInt64, n))
end

# index function that renumbers sum index k as λ[k]
struct _renumber_suminds <: IndexFunction
    λ::Vector{IndexInt}
end
(f::_renumber_suminds)(ind::QuIndex) = issumindex(ind) ? sumindex(f.λ[ind.num]) : ind

# sort a vector of operators by exchanging neighbors (which have to be of the same type), and return the sign
function _sort_with_sign!(A::Vector{BaseOperator})
    sgn = 1
    for i = 2:length(A)
        j = i
        while j>1 && A[j]<A[j-1]
            @assert optype(A[j]) == optype(A[j-1])
            sgn *= exchange!(_NoEmit(),A[j-1],A[j])
            A[j-1], A[j] = A[j], A[j-1]
            j -= 1
        end
    end
    sgn
end

# renumber the sum indices of the normal-ordered term t (sum index k becomes λ[k]) and sort it again.
# returns the new term and the sign from reordering operators
function _renumbered_term(t::QuTerm, λ::Vector{IndexInt})
    W = WorkTerm(_map_term_indices(_renumber_suminds(λ), t, nsuminds(t)))
    sgn = 1
    sort!(W.params)
    for EC in (W.expvals, W.corrs)
        for P in EC
            sgn *= _sort_with_sign!(P.v)
        end
        sort!(EC)
    end
    sgn *= _sort_with_sign!(W.bares.v)
    QuTerm(W), sgn
end

# terms are compared by their storage (any total order works, as long as it is consistent)
function _words_isless(A::QuTerm, B::QuTerm)
    length(A.w) == length(B.w) || return length(A.w) < length(B.w)
    @inbounds for k in eachindex(A.w)
        A.w[k] == B.w[k] || return A.w[k] < B.w[k]
    end
    false
end

# operators in the term are in normal order (without considering sum indices)
_ops_in_normal_form(t::QuTerm) = is_normal_form(bares(t)) && all(is_normal_form, expvals(t)) && all(is_normal_form, corrs(t))

# numbering of the sum indices in the order of their colors (with ties broken by their current numbering)
_numbering(colors::Vector{UInt64}) = (λ = zeros(IndexInt, length(colors)); λ[sortperm(colors)] .= 1:length(colors); λ)

# maximal number of numberings that are tried to find the canonical one. If there are more (for terms with many
# equivalent sum indices), sum indices are only numbered in the order of their signatures (and terms that only differ
# in the numbering of such sum indices might not be identified)
const _MAX_CANONICAL_NUMBERINGS = 1000

mutable struct _CanonicalSearch
    const t::QuTerm
    best::QuTerm
    bestsgn::Int
    vanishes::Bool
    first::QuTerm
    firstsgn::Int
    nleaves::Int
end

# the first group of sum indices with the same color (in the order of colors)
function _first_tied_color(colors::Vector{UInt64})
    for c = UInt64(1):UInt64(_ncolors(colors))
        count(==(c), colors) > 1 && return c
    end
    error("no tied colors")
end
# give sum index v a color that comes directly before the others with the same color, and refine
_individualize(t::QuTerm, colors::Vector{UInt64}, v) = (c2 = 2 .* colors; c2[v] -= 1; _refine_colors(t, _ranks(c2)))

# the term (and sign) for the leaf of the search tree that is reached by always individualizing the first sum index
function _greedy_leaf(t::QuTerm, colors::Vector{UInt64})
    while !_allunique_colors(colors)
        colors = _individualize(t, colors, findfirst(==(_first_tied_color(colors)), colors))
    end
    _renumbered_term(t, _numbering(colors))
end

# Search for the canonical numbering by individualization and refinement: if several sum indices have the same color,
# try giving each of them (in the first such group) a unique color in turn, refine the colors and continue recursively.
# Each numbering in which all colors are different is a candidate, and the one that gives the smallest term is chosen.
# Returns the term and sign of the first leaf in the search tree.
function _canonical_search!(S::_CanonicalSearch, colors::Vector{UInt64})
    if _allunique_colors(colors)
        S.nleaves += 1
        c, sgn = _renumbered_term(S.t, _numbering(colors))
        S.firstsgn == 0 && (S.first = c; S.firstsgn = sgn)
        # renumbering can make contractions possible for operators that are not neighbors (e.g., for Pauli matrices)
        if _ops_in_normal_form(c)
            if S.bestsgn == 0 || _words_isless(c, S.best)
                S.best, S.bestsgn = c, sgn
            elseif c == S.best && sgn != S.bestsgn
                S.vanishes = true
            end
        end
        return (c, sgn)
    end
    vs = findall(==(_first_tied_color(colors)), colors)
    firstleaf = _canonical_search!(S, _individualize(S.t, colors, vs[1]))
    for v in @view vs[2:end]
        S.nleaves > _MAX_CANONICAL_NUMBERINGS && break
        c2 = _individualize(S.t, colors, v)
        # if the first leaf below v gives the same term as the first leaf below vs[1], there is a symmetry
        # (renumbering of sum indices) of the term that maps v to vs[1], and the subtrees are equivalent
        leaf = _greedy_leaf(S.t, c2)
        if leaf[1] == firstleaf[1]
            # if the symmetry changes the sign, the term vanishes
            leaf[2] != firstleaf[2] && (S.vanishes = true)
            continue
        end
        _canonical_search!(S, c2)
    end
    firstleaf
end

function _canonical_search(t::QuTerm, sig::Vector{UInt64})
    S = _CanonicalSearch(t, t, 0, false, t, 0, 0)
    _canonical_search!(S, sig)
    S
end
_search_exceeded(S::_CanonicalSearch) = S.nleaves > _MAX_CANONICAL_NUMBERINGS

"""
    _canonicalize_suminds(t::QuTerm)

Renumber the sum indices of the normal-ordered term `t` canonically. Returns `(tc, sign, isnormal)` such that
`t = sign*tc`, where `isnormal` is true if `tc` is in normal form, or `nothing` if the term vanishes (because
renumbering gives the same term with opposite sign).
"""
function _canonicalize_suminds(t::QuTerm)
    nsuminds(t) ≤ 1 && return (_reorder_suminds(t), 1, false)
    sig = _suminds_signatures(t)
    if _allunique_colors(sig)
        # if the numbering is already canonical, nothing needs to be done
        issorted(sig) && return (t, 1, _ops_in_normal_form(t))
        c, sgn = _renumbered_term(t, _numbering(sig))
        return (c, sgn, _ops_in_normal_form(c))
    end
    S = _canonical_search(t, sig)
    if _search_exceeded(S)
        # (this is accepted as canonical by _suminds_canonical)
        c, sgn = _renumbered_term(t, _numbering(sig))
        return (c, sgn, _ops_in_normal_form(c))
    end
    # if no renumbering gives operators in normal order, return the first one (which then has to be normal ordered)
    S.bestsgn == 0 && return (S.first, S.firstsgn, false)
    S.vanishes ? nothing : (S.best, S.bestsgn, true)
end

"`_suminds_canonical(t::QuTerm)`: whether the sum indices of the normal-ordered term `t` are numbered canonically"
function _suminds_canonical(t::QuTerm)
    n = nsuminds(t)
    if n ≤ 1
        # the sum index has to be numbered 1
        _, ordered = _foldl_indices(_check_sumind_order, (sumindex(0).num, true), t, false)
        return ordered
    end
    sig = _suminds_signatures(t)
    # the signatures have to be in increasing order
    issorted(sig) || return false
    # if all signatures are different, the numbering is unique
    _allunique_colors(sig) && return true
    S = _canonical_search(t, sig)
    # if there are too many possible numberings, only the order of signatures is required
    _search_exceeded(S) && return true
    # otherwise, t has to be the smallest of the possible terms, and must not vanish
    S.bestsgn == 1 && !S.vanishes && S.best == t
end
