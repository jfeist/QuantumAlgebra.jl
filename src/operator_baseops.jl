# we will want to overload these operators and functions for our custom types
import Base: ==, ≈, *, +, -

export normal_form, comm

Base.length(A::Union{ExpVal,Corr}) = _groupnops(_head(A))

Base.zero(::Type{QuExpr}) = QuExpr()
Base.zero(::QuExpr) = zero(QuExpr)
Base.one(::Type{QuExpr}) = QuExpr(QuTerm())
Base.one(::Type{QuTerm}) = QuExpr(QuTerm())
Base.one(::T) where T<:Union{QuTerm,QuExpr} = one(T)

function Base.:^(A::Union{QuTerm,QuExpr},n::Int)
    if n > 0
        prod(A for _=1:n)
    elseif n == 0
        one(A)
    else
        throw(ArgumentError("Only positive exponents supported for QuantumAlgebra objects."))
    end
end
Base.literal_pow(::typeof(^),A::Union{QuTerm,QuExpr},::Val{n}) where n = A^n

==(A::QuExpr,B::QuExpr) = isequal(A.terms, B.terms)

function ≈(A::QuExpr,B::QuExpr)
    length(A.terms) != length(B.terms) && return false
    for (tA,sA) in A.terms
        sB = get(B.terms,tA) do
            return false
        end
        sA ≈ sB || return false
    end
    return true
end

Base.hash(ind::QuIndex,h::UInt) = hash(idxword(ind),h)
Base.hash(name::QuOpName,h::UInt) = hash(name.i,h)
function Base.hash(A::QuExpr,h::UInt)
    h = hash(A.terms,h)
    h
end

# Base.cmp(::AbstractArray,::AbstractArray) uses isequal and isless, so doesn't shortcut if a,b are themselves arrays
# https://github.com/JuliaLang/julia/blob/539f3ce943f59dec8aff3f2238b083f1b27f41e5/base/abstractarray.jl
function recursive_cmp(A, B)
    for (a, b) in zip(A, B)
        cc = recursive_cmp(a,b)
        cc == 0 || return cc
    end
    return cmp(length(A), length(B))
end

@inline recursive_cmp(A::BaseOperator,B::BaseOperator) = cmp(A,B)
@inline recursive_cmp(A::Param,B::Param) = cmp(A,B)
@inline recursive_cmp(A::δ,B::δ) = cmp(A,B)
# operator products are ordered first by their length
@inline function recursive_cmp(A::OpProduct,B::OpProduct)
    cc = cmp(length(A),length(B))
    cc == 0 || return cc
    for (a, b) in zip(A, B)
        cc = cmp(a,b)
        cc == 0 || return cc
    end
    return 0
end
@inline recursive_cmp(A::T,B::T) where T<:Union{ExpVal,Corr}  = recursive_cmp(ops(A),ops(B))

# compare f(A) and f(B), and return from the calling function if they are not equal
macro _cmpAB(f,rec=true)
    fun = rec ? :recursive_cmp : :cmp
    esc(quote
        cc = $fun($f(A),$f(B))
        cc == 0 || return cc
    end)
end

# compare the indices stored in words wA[offA+1:offA+nA] and wB[offB+1:offB+nB] lexicographically
# (the encoding of the indices is chosen such that comparing the words gives the order of the indices)
@inline function _cmp_index_words(wA::Words, offA::Int, nA::Int, wB::Words, offB::Int, nB::Int)
    @inbounds for k = 1:min(nA,nB)
        a, b = wA[offA+k], wB[offB+k]
        a == b || return a < b ? -1 : 1
    end
    cmp(nA,nB)
end

# order by type, name, extra information, indices
@inline function Base.cmp(A::BaseOperator,B::BaseOperator)
    hA, hB = _head(A), _head(B)
    cc = cmp(hA >> 56, hB >> 56)
    cc == 0 || return cc
    cc = cmp(_headname(hA), _headname(hB))
    cc == 0 || return cc
    cc = cmp(_headextra(hA), _headextra(hB))
    cc == 0 || return cc
    _cmp_index_words(A.w, A.off, _headnind(hA), B.w, B.off, _headnind(hB))
end

# order by name, state, indices
@inline function Base.cmp(A::Param,B::Param)
    hA, hB = _head(A), _head(B)
    cc = cmp(_headname(hA), _headname(hB))
    cc == 0 || return cc
    cc = cmp(parstate(A), parstate(B))
    cc == 0 || return cc
    _cmp_index_words(A.w, A.off, _headnind(hA), B.w, B.off, _headnind(hB))
end

@inline function Base.cmp(A::δ,B::δ)
    cc = cmp(A.iA,B.iA)
    cc == 0 || return cc
    return cmp(A.iB,B.iB)
end

function Base.cmp(A::QuTerm,B::QuTerm)
    # only evaluate each part that we need for each step of the comparison to avoid unnecessary work
    # order first by number of operators (also within expectation values)
    cc = cmp(length(A), length(B))
    cc == 0 || return cc
    # then order by operators, expectation values, correlations, params, deltas, nsuminds
    @_cmpAB bares
    @_cmpAB expvals
    @_cmpAB corrs
    @_cmpAB params
    @_cmpAB deltas
    @_cmpAB nsuminds false
    # they are equal
    return 0
end

@inline Base.isless(A::BaseOperator,B::BaseOperator) = cmp(A,B) < 0
@inline Base.isless(A::Param,B::Param) = cmp(A,B) < 0
@inline Base.isless(A::δ,B::δ) = cmp(A,B) < 0
@inline Base.isless(A::OpProduct,B::OpProduct) = recursive_cmp(A,B) < 0
@inline Base.isless(A::T,B::T) where T<:Union{ExpVal,Corr} = recursive_cmp(A,B) < 0
@inline Base.isless(A::QuTerm,B::QuTerm) = cmp(A,B) < 0

"""
    comm(A, B)

Compute the commutator ``[A,B] = A*B - B*A``.
"""
comm(A,B) = A*B - B*A

Base.adjoint(A::BaseOperator) = BaseOperator(BaseOpType_adj[Int(optype(A))],opname(A),inds(A),opextra(A))
Base.adjoint(A::OpProduct) = BaseOpProduct(reverse!(BaseOperator[adjoint(O) for O in A]))
Base.adjoint(A::T) where {T<:Union{ExpVal,Corr}} = T(adjoint(ops(A)))
Base.adjoint(A::Param) = parstate(A)=='r' ? A : Param(opname(A),parstate(A)=='n' ? 'c' : 'n',inds(A))
function Base.adjoint(A::QuTerm)
    B = QuTerm(nsuminds(A),deltas(A),adjoint.(params(A)),adjoint.(expvals(A)),adjoint.(corrs(A)),adjoint(bares(A)))
    nsuminds(B)>0 ? _reorder_suminds(B) : B
end
Base.adjoint(A::QuExpr) = QuExpr((adjoint(t),adjoint(s)) for (t,s) in A.terms)

is_normal_form(A::Union{ExpVal,Corr}) = is_normal_form(ops(A))
is_normal_form(ops::BaseOpProduct) = is_normal_form(ops.v)
# check that the operators are sorted and that no contractions occur
function is_normal_form(A::Union{Vector{BaseOperator},OpSeq})
    pend = _endpos(A)
    p = _firstpos(A)
    p == pend && return true
    while (pn = _nextpos(A,p)) < pend
        O, On = _opat(A,p), _opat(A,pn)
        isless(On,O) && return false
        contract!(_NoEmit(),O,On) === nothing || return false
        _contraction_partner(A,p) === nothing || return false
        p = pn
    end
    return true
end

function is_normal_form(t::QuTerm)
    (issorted(params(t)) && issorted(expvals(t)) && issorted(corrs(t))) || return false
    is_normal_form(bares(t)) || return false
    for EV in expvals(t)
        is_normal_form(EV) || return false
    end
    for CO in corrs(t)
        is_normal_form(CO) || return false
    end

    # ensure that sum indices are ordered, i.e., numbered in order of appearance (outside of the δs)
    _, ordered = _foldl_indices((sumind_ordered, ind) -> _check_sumind_order(sumind_ordered, ind), (sumindex(0).num, true), t, false)
    ordered || return false

    isempty(deltas(t)) && return true

    # all indices that are not in δs
    changeable_indices_set = Set{QuIndex}(indices(t,false))
    # accumulate indices that δs want to change (to compare with later δs)
    replace_inds = Set{QuIndex}()
    for (iδ,dd) in enumerate(deltas(t))
        iA,iB = dd.iA, dd.iB
        # δ is not ordered (or disappears if iA == iB)
        iA ≥ iB && return false
        # the whole term disappears in cleanup
        isintindex(iA) && isintindex(iB) && return false
        (issumindex(iA) || issumindex(iB)) && return false
        (iA in replace_inds || iB in replace_inds) && return false
        if iδ > 1
            (dd > @inbounds deltas(t)[iδ-1]) || return false
        end
        iB in changeable_indices_set && return false
        push!(replace_inds, iB)
    end
    return true
end
is_normal_form(A::QuExpr) = all(is_normal_form, keys(A.terms))

# helper for is_normal_form: sum indices have to be increasing without jumps
@inline function _check_sumind_order((last_sumind, ordered), ind)
    (ordered && issumindex(ind)) || return (last_sumind, ordered)
    ind.num ≤ last_sumind + oneunit(last_sumind) || return (last_sumind, false)
    (max(ind.num,last_sumind), true)
end

# returns a WorkTerm (that can be modified freely) or nothing if the term is zero
function _normalize_without_commutation(A::QuTerm)::Union{WorkTerm,Nothing}
    # first, clean up the δs
    if isempty(deltas(A))
        W = WorkTerm(A)
    else
        δs = sort(deltas(A))
        #println("in _normalize_without_commutation for A = $A, starting with δs = $δs")
        replacements = Dict{QuIndex,QuIndex}()
        delsuminds = IndexInt[]
        iwrite = 1
        for dd in δs
            iA = get(replacements,dd.iA,dd.iA)
            iB = get(replacements,dd.iB,dd.iB)
            if iA == iB
                # this delta just gives one
                continue
            elseif isintindex(iA) && isintindex(iB)
                # the term is zero
                return nothing
            elseif iA > iB
                iA, iB = iB, iA
            end
            if issumindex(iB) # one sum disappears, and the delta does as well
                replacements[iB] = iA
                # we will later need to shift all larger sumindices down by one
                iB.num ∈ delsuminds || push!(delsuminds,iB.num)
            elseif issumindex(iA) # one sum disappears, and the delta does as well
                replacements[iA] = iB
                # we will later need to shift all larger sumindices down by one
                iA.num ∈ delsuminds || push!(delsuminds,iA.num)
            else
                # otherwise, replace the larger index by the smaller one
                replacements[iB] = iA
                # we keep the delta, but include any possible replacements made
                δs[iwrite] = δ(iA,iB)
                iwrite += 1
            end
        end
        resize!(δs,iwrite-1)
        # after replacements, δs can be out of order
        sort!(δs)
        if !isempty(delsuminds)
            sumrepls = Dict{QuIndex,QuIndex}()
            sort!(delsuminds)
            for (shift,startind) in enumerate(delsuminds)
                endind = shift == length(delsuminds) ? nsuminds(A) : delsuminds[shift+1]-1
                for num = startind+1:endind
                    sumrepls[sumindex(num)] = sumindex(num-shift)
                end
            end
            # we want to first apply replacements and then sumrepls
            # this is the same as replacing results from replacements with sumrepls,
            # and then applying sumrepls on indices that have not been replaced yet
            for (iold,inew) in replacements
                x = get(sumrepls,inew,inew)
                x != inew && (replacements[iold] = x)
            end
            # sumrepl only for indices not in replacements
            merge!(sumrepls,replacements)
            replacements = sumrepls
        end
        # apply the replacements to all indices apart from those of the δs (which have been updated above)
        f = replace_inds(replacements)
        W = WorkTerm(_encode_term(nsuminds(A)-length(delsuminds), (δs,), (params(A),), (expvals(A),), (corrs(A),), (bares(A),), nothing, f))
    end
    # also sort all the commuting terms
    sort!(W.params)
    sort!(W.expvals)
    sort!(W.corrs)
    W
end

# Positions of operators in a sequence of operators: indices for a vector of operators,
# and word offsets for a stored product (so that both can be traversed without allocating)
_firstpos(A::Vector{BaseOperator}) = 1
_endpos(A::Vector{BaseOperator}) = length(A)+1
_nextpos(A::Vector{BaseOperator},p) = p+1
_opat(A::Vector{BaseOperator},p) = @inbounds A[p]
_firstpos(s::OpSeq) = s.off
_endpos(s::OpSeq) = s.off + s.nw
_nextpos(s::OpSeq,p) = p + _nwords(BaseOperator,s.w,p)
_opat(s::OpSeq,p) = BaseOperator(s.w,p)

# For operators where some contractions are not between neighbors in the normal order (e.g., σˣ₁ σˣ₂ σʸ₁),
# find the first operator A[kp] (after A[k+1]) that contracts with A[k] and can be moved next to it,
# i.e., that commutes with all operators in between. Returns the position of A[kp] (given that of A[k]),
# or nothing if there is no such operator.
function _contraction_partner(A,pk)
    Ak = _opat(A,pk)
    lookahead_contraction(Ak) || return nothing
    fam = opfamily(Ak)
    pend = _endpos(A)
    pk1 = _nextpos(A,pk)
    pk1 < pend || return nothing
    pprev, pkp = pk1, _nextpos(A,pk1)
    while pkp < pend
        Akp = _opat(A,pkp)
        (opfamily(Akp) == fam && lookahead_contraction(Akp)) || break
        commutes(_opat(A,pprev),Akp) || break
        # A[kp] has to commute with all operators it is moved past, not just its neighbor
        if _commutes_with_all(A,pk1,pprev,Akp) && contract!(_NoEmit(),Ak,Akp) !== nothing
            return pkp
        end
        pprev, pkp = pkp, _nextpos(A,pkp)
    end
    return nothing
end
# whether B commutes with all operators at positions p with pstart ≤ p < pstop
function _commutes_with_all(A,pstart,pstop,B)
    p = pstart
    while p < pstop
        commutes(_opat(A,p),B) || return false
        p = _nextpos(A,p)
    end
    true
end

# emitter that ignores the passed terms
struct _NoEmit end
(::_NoEmit)(args...) = nothing

# emitter used in normal_order!: e(c, δs, ops...) adds the term obtained from replacing A[i1] and A[i2]
# (with i1 < i2) in the operator product A by ops (at position i1), with prefactor c*e.prefactor
# and the additional Kronecker deltas δs, to e.collector
struct _ReplaceEmitter{T}
    A::Vector{BaseOperator}
    i1::Int
    i2::Int
    collector::QuExpr
    prefactor::T
end
function (e::_ReplaceEmitter)(c, δs::Vector{δ}, newops::BaseOperator...)
    iszero(c) && return nothing
    A, i1, i2 = e.A, e.i1, e.i2
    v = BaseOperator[]
    sizehint!(v, length(A) - 2 + length(newops))
    append!(v, view(A, 1:i1-1))
    append!(v, newops)
    append!(v, view(A, i1+1:i2-1))
    append!(v, view(A, i2+1:length(A)))
    _add_sum_term!(e.collector, QuTerm(δs, BaseOpProduct(v)), c*e.prefactor)
    nothing
end

"""
    normal_order!(ops::BaseOpProduct, term_collector, shortcut_vacA_zero=false, prefactor=1)

Bring the operator product `ops` into normal order in place. Returns the prefactor `p` of the
reordered product, such that the original `prefactor*ops` is equal to `p*ops` (after reordering)
plus the terms that are added to `term_collector` (which still have to be normal ordered).
"""
function normal_order!(ops::BaseOpProduct,term_collector,shortcut_vacA_zero=false,prefactor=one(ComplexInt))
    # do an insertion sort to get to normal ordering
    # reference: https://en.wikipedia.org/wiki/Insertion_sort
    A = ops.v
    for i = 2:length(A)
        j = i
        while j>1 && A[j]<A[j-1]
            # rewrite A[j-1] A[j] = x A[j] A[j-1] + R, where R is added to term_collector
            # (with the prefactor before the exchange)
            x = exchange!(_ReplaceEmitter(A,j-1,j,term_collector,prefactor),A[j-1],A[j])
            prefactor *= x
            iszero(prefactor) && return prefactor
            A[j-1], A[j] = A[j], A[j-1]
            j -= 1
        end
        # ⟨0| A[1] ... = 0
        shortcut_vacA_zero && annihilates_vacuum_left(A[1]) && return zero(prefactor)
    end
    # check if we have any products that simplify
    did_contractions = false
    k = 1
    while k < length(A)
        kp = k+1
        res = contract!(_ReplaceEmitter(A,k,kp,term_collector,prefactor),A[k],A[kp])
        if res === nothing
            kp = something(_contraction_partner(A,k), 0)
            if kp > 0
                res = contract!(_ReplaceEmitter(A,k,kp,term_collector,prefactor),A[k],A[kp])
            end
        end
        if res === nothing
            k += 1
            continue
        end
        did_contractions = true
        fac, op = res
        prefactor *= fac
        iszero(prefactor) && return prefactor
        # replace A[k] A[kp] by op (where op === nothing is the identity)
        if op === nothing
            deleteat!(A,(k,kp))
        else
            A[k] = op
            deleteat!(A,kp)
        end
        # we replaced [...,A[k-1],A[k],A[k+1],...,A[kp],...]
        # with        [...,A[k-1],A[k+1],...] or [...,A[k-1],Anew,A[k+1]...]
        # so reduce k by one to compare the new operator at k with A[k-1] as well
        k > 1 && (k -= 1)
    end
    if did_contractions
        # do another round of normal ordering since this could have changed after contractions
        # the prefactor of the term we were working on needs to be passed to normal_order!
        # to include in the term_collector terms
        prefactor = normal_order!(ops,term_collector,false,prefactor)
    end
    prefactor
end

# normal order (in place) the operator products in v, which are the contents of expectation values
# or correlators (of type T) in a term. Products that become empty are removed, and v is sorted.
function normal_order!(v::Vector{BaseOpProduct},::Type{T},commterms,pref=1) where {T<:Union{Corr,ExpVal}}
    for (ii,EC) in enumerate(v)
        newterms = QuExpr()
        pp = normal_order!(EC,newterms)
        for (t,s) in newterms.terms
            #@assert iszero(nsuminds(t)) && isempty(params(t)) && isempty(expvals(t)) && isempty(corrs(t))
            nv = T[T(v[jj]) for jj in eachindex(v) if jj != ii]
            isempty(bares(t)) || pushfirst!(nv, T(bares(t)))
            nt = QuTerm(deltas(t), nv)
            # include current prefactor here
            _add_sum_term!(commterms,nt,pref*s)
        end
        # only modify prefactor after the exchange
        pref *= pp
    end
    filter!(!isempty, v)
    sort!(v)
    pref
end

# iterate over a single term t with prefactor s and a sum as if they were in the same sum
termsumiter(t::QuTerm,s,sum::QuExpr) = Base.Iterators.flatten((((t,s),), sum.terms))
termsumiter(t,s,sum::QuExpr) = termsumiter(QuTerm(t),s,sum)

# empty collector that is used (and never modified) when no new terms can be generated
const _NO_NEW_TERMS = QuExpr()

function _add_with_normal_order!(A::QuExpr,t::QuTerm,s,shortcut_vacA_zero=false)
    # no need to do anything since everything is in normal order
    # in particular, no copy necessary!
    if is_normal_form(t)
       _add_sum_term!(A,t,s)
       return
    end

    W = _normalize_without_commutation(t)
    W === nothing && return
    # parts with less than two operators cannot produce new terms, so do not create collectors for them
    newbareterms = length(W.bares) < 2 ? _NO_NEW_TERMS : QuExpr()
    newexpvterms = all(E -> length(E) < 2, W.expvals) ? _NO_NEW_TERMS : QuExpr()
    newcorrterms = all(C -> length(C) < 2, W.corrs) ? _NO_NEW_TERMS : QuExpr()
    prefbare = normal_order!(W.bares,   newbareterms, shortcut_vacA_zero)
    prefexpv = normal_order!(W.expvals, ExpVal, newexpvterms)
    prefcorr = normal_order!(W.corrs,   Corr,   newcorrterms)

    # normal_form(t) = t.prefs * (prefexpv*tN.ev + ev_new) * (prefcorr*tN.co + co_new) * (prefbare*tN.bare + bare_new)
    # W has been transformed to tN in-place
    t = QuTerm(W)
    pref = prefbare*prefexpv*prefcorr
    if !iszero(pref)
        if iszero(nsuminds(t))
            # should be guaranteed to be in normal order
            _add_sum_term!(A,t,s*pref)
        else
            # might need to reorder sum indices, which can "break" our normal order
            tn = _reorder_suminds(t)
            # this will check again if it is in normal form to make sure
            _add_with_normal_order!(A,tn,s*pref,shortcut_vacA_zero)
        end
    end

    isempty(newexpvterms) && isempty(newcorrterms) && isempty(newbareterms) && return

    firstterm = true
    for (tev,sev) in termsumiter(expvals(t), prefexpv, newexpvterms)
        #@assert iszero(nsuminds(tev)) && isempty(params(tev)) && isempty(corrs(tev)) && isempty(bares(tev))
        for (tco,sco) in termsumiter(corrs(t), prefcorr, newcorrterms)
            #@assert iszero(nsuminds(tco)) && isempty(params(tco)) && isempty(expvals(tco)) && isempty(bares(tco))
            for (tba,sba) in termsumiter(bares(t), prefbare, newbareterms)
                #@assert iszero(nsuminds(tba)) && isempty(params(tba)) && isempty(expvals(tba)) && isempty(corrs(tba))

                if firstterm
                    # the first term has been added above in _add_sum_term!
                    firstterm = false
                    continue
                end

                pref = sev*sco*sba
                iszero(pref) && continue
                nt = _encode_term(nsuminds(t), (deltas(t),deltas(tev),deltas(tco),deltas(tba)), (params(t),), (expvals(tev),), (corrs(tco),), (bares(tba),), nothing, nothing)
                _add_with_normal_order!(A,nt,s*pref,shortcut_vacA_zero)
            end
        end
    end
end

"""
    normal_form(A::QuExpr; shortcut_vacA_zero=false)

Return an equivalent expression in normal-ordered form (all commutations performed).

`shortcut_vacA_zero` can be set to `true` when the expression will be multiplied
by a vacuum state from the right, allowing early elimination of annihilation terms.
"""
function normal_form(A::QuExpr,shortcut_vacA_zero=false)
    An = QuExpr()
    for (t,s) in A.terms
        _add_with_normal_order!(An,t,s,shortcut_vacA_zero)
    end
    An
end
normal_form(A::Number) = QuExpr(A)

function *(A::QuTerm,B::QuTerm)
    nA, nB = nsuminds(A), nsuminds(B)
    if nA > 0 && nB > 0
        if nA > nB
            A = shift_sumind(nB)(A)
        else
            B = shift_sumind(nA)(B)
        end
    end
    parts = ((deltas(A),deltas(B)), (params(A),params(B)), (expvals(A),expvals(B)), (corrs(A),corrs(B)), (bares(A),bares(B)))
    if nA + nB > 0
        # renumber the sum indices in the order in which they appear
        f = reorder_suminds()
        _encode_term(nA+nB, parts..., f, f)
    else
        _encode_term(0, parts..., nothing, nothing)
    end
end

*(A::Union{QuExpr,QuTerm}) = A
function *(A::QuExpr,B::QuExpr)
    QuExpr((tA*tB,sA*sB) for ((tA,sA),(tB,sB)) in Iterators.product(A.terms,B.terms))
end
*(A::Number,B::QuExpr) = QuExpr((tB,A*sB) for (tB,sB) in B.terms)
*(B::QuExpr,A::Number) = A*B
Base.:/(A::QuExpr,B::Number) = inv(B)*A

+(A::QuExpr) = A
function +(A::QuExpr,B::QuExpr)
    S = copy(A)
    for (t,s) in B.terms
        _add_with_auto_order!(S,t,s)
    end
    S
end
function +(A::QuExpr,B::Number)
    S = copy(A)
    _add_with_auto_order!(S,QuTerm(),B)
    S
end
+(B::Number,A::QuExpr) = A+B

-(A::QuExpr) = -1 * A
function -(A::QuExpr,B::QuExpr)
    S = copy(A)
    for (t,s) in B.terms
        _add_with_auto_order!(S,t,-s)
    end
    S
end
-(B::Number,A::QuExpr) = B + (-A)
-(A::QuExpr,B::Number) = A + (-B)
