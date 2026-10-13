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
        haskey(B.terms,tA) || return false
        sA ≈ B.terms[tA] || return false
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

Return the commutator ``[A,B] = AB - BA`` in normal form. It is computed directly from the
commutators of pairs of single operators (with the Leibniz rule), so that the terms of ``AB`` and
``BA`` that cancel are never generated. E.g., `comm(b[:i], b'[:j]*b[:j])` gives `δ(ij) b(i)` (after
`@boson b[:]`).
"""
function comm(A::QuExpr,B::QuExpr)
    An, Bn = _normal_form_if_needed(A), _normal_form_if_needed(B)
    C = QuExpr()
    for (tA,sA) in An.terms
        _withnum(sA) do sA
            for (tB,sB) in Bn.terms
                _withnum(sB -> _add_comm!(C,tA,tB,sA*sB), sB)
            end
        end
    end
    C
end
comm(A::QuExpr,B::Number) = QuExpr()
comm(A::Number,B::QuExpr) = QuExpr()
comm(A,B) = A*B - B*A

Base.adjoint(A::BaseOperator) = op_adjoint(A)
Base.adjoint(A::OpProduct) = BaseOpProduct(reverse!(BaseOperator[adjoint(O) for O in A]))
Base.adjoint(A::T) where {T<:Union{ExpVal,Corr}} = T(adjoint(ops(A)))
Base.adjoint(A::Param) = parstate(A)=='r' ? A : Param(opname(A),parstate(A)=='n' ? 'c' : 'n',inds(A))
function Base.adjoint(A::QuTerm)
    B = QuTerm(nsuminds(A),deltas(A),adjoint.(params(A)),adjoint.(expvals(A)),adjoint.(corrs(A)),adjoint(bares(A)))
    nsuminds(B)>0 ? _reorder_suminds(B) : B
end
Base.adjoint(A::QuExpr) = QuExpr((adjoint(t),adjoint(s)) for (t,s) in A.terms)

"""
    is_normal_form(A)

Return whether `A` (a `QuExpr`, a term, or a product of operators) is in normal form, i.e., whether
`normal_form(A) == A` (see [`normal_form`](@ref)).
"""
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

is_normal_form(t::QuTerm) = _is_normal_form_except_suminds(t) && _suminds_canonical(t)

# whether the term is in normal form, apart from the numbering of the sum indices
function _is_normal_form_except_suminds(t::QuTerm)
    (issorted(params(t)) && issorted(expvals(t)) && issorted(corrs(t))) || return false
    is_normal_form(bares(t)) || return false
    for EV in expvals(t)
        is_normal_form(EV) || return false
    end
    for CO in corrs(t)
        is_normal_form(CO) || return false
    end

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
# i.e., that commutes with all operators in between. If lookahead_meet(A[k]), A[k] can also be moved to
# the right past operators it commutes with, so that the two meet in between. Returns the positions of
# A[kp] and of the last operator before the meeting point (i.e., the position of A[k] if A[kp] is moved
# next to it), given the position of A[k], or nothing if there is no such operator.
function _contraction_partner(A,pk)
    Ak = _opat(A,pk)
    lookahead_contraction(Ak) || return nothing
    fam = opfamily(Ak)
    pend = _endpos(A)
    pk1 = _nextpos(A,pk)
    pk1 < pend || return nothing
    meet = lookahead_meet(Ak)
    # if they can meet, A[k] can be moved to the right past all following operators it commutes with (up to pr)
    pr = pk
    if meet
        while (pn = _nextpos(A,pr)) < pend && commutes(Ak,_opat(A,pn))
            pr = pn
        end
    end
    pprev, pkp = pk1, _nextpos(A,pk1)
    while pkp < pend
        Akp = _opat(A,pkp)
        if meet
            # meeting point: as far right as A[k] can be moved (but before A[kp])
            pm = pr < pkp ? pr : pprev
            if opfamily(Akp) == fam && _commutes_with_all(A,_nextpos(A,pm),pkp,Akp) && contract!(_NoEmit(),Ak,Akp) !== nothing
                return (pkp, pm)
            end
        else
            (opfamily(Akp) == fam && lookahead_contraction(Akp)) || break
            commutes(_opat(A,pprev),Akp) || break
            # A[kp] has to commute with all operators it is moved past, not just its neighbor
            if _commutes_with_all(A,pk1,pprev,Akp) && contract!(_NoEmit(),Ak,Akp) !== nothing
                return (pkp, pk)
            end
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
# (with i1 < i2) in the operator product A by ops, placed after A[ins] (with i1 ≤ ins < i2, where
# ins == i1 means at the position of A[i1]), with prefactor c*e.prefactor and the additional
# Kronecker deltas δs, to e.collector
struct _ReplaceEmitter{T}
    A::Vector{BaseOperator}
    i1::Int
    i2::Int
    ins::Int
    collector::QuExpr
    prefactor::T
end
function (e::_ReplaceEmitter)(c, δs::Vector{δ}, newops::BaseOperator...)
    iszero(c) && return nothing
    A, i1, i2, ins = e.A, e.i1, e.i2, e.ins
    v = BaseOperator[]
    sizehint!(v, length(A) - 2 + length(newops))
    append!(v, view(A, 1:i1-1))
    append!(v, view(A, i1+1:ins))
    append!(v, newops)
    append!(v, view(A, ins+1:i2-1))
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
            x = exchange!(_ReplaceEmitter(A,j-1,j,j-1,term_collector,prefactor),A[j-1],A[j])
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
        kp, ins = k+1, k
        res = contract!(_ReplaceEmitter(A,k,kp,ins,term_collector,prefactor),A[k],A[kp])
        if res === nothing
            partner = _contraction_partner(A,k)
            if partner !== nothing
                kp, ins = partner
                res = contract!(_ReplaceEmitter(A,k,kp,ins,term_collector,prefactor),A[k],A[kp])
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
        # replace A[k] A[kp] by op (where op === nothing is the identity), placed after the original A[ins]
        if op === nothing
            deleteat!(A,(k,kp))
        elseif ins == k
            A[k] = op
            deleteat!(A,kp)
        else
            deleteat!(A,(k,kp))
            insert!(A,ins,op)
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

# (inlined into the callers, and calls the implementation with the type of s known, see `_withnum`)
@inline _add_with_normal_order!(A::QuExpr,t::QuTerm,s,shortcut_vacA_zero=false) =
    _withnum(s -> _add_with_normal_order_impl!(A,t,s,shortcut_vacA_zero), s)
function _add_with_normal_order_impl!(A::QuExpr,t::QuTerm,s,shortcut_vacA_zero)
    if _is_normal_form_except_suminds(t)
        if nsuminds(t) ≤ 1
            # no need to do anything since everything is in normal order
            # in particular, no copy necessary!
            if _suminds_canonical(t)
                _add_sum_term!(A,t,s)
                return
            end
        else
            # only the sum indices might have to be renumbered
            res = _canonicalize_suminds(t)
            res === nothing && return
            tc, sgn, isnormal = res
            if isnormal
                _add_sum_term!(A,tc,s*sgn)
                return
            end
            t, s = tc, s*sgn
        end
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
            # number the sum indices canonically (which can change the order of operators)
            res = _canonicalize_suminds(t)
            if res !== nothing
                tc, sgn, isnormal = res
                # if tc is not guaranteed to be in normal form, this will check again
                isnormal ? _add_sum_term!(A,tc,s*pref*sgn) : _add_with_normal_order!(A,tc,s*pref*sgn,shortcut_vacA_zero)
            end
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
    normal_form(A::QuExpr, shortcut_vacA_zero=false)

Return the normal form of `A`, an equal expression in which the operators of each term are normal
ordered (creation operators to the left of annihilation operators, using the algebra of the
operators, e.g., ``a a^† = 1 + a^† a``), products of operators that simplify are simplified, and the
terms are in a canonical form (e.g., summation indices are numbered canonically), so that equal
expressions have the same normal form. E.g., `normal_form(a*a')` gives `1 + a†() a()`.

If the vacuum ``⟨0|`` will be applied from the left (see [`vacA`](@ref)), `shortcut_vacA_zero=true`
drops the terms that it annihilates during normal ordering already.

See also [`auto_normal_form`](@ref), [`is_normal_form`](@ref).
"""
function normal_form(A::QuExpr,shortcut_vacA_zero=false)
    An = QuExpr()
    for (t,s) in A.terms
        _add_with_normal_order!(An,t,s,shortcut_vacA_zero)
    end
    An
end
normal_form(A::Number) = QuExpr(A)

_normal_form_if_needed(A::QuExpr) = all(is_normal_form, keys(A.terms)) ? A : normal_form(A)

## commutators

# The commutator of operator products X = A₁⋯Aₙ and Y = B₁⋯Bₘ is computed with the graded
# Leibniz rule
#   [X,Y}_ε = ∑_{i,j} ε(A_{>i},Y) ε(A_i,B_{<j}) A_{<i} B_{<j} [A_i,B_j}_ε B_{>j} A_{>i},
# where [P,Q}_ε = PQ - ε(P,Q) QP, and ε(P,Q) = ±1 is the sign obtained when exchanging P and Q
# (-1 for an odd number of exchanges of fermions with the same name). [A_i,B_j}_ε is thus the
# anticommutator for such fermions, and the commutator otherwise. If ε(X,Y) = 1, [X,Y] = [X,Y}_ε,
# and only the pairs A_i,B_j that do not (anti)commute contribute. Otherwise, XY - YX is computed.

# ε(A,B) for single operators: -1 for fermions with the same name (which anticommute)
_exchange_sign(A::BaseOperator,B::BaseOperator) =
    (opfamily(A) == FermionFamily && opfamily(B) == FermionFamily && opname(A) == opname(B)) ? -1 : 1
# whether [A,B}_ε = 0 (for fermions with the same name, whether they anticommute)
_graded_commutes(A::BaseOperator,B::BaseOperator,ε::Int) =
    ε == 1 ? commutes(A,B) : (optype(A) == optype(B) || !maybe_same_site(A,B))

# add s [tA,tB] (in normal form) to C, where tA and tB are terms in normal form
function _add_comm!(C::QuExpr,tA::QuTerm,tB::QuTerm,s)
    X, Y = bares(tA), bares(tB)
    (isempty(X) || isempty(Y)) && return C
    # check whether any pair contributes, and get the total sign ε(X,Y)
    # (sum indices are not renumbered yet, but this does not change whether operators (anti)commute)
    εXY = 1
    contributes = false
    for A in X, B in Y
        ε = _exchange_sign(A,B)
        εXY *= ε
        contributes || (contributes = !_graded_commutes(A,B,ε))
    end
    if εXY == -1
        # (here, [X,Y] = [X,Y}_ε - 2YX does not vanish even if all pairs anticommute)
        _add_with_normal_order!(C,tA*tB,s)
        _add_with_normal_order!(C,tB*tA,-s)
        return C
    end
    contributes || return C
    # the product of the terms combines the c-number parts and numbers the sum indices consistently
    t = tA*tB
    ops = collect(bares(t))
    n, m = length(X), length(Y)
    # ε(A_{≥k},Y) for k = 1…n+1
    εrest = ones(Int,n+1)
    for k = n:-1:1
        ε = 1
        for j = n+1:n+m
            ε *= _exchange_sign(ops[k],ops[j])
        end
        εrest[k] = ε*εrest[k+1]
    end
    for i = 1:n
        sgn = εrest[i+1]
        for j = n+1:n+m
            ε = _exchange_sign(ops[i],ops[j])
            _graded_commutes(ops[i],ops[j],ε) || _add_pair_comm!(C,t,ops,n,i,j,ε,s*sgn)
            sgn *= ε
        end
    end
    C
end

# add s A_{<i} B_{<j} [A_i,B_j}_ε B_{>j} A_{>i} to C, where ops = [A₁,…,Aₙ,B₁,…,Bₘ] are the operators of t
# (and j is the position of B_j in ops)
function _add_pair_comm!(C::QuExpr,t::QuTerm,ops::Vector{BaseOperator},n,i,j,ε,s)
    A, B = ops[i], ops[j]
    if A < B
        # B A = x A B + R, so [A,B}_ε = (1 - εx) A B - ε R
        x = exchange!(_CommEmitter(C,t,ops,n,i,j,-ε*s),B,A)
        _CommEmitter(C,t,ops,n,i,j,s)(1-ε*x,_NO_DELTAS,A,B)
    else
        # A B = x B A + R, so [A,B}_ε = (x - ε) B A + R
        emit = _CommEmitter(C,t,ops,n,i,j,s)
        x = exchange!(emit,A,B)
        emit(x-ε,_NO_DELTAS,B,A)
    end
    nothing
end

# emitter for the terms of [A_i,B_j}_ε: e(c, δs, newops...) adds c*e.s A_{<i} B_{<j} newops B_{>j} A_{>i}
# with the additional deltas δs (and the c-number parts of t) to e.C
struct _CommEmitter{T}
    C::QuExpr
    t::QuTerm
    ops::Vector{BaseOperator}
    n::Int
    i::Int
    j::Int
    s::T
end
function (e::_CommEmitter)(c, δs::Vector{δ}, newops::BaseOperator...)
    iszero(c) && return nothing
    t, ops, n, i, j = e.t, e.ops, e.n, e.i, e.j
    newbares = (view(ops,1:i-1), view(ops,n+1:j-1), newops, view(ops,j+1:length(ops)), view(ops,i+1:n))
    nt = _encode_term(nsuminds(t), (deltas(t),δs), (params(t),), (expvals(t),), (corrs(t),), newbares, nothing, nothing)
    _add_with_normal_order!(e.C,nt,c*e.s)
    nothing
end

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
    C = QuExpr()
    for (tB,sB) in B.terms
        _withnum(sB) do sB
            for (tA,sA) in A.terms
                _withnum(sA -> _add_with_auto_order!(C,tA*tB,simplify_number(sA*sB)), sA)
            end
        end
    end
    C
end
function *(A::Number,B::QuExpr)
    C = QuExpr()
    for (tB,sB) in B.terms
        _withnum(sB -> _add_with_auto_order!(C,tB,simplify_number(A*sB)), sB)
    end
    C
end
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
