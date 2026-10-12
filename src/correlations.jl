using Combinatorics

export expval_as_corrs, corr_as_expvals

"`OpPartition`: a partition of operators 1:N into groups, each given by the (sorted) indices of the operators in it"
const OpPartition = Vector{Vector{Int}}
_partition(C) = OpPartition([collect(Int,grp) for grp in C])

# the existing correlators (T=Corr) or expectation values (T=ExpVal) of A, together with new ones
# made from groups of the operators in barevec (as given by C), as a sorted vector of operator products
function _merge_groups(old, barevec::Vector{BaseOperator}, C::OpPartition)
    groups = OpSeq[ops(x) for x in old]
    sizehint!(groups, length(groups)+length(C))
    # store the new groups in a common buffer
    buf = _words(sum(grp -> sum(i -> _nwords(barevec[i]), grp; init=0), C; init=0))
    pos = 1
    for grp in C
        start = pos
        for i in grp
            pos = _write!(buf, pos, barevec[i], nothing)
        end
        push!(groups, OpSeq(buf, start, length(grp), pos - start))
    end
    sort!(groups)
end

# replace the bare operators `barevec` of A by a product of correlators of groups of them, as given by C
term2corr(A::QuTerm,barevec::Vector{BaseOperator},C::OpPartition) =
    QuTerm(nsuminds(A),deltas(A),params(A),expvals(A),_merge_groups(corrs(A),barevec,C),BaseOpProduct())

"""
    expval_as_corrs(expr::QuExpr)

Take an expression `expr=A B C + D E...` and write its expectation value in
terms of correlations ``⟨A⟩_c, ⟨B⟩_c, ⟨AB⟩_c, ⟨ABC⟩_c, \\ldots``. Note that
``⟨A⟩_c = ⟨A⟩``.

E.g., `expval_as_corrs(a'(:n)*a(:n))` returns ``⟨a^\\dagger_n a_n⟩_c +
⟨a^\\dagger_n⟩_c ⟨a_n⟩_c`` (which is equal to ``⟨a^\\dagger_n a_n⟩``), while
`expval_as_corrs(a'(:n)*a(:m)*a(:n))` returns ``\\langle a_{n}^\\dagger a_{m}
a_{n} \\rangle_{c} + \\langle a_{n}^\\dagger \\rangle_{c} \\langle a_{m}
\\rangle_{c} \\langle a_{n} \\rangle_{c} + \\langle a_{n}^\\dagger \\rangle_{c}
\\langle a_{m} a_{n} \\rangle_{c} + \\langle a_{m} \\rangle_{c} \\langle
a_{n}^\\dagger a_{n} \\rangle_{c} + \\langle a_{n} \\rangle_{c} \\langle
a_{n}^\\dagger a_{m} \\rangle_{c}``.

See also: [`expval`](@ref), [`corr`](@ref)"""
function expval_as_corrs(A::QuExpr)
    _expand_bare_products(_add_corrs!, A)
end
expval_as_corrs(A::Number) = expval_as_corrs(QuExpr(A))

# Rewrite (the expectation value or correlator of) the bare operators in each term of A as a sum of
# products of correlators (add! = _add_corrs!) or of expectation values (add! = _add_expvs!).
function _expand_bare_products(add!, A::QuExpr)
    newA = QuExpr()
    for (t,s) in A.terms
        # first calculate the expansion for the term with the "bare" indices, which means that the sum indices
        # are assumed to be distinct from the other indices of the term
        add!(newA,t,s)
        nsuminds(t)>0 && !isempty(bares(t)) || continue
        # then, for terms where a sum index is equal to one of the non-summed indices, the operator product can
        # change (by reordering or simplification). in that case, subtract the "wrong" expansion where index
        # identity is not taken into account and add the correct one
        extinds = unique(extindices(bares(t)))
        for ii in 1:nsuminds(t)
            sumind = sumindex(ii)
            for ind in extinds
                f = replace_inds(sumind=>ind,(sumindex(jj)=>sumindex(jj-1) for jj=ii+1:nsuminds(t))...)
                tb = f(bares(t))
                is_normal_form(tb) && continue
                W = _normalize_without_commutation(QuTerm(nsuminds(t)-1,f.(deltas(t)),f.(params(t)),f.(expvals(t)),f.(corrs(t)),tb))
                # if the δs make the term vanish, both the "right" and "wrong" terms are zero
                W === nothing && continue
                twrong = QuTerm(W)
                # the correct expansion is that of the normal-ordered term (which can consist of several terms)
                right = QuExpr()
                _add_with_normal_order!(right,twrong,s)
                for (tr,sr) in right.terms
                    add!(newA,tr,sr)
                end
                add!(newA,twrong,-s,_add_with_normal_order!)
            end
        end
    end
    newA
end

function _add_corrs!(A,t,s,addfun! = _add_sum_term!)
    if isempty(bares(t))
        addfun!(A,t,s)
    else
        barevec = collect(bares(t))
        for C in expval2corrs_inds(length(barevec))
            addfun!(A,term2corr(t,barevec,C),s)
        end
    end
end

function _add_expvs!(A,t,s,addfun! = _add_sum_term!)
    if isempty(bares(t))
        addfun!(A,t,s)
    else
        barevec = collect(bares(t))
        for (C,sC) in corr2expvals_inds(length(barevec))
            addfun!(A,term2expvals(t,barevec,C),s*sC)
        end
    end
end

"""
    CorrTup_isless(a,b)

isless for Tuples of integers that represent Corr of sorted operators (with n representing An such that n<m == An<Am).
(n,m,...) ≡ Corr(An*Am*...). Defined in such a way that the same order is obtained as with `BaseOperator` objects
"""
CorrTup_isless(a::NTuple{N,Int},b::NTuple{M,Int}) where {N,M} = N==M ? a<b : N<M

"""
    CorrPerm_isless(a,b)

isless for Tuples of Tuples representing products of Corr (see above) of a permutation of operators.
E.g., ((1,3),(2,)) represents <A1 A3>_C <A2>_C.
It is assumed that the total number of operators in a and b is equal, i.e., that `sum(length.(a)) == sum(length.(b))`."""
function CorrPerm_isless(a::Tuple,b::Tuple)
    for (ca,cb) in zip(a,b)
        ca == cb || return CorrTup_isless(ca,cb)
    end
    # if we reach here, they are equal
    return false
end

"""
    expval2corrs_inds(N::Int)

for N operators, create an array of partitions (see [`OpPartition`](@ref)) that
represents the terms in a sum of products of correlators. Each partition
corresponds to a sum term. Internally, they are computed as tuples of tuples,
see [`CorrTup_isless`](@ref) and [`CorrPerm_isless`](@ref) for details. The
returned array and terms are sorted such that if the N operators are sorted,
the represented expression is also sorted with the conventions of the
QuantumAlgebra package. This allows to directly return a normal-ordered form.
"""
function expval2corrs_inds(N::Int)
    # get return value for N from cache if present, otherwise calculate and cache it
    get!(_EXPVAL2CORRS_CACHE,N) do
        terms = Any[Tuple(tuple.(1:N))]
        for n = 2:N-1
            append!(terms,ncomb_inds(n,1:N))
        end
        _partition.(Any[(tuple(1:N...),),sort!(terms,lt=CorrPerm_isless)...])
    end
end
# preload the cache for N=0 and N=1 so we do not have to special-case those above
const _EXPVAL2CORRS_CACHE = _COWDict{Int,Vector{OpPartition}}(Dict(0 => [_partition(((),))], 1 => [_partition(((1,),))]))

function ncomb_inds(n,inds,used_combs=Set())
    terms = []
    for c in combinations(inds,n)
        C = Tuple(c)
        C in used_combs && continue
        push!(used_combs,C)
        notCinds = setdiff(inds,c)
        S = tuple.(notCinds)
        push!(terms,Tuple(sort!([C,S...],lt=CorrTup_isless)))
        if n<length(inds)
            # pass a copy of used_combs here so getting, e.g., <jk>*C here does not prevent a term <jk> <a><b> later
            for CC in ncomb_inds(n,notCinds,copy(used_combs))
                push!(terms,Tuple(sort!([C,CC...],lt=CorrTup_isless)))
            end
        end
    end
    terms
end

# replace the bare operators `barevec` of A by a product of expectation values of groups of them, as given by C
term2expvals(A::QuTerm,barevec::Vector{BaseOperator},C::OpPartition) =
    QuTerm(nsuminds(A),deltas(A),params(A),_merge_groups(expvals(A),barevec,C),corrs(A),BaseOpProduct())

_canon_op_to_num(A::BaseOperator) = _canon_ind_to_num(only(assignedinds(inds(A))))
_expval2tuple(A::ExpVal) = Tuple(_canon_op_to_num(O) for O in ops(A))
const _CORR2EXPVALS_CACHE = _COWDict{Int,Vector{Pair{OpPartition,Number}}}()
function corr2expvals_inds(N)
    get!(_CORR2EXPVALS_CACHE,N) do
        A = QuExpr(QuTerm(BaseOpProduct(BosonDestroy.(:a,_canon_ind.(1:N)))))
        # expval_as_corrs(As...) = <As...> in terms of correlators <A1 A2...>c
        # <As...>c = <As...> - <As...> + <As...>c, where RHS expression only contains <As...> and lower-order correlators
        cA = expval(A) - expval_as_corrs(A) + corr(A)
        # rewrite the lower-order correlators in terms of expectation values
        ncA = QuExpr()
        for (t,s) in cA.terms
            nts = QuExpr(QuTerm(nsuminds(t),deltas(t),params(t),expvals(t),Corr[],bares(t)))
            for co in corrs(t)
                cot = QuTerm(ops(co))
                cov = collect(bares(cot))
                nts *= QuExpr((term2expvals(cot,cov,C)=>sC for (C,sC) in corr2expvals_inds(length(co))))
            end
            for (ntt,ns) in nts.terms
                _add_with_normal_order!(ncA,ntt,ns*s)
            end
        end
        Pair{OpPartition,Number}[_partition(_expval2tuple.(expvals(t)))=>s for (t,s) in sort!(collect(ncA.terms),by=x->x[1])]
    end
end

"""
    corr_as_expvals(A::QuExpr)

Rewrite connected correlators in `A` into ordinary expectation values.
"""
function corr_as_expvals(A::QuExpr)
    _expand_bare_products(_add_expvs!, A)
end

corr_as_expvals(A::Number) = corr_as_expvals(QuExpr(A))
