using Combinatorics

export expval_as_corrs

function term2corr(A::OpTerm,C::Tuple)
    corrs = deepcopy(A.corrs)
    sizehint!(corrs,length(corrs)+length(C))
    for grp in C
        cc = Corr(BaseOpProduct(collect(map(i -> A.bares.v[i], grp))))
        push!(corrs,cc)
    end
    sort!(corrs)
    OpTerm(A.nsuminds,A.δs,A.params,A.expvals,corrs,BaseOpProduct())
end

"""
    expval_as_corrs(expr::Operator)

Take an expression `expr=A B C + D E...` and write its expectation value in
terms of correlations ``⟨A⟩_c, ⟨B⟩_c, ⟨AB⟩_c, ⟨ABC⟩_c, \\ldots``. Note that
``⟨A⟩_c = ⟨A⟩``.

E.g., `expval_as_corrs(adag(:n)*a(:n))` returns ``⟨a^\\dagger_n a_n⟩_c +
⟨a^\\dagger_n⟩_c ⟨a_n⟩_c`` (which is equal to ``⟨a^\\dagger_n a_n⟩``), while
`expval_as_corrs(adag(:n)*a(:m)*a(:n))` returns ``\\langle a_{n}^\\dagger a_{m}
a_{n} \\rangle_{c} + \\langle a_{n}^\\dagger \\rangle_{c} \\langle a_{m}
\\rangle_{c} \\langle a_{n} \\rangle_{c} + \\langle a_{n}^\\dagger \\rangle_{c}
\\langle a_{m} a_{n} \\rangle_{c} + \\langle a_{m} \\rangle_{c} \\langle
a_{n}^\\dagger a_{n} \\rangle_{c} + \\langle a_{n} \\rangle_{c} \\langle
a_{n}^\\dagger a_{m} \\rangle_{c}``.

See also: [`expval`](@ref), [`corr`](@ref)"""
function expval_as_corrs(A::OpSum)
    newA = OpSum()
    for (t,s) in A.terms
        # first calculate the correlation for the term in the sum with the "bare" indices, which means that the sum index
        # is assumed to be distinct from the indices of the expressions
        # then, calculate the correction expression for the term where the sum index is the same as any of the non-summed indices,
        # and add the correction between that and the "incorrect" term where index identity is not taken into account
        _add_corrs!(newA,t,s)

        if t.nsuminds>0 && !isempty(t.bares)
            extinds = unique(extindices(t.bares))
            for ii in 1:t.nsuminds
                sumind = sumindex(ii)
                for ind in extinds
                    f = replace_inds(sumind=>ind,(sumindex(jj)=>sumindex(jj-1) for jj=ii+1:t.nsuminds)...)
                    tb = f(t.bares)
                    tbc = deepcopy(tb)
                    normal_order!(tb,NotASum())
                    if tb != tbc
                        tright = _normalize_without_commutation(OpTerm(t.nsuminds-1,f.(t.δs),f.(t.params),f.(t.expvals),f.(t.corrs),tb))
                        _add_corrs!(newA,tright,s)
                        twrong = OpTerm(tright.nsuminds,tright.δs,tright.params,tright.expvals,tright.corrs,tbc)
                        _add_corrs!(newA,twrong,-s,_add_with_normal_order!)
                    end
                end
            end
        end
    end
    newA
end

function _add_corrs!(A,t,s,addfun! = _add_sum_term!)
    if isempty(t.bares)
        addfun!(A,t,s)
    else
        for C in expval2corrs_inds(length(t.bares))
            addfun!(A,term2corr(t,C),s)
        end
    end
end

"""
    CorrTup_isless(a,b)

isless for Tuples of integers that represent Corr of sorted Operators (with n representing An such that n<m == An<Am).
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

for N operators, create an array of tuples of tuples that represents the terms
in a sum of products of correlators. Each tuple corresponds to a sum term, see
[`CorrTup_isless`](@ref) and [`CorrPerm_isless`](@ref) for details of the
format. The returned array and terms are sorted such that if the N operators are
sorted, the represented expression is also sorted with the conventions of the
QuantumAlgebra package. This allows to directly return a normal-ordered form.
"""
function expval2corrs_inds(N::Int)
    # get return value for N from cache if present, otherwise calculate and cache it
    get!(_EXPVAL2CORRS_CACHE,N) do
        terms = Any[Tuple(tuple.(1:N))]
        for n = 2:N-1
            append!(terms,ncomb_inds(n,1:N))
        end
        [(tuple(1:N...),),sort!(terms,lt=CorrPerm_isless)...]
    end
end
# preload the cache for N=0 and N=1 so we do not have to special-case those above
const _EXPVAL2CORRS_CACHE = Dict{Int,Vector{Any}}(0 => [((),)], 1 => [((1,),)])

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

function term2expvals(A::OpTerm,C::Tuple)
    expvals = deepcopy(A.expvals)
    sizehint!(expvals,length(expvals)+length(C))
    for grp in C
        cc = ExpVal(BaseOpProduct(collect(map(i -> A.bares.v[i], grp))))
        push!(expvals,cc)
    end
    sort!(expvals)
    OpTerm(A.nsuminds,A.δs,A.params,expvals,A.corrs,BaseOpProduct())
end

_canon_op_to_num(A::BaseOperator) = _canon_ind_to_num(only(A.inds))
_expval2tuple(A::ExpVal) = Tuple(_canon_op_to_num.(A.ops.v))
const _CORR2EXPVALS_CACHE = Dict{Int,Vector{Any}}()
function corr2expvals_inds(N)
    get!(_CORR2EXPVALS_CACHE,N) do
        A = OpSum(OpTerm(BaseOpProduct(BosonDestroy.(:a,_canon_ind.(1:N)))))
        # expval_as_corrs(As...) = <As...> in terms of correlators <A1 A2...>c
        # <As...>c = <As...> - <As...> + <As...>c, where RHS expression only contains <As...> and lower-order correlators
        cA = expval(A) - expval_as_corrs(A) + corr(A)
        # rewrite the lower-order correlators in terms of expectation values
        ncA = OpSum()
        for (t,s) in cA.terms
            nt = deepcopy(t)
            empty!(nt.corrs)
            nts = OpSum(nt)
            for co in t.corrs
                nts *= OpSum((term2expvals(OpTerm(co.ops),C)=>sC for (C,sC) in corr2expvals_inds(length(co))))
            end
            for (ntt,ns) in nts.terms
                _add_with_normal_order!(ncA,ntt,ns*s)
            end
        end
        [Tuple(_expval2tuple.(t.expvals))=>s for (t,s) in sort!(collect(ncA.terms),by=x->x[1])]
    end
end
