using Base: tail
using Combinatorics

export ascorr

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
    ascorr(expr::Operator)

Take an expression `expr=A B C + D E...` and write its expectation value in
terms of single-body expectation values ``⟨A⟩, ⟨B⟩, \\ldots``, and many-body
correlations ``⟨AB⟩_c, ⟨ABC⟩_c``, etc.

E.g., `ascorr(adag(:n)*a(:n))` returns ``⟨a^\\dagger_n a_n⟩_c + ⟨a^\\dagger_n⟩
⟨a_n⟩`` (which is equal to ``⟨a^\\dagger_n a_n⟩``), while
`ascorr(adag(:n)*a(:m)*a(:n))` returns ``\\langle a_{n}^\\dagger a_{m} a_{n}
\\rangle_{c} + \\langle a_{n}^\\dagger \\rangle \\langle a_{m} \\rangle \\langle
a_{n} \\rangle + \\langle a_{n}^\\dagger \\rangle \\langle a_{m} a_{n}
\\rangle_{c} + \\langle a_{m} \\rangle \\langle a_{n}^\\dagger a_{n}
\\rangle_{c} + \\langle a_{n} \\rangle \\langle a_{n}^\\dagger a_{m}
\\rangle_{c}``.

See also: [`ExpVal`](@ref), [`Corr`](@ref)"""
function ascorr(A::OpSum)
    newA = OpSum()
    for (t,s) in A.terms
        lt = length(t.bares)
        if lt==0
            _add_sum_term!(newA,t,s)
        else
            for C in prodcorr_inds(lt)
                _add_sum_term!(newA,term2corr(t,C),s)
            end
        end
    end
    newA
end

##################################################################################################################################
#### MISSING: ####################################################################################################################
##################################################################################################################################
#     # first calculate the correlation for the term in the sum with the "bare" indices, which means that the sum index
#     # is assumed to be distinct from the indices of the expressions
#     # then, calculate the correction expression for the term where the sum index is the same as any of the non-summed indices,
#     # and add the correction between that and the "incorrect" term where index identity is not taken into account
#     tmp = ascorr(A.A)
#     res = OpSumAnalytic(A.ind,tmp)
#     for ii in setdiff(indexset(A),sumindexset(A))
#         res += ascorr(replace_index(A.A,A.ind,ii)) - replace_index(tmp,A.ind,ii)
#     end
#     res
# end

"""
    CorrExpTup_isless(a,b)

isless for Tuples of integers that represent ExpVals and Corrs of sorted Operators (with n representing An such that n<m == An<Am).
(n,) ⧋ ExpVal(An)
(n,m,...) ⧋ Corr(An*Am*...)
Defined in such a way that the same order is obtained as with `Operator` objects
"""
CorrExpTup_isless(a::Tuple{Int},b::Tuple{Int}) = a[1] < b[1]
CorrExpTup_isless(a::Tuple{Int},b::NTuple{N,Int}) where {N} = true  # ExpVal < Corr
CorrExpTup_isless(a::NTuple{N,Int},b::Tuple{Int}) where {N} = false # Corr > ExpVal
CorrExpTup_isless(a::NTuple{N,Int},b::NTuple{M,Int}) where {N,M} = N==M ? a<b : N<M

"""
    CorrExpPerm_isless(a,b)

isless for Tuples of Tuples representing products of Corr and ExpVal (see above) of a permutation of operators.
E.g., ((1,3),(2,)) represents <A1 A3>_C <A2>.
It is assumed that the total number of operators in a and b is equal, i.e., that `sum(length.(a)) == sum(length.(b))`."""
function CorrExpPerm_isless(a::Tuple,b::Tuple)
    for (ca,cb) in zip(a,b)
        ca == cb || return CorrExpTup_isless(ca,cb)
    end
    # if we reach here, they are equal
    return false
end

"""
    prodcorr_inds(N::Int)

for N operators, create an array of tuples of tuples that represents the terms in a sum of products
of expectation values and correlators. Each tuple corresponds to a sum term, see
[`CorrExpTup_isless`](@ref) and [`CorrExpPerm_isless`](@ref) for details of the format.
The returned array and terms are sorted such that if the N operators are sorted, the represented
expression is also sorted with the conventions of the QuantumAlgebra package.
This allows to directly construct the nested OpSums and OpProds without having to go through simplification.
"""
function prodcorr_inds(N::Int)
    # get return value for N from cache if present, otherwise calculate and cache it
    get!(_PRODCORR_INDS_CACHE,N) do
        terms = Any[Tuple(tuple.(1:N))]
        for n = 2:N-1
            append!(terms,ncomb_inds(n,1:N))
        end
        [(tuple(1:N...),),sort!(terms,lt=CorrExpPerm_isless)...]
    end
end
# preload the cache for N=0 and N=1 so we do not have to special-case those above
const _PRODCORR_INDS_CACHE = Dict{Int,Vector{Any}}(0 => [((),)], 1 => [((1,),)])

function ncomb_inds(n,inds,used_combs=Set())
    terms = []
    for c in combinations(inds,n)
        C = Tuple(c)
        C in used_combs && continue
        push!(used_combs,C)
        notCinds = setdiff(inds,c)
        S = tuple.(notCinds)
        push!(terms,Tuple(sort!([C,S...],lt=CorrExpTup_isless)))
        if n<length(inds)
            # pass a copy of used_combs here so getting, e.g., <jk>*C here does not prevent a term <jk> <a><b> later
            for CC in ncomb_inds(n,notCinds,copy(used_combs))
                push!(terms,Tuple(sort!([C,CC...],lt=CorrExpTup_isless)))
            end
        end
    end
    terms
end
