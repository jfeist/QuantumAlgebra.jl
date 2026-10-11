const _GROUPALIASES = Dict{QuOpName,Vector{QuOpName}}()

add_groupaliases(groupname,names) = (_GROUPALIASES[QuOpName(groupname)] = QuOpName.(collect(names)))

function unalias(A::BaseOperator)::BaseOperator
    if opname(A) ∈ keys(_GROUPALIASES)
        names = _GROUPALIASES[opname(A)]
        ind = first(inds(A))
        @assert isintindex(ind)
        BaseOperator(optype(A),names[ind.num],Base.tail(inds(A)),opextra(A))
    else
        A
    end
end

unalias(A::T) where T<:Union{ExpVal,Corr} = T(BaseOpProduct(BaseOperator[unalias(O) for O in ops(A)]))
