const _GROUPALIASES = Dict{QuOpName,Vector{QuOpName}}()

add_groupaliases(groupname,names) = (_GROUPALIASES[QuOpName(groupname)] = QuOpName.(collect(names)))

function unalias(A::BaseOperator)::BaseOperator
    if A.name ∈ keys(_GROUPALIASES)
        names = _GROUPALIASES[A.name]
        ind = first(A.inds)
        @assert isintindex(ind)
        BaseOperator(A.t,names[ind.num],Base.tail(A.inds))
    else
        A
    end
end

unalias(A::T) where T<:Union{ExpVal,Corr} = T(BaseOpProduct(unalias.(A.ops.v)))
