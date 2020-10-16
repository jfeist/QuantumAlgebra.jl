abstract type IndexFunction end

(f::IndexFunction)(A::BaseOperator) = BaseOperator(A.t,A.a,A.name,f.(A.inds))
(f::IndexFunction)(d::δ) = δ(f(d.iA),f(d.iB))
(f::IndexFunction)(par::Param) = Param(par.name,par.state,f.(par.inds))
(f::IndexFunction)(ops::BaseOpProduct) = BaseOpProduct(f.(ops.v))
(f::IndexFunction)(ev::ExpVal) = ExpVal(f(ev.ops))
(f::IndexFunction)(ev::Corr) = Corr(f(ev.ops))
(f::IndexFunction)(A::OpTerm,nsuminds=A.nsuminds) = OpTerm(nsuminds,f.(A.δs),f.(A.params),f.(A.expvals),f.(A.corrs),f(A.bares))
(f::IndexFunction)(A::OpSum) = OpSum((f(t),s) for (t,s) in A.terms)

mutable struct canon_inds <: IndexFunction
    icurr::IndexInt
    canon_inds() = new(0)
end
(f::canon_inds)(ind::OpIndex) = isnoindex(ind) ? ind : OpIndex('c',f.icurr+=one(f.icurr))

mutable struct canon_inds_remember <: IndexFunction
    icurr::IndexInt
    replacements::Dict{OpIndex,OpIndex}
    canon_inds_remember() = new(0,Dict{OpIndex,OpIndex}())
end
function (f::canon_inds_remember)(ind::OpIndex)
    isnoindex(ind) && return ind
    i = OpIndex('c',f.icurr+=one(f.icurr))
    f.replacements[i] = ind
    i
end

struct shift_sumind <: IndexFunction
    n::IndexInt
end
(f::shift_sumind)(ind::OpIndex) = issumindex(ind) ? sumindex(ind.num+f.n) : ind

struct replace_inds <: IndexFunction
    replacements::Dict{OpIndex,OpIndex}
end
(f::replace_inds)(ind::OpIndex) = get(f.replacements,ind,ind)

mutable struct reorder_suminds <: IndexFunction
    icurr::IndexInt
    replacements::Dict{OpIndex,OpIndex}
    reorder_suminds() = new(0,Dict{OpIndex,OpIndex}())
end
(f::reorder_suminds)(ind::OpIndex) = issumindex(ind) ? get!(()->sumindex(f.icurr += one(f.icurr)), f.replacements, ind) : ind

mutable struct popped_out_suminds <: IndexFunction
    conserved_ind::OpIndex
    icurr::IndexInt
    suminds::Vector{OpIndex}
    popped_out_suminds(n) = new(n,0,OpIndex[])
end
function (f::popped_out_suminds)(ind::OpIndex)
    if isnoindex(ind)
        ind
    elseif ind==f.conserved_ind
        sumindex(1)
    else
        i = OpIndex('c',f.icurr+=one(f.icurr))
        push!(f.suminds,ind)
        i
    end
end

hasind(ind::OpIndex,d::δ) = ind in (d.iA,d.iB)
hasind(ind::OpIndex,A::Union{BaseOperator,Param}) = ind in A.inds
hasind(ind::OpIndex,ev::Union{ExpVal,Corr}) = hasind(ind,ev.ops)
hasind(ind::OpIndex,A::BaseOpProduct) = any(hasind.((ind,),A.v))

split_by_ind(ind::OpIndex,v) = (iwiths = hasind.((ind,),v); (v[.!iwiths], v[iwiths]))