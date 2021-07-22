abstract type IndexFunction end

(f::IndexFunction)(A::BaseOperator) = BaseOperator(A.t,A.name,f.(A.inds))
(f::IndexFunction)(d::δ) = δ(f(d.iA),f(d.iB))
(f::IndexFunction)(par::Param) = Param(par.name,par.state,f.(par.inds))
(f::IndexFunction)(ops::BaseOpProduct) = BaseOpProduct(f.(ops.v))
(f::IndexFunction)(ev::T) where T<:Union{ExpVal,Corr} = T(f(ev.ops))
(f::IndexFunction)(A::QuTerm,nsuminds=A.nsuminds) = QuTerm(nsuminds,f.(A.δs),f.(A.params),f.(A.expvals),f.(A.corrs),f(A.bares))
# some index functions have to be reset for each new term in a sum
(f::IndexFunction)(A::QuExpr) = _map_quexpr_ops(x->(adapt_for_new_term!(f); f(x)),A)
# by default, nothing to do
adapt_for_new_term!(f::IndexFunction) = f

function _canon_ind(n)
    i,c = divrem(n-1,15)
    i == 0 ? QuIndex('i'+c) : QuIndex('i'+c,i+1)
end
_canon_ind_to_num(ind::QuIndex) = ind.sym-'i'+1 + (ind.num==typemin(ind.num) ? 0 : (ind.num-1)*15)

mutable struct canon_inds <: IndexFunction
    icurr::IndexInt
    canon_inds() = new(0)
end
(f::canon_inds)(ind::QuIndex) = isnoindex(ind) ? ind : _canon_ind(f.icurr+=one(f.icurr))
adapt_for_new_term!(f::canon_inds) = (f.icurr = 0; f)

mutable struct canon_inds_remember <: IndexFunction
    icurr::IndexInt
    replacements::Dict{QuIndex,QuIndex}
    canon_inds_remember() = new(0,Dict{QuIndex,QuIndex}())
end
function (f::canon_inds_remember)(ind::QuIndex)
    isnoindex(ind) && return ind
    i = _canon_ind(f.icurr+=one(f.icurr))
    f.replacements[i] = ind
    i
end
adapt_for_new_term!(f::canon_inds_remember) = error("canon_inds_remember should not be applied on QuExpr.")

struct shift_sumind <: IndexFunction
    n::IndexInt
end
(f::shift_sumind)(ind::QuIndex) = issumindex(ind) ? sumindex(ind.num+f.n) : ind

struct replace_inds <: IndexFunction
    replacements::Dict{QuIndex,QuIndex}
    replace_inds(args...) = new(Dict{QuIndex,QuIndex}(args...))
end
(f::replace_inds)(ind::QuIndex) = get(f.replacements,ind,ind)

mutable struct reorder_suminds <: IndexFunction
    icurr::IndexInt
    replacements::Dict{QuIndex,QuIndex}
    reorder_suminds() = new(0,Dict{QuIndex,QuIndex}())
end
(f::reorder_suminds)(ind::QuIndex) = issumindex(ind) ? get!(()->sumindex(f.icurr += one(f.icurr)), f.replacements, ind) : ind
adapt_for_new_term!(f::reorder_suminds) = (f.icurr=0; empty!(f.replacements); f)

mutable struct popped_out_suminds <: IndexFunction
    conserved_ind::QuIndex
    icurr::IndexInt
    suminds::Vector{QuIndex}
    popped_out_suminds(n) = new(n,0,QuIndex[])
end
function (f::popped_out_suminds)(ind::QuIndex)
    if isnoindex(ind)
        ind
    elseif ind==f.conserved_ind
        sumindex(1)
    else
        push!(f.suminds,ind)
        _canon_ind(f.icurr+=one(f.icurr))
    end
end
adapt_for_new_term!(f::popped_out_suminds) = error("popped_out_suminds should not be applied on QuExpr.")

hasind(ind::QuIndex,d::δ) = ind in (d.iA,d.iB)
hasind(ind::QuIndex,A::Union{BaseOperator,Param}) = ind in A.inds
hasind(ind::QuIndex,ev::Union{ExpVal,Corr}) = hasind(ind,ev.ops)
hasind(ind::QuIndex,A::BaseOpProduct) = any(hasind.((ind,),A.v))

split_by_ind(ind::QuIndex,v) = (iwiths = hasind.((ind,),v); (v[.!iwiths], v[iwiths]))

indices(d::δ)::Vector{QuIndex} = [d.iA,d.iB]
indices(A::Union{BaseOperator,Param})::Vector{QuIndex} = collect(assignedinds(A.inds)) # collect to guarantee copy
indices(v::AbstractVector)::Vector{QuIndex} = mapreduce(indices, vcat, v; init=QuIndex[])
indices(A::BaseOpProduct)::Vector{QuIndex} = indices(A.v)
indices(A::Union{ExpVal,Corr})::Vector{QuIndex} = indices(A.ops)
indices(A::QuTerm, withδs::Bool=true)::Vector{QuIndex} = begin
    args = withδs ? (A.δs,A.params,A.expvals,A.corrs,A.bares) : (A.params,A.expvals,A.corrs,A.bares)
    mapreduce(indices, vcat, args; init=QuIndex[])
end
indices(A::QuExpr)::Vector{QuIndex} = mapreduce(indices, vcat, sort!(collect(keys(A.terms))); init=QuIndex[])

"`extindices(A)` return externally visible indices of an expression"
extindices(A) = filter(!issumindex,indices(A))

"`symmetric_index_nums(A)` return sequence of numbers of exchange-symmetric indices"
symmetric_index_nums(A) = symmetric_index_nums(QuExpr(A))
function symmetric_index_nums(A::QuExpr)
    inds = extindices(A)
    Nsyms = [1]
    for ii=2:length(inds)
        i1,i2 = inds[ii-1], inds[ii]
        Aexc = normal_form(replace_inds(i1=>i2,i2=>i1)(A))
        if A == Aexc
            Nsyms[end] += 1
        else
            push!(Nsyms,1)
        end
    end
    Nsyms
end
