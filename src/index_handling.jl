abstract type IndexFunction end

(f::IndexFunction)(A::BaseOperator) = BaseOperator(A.t,A.name,f.(A.inds))
(f::IndexFunction)(d::δ) = δ(f(d.iA),f(d.iB))
(f::IndexFunction)(par::Param) = Param(par.name,par.state,f.(par.inds))
(f::IndexFunction)(ops::BaseOpProduct) = BaseOpProduct(f.(ops.v))
(f::IndexFunction)(ev::T) where T<:Union{ExpVal,Corr} = T(f(ev.ops))
(f::IndexFunction)(A::OpTerm,nsuminds=A.nsuminds) = OpTerm(nsuminds,f.(A.δs),f.(A.params),f.(A.expvals),f.(A.corrs),f(A.bares))
# some index functions have to be reset for each new term in a sum
(f::IndexFunction)(A::OpSum) = _map_opsum_ops(x->(adapt_for_new_term!(f); f(x)),A)
# by default, nothing to do
adapt_for_new_term!(f::IndexFunction) = f

function _canon_ind(n)
    i,c = divrem(n-1,15)
    i == 0 ? OpIndex('i'+c) : OpIndex('i'+c,i+1)
end
_canon_ind_to_num(ind::OpIndex) = ind.sym-'i'+1 + (ind.num==typemin(ind.num) ? 0 : (ind.num-1)*15)

mutable struct canon_inds <: IndexFunction
    icurr::IndexInt
    canon_inds() = new(0)
end
(f::canon_inds)(ind::OpIndex) = isnoindex(ind) ? ind : _canon_ind(f.icurr+=one(f.icurr))
adapt_for_new_term!(f::canon_inds) = (f.icurr = 0; f)

mutable struct canon_inds_remember <: IndexFunction
    icurr::IndexInt
    replacements::Dict{OpIndex,OpIndex}
    canon_inds_remember() = new(0,Dict{OpIndex,OpIndex}())
end
function (f::canon_inds_remember)(ind::OpIndex)
    isnoindex(ind) && return ind
    i = _canon_ind(f.icurr+=one(f.icurr))
    f.replacements[i] = ind
    i
end
adapt_for_new_term!(f::canon_inds_remember) = error("canon_inds_remember should not be applied on OpSum.")

struct shift_sumind <: IndexFunction
    n::IndexInt
end
(f::shift_sumind)(ind::OpIndex) = issumindex(ind) ? sumindex(ind.num+f.n) : ind

struct replace_inds <: IndexFunction
    replacements::Dict{OpIndex,OpIndex}
    replace_inds(args...) = new(Dict{OpIndex,OpIndex}(args...))
end
(f::replace_inds)(ind::OpIndex) = get(f.replacements,ind,ind)

mutable struct reorder_suminds <: IndexFunction
    icurr::IndexInt
    replacements::Dict{OpIndex,OpIndex}
    reorder_suminds() = new(0,Dict{OpIndex,OpIndex}())
end
(f::reorder_suminds)(ind::OpIndex) = issumindex(ind) ? get!(()->sumindex(f.icurr += one(f.icurr)), f.replacements, ind) : ind
adapt_for_new_term!(f::reorder_suminds) = (f.icurr=0; empty!(f.replacements); f)

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
        push!(f.suminds,ind)
        _canon_ind(f.icurr+=one(f.icurr))
    end
end
adapt_for_new_term!(f::popped_out_suminds) = error("popped_out_suminds should not be applied on OpSum.")

hasind(ind::OpIndex,d::δ) = ind in (d.iA,d.iB)
hasind(ind::OpIndex,A::Union{BaseOperator,Param}) = ind in A.inds
hasind(ind::OpIndex,ev::Union{ExpVal,Corr}) = hasind(ind,ev.ops)
hasind(ind::OpIndex,A::BaseOpProduct) = any(hasind.((ind,),A.v))

split_by_ind(ind::OpIndex,v) = (iwiths = hasind.((ind,),v); (v[.!iwiths], v[iwiths]))

indices(d::δ)::Vector{OpIndex} = [d.iA,d.iB]
indices(A::Union{BaseOperator,Param})::Vector{OpIndex} = collect(assignedinds(A.inds))
indices(v::AbstractVector)::Vector{OpIndex} = vcat(indices.(v)...)
indices(A::BaseOpProduct)::Vector{OpIndex} = indices(A.v)
indices(A::Union{ExpVal,Corr})::Vector{OpIndex} = indices(A.ops)
indices(A::OpTerm)::Vector{OpIndex} = vcat(indices.((A.δs,A.params,A.expvals,A.corrs,A.bares))...)
indices(A::OpSum)::Vector{OpIndex} = vcat(indices.(sort!(collect(keys(A.terms))))...)

"`extindices(A)` return externally visible indices of an expression"
extindices(A) = filter(!issumindex,indices(A))

"`symmetric_index_nums(A)` return sequence of numbers of exchange-symmetric indices"
symmetric_index_nums(A) = symmetric_index_nums(OpSum(A))
function symmetric_index_nums(A::OpSum)
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
