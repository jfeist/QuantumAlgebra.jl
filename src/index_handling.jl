abstract type IndexFunction end

(f::IndexFunction)(A::BaseOperator) = BaseOperator(A.t,A.a,A.name,f.(A.inds))
(f::IndexFunction)(d::δ) = δ(f(d.iA),f(d.iB))
(f::IndexFunction)(par::Param) = Param(par.name,par.state,f.(par.inds))
(f::IndexFunction)(ops::BaseOpProduct) = BaseOpProduct(f.(ops.v))
(f::IndexFunction)(ev::ExpVal) = ExpVal(f(ev.ops))
(f::IndexFunction)(ev::Corr) = Corr(f(ev.ops))
(f::IndexFunction)(A::OpTerm) = OpTerm(A.nsuminds,f.(A.δs),f.(A.params),f.(A.expvals),f.(A.corrs),f(A.bares))

struct shift_sumind <: IndexFunction
    n::IndexInt
end
(f::shift_sumind)(ind::OpIndex) = issumindex(ind) ? sumindex(ind.num+f.n) : ind

struct replace_inds <: IndexFunction
    replacements::Dict{OpIndex,OpIndex}
end
(f::replace_inds)(ind::OpIndex) = get(f.replacements,ind,ind)

struct reorder_suminds <: IndexFunction
    replacements::Dict{OpIndex,OpIndex}
    icurr::Ref{IndexInt}
    reorder_suminds() = new(Dict{OpIndex,OpIndex}(),Ref(zero(IndexInt)))
end
(f::reorder_suminds)(ind::OpIndex) = issumindex(ind) ? get!(()->sumindex(f.icurr[] += 1), f.replacements, ind) : ind
