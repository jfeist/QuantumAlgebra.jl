# we will want to overload these operators and functions for our custom types
import Base: ==, ≈, *, +, -

export normal_form, comm

Base.length(A::BaseOpProduct) = length(A.v)
Base.length(A::Union{ExpVal,Corr}) = length(A.ops)
Base.length(A::OpTerm) = length(A.bares) + mapreduce(length,+,A.expvals;init=0) + mapreduce(length,+,A.corrs;init=0)

noaliascopy(A::Param) = Param(A.name,A.state,copy(A.inds))
# δ is isbitstype
noaliascopy(A::δ) = A
noaliascopy(A::BaseOperator) = BaseOperator(A.t,A.name,copy(A.inds))
noaliascopy(A::T) where T<:Union{ExpVal,Corr} = T(noaliascopy(A.ops))
noaliascopy(A::BaseOpProduct) = BaseOpProduct(noaliascopy(A.v))
noaliascopy(v::Vector{T}) where T<:Union{δ,Param,BaseOperator,ExpVal,Corr} = noaliascopy.(v)
noaliascopy(A::OpTerm) = OpTerm(A.nsuminds, noaliascopy(A.δs), noaliascopy(A.params),
                                noaliascopy(A.expvals), noaliascopy(A.corrs), noaliascopy(A.bares))
noaliascopy(A::OpSum) = begin
    B = OpSum()
    for (t,s) in A.terms
        B.terms[noaliascopy(t)] = s
    end
    B
end

==(A::BaseOperator,B::BaseOperator) = A.t == B.t && A.name == B.name && A.inds == B.inds
==(A::BaseOpProduct,B::BaseOpProduct) = A.v == B.v
==(A::Param,B::Param) = A.name == B.name && A.state == B.state && A.inds == B.inds
==(A::T,B::T) where T<:Union{ExpVal,Corr} = A.ops == B.ops
==(A::OpTerm,B::OpTerm) = (A.nsuminds == B.nsuminds && A.δs == B.δs && A.params == B.params &&
                           A.expvals == B.expvals && A.corrs == B.corrs && A.bares == B.bares)
==(A::OpSum,B::OpSum) = A.terms == B.terms

function ≈(A::OpSum,B::OpSum)
    length(A.terms) != length(B.terms) && return false
    for (tA,sA) in A.terms
        sB = get(B.terms,tA) do
            return false
        end
        sA ≈ sB || return false
    end
    return true
end

function Base.hash(A::BaseOperator,h::UInt)
    h = hash(A.t,h)
    h = hash(A.name,h)
    h = hash(A.inds,h)
    h
end
function Base.hash(A::Param,h::UInt)
    h = hash(A.name,h)
    h = hash(A.state,h)
    h = hash(A.inds,h)
    h
end
function Base.hash(A::BaseOpProduct,h::UInt)
    h = hash(A.v,h)
    h
end
function Base.hash(A::T,h::UInt) where T<:Union{ExpVal,Corr}
    h = hash(T,h)
    h = hash(A.ops,h)
    h
end
function Base.hash(A::OpTerm,h::UInt)
    h = hash(A.nsuminds,h)
    h = hash(A.δs,h)
    h = hash(A.expvals,h)
    h = hash(A.corrs,h)
    h = hash(A.bares,h)
    h
end
function Base.hash(A::OpSum,h::UInt)
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

@inline recursive_cmp(A::BaseOperator,B::BaseOperator) = cmp(A,B)
@inline recursive_cmp(A::Param,B::Param) = cmp(A,B)
@inline recursive_cmp(A::δ,B::δ) = cmp(A,B)
@inline recursive_cmp(A::BaseOpProduct,B::BaseOpProduct) = (cc = cmp(length(A),length(B)); cc == 0 ? recursive_cmp(A.v,B.v) : cc)
@inline recursive_cmp(A::T,B::T) where T<:Union{ExpVal,Corr}  = recursive_cmp(A.ops,B.ops)

macro _cmpAB(member,rec=true)
    fun = rec ? :recursive_cmp : :cmp
    esc(quote
        cc = $fun(A.$member,B.$member)
        cc == 0 || return cc
    end)
end

@inline function Base.cmp(A::BaseOperator,B::BaseOperator)
    @_cmpAB t false
    @_cmpAB name false
    @_cmpAB inds false
    # they are equal
    return 0
end

@inline function Base.cmp(A::Param,B::Param)
    @_cmpAB name false
    @_cmpAB state false
    @_cmpAB inds false
    # they are equal
    return 0
end

@inline function Base.cmp(A::δ,B::δ)
    @_cmpAB iA false
    @_cmpAB iB false
    # they are equal
    return 0
end

function Base.cmp(A::OpTerm,B::OpTerm)
    # only evaluate each part that we need for each step of the comparison to avoid unnecessary work
    # order first by number of operators (also within expectation values)
    cc = cmp(length(A), length(B))
    cc == 0 || return cc
    # then order by operators, expectation values, correlations, params, deltas, nsuminds
    @_cmpAB bares
    @_cmpAB expvals
    @_cmpAB corrs
    @_cmpAB params
    @_cmpAB δs
    @_cmpAB nsuminds false
    # they are equal
    return 0
end

@inline Base.isless(A::BaseOperator,B::BaseOperator) = cmp(A,B) < 0
@inline Base.isless(A::Param,B::Param) = cmp(A,B) < 0
@inline Base.isless(A::δ,B::δ) = cmp(A,B) < 0
@inline Base.isless(A::BaseOpProduct,B::BaseOpProduct) = recursive_cmp(A,B) < 0
@inline Base.isless(A::T,B::T) where T<:Union{ExpVal,Corr} = recursive_cmp(A,B) < 0
@inline Base.isless(A::OpTerm,B::OpTerm) = cmp(A,B) < 0

comm(A,B) = A*B - B*A

Base.adjoint(A::BaseOperator) = BaseOperator(OpType_adj[Int(A.t)],A.name,A.inds)
Base.adjoint(A::BaseOpProduct) = BaseOpProduct(adjoint.(view(A.v, lastindex(A.v):-1:1)))
Base.adjoint(A::T) where {T<:Union{ExpVal,Corr}} = T(adjoint(A.ops))
Base.adjoint(A::Param) = A.state=='r' ? A : Param(A.name,A.state=='n' ? 'c' : 'n',A.inds)
function Base.adjoint(A::OpTerm)
    B = OpTerm(A.nsuminds,A.δs,adjoint.(A.params),adjoint.(A.expvals),adjoint.(A.corrs),adjoint(A.bares))
    B.nsuminds>0 ? reorder_suminds()(B) : B
end
Base.adjoint(A::OpSum) = OpSum((adjoint(t),adjoint(s)) for (t,s) in A.terms)

is_normal_form(A::Union{ExpVal,Corr}) = is_normal_form(A.ops)
is_normal_form(ops::BaseOpProduct) = is_normal_form(ops.v)
function is_normal_form(v::Vector{BaseOperator})
    issorted(v) || return false
    # now check that no contractions occur
    for k in 1:length(v)-1
        O = v[k]
        dosimplify, fac, op = _contract(O,v[k+1])
        dosimplify && return false
        if O.t in (TLSx_,TLSy_,TLSz_)
            # lookahead while we can commute through
            for kp = k+2:length(v)
                v[kp].t in (TLSx_,TLSy_,TLSz_) || break
                _exchange(v[kp-1],v[kp]) == (1,nothing) || break
                dosimplify, fac, op = _contract(O,v[kp])
                dosimplify && return false
            end
        end
    end
    return true
end
function is_normal_form(t::OpTerm)
    (issorted(t.params) && issorted(t.expvals) && issorted(t.corrs)) || return false
    is_normal_form(t.bares) || return false
    for EV in t.expvals
        is_normal_form(EV) || return false
    end
    for CO in t.corrs
        is_normal_form(CO) || return false
    end

    # ensure that sum indices are ordered
    changeable_indices = vcat(indices.((t.params,t.expvals,t.corrs,t.bares))...)
    last_sumind = sumindex(0).num
    for ind in changeable_indices
        #@show ind, last_sumind
        issumindex(ind) || continue
        # sumindex have to be increasing without jumps
        ind.num ≤ last_sumind + oneunit(last_sumind) || return false
        last_sumind = max(ind.num,last_sumind)
    end

    isempty(t.δs) && return true

    # all indices that are not in δs
    changeable_indices_set = Set{OpIndex}(changeable_indices)
    # accumulate indices that δs want to change (to compare with later δs)
    replace_inds = Set{OpIndex}()
    for (iδ,dd) in enumerate(t.δs)
        iA,iB = dd.iA, dd.iB
        # δ is not ordered (or disappears if iA == iB)
        iA ≥ iB && return false
        # the whole term disappears in cleanup
        isintindex(iA) && isintindex(iB) && return false
        (issumindex(iA) || issumindex(iB)) && return false
        (iA in replace_inds || iB in replace_inds) && return false
        if iδ > 1
            (dd > @inbounds t.δs[iδ-1]) || return false
        end
        iB in changeable_indices_set && return false
        push!(replace_inds, iB)
    end
    return true
end
function is_normal_form(A::OpSum)
    for t in keys(A.terms)
        is_normal_form(t) || return false
    end
    return true
end

function _normalize_without_commutation(A::OpTerm)::Union{OpTerm,Nothing}
    # first, clean up the δs
    if isempty(A.δs)
        # we always want _normalize_without_commutation to return a new object so we can modify it later
        A = noaliascopy(A)
    else
        δs = sort(A.δs)
        #println("in _normalize_without_commutation for A = $A, starting with δs = $δs")
        replacements = Dict{OpIndex,OpIndex}()
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
                push!(delsuminds,iB.num)
            elseif issumindex(iA) # one sum disappears, and the delta does as well
                replacements[iA] = iB
                # we will later need to shift all larger sumindices down by one
                push!(delsuminds,iA.num)
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
            sumrepls = Dict{OpIndex,OpIndex}()
            sort!(delsuminds)
            for (shift,startind) in enumerate(delsuminds)
                endind = shift == length(delsuminds) ? A.nsuminds : delsuminds[shift+1]-1
                for num = startind+1:endind
                    sumrepls[sumindex(num)] = sumindex(num-shift)
                end
            end
            # we want to first apply replacements and then sumrepls
            # this is the same as replacing results from replacements with sumrepls,
            # and then applying sumrepls on indices that have not been replaced yet
            for (iold,inew) in replacements
                replacements[iold] = get(sumrepls,inew,inew)
            end
            for (iold,inew) in sumrepls
                if iold ∉ keys(replacements)
                    replacements[iold] = inew
                end
            end
        end
        nsuminds = A.nsuminds-length(delsuminds)
        f = replace_inds(replacements)
        #println("in _normalize_without_commutation, deleting sum indices $delsuminds and doing replacements $replacements.")
        A = OpTerm(nsuminds,δs,f.(A.params),f.(A.expvals),f.(A.corrs),f(A.bares))
    end
    # also sort all the commuting terms
    sort!(A.params)
    sort!(A.expvals)
    sort!(A.corrs)
    A
end

# levicivita_lut[a,b] contains the Levi-Cevita symbol ϵ_abc
# for c=6-a-b, i.e, when a,b,c is a permutation of 1,2,3
const levicivita_lut = [0 1 -1; -1 0 1; 1 -1 0]
function ϵ_ab(A::BaseOperator,B::BaseOperator)
    # a+b+c == 6 (since a,b,c is a permutation of 1,2,3)
    a = Int(A.t) - Int(TLSx_) + 1
    b = Int(B.t) - Int(TLSx_) + 1
    c = OpType(Int(TLSx_) - 1 + (6 - a - b))
    s = @inbounds levicivita_lut[a,b]
    c, s
end

struct ExchangeResult
    pref::Complex{Rational{Int}}
    δs::Vector{δ}
    op::Union{Nothing,BaseOperator}
end

# rewrite B A as x A B + y, with A < B in our ordering
# returns (x::Int,y::Union{ExchangeResult,Nothing})
function _exchange(A::BaseOperator,B::BaseOperator)::Tuple{Int,Union{ExchangeResult,Nothing}}
    if A.t == B.t
        # these operators always commute with the same type
        A.t in (BosonDestroy_,BosonCreate_,TLSCreate_,TLSDestroy_,TLSx_,TLSy_,TLSz_) && return (1,nothing)
        # these operators anticommute if they refer to the same species (name), and commute otherwise
        A.t in (FermionDestroy_,FermionCreate_) && return (A.name == B.name ? -1 : 1,nothing)
    end

    # different types of operators commute
    if A.t in (BosonDestroy_,BosonCreate_) && B.t in (FermionDestroy_,FermionCreate_,TLSx_,TLSy_,TLSz_,TLSDestroy_,TLSCreate_)
        return (1,nothing)
    elseif A.t in (FermionDestroy_,FermionCreate_) && B.t in (BosonDestroy_,BosonCreate_,TLSx_,TLSy_,TLSz_,TLSDestroy_,TLSCreate_)
        return (1,nothing)
    elseif A.t in (TLSx_,TLSy_,TLSz_,TLSDestroy_,TLSCreate_) && B.t in (BosonDestroy_,BosonCreate_,FermionDestroy_,FermionCreate_)
        return (1,nothing)
    end

    # a(i) adag(j) = adag(j) a(i) + δij
    if A.t == BosonCreate_ && B.t == BosonDestroy_
        if A.name == B.name && (dd = δ(A.inds,B.inds)) !== nothing
            return (1, ExchangeResult(1,dd,nothing))
        else
            return (1, nothing)
        end
    end

    # f(i) fdag(j) = -fdag(j) f(i) + δij
    if A.t == FermionCreate_ && B.t == FermionDestroy_
        if A.name == B.name
            dd = δ(A.inds,B.inds)
            return (-1, dd === nothing ? nothing : ExchangeResult(1,dd,nothing))
        else
            return (1, nothing)
        end
    end

    if A.t == TLSCreate_ && B.t == TLSDestroy_
        if A.name != B.name || (dd = δ(A.inds,B.inds)) === nothing
            return (1,nothing)
        elseif isempty(dd)
            # indices were all the same, σ- σ+ = 1 - σ+ σ-
            return (-1, ExchangeResult(1, dd, nothing))
        else
            # σ-_i σ+_j = σ+_j σ-_i - δij σz_i = σ+_j σ-_i + δij (1 - 2 σ+_i σ-_i)
            # note that we return a TLSz to "fit" in ExchangeResult, but need to undo that later on
            # also note that we pass "+TLSz" instead of "-TLSz" so we do not have to flip the sign of the "1" term in the sort algorithm
            return (1, ExchangeResult(1, dd, TLSz(A.name,A.inds)))
        end
    end

    # B A as x A B + y, with A < B
    if A.t in (TLSx_,TLSy_,TLSz_) && B.t in (TLSx_,TLSy_,TLSz_)
        if A.name != B.name || (dd = δ(A.inds,B.inds)) === nothing
            return (1,nothing)
        else
            # we need ϵbac
            c, s = ϵ_ab(B,A)
            if isempty(dd)
                # indices were all the same,
                # σb σa = δab + i ϵbac σc
                # since A < B, we know that a != b
                return (0, ExchangeResult(im*s,dd,BaseOperator(c,A.name,A.inds)))
            else
                # [σb_i,σa_j] = 2i εbac δij σc_i =>
                # σb_j σa_i = σa_i σb_j + 2i εbac δij σc_i
                return (1, ExchangeResult(2im*s,dd,BaseOperator(c,A.name,A.inds)))
            end
        end
    end

    if (A.t == TLSCreate_ && B.t in (TLSx_,TLSy_,TLSz_)) || (A.t in (TLSx_,TLSy_,TLSz_) && B.t == TLSDestroy_)
        throw(ArgumentError("QuantumAlgebra currently does not support mixing 'normal' (x,y,z) Pauli matrices with jump operators (+,-)."))
    end

    error("_exchange should never reach this! A=$A, B=$B.")
end

const ComplexInt = Complex{Int}

function _contract(A::BaseOperator,B::BaseOperator)::Tuple{Bool,ComplexInt,Union{BaseOperator,Nothing}}
    if A.t in (TLSCreate_,TLSDestroy_,FermionDestroy_,FermionCreate_) && A == B
        return (true,zero(ComplexInt),nothing)
    elseif A.t in (TLSx_,TLSy_,TLSz_) && B.t in (TLSx_,TLSy_,TLSz_) && A.name == B.name && A.inds == B.inds
        # σa σb = δab + i ϵabc σc
        if B.t == A.t
            return (true,one(ComplexInt),nothing)
        else
            c, s = ϵ_ab(A,B)
            return (true, ComplexInt(0,s), BaseOperator(c,A.name,A.inds))
        end
    else
        return (false,zero(ComplexInt),nothing)
    end
end

function normal_order!(ops::BaseOpProduct,term_collector)
    # do an insertion sort to get to normal ordering
    # reference: https://en.wikipedia.org/wiki/Insertion_sort
    A = ops.v
    prefactor::ComplexInt = one(ComplexInt)
    for i = 2:length(A)
        j = i
        while j>1 && A[j]<A[j-1]
            # need to commute A[j-1] and A[j]
            pp, exc_res = _exchange(A[j],A[j-1])
            #println("exchanging A[$j] = $(A[j]) and A[$kk] = $(A[kk]) gave result: $pp, $exc_res.")
            if exc_res !== nothing
                if exc_res.op === nothing
                    onew = BaseOpProduct([A[1:j-2]; A[j+1:end]])
                elseif A[j].t == TLSCreate_
                    @assert exc_res.op.t == TLSz_
                    # we got σz_i, have to replace it by (1 - 2 σ+_i σ-_i) (which is really -σz, see explanation in _exchange)
                    onew = BaseOpProduct([A[1:j-2]; TLSCreate(exc_res.op.name,exc_res.op.inds); TLSDestroy(exc_res.op.name,exc_res.op.inds); A[j+1:end]])
                    t = OpTerm(exc_res.δs, onew)
                    _add_sum_term!(term_collector,t,-2exc_res.pref*prefactor)
                    # this one gives the "1" term that is used below
                    onew = BaseOpProduct([A[1:j-2]; A[j+1:end]])
                else
                    onew = BaseOpProduct([A[1:j-2]; exc_res.op; A[j+1:end]])
                end
                t = OpTerm(exc_res.δs, onew)
                _add_sum_term!(term_collector,t,exc_res.pref*prefactor)
                #println("adding term = $(exc_res.pref*prefactor) * $t, term_collector = $(term_collector)")
            end
            # only modify prefactor after the exchange
            prefactor *= pp
            iszero(prefactor) && return prefactor

            # now finally exchange the two
            A[j-1], A[j] = A[j], A[j-1]
            j -= 1
        end
    end
    # check if we have any products that simplify
    did_contractions = false
    k = 1
    while k < length(A)
        dosimplify, fac, op = _contract(A[k],A[k+1])
        if dosimplify
            did_contractions = true
            prefactor *= fac
            iszero(prefactor) && return prefactor
            if op === nothing
                deleteat!(A,(k,k+1))
                #println("deleted indices $(k-1) and $k, A is now: $A, prefactor is now: $prefactor")
            else
                A[k] = op
                deleteat!(A,k+1)
                #println("deleted index $k and replaced $(k-1) by $op. A is now: $A, prefactor is now: $prefactor")
            end
            # we replaced [...,A[k-1],A[k],A[k+1],A[k+2],...]
            # with        [...,A[k-1],A[k+2],...]
            # or with     [...,A[k-1],Anew,A[k+2]...]
            # so reduce k by one to compare the new operator at k with A[k-1] as well
            # (with the ordering we have, it's probably impossible for A[k-1] to contract
            # with Anew if it didn't contract with A[k], but let's be safe)
            k > 1 && (k -= 1)
        else
            if A[k].t in (TLSx_,TLSy_,TLSz_)
                # lookahead while we can commute through
                # we already checked with k+1, so start at k+2
                for kp = k+2:length(A)
                    A[kp].t in (TLSx_,TLSy_,TLSz_) || break
                    _exchange(A[kp-1],A[kp]) == (1,nothing) || break
                    dosimplify, fac, op = _contract(A[k],A[kp])
                    if dosimplify
                        did_contractions = true
                        prefactor *= fac
                        iszero(prefactor) && return prefactor
                        if op === nothing
                            deleteat!(A,(k,kp))
                            #println("deleted indices $(k-1) and $k, A is now: $A, prefactor is now: $prefactor")
                        else
                            A[k] = op
                            deleteat!(A,kp)
                            #println("deleted index $k and replaced $(k-1) by $op. A is now: $A, prefactor is now: $prefactor")
                        end
                        # k -= 2 here because we have k += 1 just below
                        k > 1 && (k -= 2)
                        break
                    end
                end
            end
            k += 1
        end
    end
    if did_contractions
        # do another round of normal ordering since this could have changed after contractions
        prefactor *= normal_order!(ops,term_collector)
    end
    prefactor
end

normal_order!(A::Union{Corr,ExpVal},term_collector) = normal_order!(A.ops,term_collector)
function normal_order!(v::Vector{T},commterms,pref=1) where {T<:Union{Corr,ExpVal}}
    for (ii,EC) in enumerate(v)
        newterms = OpSum()
        pp = normal_order!(EC,newterms)
        for (t,s) in newterms.terms
            #@assert iszero(t.nsuminds) && isempty(t.params) && isempty(t.expvals) && isempty(t.corrs)
            nv = [isempty(t.bares) ? T[] : T(t.bares); v[1:ii-1]; v[ii+1:end]]
            nt = OpTerm(t.δs, noaliascopy(nv))
            # include current prefactor here
            _add_sum_term!(commterms,nt,pref*s)
        end
        # only modify prefactor after the exchange
        pref *= pp
    end
    filter!(EC -> !isempty(EC.ops), v)
    sort!(v)
    pref
end

# iterate over a single term t with prefactor s and a sum as if they were in the same sum
termsumiter(t::OpTerm,s,sum::OpSum) = Base.Iterators.flatten((((t,s),), sum.terms))
termsumiter(t,s,sum::OpSum) = termsumiter(OpTerm(t),s,sum)

function _add_with_normal_order!(A::OpSum,t::OpTerm,s)
    # no need to do anything since everything is in normal order
    # in particular, no copy necessary!
    if is_normal_form(t)
       _add_sum_term!(A,t,s)
       return
    end

    t = _normalize_without_commutation(t)
    t === nothing && return
    newbareterms = OpSum()
    newexpvterms = OpSum()
    newcorrterms = OpSum()
    prefbare = normal_order!(t.bares,   newbareterms)
    prefexpv = normal_order!(t.expvals, newexpvterms)
    prefcorr = normal_order!(t.corrs,   newcorrterms)

    # normal_form(t) = t.prefs * (prefexpv*tN.ev + ev_new) * (prefcorr*tN.co + co_new) * (prefbare*tN.bare + bare_new)
    # t has already been transformed to tN in-place and is guaranteed to be in normal order
    pref = prefbare*prefexpv*prefcorr
    iszero(t.nsuminds) || (t = reorder_suminds()(t))
    iszero(pref) || _add_sum_term!(A,t,s*pref)
    isempty(newexpvterms) && isempty(newcorrterms) && isempty(newbareterms) && return

    firstterm = true
    for (tev,sev) in termsumiter(t.expvals, prefexpv, newexpvterms)
        #@assert iszero(tev.nsuminds) && isempty(tev.params) && isempty(tev.corrs) && isempty(tev.bares)
        for (tco,sco) in termsumiter(t.corrs, prefcorr, newcorrterms)
            #@assert iszero(tco.nsuminds) && isempty(tco.params) && isempty(tco.expvals) && isempty(tco.bares)
            for (tba,sba) in termsumiter(t.bares, prefbare, newbareterms)
                #@assert iszero(tba.nsuminds) && isempty(tba.params) && isempty(tba.expvals) && isempty(tba.corrs)

                if firstterm
                    # the first term has been added above in _add_sum_term!
                    firstterm = false
                    continue
                end

                pref = sev*sco*sba
                iszero(pref) && continue
                # IMPTE: nt is cleaned afterwards with _normalize_without_commutation,
                # which creates copies of everything, so we can reuse arrays here
                nt = OpTerm(t.nsuminds,[t.δs;tev.δs;tco.δs;tba.δs],t.params,tev.expvals,tco.corrs,tba.bares)
                _add_with_normal_order!(A,nt,s*pref)
            end
        end
    end
end

function normal_form(A::OpSum)
    An = OpSum()
    for (t,s) in A.terms
        _add_with_normal_order!(An,t,s)
    end
    An
end

function *(A::OpTerm,B::OpTerm)
    nsuminds = A.nsuminds + B.nsuminds
    if A.nsuminds > 0 && B.nsuminds > 0
        if A.nsuminds > B.nsuminds
            A = shift_sumind(B.nsuminds)(A)
        else
            B = shift_sumind(A.nsuminds)(B)
        end
    end
    if nsuminds > 0
        f = reorder_suminds()
        apflat = (args...) -> f.(Iterators.flatten(args))
    else
        apflat = (args...) -> vcat(args...)
    end
    δs = apflat(A.δs, B.δs)
    params = apflat(A.params, B.params)
    expvals = apflat(A.expvals, B.expvals)
    corrs = apflat(A.corrs, B.corrs)
    barevec = apflat(A.bares.v, B.bares.v)
    OpTerm(nsuminds,δs,params,expvals,corrs,BaseOpProduct(barevec))
end

*(A::Union{OpSum,OpTerm}) = A
function *(A::OpSum,B::OpSum)
    OpSum((tA*tB,sA*sB) for ((tA,sA),(tB,sB)) in Iterators.product(A.terms,B.terms))
end
*(A::Number,B::OpSum) = OpSum((tB,A*sB) for (tB,sB) in B.terms)
*(B::OpSum,A::Number) = A*B

+(A::OpSum) = A
function +(A::OpSum,B::OpSum)
    S = noaliascopy(A)
    for (t,s) in B.terms
        _add_with_auto_order!(S,t,s)
    end
    S
end
function +(A::OpSum,B::Number)
    S = noaliascopy(A)
    _add_with_auto_order!(S,OpTerm(),B)
    S
end
+(B::Number,A::OpSum) = A+B

-(A::OpSum) = -1 * A
function -(A::OpSum,B::OpSum)
    S = noaliascopy(A)
    for (t,s) in B.terms
        _add_with_auto_order!(S,t,-s)
    end
    S
end
-(B::Number,A::OpSum) = B + (-A)
-(A::OpSum,B::Number) = A + (-B)
