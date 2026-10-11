# we will want to overload these operators and functions for our custom types
import Base: ==, ≈, *, +, -

export normal_form, comm

Base.length(A::Union{ExpVal,Corr}) = length(ops(A))
Base.length(A::QuTerm) = length(bares(A)) + mapreduce(length,+,expvals(A);init=0) + mapreduce(length,+,corrs(A);init=0)

Base.zero(::Type{QuExpr}) = QuExpr()
Base.zero(::QuExpr) = zero(QuExpr)
Base.one(::Type{QuExpr}) = QuExpr(QuTerm())
Base.one(::Type{QuTerm}) = QuExpr(QuTerm())
Base.one(::T) where T<:Union{QuTerm,QuExpr} = one(T)

function Base.:^(A::Union{QuTerm,QuExpr},n::Int)
    if n > 0
        prod(A for _=1:n)
    elseif n == 0
        one(A)
    else
        throw(ArgumentError("Only positive exponents supported for QuantumAlgebra objects."))
    end
end
Base.literal_pow(::typeof(^),A::Union{QuTerm,QuExpr},::Val{n}) where n = A^n

==(A::BaseOperator,B::BaseOperator) = optype(A) == optype(B) && opname(A) == opname(B) && inds(A) == inds(B)
==(A::BaseOpProduct,B::BaseOpProduct) = A.v == B.v
==(A::Param,B::Param) = opname(A) == opname(B) && parstate(A) == parstate(B) && inds(A) == inds(B)
==(A::T,B::T) where T<:Union{ExpVal,Corr} = ops(A) == ops(B)
==(A::QuTerm,B::QuTerm) = (nsuminds(A) == nsuminds(B) && deltas(A) == deltas(B) && params(A) == params(B) &&
                           expvals(A) == expvals(B) && corrs(A) == corrs(B) && bares(A) == bares(B))
==(A::QuExpr,B::QuExpr) = isequal(A.terms, B.terms)

function ≈(A::QuExpr,B::QuExpr)
    length(A.terms) != length(B.terms) && return false
    for (tA,sA) in A.terms
        sB = get(B.terms,tA) do
            return false
        end
        sA ≈ sB || return false
    end
    return true
end

function Base.hash(ind::QuIndex,h::UInt)
    ### CAREFUL: This relies on the fact that QuIndex has two 32-bit components
    ### and is the same size as a UInt64
    GC.@preserve ind begin
        hash(Base.unsafe_load(Ptr{UInt64}(Base.unsafe_convert(Ptr{Cvoid},Ref(ind)))), h)
    end
end
function Base.hash(name::QuOpName,h::UInt)
    h = hash(name.i,h)
    h
end
function Base.hash(A::BaseOperator,h::UInt)
    h = hash(UInt(optype(A)),h)
    h = hash(opname(A),h)
    h = hash(inds(A),h)
    h
end
function Base.hash(A::Param,h::UInt)
    h = hash(opname(A),h)
    h = hash(parstate(A),h)
    h = hash(inds(A),h)
    h
end
function Base.hash(A::BaseOpProduct,h::UInt)
    h = hash(A.v,h)
    h
end
function Base.hash(A::T,h::UInt) where T<:Union{ExpVal,Corr}
    h = hash(T,h)
    h = hash(ops(A),h)
    h
end
function Base.hash(A::QuTerm,h::UInt)
    h = hash(nsuminds(A),h)
    h = hash(deltas(A),h)
    h = hash(expvals(A),h)
    h = hash(corrs(A),h)
    h = hash(bares(A),h)
    h
end
function Base.hash(A::QuExpr,h::UInt)
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
@inline recursive_cmp(A::T,B::T) where T<:Union{ExpVal,Corr}  = recursive_cmp(ops(A),ops(B))

# compare f(A) and f(B), and return from the calling function if they are not equal
macro _cmpAB(f,rec=true)
    fun = rec ? :recursive_cmp : :cmp
    esc(quote
        cc = $fun($f(A),$f(B))
        cc == 0 || return cc
    end)
end

@inline function Base.cmp(A::BaseOperator,B::BaseOperator)
    @_cmpAB optype false
    @_cmpAB opname false
    @_cmpAB inds false
    # they are equal
    return 0
end

@inline function Base.cmp(A::Param,B::Param)
    @_cmpAB opname false
    @_cmpAB parstate false
    @_cmpAB inds false
    # they are equal
    return 0
end

@inline function Base.cmp(A::δ,B::δ)
    cc = cmp(A.iA,B.iA)
    cc == 0 || return cc
    return cmp(A.iB,B.iB)
end

function Base.cmp(A::QuTerm,B::QuTerm)
    # only evaluate each part that we need for each step of the comparison to avoid unnecessary work
    # order first by number of operators (also within expectation values)
    cc = cmp(length(A), length(B))
    cc == 0 || return cc
    # then order by operators, expectation values, correlations, params, deltas, nsuminds
    @_cmpAB bares
    @_cmpAB expvals
    @_cmpAB corrs
    @_cmpAB params
    @_cmpAB deltas
    @_cmpAB nsuminds false
    # they are equal
    return 0
end

@inline Base.isless(A::BaseOperator,B::BaseOperator) = cmp(A,B) < 0
@inline Base.isless(A::Param,B::Param) = cmp(A,B) < 0
@inline Base.isless(A::δ,B::δ) = cmp(A,B) < 0
@inline Base.isless(A::BaseOpProduct,B::BaseOpProduct) = recursive_cmp(A,B) < 0
@inline Base.isless(A::T,B::T) where T<:Union{ExpVal,Corr} = recursive_cmp(A,B) < 0
@inline Base.isless(A::QuTerm,B::QuTerm) = cmp(A,B) < 0

"""
    comm(A, B)

Compute the commutator ``[A,B] = A*B - B*A``.
"""
comm(A,B) = A*B - B*A

Base.adjoint(A::BaseOperator) = BaseOperator(BaseOpType_adj[Int(optype(A))],opname(A),inds(A))
Base.adjoint(A::BaseOpProduct) = BaseOpProduct(adjoint.(view(A.v, lastindex(A.v):-1:1)))
Base.adjoint(A::T) where {T<:Union{ExpVal,Corr}} = T(adjoint(ops(A)))
Base.adjoint(A::Param) = parstate(A)=='r' ? A : Param(opname(A),parstate(A)=='n' ? 'c' : 'n',inds(A))
function Base.adjoint(A::QuTerm)
    B = QuTerm(nsuminds(A),deltas(A),adjoint.(params(A)),adjoint.(expvals(A)),adjoint.(corrs(A)),adjoint(bares(A)))
    nsuminds(B)>0 ? reorder_suminds()(B) : B
end
Base.adjoint(A::QuExpr) = QuExpr((adjoint(t),adjoint(s)) for (t,s) in A.terms)

is_normal_form(A::Union{ExpVal,Corr}) = is_normal_form(ops(A))
is_normal_form(ops::BaseOpProduct) = is_normal_form(ops.v)
function is_normal_form(v::Vector{BaseOperator})
    issorted(v) || return false
    # now check that no contractions occur
    for k in 1:length(v)-1
        O = v[k]
        dosimplify, fac, op = _contract(O,v[k+1])
        dosimplify && return false
        if optype(O) in (TLSx_,TLSy_,TLSz_)
            # lookahead while we can commute through
            for kp = k+2:length(v)
                optype(v[kp]) in (TLSx_,TLSy_,TLSz_) || break
                _exchange(v[kp-1],v[kp]) == (1,nothing) || break
                all(m -> _exchange(v[m],v[kp]) == (1,nothing), k+1:kp-2) || continue
                dosimplify, fac, op = _contract(O,v[kp])
                dosimplify && return false
            end
        end
    end
    return true
end
function is_normal_form(t::QuTerm)
    (issorted(params(t)) && issorted(expvals(t)) && issorted(corrs(t))) || return false
    is_normal_form(bares(t)) || return false
    for EV in expvals(t)
        is_normal_form(EV) || return false
    end
    for CO in corrs(t)
        is_normal_form(CO) || return false
    end

    # ensure that sum indices are ordered
    changeable_indices = indices(t,false)
    last_sumind = sumindex(0).num
    for ind in Base.Iterators.filter(issumindex,changeable_indices)
        # sumindex have to be increasing without jumps
        ind.num ≤ last_sumind + oneunit(last_sumind) || return false
        last_sumind = max(ind.num,last_sumind)
    end

    isempty(deltas(t)) && return true

    # all indices that are not in δs
    changeable_indices_set = Set{QuIndex}(changeable_indices)
    # accumulate indices that δs want to change (to compare with later δs)
    replace_inds = Set{QuIndex}()
    for (iδ,dd) in enumerate(deltas(t))
        iA,iB = dd.iA, dd.iB
        # δ is not ordered (or disappears if iA == iB)
        iA ≥ iB && return false
        # the whole term disappears in cleanup
        isintindex(iA) && isintindex(iB) && return false
        (issumindex(iA) || issumindex(iB)) && return false
        (iA in replace_inds || iB in replace_inds) && return false
        if iδ > 1
            (dd > @inbounds deltas(t)[iδ-1]) || return false
        end
        iB in changeable_indices_set && return false
        push!(replace_inds, iB)
    end
    return true
end
is_normal_form(A::QuExpr) = all(is_normal_form, keys(A.terms))

# returns a WorkTerm (that can be modified freely) or nothing if the term is zero
function _normalize_without_commutation(A::QuTerm)::Union{WorkTerm,Nothing}
    # first, clean up the δs
    if isempty(deltas(A))
        W = WorkTerm(A)
    else
        δs = sort(deltas(A))
        #println("in _normalize_without_commutation for A = $A, starting with δs = $δs")
        replacements = Dict{QuIndex,QuIndex}()
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
                iB.num ∈ delsuminds || push!(delsuminds,iB.num)
            elseif issumindex(iA) # one sum disappears, and the delta does as well
                replacements[iA] = iB
                # we will later need to shift all larger sumindices down by one
                iA.num ∈ delsuminds || push!(delsuminds,iA.num)
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
            sumrepls = Dict{QuIndex,QuIndex}()
            sort!(delsuminds)
            for (shift,startind) in enumerate(delsuminds)
                endind = shift == length(delsuminds) ? nsuminds(A) : delsuminds[shift+1]-1
                for num = startind+1:endind
                    sumrepls[sumindex(num)] = sumindex(num-shift)
                end
            end
            # we want to first apply replacements and then sumrepls
            # this is the same as replacing results from replacements with sumrepls,
            # and then applying sumrepls on indices that have not been replaced yet
            for (iold,inew) in replacements
                x = get(sumrepls,inew,inew)
                x != inew && (replacements[iold] = x)
            end
            # sumrepl only for indices not in replacements
            merge!(sumrepls,replacements)
            replacements = sumrepls
        end
        f = replace_inds(replacements)
        W = WorkTerm(nsuminds(A)-length(delsuminds), δs, f.(params(A)), [f(ops(E)) for E in expvals(A)],
                     [f(ops(C)) for C in corrs(A)], f(bares(A)))
    end
    # also sort all the commuting terms
    sort!(W.params)
    sort!(W.expvals)
    sort!(W.corrs)
    W
end

# levicivita_lut[a,b] contains the Levi-Cevita symbol ϵ_abc
# for c=6-a-b, i.e, when a,b,c is a permutation of 1,2,3
const levicivita_lut = ((0,1,-1), (-1,0,1), (1,-1,0))
function ϵ_ab(A::BaseOperator,B::BaseOperator)
    # a+b+c == 6 (since a,b,c is a permutation of 1,2,3)
    a = Int(optype(A)) - Int(TLSx_) + 1
    b = Int(optype(B)) - Int(TLSx_) + 1
    c = BaseOpType(Int(TLSx_) - 1 + (6 - a - b))
    s = @inbounds levicivita_lut[a][b]
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
    if optype(A) == optype(B)
        # these operators always commute with the same type
        optype(A) in (BosonDestroy_,BosonCreate_,TLSCreate_,TLSDestroy_,TLSx_,TLSy_,TLSz_) && return (1,nothing)
        # these operators anticommute if they refer to the same species (name), and commute otherwise
        optype(A) in (FermionDestroy_,FermionCreate_) && return (opname(A) == opname(B) ? -1 : 1,nothing)
    end

    # different types of operators commute
    if optype(A) in (BosonDestroy_,BosonCreate_) && optype(B) in (FermionDestroy_,FermionCreate_,TLSx_,TLSy_,TLSz_,TLSDestroy_,TLSCreate_)
        return (1,nothing)
    elseif optype(A) in (FermionDestroy_,FermionCreate_) && optype(B) in (BosonDestroy_,BosonCreate_,TLSx_,TLSy_,TLSz_,TLSDestroy_,TLSCreate_)
        return (1,nothing)
    elseif optype(A) in (TLSx_,TLSy_,TLSz_,TLSDestroy_,TLSCreate_) && optype(B) in (BosonDestroy_,BosonCreate_,FermionDestroy_,FermionCreate_)
        return (1,nothing)
    end

    # a(i) a'(j) = a'(j) a(i) + δij
    if optype(A) == BosonCreate_ && optype(B) == BosonDestroy_
        if opname(A) == opname(B) && (dd = δ(inds(A),inds(B))) !== nothing
            return (1, ExchangeResult(1,dd,nothing))
        else
            return (1, nothing)
        end
    end

    # f(i) f'(j) = -f'(j) f(i) + δij
    if optype(A) == FermionCreate_ && optype(B) == FermionDestroy_
        if opname(A) == opname(B)
            dd = δ(inds(A),inds(B))
            return (-1, dd === nothing ? nothing : ExchangeResult(1,dd,nothing))
        else
            return (1, nothing)
        end
    end

    if optype(A) == TLSCreate_ && optype(B) == TLSDestroy_
        if opname(A) != opname(B) || (dd = δ(inds(A),inds(B))) === nothing
            return (1,nothing)
        elseif isempty(dd)
            # indices were all the same, σ- σ+ = 1 - σ+ σ-
            return (-1, ExchangeResult(1, dd, nothing))
        else
            # σ-_i σ+_j = σ+_j σ-_i - δij σz_i = σ+_j σ-_i + δij (1 - 2 σ+_i σ-_i)
            # note that we return a TLSz to "fit" in ExchangeResult, but need to undo that later on
            # also note that we pass "+TLSz" instead of "-TLSz" so we do not have to flip the sign of the "1" term in the sort algorithm
            return (1, ExchangeResult(1, dd, TLSz(opname(A),inds(A))))
        end
    end

    # B A as x A B + y, with A < B
    if optype(A) in (TLSx_,TLSy_,TLSz_) && optype(B) in (TLSx_,TLSy_,TLSz_)
        if opname(A) != opname(B) || (dd = δ(inds(A),inds(B))) === nothing
            return (1,nothing)
        else
            # we need ϵbac
            c, s = ϵ_ab(B,A)
            if isempty(dd)
                # indices were all the same,
                # σb σa = δab + i ϵbac σc
                # since A < B, we know that a != b
                return (0, ExchangeResult(im*s,dd,BaseOperator(c,opname(A),inds(A))))
            else
                # [σb_i,σa_j] = 2i εbac δij σc_i =>
                # σb_j σa_i = σa_i σb_j + 2i εbac δij σc_i
                return (1, ExchangeResult(2im*s,dd,BaseOperator(c,opname(A),inds(A))))
            end
        end
    end

    if (optype(A) == TLSCreate_ && optype(B) in (TLSx_,TLSy_,TLSz_)) || (optype(A) in (TLSx_,TLSy_,TLSz_) && optype(B) == TLSDestroy_)
        throw(ArgumentError("QuantumAlgebra currently does not support mixing 'normal' (x,y,z) Pauli matrices with jump operators (+,-)."))
    end

    error("_exchange should never reach this! A=$A, B=$B.")
end

const ComplexInt = Complex{Int}

function _contract(A::BaseOperator,B::BaseOperator)::Tuple{Bool,ComplexInt,Union{BaseOperator,Nothing}}
    if optype(A) in (TLSCreate_,TLSDestroy_,FermionDestroy_,FermionCreate_) && A == B
        return (true,zero(ComplexInt),nothing)
    elseif optype(A) in (TLSx_,TLSy_,TLSz_) && optype(B) in (TLSx_,TLSy_,TLSz_) && opname(A) == opname(B) && inds(A) == inds(B)
        # σa σb = δab + i ϵabc σc
        if optype(B) == optype(A)
            return (true,one(ComplexInt),nothing)
        else
            c, s = ϵ_ab(A,B)
            return (true, ComplexInt(0,s), BaseOperator(c,opname(A),inds(A)))
        end
    else
        return (false,zero(ComplexInt),nothing)
    end
end

function normal_order!(ops::BaseOpProduct,term_collector,shortcut_vacA_zero=false,prefactor=one(ComplexInt))
    # do an insertion sort to get to normal ordering
    # reference: https://en.wikipedia.org/wiki/Insertion_sort
    A = ops.v
    for i = 2:length(A)
        j = i
        while j>1 && A[j]<A[j-1]
            # need to commute A[j-1] and A[j]
            pp, exc_res = _exchange(A[j],A[j-1])
            #println("exchanging A[$j] = $(A[j]) and A[$(j-1)] = $(A[j-1]) gave result: $pp, $exc_res.")
            if exc_res !== nothing
                if exc_res.op === nothing
                    onew = BaseOpProduct([A[1:j-2]; A[j+1:end]])
                elseif optype(A[j]) == TLSCreate_
                    @assert optype(exc_res.op) == TLSz_
                    # we got σz_i, have to replace it by (1 - 2 σ+_i σ-_i) (which is really -σz, see explanation in _exchange)
                    onew = BaseOpProduct([A[1:j-2]; TLSCreate(opname(exc_res.op),inds(exc_res.op)); TLSDestroy(opname(exc_res.op),inds(exc_res.op)); A[j+1:end]])
                    t = QuTerm(exc_res.δs, onew)
                    _add_sum_term!(term_collector,t,-2exc_res.pref*prefactor)
                    # this one gives the "1" term that is used below
                    onew = BaseOpProduct([A[1:j-2]; A[j+1:end]])
                else
                    onew = BaseOpProduct([A[1:j-2]; exc_res.op; A[j+1:end]])
                end
                t = QuTerm(exc_res.δs, onew)
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
        if shortcut_vacA_zero && optype(A[1]) in (BosonCreate_,FermionCreate_,TLSCreate_)
            return zero(prefactor)
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
            if optype(A[k]) in (TLSx_,TLSy_,TLSz_)
                # lookahead while we can commute through
                # we already checked with k+1, so start at k+2
                for kp = k+2:length(A)
                    optype(A[kp]) in (TLSx_,TLSy_,TLSz_) || break
                    _exchange(A[kp-1],A[kp]) == (1,nothing) || break
                    # A[kp] has to commute with all operators it is moved past, not just its neighbor
                    all(m -> _exchange(A[m],A[kp]) == (1,nothing), k+1:kp-2) || continue
                    dosimplify, fac, op = _contract(A[k],A[kp])
                    if dosimplify
                        did_contractions = true
                        prefactor *= fac
                        iszero(prefactor) && return prefactor
                        if op === nothing
                            deleteat!(A,(k,kp))
                            #println("deleted indices $k and $kp, A is now: $A, prefactor is now: $prefactor")
                        else
                            A[k] = op
                            deleteat!(A,kp)
                            #println("deleted index $kp and replaced $k by $op. A is now: $A, prefactor is now: $prefactor")
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
        #println("after contractions, A is: $A, prefactor is: $prefactor, term_collector is: $term_collector")
        # do another round of normal ordering since this could have changed after contractions
        # the prefactor of the term we were working on needs to be passed to normal_order!
        # to include in the term_collector terms
        prefactor = normal_order!(ops,term_collector,false,prefactor)
    end
    prefactor
end

# normal order (in place) the operator products in v, which are the contents of expectation values
# or correlators (of type T) in a term. Products that become empty are removed, and v is sorted.
function normal_order!(v::Vector{BaseOpProduct},::Type{T},commterms,pref=1) where {T<:Union{Corr,ExpVal}}
    for (ii,EC) in enumerate(v)
        newterms = QuExpr()
        pp = normal_order!(EC,newterms)
        for (t,s) in newterms.terms
            #@assert iszero(nsuminds(t)) && isempty(params(t)) && isempty(expvals(t)) && isempty(corrs(t))
            # copy the other products since they can still be modified afterwards
            nv = T[T(copy(v[jj])) for jj in eachindex(v) if jj != ii]
            isempty(bares(t)) || pushfirst!(nv, T(bares(t)))
            nt = QuTerm(deltas(t), nv)
            # include current prefactor here
            _add_sum_term!(commterms,nt,pref*s)
        end
        # only modify prefactor after the exchange
        pref *= pp
    end
    filter!(!isempty, v)
    sort!(v)
    pref
end

# iterate over a single term t with prefactor s and a sum as if they were in the same sum
termsumiter(t::QuTerm,s,sum::QuExpr) = Base.Iterators.flatten((((t,s),), sum.terms))
termsumiter(t,s,sum::QuExpr) = termsumiter(QuTerm(t),s,sum)

function _add_with_normal_order!(A::QuExpr,t::QuTerm,s,shortcut_vacA_zero=false)
    # no need to do anything since everything is in normal order
    # in particular, no copy necessary!
    if is_normal_form(t)
       _add_sum_term!(A,t,s)
       return
    end

    W = _normalize_without_commutation(t)
    W === nothing && return
    newbareterms = QuExpr()
    newexpvterms = QuExpr()
    newcorrterms = QuExpr()
    prefbare = normal_order!(W.bares,   newbareterms, shortcut_vacA_zero)
    prefexpv = normal_order!(W.expvals, ExpVal, newexpvterms)
    prefcorr = normal_order!(W.corrs,   Corr,   newcorrterms)

    # normal_form(t) = t.prefs * (prefexpv*tN.ev + ev_new) * (prefcorr*tN.co + co_new) * (prefbare*tN.bare + bare_new)
    # W has been transformed to tN in-place
    t = QuTerm(W)
    pref = prefbare*prefexpv*prefcorr
    if !iszero(pref)
        if iszero(nsuminds(t))
            # should be guaranteed to be in normal order
            _add_sum_term!(A,t,s*pref)
        else
            # might need to reorder sum indices, which can "break" our normal order
            tn = reorder_suminds()(t)
            # this will check again if it is in normal form to make sure
            _add_with_normal_order!(A,tn,s*pref,shortcut_vacA_zero)
        end
    end

    isempty(newexpvterms) && isempty(newcorrterms) && isempty(newbareterms) && return

    firstterm = true
    for (tev,sev) in termsumiter(expvals(t), prefexpv, newexpvterms)
        #@assert iszero(nsuminds(tev)) && isempty(params(tev)) && isempty(corrs(tev)) && isempty(bares(tev))
        for (tco,sco) in termsumiter(corrs(t), prefcorr, newcorrterms)
            #@assert iszero(nsuminds(tco)) && isempty(params(tco)) && isempty(expvals(tco)) && isempty(bares(tco))
            for (tba,sba) in termsumiter(bares(t), prefbare, newbareterms)
                #@assert iszero(nsuminds(tba)) && isempty(params(tba)) && isempty(expvals(tba)) && isempty(corrs(tba))

                if firstterm
                    # the first term has been added above in _add_sum_term!
                    firstterm = false
                    continue
                end

                pref = sev*sco*sba
                iszero(pref) && continue
                nt = QuTerm(nsuminds(t),[deltas(t);deltas(tev);deltas(tco);deltas(tba)],params(t),expvals(tev),corrs(tco),bares(tba))
                _add_with_normal_order!(A,nt,s*pref,shortcut_vacA_zero)
            end
        end
    end
end

"""
    normal_form(A::QuExpr; shortcut_vacA_zero=false)

Return an equivalent expression in normal-ordered form (all commutations performed).

`shortcut_vacA_zero` can be set to `true` when the expression will be multiplied
by a vacuum state from the right, allowing early elimination of annihilation terms.
"""
function normal_form(A::QuExpr,shortcut_vacA_zero=false)
    An = QuExpr()
    for (t,s) in A.terms
        _add_with_normal_order!(An,t,s,shortcut_vacA_zero)
    end
    An
end
normal_form(A::Number) = QuExpr(A)

*(A::BaseOpProduct,B::BaseOpProduct) = BaseOpProduct(vcat(A.v,B.v))

function *(A::QuTerm,B::QuTerm)
    nA, nB = nsuminds(A), nsuminds(B)
    if nA > 0 && nB > 0
        if nA > nB
            A = shift_sumind(nB)(A)
        else
            B = shift_sumind(nA)(B)
        end
    end
    if nA + nB > 0
        # NB: reorder_suminds numbers the sum indices in the order they are encountered,
        # so the order of the calls below matters (has to be the order of the fields of QuTerm)
        f = reorder_suminds()
        apflat = (args...) -> f.(Iterators.flatten(args))
        δs = apflat(deltas(A), deltas(B))
        ps = apflat(params(A), params(B))
        evs = apflat(expvals(A), expvals(B))
        cos = apflat(corrs(A), corrs(B))
        barevec = BaseOpProduct(apflat(bares(A), bares(B)))
        QuTerm(nA+nB, δs, ps, evs, cos, barevec)
    else
        QuTerm(0, vcat(deltas(A), deltas(B)), vcat(params(A), params(B)),
               vcat(expvals(A), expvals(B)), vcat(corrs(A), corrs(B)), bares(A) * bares(B))
    end
end

*(A::Union{QuExpr,QuTerm}) = A
function *(A::QuExpr,B::QuExpr)
    QuExpr((tA*tB,sA*sB) for ((tA,sA),(tB,sB)) in Iterators.product(A.terms,B.terms))
end
*(A::Number,B::QuExpr) = QuExpr((tB,A*sB) for (tB,sB) in B.terms)
*(B::QuExpr,A::Number) = A*B
Base.:/(A::QuExpr,B::Number) = inv(B)*A

+(A::QuExpr) = A
function +(A::QuExpr,B::QuExpr)
    S = copy(A)
    for (t,s) in B.terms
        _add_with_auto_order!(S,t,s)
    end
    S
end
function +(A::QuExpr,B::Number)
    S = copy(A)
    _add_with_auto_order!(S,QuTerm(),B)
    S
end
+(B::Number,A::QuExpr) = A+B

-(A::QuExpr) = -1 * A
function -(A::QuExpr,B::QuExpr)
    S = copy(A)
    for (t,s) in B.terms
        _add_with_auto_order!(S,t,-s)
    end
    S
end
-(B::Number,A::QuExpr) = B + (-A)
-(A::QuExpr,B::Number) = A + (-B)
