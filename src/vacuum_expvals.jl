export Avac, vacA, vacExpVal

# The modes that are assumed to be in the vacuum state are given by operators without indices (e.g., `a`
# for `@boson a`, or the variables of operators with indices, e.g., `a` for `@boson a[:]`, which stand for
# all modes of the operator), and are identified by their names.
function _mode_names(ex::QuExpr, allowinds::Bool=false)
    names = Set{QuOpName}()
    for t in keys(ex.terms), O in bares(t)
        O = unalias(O)
        allowinds || isempty(assignedinds(inds(O))) ||
            throw(ArgumentError("modes in vacuum have to be given without indices (e.g., `a` for an operator defined with `@boson a[:]`), got $ex"))
        push!(names, opname(O))
    end
    isempty(names) && throw(ArgumentError("modes in vacuum have to be given as operators, got $ex"))
    names
end
function _mode_names(q::QuExprCstrctArrAlias)
    # (the operator with arbitrary indices gives its name)
    ex = q.constructor(ntuple(_ -> 1, q.numdims)...)
    ex isa QuExpr || throw(ArgumentError("modes in vacuum have to be given as operators, got $q"))
    _mode_names(ex, true)
end
_mode_names(v::_PairVariable) = Set([v.ops.name])
_mode_names(x) = throw(ArgumentError("modes in vacuum have to be given as operators, got $x"))

_parse_modes_in_vacuum(modes::Nothing) = nothing
_parse_modes_in_vacuum(modes::Set{QuOpName}) = modes
_parse_modes_in_vacuum(modes::Union{Tuple,AbstractArray,AbstractSet}) = mapreduce(_mode_names, union, modes; init=Set{QuOpName}())
_parse_modes_in_vacuum(modes) = _mode_names(modes)

# helper function that is only called on normal-ordered terms that only contain
# operators that are in the vacuum state: act with the operators on the vacuum
# to the right as long as the vacuum is an eigenstate of them
function _Avac(BO::BaseOpProduct,fac)
    ii = length(BO)
    # eigenvalues are usually integers (and multiplied into fac directly), but can be other numbers
    otherfac = 1
    while ii >= 1
        λ = vacuum_eigenvalue_right(BO[ii])
        λ === nothing && break
        if λ isa Int
            iszero(λ) && return (BO,zero(fac))
            fac *= λ
        else
            otherfac *= _vacuum_value(λ)
        end
        ii -= 1
    end
    ii == length(BO) && return (BO,fac)
    return (BaseOpProduct(BO[1:ii]), isone(otherfac) ? fac : fac*otherfac)
end

# helper function that is only called on normal-ordered terms that only contain
# operators that are in the vacuum state: act with the operators on the vacuum
# to the left as long as the vacuum is an eigenstate of them
function _vacA(BO::BaseOpProduct,fac)
    ii = 1
    # eigenvalues are usually integers (and multiplied into fac directly), but can be other numbers
    otherfac = 1
    while ii <= length(BO)
        λ = vacuum_eigenvalue_left(BO[ii])
        λ === nothing && break
        if λ isa Int
            iszero(λ) && return (BO,zero(fac))
            fac *= λ
        else
            otherfac *= _vacuum_value(λ)
        end
        ii += 1
    end
    ii == 1 && return (BO,fac)
    return (BaseOpProduct(BO[ii:end]), isone(otherfac) ? fac : fac*otherfac)
end

function _apply_on_vacuum_modes(func,A::QuTerm,fac,modes_in_vacuum::Nothing)
    BA = BaseOpProduct(bares(A))
    BO, fac = func(BA,fac)
    Anew = BO === BA ? A : QuTerm(nsuminds(A),deltas(A),params(A),expvals(A),corrs(A),BO)
    return (Anew, fac)
end

function _apply_on_vacuum_modes(func,A::QuTerm,fac,modes_in_vacuum::Set{QuOpName})
    vacmodes = BaseOperator[]
    novacmodes = BaseOperator[]
    for O in bares(A)
        push!(opname(unalias(O)) ∈ modes_in_vacuum ? vacmodes : novacmodes, O)
    end
    BO_vacmodes = BaseOpProduct(vacmodes)
    BO_vacnew, fac = func(BO_vacmodes,fac)
    BO_vacnew === BO_vacmodes && return (A,fac)
    # we know that operators were normal ordered and our operations did not
    # change this, so just sort without commutations
    newv = sort!([novacmodes; collect(BO_vacnew)])
    Anew = QuTerm(nsuminds(A),deltas(A),params(A),expvals(A),corrs(A),BaseOpProduct(newv))
    return (Anew, fac)
end

"""
    Avac(A, modes_in_vacuum=nothing)
    vacA(A, modes_in_vacuum=nothing)

Simplify the expression `A` acting on the vacuum from the right (`Avac`, for ``A|0⟩``) or from the
left (`vacA`, for ``⟨0|A``): `Avac(A)` returns an expression ``A'`` with ``A'|0⟩ = A|0⟩`` (and
similarly for `vacA`), where the terms that annihilate the vacuum are removed and operators for which
the vacuum is an eigenstate (e.g., ``σ^z``) are replaced by their eigenvalues. E.g., `Avac(a*a')`
gives `1`, and `Avac(a*a'*b[1])` gives `0`.

`modes_in_vacuum` restricts the vacuum to the given modes (the others are left as they are): an
operator or a collection of operators, given as the variables defined by the macros (e.g., `a` after
`@boson a`, or `b` after `@boson b[:]`, which stands for all modes of `b`).
"""
Avac(A::QuTerm,fac,modes_in_vacuum) = _apply_on_vacuum_modes(_Avac, A, fac, modes_in_vacuum)
@doc (@doc Avac) vacA
vacA(A::QuTerm,fac,modes_in_vacuum) = _apply_on_vacuum_modes(_vacA, A, fac, modes_in_vacuum)

Avac(A::QuExpr,modes_in_vacuum=nothing) = _apply_on_vacuum(Avac,A,modes_in_vacuum)
vacA(A::QuExpr,modes_in_vacuum=nothing) = _apply_on_vacuum(vacA,A,modes_in_vacuum)
function _apply_on_vacuum(func::F,A::QuExpr,modes_in_vacuum) where F
    mv = _parse_modes_in_vacuum(modes_in_vacuum)
    An = QuExpr()
    for (t,s) in normal_form(A).terms
        _withnum(s) do s
            tn, sn = func(t,s,mv)
            _add_with_auto_order!(An,tn,simplify_number(sn))
        end
    end
    An
end
Avac(s::Number,modes_in_vacuum=nothing) = QuExpr(s)
vacA(s::Number,modes_in_vacuum=nothing) = QuExpr(s)

function _TLS_to_pm_normal(A::QuExpr,shortcut_vacA_zero)
    An = QuExpr()
    terms_to_clean = collect(A.terms)
    while !isempty(terms_to_clean)
        (t,s) = pop!(terms_to_clean)
        _withnum(s -> _TLS_to_pm_normal!(An,terms_to_clean,t,s,shortcut_vacA_zero), s)
    end
    An
end
function _TLS_to_pm_normal!(An::QuExpr,terms_to_clean,t::QuTerm,s,shortcut_vacA_zero)
    v = collect(bares(t))
    iTLS = findfirst(A->optype(A) ∈ (TLSx_,TLSy_,TLSz_), v)
    if iTLS === nothing
        _add_with_normal_order!(An,t,s,shortcut_vacA_zero) # last argument is shortcut_vacA_zero
    else
        O = v[iTLS]
        if optype(O) == TLSx_
            # σˣ = σ⁺ + σ⁻
            for On in (TLSCreate(opname(O),inds(O)),TLSDestroy(opname(O),inds(O)))
                vn = [v[1:iTLS-1]; On; v[iTLS+1:end]]
                tn = QuTerm(nsuminds(t),deltas(t),params(t),expvals(t),corrs(t),BaseOpProduct(vn))
                push!(terms_to_clean,tn=>s)
            end
        elseif optype(O) == TLSy_
            # σʸ = -i σ⁺ + i σ⁻
            for (On,sn) in ((TLSCreate(opname(O),inds(O)),-1im),(TLSDestroy(opname(O),inds(O)),1im))
                vn = [v[1:iTLS-1]; On; v[iTLS+1:end]]
                tn = QuTerm(nsuminds(t),deltas(t),params(t),expvals(t),corrs(t),BaseOpProduct(vn))
                push!(terms_to_clean,tn=>s*sn)
            end
        elseif optype(O) == TLSz_
            # σᶻ = 2σ⁺σ⁻ - 1
            for (On,sn) in (([TLSCreate(opname(O),inds(O)),TLSDestroy(opname(O),inds(O))],2),(BaseOperator[],-1))
                vn = [v[1:iTLS-1]; On; v[iTLS+1:end]]
                tn = QuTerm(nsuminds(t),deltas(t),params(t),expvals(t),corrs(t),BaseOpProduct(vn))
                push!(terms_to_clean,tn=>s*sn)
            end
        end
    end
    nothing
end

"""
    vacExpVal(A, S=1, modes_in_vacuum=nothing)

Return the vacuum expectation value ``⟨0|S^\\dagger A S|0⟩``, i.e., the expectation value
``⟨ψ|A|ψ⟩`` in the state ``|ψ⟩ = S|0⟩`` (which is not normalized automatically). E.g.,
`vacExpVal(a'*a, 1 + a')` gives `1`. With `modes_in_vacuum`, only the given modes are assumed to be
in the vacuum, and the result can contain operators of the other modes (see [`Avac`](@ref)).
"""
function vacExpVal(A::QuExpr,stateop::QuExpr=QuExpr(QuTerm()),modes_in_vacuum=nothing)
    # simplify down as much as possible by applying vacuum from left and right
    # convert TLSx/y/z operators to TLSCreate/TLSDestroy to ensure that no bare operators survive

    mv = _parse_modes_in_vacuum(modes_in_vacuum)

    # since we will have <vac|stateop' * A * stateop|vac>, we can simplify
    # stateop by applying vacuum from right
    stateop = Avac(stateop,mv)
    # second argument is shortcut_vacA_zero, if mv is nothing, we can set it to true
    x = normal_form(stateop' * A, isnothing(mv))
    # same here, simplify since we will have <vac|stateop' * A
    x = vacA(x,mv)
    x = _TLS_to_pm_normal(x*stateop, isnothing(mv)) # second argument is shortcut_vacA_zero

    vAv = QuExpr()
    for (t,s) in x.terms
        _withnum(s) do s
            tn,sn = Avac(vacA(t,s,mv)...,mv)
            isnothing(mv) && @assert sn==0 || isempty(bares(tn))
            _add_sum_term!(vAv,tn,simplify_number(sn))
        end
    end
    vAv
end

vacExpVal(A::QuExpr,stateop::Number,modes_in_vacuum=nothing) = vacExpVal(A,QuExpr(stateop),modes_in_vacuum)
vacExpVal(A::Number,stateop=1,modes_in_vacuum=nothing) = vacExpVal(QuExpr(A),stateop,modes_in_vacuum)
