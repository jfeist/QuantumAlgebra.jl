export Avac, vacA, vacExpVal

function _get_mode_name(ex::QuExpr)
    opprod = to_opprod(ex)
    op = unalias(only(opprod))
    isempty(assignedinds(inds(op))) || throw(ArgumentError("Modes in vacuum should be given without indices"))
    opname(op)
end

_parse_modes_in_vacuum(modes::Nothing) = nothing
_parse_modes_in_vacuum(modes::Set{QuOpName}) = modes
_parse_modes_in_vacuum(modes::QuExpr) = _parse_modes_in_vacuum((modes,))
_parse_modes_in_vacuum(modes) = Set{QuOpName}(_get_mode_name.(modes))

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
    Avac(A; modes_in_vacuum=nothing)
    vacA(A; modes_in_vacuum=nothing)

Simplify an operator expression by acting on the vacuum from the right (`Avac`)
or from the left (`vacA`). Returns a new `QuExpr` with terms that annihilate the
vacuum removed and with any TLS `σz` eigenvalue simplifications applied. Set
`modes_in_vacuum` to a set or iterable of modes to restrict the vacuum
assumption; defaults to all modes.
"""
Avac(A::QuTerm,fac,modes_in_vacuum) = _apply_on_vacuum_modes(_Avac, A, fac, modes_in_vacuum)
vacA(A::QuTerm,fac,modes_in_vacuum) = _apply_on_vacuum_modes(_vacA, A, fac, modes_in_vacuum)

Avac(A::QuExpr,modes_in_vacuum=nothing) = (mv = _parse_modes_in_vacuum(modes_in_vacuum); QuExpr(Avac(t,s,mv) for (t,s) in normal_form(A).terms))
vacA(A::QuExpr,modes_in_vacuum=nothing) = (mv = _parse_modes_in_vacuum(modes_in_vacuum); QuExpr(vacA(t,s,mv) for (t,s) in normal_form(A).terms))
Avac(s::Number,modes_in_vacuum=nothing) = QuExpr(s)
vacA(s::Number,modes_in_vacuum=nothing) = QuExpr(s)

function _TLS_to_pm_normal(A::QuExpr,shortcut_vacA_zero)
    An = QuExpr()
    terms_to_clean = collect(A.terms)
    while !isempty(terms_to_clean)
        (t,s) = pop!(terms_to_clean)
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
    end
    An
end

"""
    vacExpVal(A::QuExpr, S::QuExpr=1; modes_in_vacuum=nothing)

Calculate the vacuum expectation value ``⟨0|S^\\\\dagger A S|0⟩``, i.e. the
expectation value ``⟨ψ|A|ψ⟩`` for the state ``|ψ⟩ = S|0⟩``. Optionally restrict
the vacuum assumption to specific modes via `modes_in_vacuum`.
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
        tn,sn = Avac(vacA(t,s,mv)...,mv)
        isnothing(mv) && @assert sn==0 || isempty(bares(tn))
        _add_sum_term!(vAv,tn,simplify_number(sn))
    end
    vAv
end

vacExpVal(A::QuExpr,stateop::Number,modes_in_vacuum=nothing) = vacExpVal(A,QuExpr(stateop),modes_in_vacuum)
vacExpVal(A::Number,stateop=1,modes_in_vacuum=nothing) = vacExpVal(QuExpr(A),stateop,modes_in_vacuum)
