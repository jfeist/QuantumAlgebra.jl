export Avac, vacA, vacExpVal

function _get_mode_name(ex::QuExpr)
    opprod = to_opprod(ex)
    op = unalias(only(opprod.v))
    isempty(op.inds) || throw(ArgumentError("Modes in vacuum should be given without indices"))
    op.name
end

_parse_modes_in_vacuum(modes::Nothing) = nothing
_parse_modes_in_vacuum(modes::QuExpr) = _parse_modes_in_vacuum((modes,))
_parse_modes_in_vacuum(modes) = Set{QuOpName}(_get_mode_name.(modes))

# helper function that is only called on normal-ordered terms that only contain
# operators that are in the vacuum state
function _Avac(BO::BaseOpProduct,fac)
    ii = length(BO.v)
    while ii >= 1
        O = BO.v[ii]
        if O.t in (BosonCreate_,FermionCreate_,TLSCreate_)
            ii==length(BO.v) && return (BO,fac)
            break
        elseif O.t in (BosonDestroy_,FermionDestroy_,TLSDestroy_)
            return (BO,zero(fac))
        elseif O.t == TLSz_
            # vacuum is eigenvalue of σz - since the term is in normal order and
            # we break at TLSx_,TLSy_, we know it acts on the vacuum
            fac = -fac
        elseif O.t in (TLSx_,TLSy_)
            break
        else
            error("should not be reached")
        end
        ii -= 1
    end
    return (BaseOpProduct(BO.v[1:ii]), fac)
end

# helper function that is only called on normal-ordered terms that only contain
# operators that are in the vacuum state
function _vacA(BO::BaseOpProduct,fac)
    ii = 1
    while ii <= length(BO.v)
        O = BO.v[ii]
        if O.t in (BosonDestroy_,FermionDestroy_,TLSDestroy_)
            ii==1 && return (BO,fac)
            break
        elseif O.t in (BosonCreate_,FermionCreate_,TLSCreate_)
            return (BO,zero(fac))
        elseif O.t == TLSz_
            # vacuum is eigenvalue of σz - since the term is in normal order and
            # we break at TLSx_,TLSy_, we know it acts on the vacuum
            fac = -fac
        elseif O.t in (TLSx_,TLSy_)
            break
        else
            error("should not be reached")
        end
        ii += 1
    end
    return (BaseOpProduct(BO.v[ii:end]), fac)
end

function _apply_on_vacuum_modes(func,A::QuTerm,fac,modes_in_vacuum::Nothing)
    BO, fac = func(A.bares,fac)
    Anew = BO === A.bares ? A : QuTerm(A.nsuminds,A.δs,A.params,A.expvals,A.corrs,BO)
    return (Anew, fac)
end

function _apply_on_vacuum_modes(func,A::QuTerm,fac,modes_in_vacuum::Set{QuOpName})
    BO_vacmodes = BaseOpProduct()
    BO_novacmodes = BaseOpProduct()
    for O in A.bares.v
        Oua = unalias(O)
        BO = Oua.name ∈ modes_in_vacuum ? BO_vacmodes : BO_novacmodes
        push!(BO.v, O)
    end
    BO_vacnew, fac = func(BO_vacmodes,fac)
    BO_vacnew === BO_vacmodes && return (A,fac)
    # we know that operators were normal ordered and our operations did not
    # change this, so just sort without commutations
    newv = sort!([BO_novacmodes.v; BO_vacnew.v])
    Anew = QuTerm(A.nsuminds,A.δs,A.params,A.expvals,A.corrs,BaseOpProduct(newv))
    return (Anew, fac)
end

Avac(A::QuTerm,fac,modes_in_vacuum) = _apply_on_vacuum_modes(_Avac, A, fac, modes_in_vacuum)
vacA(A::QuTerm,fac,modes_in_vacuum) = _apply_on_vacuum_modes(_vacA, A, fac, modes_in_vacuum)

Avac(A::QuExpr,modes_in_vacuum=nothing) = (mv = _parse_modes_in_vacuum(modes_in_vacuum); QuExpr(Avac(t,s,mv) for (t,s) in normal_form(A).terms))
vacA(A::QuExpr,modes_in_vacuum=nothing) = (mv = _parse_modes_in_vacuum(modes_in_vacuum); QuExpr(vacA(t,s,mv) for (t,s) in normal_form(A).terms))
Avac(s::Number,modes_in_vacuum=nothing) = QuExpr(s)
vacA(s::Number,modes_in_vacuum=nothing) = QuExpr(s)


"""
    Avac(A::QuExpr), vacA(A::QuExpr)

Simplify operator by assuming it is applied to the vacuum from the left or
right, respectively. To be precise, `Avac(A)` returns ``A'`` such that ``A'|0⟩ =
A|0⟩``, while `vacA(A)` does the same for ``⟨0|A``."""
Avac, vacA

function _TLS_to_pm_normal(A::QuExpr)
    An = QuExpr()
    terms_to_clean = collect(A.terms)
    while !isempty(terms_to_clean)
        (t,s) = pop!(terms_to_clean)
        v = t.bares.v
        iTLS = findfirst(A->A.t ∈ (TLSx_,TLSy_,TLSz_), v)
        if iTLS === nothing
            _add_with_normal_order!(An,t,s,true) # last argument is shortcut_vacA_zero
        else
            O = v[iTLS]
            if O.t == TLSx_
                # σˣ = σ⁺ + σ⁻
                for On in (TLSCreate(O.name,O.inds),TLSDestroy(O.name,O.inds))
                    vn = [v[1:iTLS-1]; On; v[iTLS+1:end]]
                    tn = QuTerm(t.nsuminds,t.δs,t.params,t.expvals,t.corrs,BaseOpProduct(vn))
                    push!(terms_to_clean,tn=>s)
                end
            elseif O.t == TLSy_
                # σʸ = -i σ⁺ + i σ⁻
                for (On,sn) in ((TLSCreate(O.name,O.inds),-1im),(TLSDestroy(O.name,O.inds),1im))
                    vn = [v[1:iTLS-1]; On; v[iTLS+1:end]]
                    tn = QuTerm(t.nsuminds,t.δs,t.params,t.expvals,t.corrs,BaseOpProduct(vn))
                    push!(terms_to_clean,tn=>s*sn)
                end
            elseif O.t == TLSz_
                # σᶻ = 2σ⁺σ⁻ - 1
                for (On,sn) in (([TLSCreate(O.name,O.inds),TLSDestroy(O.name,O.inds)],2),(BaseOperator[],-1))
                    vn = [v[1:iTLS-1]; On; v[iTLS+1:end]]
                    tn = QuTerm(t.nsuminds,t.δs,t.params,t.expvals,t.corrs,BaseOpProduct(vn))
                    push!(terms_to_clean,tn=>s*sn)
                end
            end
        end
    end
    An
end

"""
    vacExpVal(A::QuExpr,S::QuExpr=1)

Calculate the vacuum expectation value ``⟨0|S^\\dagger A S|0⟩``, i.e., the
expectation value ``⟨ψ|A|ψ⟩`` for the state defined by ``|ψ⟩= S|0⟩```.
"""
function vacExpVal(A::QuExpr,stateop::QuExpr=QuExpr(QuTerm()),modes_in_vacuum::Nothing=nothing)
    # simplify down as much as possible by applying vacuum from left and right
    # convert TLSx/y/z operators to TLSCreate/TLSDestroy to ensure that no bare operators survive

    # currently, only nothing is supported, pending implementation of modes_in_vacuum for normal_form(A, shortcut_vacA_zero=true)
    mv = _parse_modes_in_vacuum(modes_in_vacuum)

    # since we will have <vac|stateop' * A * stateop|vac>, we can simplify
    # stateop by applying vacuum from right
    stateop = Avac(stateop,mv)
    x = normal_form(stateop' * A, true) # second argument is shortcut_vacA_zero
    # same here, simplify since we will have <vac|stateop' * A
    x = vacA(x)
    x = _TLS_to_pm_normal(x*stateop)

    vAv = QuExpr()
    for (t,s) in x.terms
        tn,sn = Avac(vacA(t,s,mv)...,mv)
        @assert sn==0 || isempty(tn.bares)
        _add_sum_term!(vAv,tn,simplify_number(sn))
    end
    vAv
end

vacExpVal(A::QuExpr,stateop::Number,modes_in_vacuum=nothing) = vacExpVal(A,QuExpr(stateop),modes_in_vacuum)
vacExpVal(A::Number,stateop=1,modes_in_vacuum=nothing) = vacExpVal(QuExpr(A),stateop,modes_in_vacuum)
