export Avac, vacA, vacExpVal

function Avac(A::OpTerm,fac)
    ii = length(A.bares.v)
    while ii >= 1
        O = A.bares.v[ii]
        if O.t in (BosonCreate_,FermionCreate_,TLSCreate_)
            ii==length(A.bares.v) && return (A,fac)
            break
        elseif O.t in (BosonDestroy_,FermionDestroy_,TLSDestroy_)
            return (A,zero(fac))
        elseif O.t == TLSz_
            # vacuum is eigenvalue of σz
            fac = -fac
        elseif O.t in (TLSx_,TLSy_)
            break
        else
            error("should not be reached")
        end
        ii -= 1
    end
    Anew = OpTerm(A.nsuminds,A.δs,A.params,A.expvals,A.corrs,BaseOpProduct(A.bares.v[1:ii]))
    return (Anew, fac)
end

function vacA(A::OpTerm,fac)
    ii = 1
    while ii <= length(A.bares.v)
        O = A.bares.v[ii]
        if O.t in (BosonDestroy_,FermionDestroy_,TLSDestroy_)
            ii==1 && return (A,fac)
            break
        elseif O.t in (BosonCreate_,FermionCreate_,TLSCreate_)
            return (A,zero(fac))
        elseif O.t == TLSz_
            # vacuum is eigenvalue of σz
            fac = -fac
        elseif O.t in (TLSx_,TLSy_)
            break
        else
            error("should not be reached")
        end
        ii += 1
    end
    Anew = OpTerm(A.nsuminds,A.δs,A.params,A.expvals,A.corrs,BaseOpProduct(A.bares.v[ii:end]))
    return (Anew, fac)
end

Avac(A::OpSum) = OpSum(Avac(t,s) for (t,s) in normal_form(A).terms)
vacA(A::OpSum) = OpSum(vacA(t,s) for (t,s) in normal_form(A).terms)

"""
    Avac(A::OpSum), vacA(A::OpSum)

Simplify operator by assuming it is applied to the vacuum from the left or
right, respectively. To be precise, `Avac(A)` returns ``A'`` such that ``A'|0⟩ =
A|0⟩``, while `vacA(A)` does the same for ``⟨0|A``."""
Avac, vacA

function _TLS_to_pm_normal(A::OpSum)
    An = OpSum()
    terms_to_clean = collect(A.terms)
    while !isempty(terms_to_clean)
        (t,s) = pop!(terms_to_clean)
        v = t.bares.v
        iTLS = findfirst(A->A.t ∈ (TLSx_,TLSy_,TLSz_), v)
        if iTLS === nothing
            _add_with_normal_order!(An,t,s)
        else
            O = v[iTLS]
            if O.t == TLSx_
                # σˣ = σ⁺ + σ⁻
                for On in (TLSCreate(O.name,O.inds),TLSDestroy(O.name,O.inds))
                    vn = [v[1:iTLS-1]; On; v[iTLS+1:end]]
                    tn = OpTerm(t.nsuminds,t.δs,t.params,t.expvals,t.corrs,BaseOpProduct(vn))
                    push!(terms_to_clean,tn=>s)
                end
            elseif O.t == TLSy_
                # σʸ = -i σ⁺ + i σ⁻
                for (On,sn) in ((TLSCreate(O.name,O.inds),-1im),(TLSDestroy(O.name,O.inds),1im))
                    vn = [v[1:iTLS-1]; On; v[iTLS+1:end]]
                    tn = OpTerm(t.nsuminds,t.δs,t.params,t.expvals,t.corrs,BaseOpProduct(vn))
                    push!(terms_to_clean,tn=>s*sn)
                end
            elseif O.t == TLSz_
                # σᶻ = 2σ⁺σ⁻ - 1
                for (On,sn) in (([TLSCreate(O.name,O.inds),TLSDestroy(O.name,O.inds)],2),(BaseOperator[],-1))
                    vn = [v[1:iTLS-1]; On; v[iTLS+1:end]]
                    tn = OpTerm(t.nsuminds,t.δs,t.params,t.expvals,t.corrs,BaseOpProduct(vn))
                    push!(terms_to_clean,tn=>s*sn)
                end
            end
        end
    end
    An
end

"""
    vacExpVal(A::OpSum,S::OpSum=1)

Calculate the vacuum expectation value ``⟨0|S^\\dagger A S|0⟩``, i.e., the
expectation value ``⟨ψ|A|ψ⟩`` for the state defined by ``|ψ⟩= S|0⟩```.
"""
function vacExpVal(A::OpSum,stateop::OpSum=OpSum(OpTerm()))
    # simplify down as much as possible by applying vacuum from left and right
    # convert TLSx/y/z operators to TLSCreate/TLSDestroy to ensure that no bare operators survive
    vsAsv = vacA(Avac(_TLS_to_pm_normal(stateop' * A * stateop)))
    _map_opsum_ops(t -> (@assert isempty(t.bares); t), vsAsv)
end
