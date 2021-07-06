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

"""
    vacExpVal(A::OpSum,S::OpSum=1)

Calculate the vacuum expectation value ``⟨0|S^\\dagger A S|0⟩``, i.e., the
expectation value ``⟨ψ|A|ψ⟩`` for the state defined by ``|ψ⟩= S|0⟩```.
"""
function vacExpVal(A::OpSum,stateop::OpSum=OpSum(OpTerm()))
    # simplify down as much as possible by applying vacuum from left and right
    vsAsv = vacA(Avac(normal_form(stateop' * A * stateop)))
    # only operators that should survive here as operators are σx or σy 
    OpSum(_vacExpVal(t,s) for (t,s) in vsAsv.terms)
end

function _vacExpVal(A::OpTerm,fac)
    # the terms that survive until here have at most a single σ of any particle
    # so it does not matter if we are inside a product
    isempty(A.bares) && return (A,fac)
    for O in A.bares.v
        @assert O.t in (TLSx_,TLSy_)
    end
    if length(A.bares) > 1
        throw(ArgumentError("vacExpVal implementation currently wrong if more than one σˣ or σʸ survives. Got $(A.bares)"))
    end
    return (A,zero(fac))
end
