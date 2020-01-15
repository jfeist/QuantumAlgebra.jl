export Avac, vacA, vacExpVal

Avac(A::Union{BosonDestroy,FermionDestroy,σminus}) = scal(0)
Avac(A::Union{BosonCreate,FermionCreate,σplus}) = A
# vacuum is an eigenstate of σz
Avac(A::σ) = (A.a == z) ? scal(-1) : A
Avac(A::OpSum) = Avac(A.A) + Avac(A.B)
Avac(A::OpSumAnalytic) = OpSumAnalytic(A.ind,Avac(A.A))
Avac(A::OpProd) = begin
    Bv = Avac(A.B)
    Bv != A.B && return Avac(A.A*Bv)
    # if A.A and A.B commute, try what happens with Avac(A.B*A.A)
    if comm(A.A,A.B) == scal(0)
        Av = Avac(A.A)
        Av != A.A && return Avac(A.B*Av)
    end
    # did not find a way to simplify
    A
end
Avac(A::Scalar) = A

vacA(A::Union{BosonDestroy,FermionDestroy,σminus}) = A
vacA(A::Union{BosonCreate,FermionCreate,σplus}) = scal(0)

# vacuum is an eigenstate of σz
vacA(A::σ) = (A.a == z) ? scal(-1) : A
vacA(A::OpSum) = vacA(A.A) + vacA(A.B)
vacA(A::OpSumAnalytic) = OpSumAnalytic(A.ind,vacA(A.A))
vacA(A::OpProd) = begin
    if A.A isa Scalar
        A.A*vacA(A.B)
    else
        # for applying ⟨0|AB, see if ⟨0|A changes A, and if so,
        # keep going with the new operator
        vA = vacA(A.A)
        vA==A.A ?  A : vacA(vA * A.B)
    end
end
vacA(A::Scalar) = A

"""
    Avac(A::Operator), vacA(A::Operator)

Simplify operator by assuming it is applied to the vacuum from the left or
right, respectively. To be precise, `Avac(A)` returns ``A'`` such that ``A'|0⟩ =
A|0⟩``, while `vacA(A)` does the same for ``⟨0|A``."""
Avac, vacA

"""
    vacExpVal(A::Operator,S::Operator=scal(1))

Calculate the vacuum expectation value ``⟨0|S^\\dagger A S|0⟩``, i.e., the
expectation value ``⟨ψ|A|ψ⟩`` for the state defined by ``|ψ⟩= S|0⟩```.
"""
function vacExpVal(A::Operator,stateop::Operator=scal(1))
    # simplify down as much as possible by applying vacuum from left and right
    vsAsv = vacA(Avac(stateop' * A * stateop))
    # only operators that should survive here as operators are σs
    _vacExpVal(vsAsv)
end
_vacExpVal(A::Scalar) = A
# the terms that survive until here have at most a single σ of any particle
# so it does not matter if we are inside a product
_vacExpVal(A::σ) = A.a == z ? scal(-1) : scal(0)
_vacExpVal(A::OpSum) = _vacExpVal(A.A) + _vacExpVal(A.B)
_vacExpVal(A::OpSumAnalytic) = OpSumAnalytic(A.ind,_vacExpVal(A.A))
# we know that the operators here commute (all a and a† have disappeared, and at most a single σ remaining for each particle)
_vacExpVal(A::OpProd) = _vacExpVal(A.A) * _vacExpVal(A.B)
