export ascorr, CorrOrExp

"""
    ascorr(expr::Operator)

Take an expression `expr=A B C + D E...` and write its expectation value in
terms of single-body expectation values ``⟨A⟩, ⟨B⟩, \\ldots``, and many-body
correlations ``⟨AB⟩_c, ⟨ABC⟩_c``, etc. Currently, up to fourth-order
correlations (i.e., products of four operators) are supported.

E.g., `ascorr(adag(:n)*a(:n))` returns ``⟨a^\\dagger_n a_n⟩_c + ⟨a^\\dagger_n⟩
⟨a_n⟩`` (which is equal to ``⟨a^\\dagger_n a_n⟩``), while
`ascorr(adag(:n)*a(:m)*a(:n))` returns ``\\langle a_{n}^\\dagger a_{m} a_{n}
\\rangle_{c} + \\langle a_{n}^\\dagger \\rangle \\langle a_{m} \\rangle \\langle
a_{n} \\rangle + \\langle a_{n}^\\dagger \\rangle \\langle a_{m} a_{n}
\\rangle_{c} + \\langle a_{m} \\rangle \\langle a_{n}^\\dagger a_{n}
\\rangle_{c} + \\langle a_{n} \\rangle \\langle a_{n}^\\dagger a_{m}
\\rangle_{c}``.

See also: [`ExpVal`](@ref), [`Corr`](@ref)"""
function ascorr end

ascorr(A::Scalar) = A
for op in (adag,a,f,fdag,σ,σminus,σplus)
    @eval ascorr(A::$op) = ExpVal(A)
end
ascorr(A::OpSum) = ascorr(A.A) + ascorr(A.B)
ascorr(A::OpSumAnalytic) = begin
    # first calculate the correlation for the term in the sum with the "bare" indices, which means that the sum index
    # is assumed to be distinct from the indices of the expressions
    # then, calculate the correction expression for the term where the sum index is the same as any of the non-summed indices,
    # and add the correction between that and the "incorrect" term where index identity is not taken into account
    tmp = ascorr(A.A)
    res = OpSumAnalytic(A.ind,tmp)
    for ii in setdiff(indexset(A),sumindexset(A))
        res += ascorr(replace_index(A.A,A.ind,ii)) - replace_index(tmp,A.ind,ii)
    end
    res
end
function ascorr(A::OpProd)::Operator
    A.A isa OpSumAnalytic && error("should not occur!")
    A.A isa scal && A.B isa OpSumAnalytic && return A.A*ascorr(A.B)
    preftup, exptup, optup = prodtuples(A)
    pref = prod((scal(1),preftup...,exptup...))
    if length(optup)==0
        pref
    elseif length(optup)==1
        pref*ExpVal(optup[1])
    elseif length(optup)==2
        A,B = optup
        pref*(Corr(A*B) + ExpVal(A)*ExpVal(B))
    elseif length(optup)==3
        A,B,C = optup
        pref*(Corr(A*B*C) + ExpVal(A)*Corr(B*C) + ExpVal(B)*Corr(A*C) + ExpVal(C)*Corr(A*B) + ExpVal(A)*ExpVal(B)*ExpVal(C))
    elseif length(optup)==4
        A,B,C,D = optup
        pref*(Corr(A*B*C*D) + ExpVal(A)*ExpVal(B)*ExpVal(C)*ExpVal(D)
            + ExpVal(A)*Corr(B*C*D) + ExpVal(B)*Corr(A*C*D) + ExpVal(C)*Corr(A*B*D) + ExpVal(D)*Corr(A*B*C)
            + ExpVal(A)*ExpVal(B)*Corr(C*D) + Corr(A*B)*ExpVal(C)*ExpVal(D) + Corr(A*B)*Corr(C*D)
            + ExpVal(A)*ExpVal(C)*Corr(B*D) + Corr(A*C)*ExpVal(B)*ExpVal(D) + Corr(A*C)*Corr(B*D)
            + ExpVal(A)*ExpVal(D)*Corr(B*C) + Corr(A*D)*ExpVal(B)*ExpVal(C) + Corr(A*D)*Corr(B*C))
    else
        throw(ArgumentError("ERROR: Only correlations up to fourth order are implemented for now"))
    end
end
CorrOrExp(A::Operator) = length(A)==1 ? ExpVal(A) : Corr(A)
