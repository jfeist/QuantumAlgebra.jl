@testset "qudits" begin
    q = qudit(:q, 3)
    # (the same qudit with an index)
    qi = qudit(:q, 3; ndims=1)
    @test string(q) == "q (3-level qudit)"
    @test string(qi) == "q[:] (3-level qudit)"
    @test qi(2,1)'[:i] == qi(1,2)[:i]
    @test_throws ArgumentError qi(2,1)[:i,:j]
    @test string(qi(2,1)[:i]) == "q²¹(i)"
    # products |i⟩⟨1|1⟩⟨j| on the same qudit are displayed as |i⟩⟨j|
    @test string(qi(2,3)[:i]) == "q²³(i)"
    @test string(q(1,1)) == "1 - q²²() - q³³()"
    @test latex(qi(2,3)[:i]) == "{q}_{i}^{23}"
    @test julia_expression(expval(qi(2,3)[:i])) == :(q_2_3[i])
    @test julia_expression(expval(qi(3,2)[:i]),[:q_2_3]) == :(conj(q_2_3[i]))
    # ...but only if this does not require moving operators past ones that they might not commute with
    @test string(qi(2,1)[:i]*qi(3,1)[:k]*qi(1,2)[:k]*qi(1,3)[:i]) == "q²¹(i) q³²(k) q¹³(i)"
    @test string(qi(2,1)[1]*qi(3,1)[2]*qi(1,2)[2]*qi(1,3)[1]) == "q²³(1) q³²(2)"
    @test latex(qi(1,3)[:i]) == "{q}_{i}^{13}"
    @test julia_expression(expval(qi(1,3)[:i])) == :(q_1_3[i])
    @test_throws ArgumentError q(4,1)
    @test_throws ArgumentError q(0,2)
    @test_throws ArgumentError qudit(:q, 1)

    # |1⟩⟨1| is eliminated through completeness
    @test q(1,1) == 1 - q(2,2) - q(3,3)
    @test normal_form(q(1,1) + q(2,2) + q(3,3)) == QuExpr(1)
    @test normal_form(q(1,2)*q(2,3)) == q(1,3)
    @test normal_form(q(1,2)*q(2,1)) == q(1,1)
    @test iszero(normal_form(q(1,2)*q(1,3)))
    @test normal_form(qi(2,1)[:i]*qi(1,3)[:j]) == qi(2,1)[:i]*qi(1,3)[:j]
    @test normal_form(qi(1,3)[:j]*qi(2,1)[:i]) == normal_form(qi(2,1)[:i]*qi(1,3)[:j] - myδ(:i,:j)*qi(2,3)[:i])
    for i=1:3, j=1:3
        @test qi(i,j)[:k]' == qi(j,i)[:k]
        # transition operators are idempotent projectors for i == j
        i == j && @test normal_form(q(i,i)*q(i,i)) == normal_form(q(i,i))
    end
    # operators with the same name but different numbers of levels cannot be combined
    @test_throws ArgumentError normal_form(q(1,2)*qudit(:q, 4)(2,1))

    # vacuum (all qudits in level 1)
    @test vacExpVal(q(1,1)) == QuExpr(1)
    @test iszero(vacExpVal(q(2,2)))
    @test vacExpVal(q(3,3), q(3,1)) == QuExpr(1)
    @test vacExpVal(q(2,3), q(2,1) + q(3,1)) == QuExpr(1)
    @test iszero(Avac(qi(1,2)[:i]))
    @test Avac(qi(2,1)[:i]) == qi(2,1)[:i]
    @test iszero(vacA(qi(2,1)[:i]))

    # a two-level qudit is equivalent to a two-level system with σ⁺ = |2⟩⟨1| and σ⁻ = |1⟩⟨2|
    QuantumAlgebra.use_σpm(true)
    @tlspm_ops s
    tls2qudit(A) = QuExpr((QuTerm(QuantumAlgebra.nsuminds(t), QuantumAlgebra.deltas(t), QuantumAlgebra.params(t), QuantumAlgebra.expvals(t), QuantumAlgebra.corrs(t),
                                  BaseOpProduct([QuantumAlgebra.qudit_op(QuantumAlgebra.opname(O), 2, (QuantumAlgebra.optype(O) == QuantumAlgebra.TLSCreate_ ? (2,1) : (1,2))..., QuantumAlgebra.inds(O))
                                                 for O in QuantumAlgebra.bares(t)])), x) for (t,x) in A.terms)
    for X in (sm(:i)*sp(:j), sm(:i)*sp(:j)*sm(:k)*sp(:i), sm(1)*sp(1)*sm(1), (sp(:i)+sm(:j))^3)
        @test normal_form(tls2qudit(normal_form(X))) == normal_form(tls2qudit(X))
    end
    QuantumAlgebra.use_σpm(false)

    # generalized Gell-Mann matrices, normalized as Tr(λ_a λ_b) = 2δ_ab
    λ = gellmann(qi)
    @test gellmann(q)[1] == q(1,2) + q(2,1)
    @test length(λ) == 8
    @test string(λ[8]) == "λ₈[:]"
    λmats = (ComplexF64[0 1 0; 1 0 0; 0 0 0], ComplexF64[0 -im 0; im 0 0; 0 0 0], ComplexF64[1 0 0; 0 -1 0; 0 0 0],
             ComplexF64[0 0 1; 0 0 0; 1 0 0], ComplexF64[0 0 -im; 0 0 0; im 0 0], ComplexF64[0 0 0; 0 0 1; 0 1 0],
             ComplexF64[0 0 0; 0 0 -im; 0 im 0], ComplexF64[1 0 0; 0 1 0; 0 0 -2]/sqrt(3))
    modes = [((:q,(QuIndex(1),)), 3)]
    for a = 1:8
        @test matrep(λ[a][1], modes) ≈ λmats[a]
        @test λ[a][1]' == λ[a][1]
    end
    @test normal_form(comm(λ[1][:i], λ[2][:i])) == normal_form(2im*λ[3][:i])
    # with an exact square root: SU(3) structure constant f_458 = √3/2, i.e., [λ₄,λ₅] = 2i (λ₃/2 + √3/2 λ₈)
    Sym = SymPyPythonCall.Sym
    λs = gellmann(qi; sqrt = x -> sqrt(Sym(x)))
    @test iszero(normal_form(comm(λs[4][:i], λs[5][:i]) - 1im*(λs[3][:i] + sqrt(Sym(3))*λs[8][:i])))
    # normalization factors that are rational are exact, e.g., for d = 2: (σˣ, -σʸ, -σᶻ)
    p2 = qudit(:p, 2; ndims=1)
    λ2 = gellmann(p2)
    @test λ2[3][:i] == p2(1,1)[:i] - p2(2,2)[:i]
    @test matrep(λ2[2][1], [((:p,(QuIndex(1),)), 2)]) == ComplexF64[0 -im; im 0]
    @test normal_form(λ[1][:i]^2) == normal_form(qi(1,1)[:i] + qi(2,2)[:i])

    # Heisenberg equations of motion for a driven three-level system
    H = Pr"Δ"*q(2,2) + Pr"Ω"*(q(2,1) + q(1,2))
    @test normal_form(heisenberg_eom(q(1,2),H)) == normal_form(-1im*Pr"Δ"*q(1,2) + 1im*Pr"Ω"*(q(2,2) - q(1,1)))
end
