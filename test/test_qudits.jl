@testset "qudits" begin
    q = qudit_ops(:q, 3)
    @test string(q) == "q (constructor for 3-level qudit)"
    @test string(q(2,1,:i)) == "q²¹(i)"
    # products |i⟩⟨1|1⟩⟨j| on the same qudit are displayed as |i⟩⟨j|
    @test string(q(2,3,:i)) == "q²³(i)"
    @test string(q(1,1)) == "1 - q²²() - q³³()"
    @test latex(q(2,3,:i)) == "{q}_{i}^{23}"
    @test julia_expression(expval(q(2,3,:i))) == :(q_2_3[i])
    @test julia_expression(expval(q(3,2,:i)),[:q_2_3]) == :(conj(q_2_3[i]))
    # ...but only if this does not require moving operators past ones that they might not commute with
    @test string(q(2,1,:i)*q(3,1,:k)*q(1,2,:k)*q(1,3,:i)) == "q²¹(i) q³²(k) q¹³(i)"
    @test string(q(2,1,1)*q(3,1,2)*q(1,2,2)*q(1,3,1)) == "q²³(1) q³²(2)"
    @test latex(q(1,3,:i)) == "{q}_{i}^{13}"
    @test julia_expression(expval(q(1,3,:i))) == :(q_1_3[i])
    @test_throws ArgumentError q(4,1)
    @test_throws ArgumentError q(0,2)
    @test_throws ArgumentError qudit_ops(:q, 1)

    # |1⟩⟨1| is eliminated through completeness
    @test q(1,1) == 1 - q(2,2) - q(3,3)
    @test normal_form(q(1,1) + q(2,2) + q(3,3)) == QuExpr(1)
    @test normal_form(q(1,2)*q(2,3)) == q(1,3)
    @test normal_form(q(1,2)*q(2,1)) == q(1,1)
    @test iszero(normal_form(q(1,2)*q(1,3)))
    @test normal_form(q(2,1,:i)*q(1,3,:j)) == q(2,1,:i)*q(1,3,:j)
    @test normal_form(q(1,3,:j)*q(2,1,:i)) == normal_form(q(2,1,:i)*q(1,3,:j) - myδ(:i,:j)*q(2,3,:i))
    for i=1:3, j=1:3
        @test q(i,j,:k)' == q(j,i,:k)
        # transition operators are idempotent projectors for i == j
        i == j && @test normal_form(q(i,i)*q(i,i)) == normal_form(q(i,i))
    end
    # operators with the same name but different numbers of levels cannot be combined
    @test_throws ArgumentError normal_form(q(1,2)*qudit_ops(:q, 4)(2,1))

    # vacuum (all qudits in level 1)
    @test vacExpVal(q(1,1)) == QuExpr(1)
    @test iszero(vacExpVal(q(2,2)))
    @test vacExpVal(q(3,3), q(3,1)) == QuExpr(1)
    @test vacExpVal(q(2,3), q(2,1) + q(3,1)) == QuExpr(1)
    @test iszero(Avac(q(1,2,:i)))
    @test Avac(q(2,1,:i)) == q(2,1,:i)
    @test iszero(vacA(q(2,1,:i)))

    # a two-level qudit is equivalent to a two-level system with σ⁺ = |2⟩⟨1| and σ⁻ = |1⟩⟨2|
    s = qudit_ops(:s, 2)
    QuantumAlgebra.use_σpm(true)
    @tlspm_ops s
    tls2qudit(A) = QuExpr((QuTerm(QuantumAlgebra.nsuminds(t), QuantumAlgebra.deltas(t), QuantumAlgebra.params(t), QuantumAlgebra.expvals(t), QuantumAlgebra.corrs(t),
                                  BaseOpProduct([QuantumAlgebra.qudit_op(QuantumAlgebra.opname(O), 2, (QuantumAlgebra.optype(O) == QuantumAlgebra.TLSCreate_ ? (2,1) : (1,2))..., QuantumAlgebra.inds(O))
                                                 for O in QuantumAlgebra.bares(t)])), x) for (t,x) in A.terms)
    for X in (sm(:i)*sp(:j), sm(:i)*sp(:j)*sm(:k)*sp(:i), sm(1)*sp(1)*sm(1), (sp(:i)+sm(:j))^3)
        @test normal_form(tls2qudit(normal_form(X))) == normal_form(tls2qudit(X))
    end
    QuantumAlgebra.use_σpm(false)

    # generalized Gell-Mann matrices (with unnormalized diagonal ones)
    λ = gellmann_ops(q)
    @test length(λ) == 8
    λmats = (ComplexF64[0 1 0; 1 0 0; 0 0 0], ComplexF64[0 -im 0; im 0 0; 0 0 0], ComplexF64[1 0 0; 0 -1 0; 0 0 0],
             ComplexF64[0 0 1; 0 0 0; 1 0 0], ComplexF64[0 0 -im; 0 0 0; im 0 0], ComplexF64[0 0 0; 0 0 1; 0 1 0],
             ComplexF64[0 0 0; 0 0 -im; 0 im 0], ComplexF64[1 0 0; 0 1 0; 0 0 -2])
    modes = [((:q,(QuIndex(1),)), 3)]
    for a = 1:8
        @test matrep(λ[a](1), modes) == λmats[a]
        @test λ[a](1)' == λ[a](1)
    end
    @test normal_form(comm(λ[1](:i), λ[2](:i))) == normal_form(2im*λ[3](:i))
    # SU(3) structure constant f_458 = √3/2, i.e., [λ₄,λ₅] = 2i (λ₃/2 + √3/2 λ₈) = i (λ₃ + D₂)
    @test normal_form(comm(λ[4](:i), λ[5](:i))) == normal_form(1im*(λ[3](:i) + λ[8](:i)))
    @test normal_form(λ[1](:i)^2) == normal_form(q(1,1,:i) + q(2,2,:i))

    # Heisenberg equations of motion for a driven three-level system
    H = Pr"Δ"*q(2,2) + Pr"Ω"*(q(2,1) + q(1,2))
    @test normal_form(heisenberg_eom(q(1,2),H)) == normal_form(-1im*Pr"Δ"*q(1,2) + 1im*Pr"Ω"*(q(2,2) - q(1,1)))
end
