@testset "u(n)" begin
    @un 3 S
    # (the same generators with an index)
    Si = un(:S, 3; ndims=1)
    @test string(S) == "S (u(3) generators)"
    @test string(Si) == "S[:] (u(3) generators)"
    @test string(Si(1,2)[:i]) == "S¹²(i)"
    @test latex(Si(1,2)[:i]) == "{S}_{i}^{12}"
    @test S(1,2)' == S(2,1) && S(3,3)' == S(3,3)
    @test Si(1,3)'[:i] == Si(3,1)[:i]
    @test_throws ArgumentError S(1,4)
    @test_throws ArgumentError un(:S, 1)

    # [S_ab, S_cd] = δ_bc S_ad - δ_ad S_cb
    for a=1:3, b=1:3, c=1:3, d=1:3
        @test normal_form(comm(S(a,b), S(c,d))) == normal_form((b==c)*S(a,d) - (a==d)*S(c,b))
    end
    # order: raising (a > b), lowering (a < b), diagonal
    @test normal_form(S(1,2)*S(2,1)) == S(2,1)*S(1,2) + S(1,1) - S(2,2)
    @test normal_form(S(3,3)*S(1,3)*S(3,2)) == normal_form(S(3,2)*S(1,3)*S(3,3) + S(1,2)*S(3,3))
    # different sites
    @test normal_form(comm(Si(1,2)[:i], Si(2,3)[:j])) == normal_form(myδ(:i,:j)*Si(1,3)[:i])
    @test iszero(normal_form(comm(Si(1,2)[1], Si(2,1)[2])))
    # n = 2 is su(2) (plus the central S_11 + S_22), with J₊ = S_21, J₋ = S_12, J_z = (S_22 - S_11)/2
    @un 2 U
    Jp, Jm, Jz = U(2,1), U(1,2), (U(2,2) - U(1,1))/2
    @test normal_form(comm(Jp, Jm)) == normal_form(2Jz)
    @test normal_form(comm(Jz, Jp)) == Jp

    # the linear and quadratic Casimir operators commute with all generators, and the Jacobi identity holds
    for n in (3, 4)
        T = un(Symbol(:T,n), n)
        gens = [T(a,b) for a=1:n for b=1:n]
        C1 = sum(T(a,a) for a=1:n)
        C2 = sum(T(a,b)*T(b,a) for a=1:n for b=1:n)
        @test all(G -> iszero(normal_form(comm(C1, G))) && iszero(normal_form(comm(C2, G))), gens)
        for _ = 1:20
            A, B, D = rand(gens, 3)
            @test iszero(normal_form(comm(A, comm(B, D)) + comm(B, comm(D, A)) + comm(D, comm(A, B))))
        end
    end

    # collective operators of qudits, S_ab = ∑_k |a⟩_k⟨b|_k, satisfy the same commutation relations
    @qudit 3 q[:]
    Q(a,b) = ∑(:k, q(a,b)[:k])
    for (a,b,c,d) in ((1,2,2,3), (2,1,1,2), (3,1,2,3), (2,2,2,1))
        @test normal_form(comm(Q(a,b), Q(c,d))) == normal_form((b==c)*Q(a,d) - (a==d)*Q(c,b))
    end

    # there is no vacuum state
    @test_throws ArgumentError vacExpVal(S(1,1))
end
