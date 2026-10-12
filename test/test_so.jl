@testset "so(n)" begin
    @so 3 L
    # (the same generators with an index)
    Li = so(:L, 3; ndims=1)
    @test string(L) == "L (so(3) generators)"
    @test string(Li) == "L[:] (so(3) generators)"
    @test string(Li(1,2)[:i]) == "L¹²(i)"
    @test latex(Li(1,2)[:i]) == "{L}_{i}^{12}"
    # antisymmetry and hermiticity
    @test Li(2,1)[:i] == -Li(1,2)[:i]
    @test iszero(L(3,3))
    @test L(1,3)' == L(1,3)
    @test_throws ArgumentError L(1,4)
    @test_throws ArgumentError so(:L, 1)

    # so(3): (L_x, L_y, L_z) = (L₂₃, L₃₁, L₁₂) with [L_a, L_b] = i ϵ_abc L_c
    Lx, Ly, Lz = L(2,3), L(3,1), L(1,2)
    @test normal_form(comm(Lx, Ly)) == im*Lz
    @test normal_form(comm(Ly, Lz)) == im*Lx
    @test normal_form(comm(Lz, Lx)) == im*Ly
    # ... and is isomorphic to su(2), with J± = L_x ± i L_y and J_z = L_z
    @test normal_form(comm(Lz, Lx + im*Ly)) == normal_form(Lx + im*Ly)
    @test normal_form(comm(Lx + im*Ly, Lx - im*Ly)) == 2Lz
    # different sites
    @test normal_form(comm(Li(1,2)[:i], Li(2,3)[:j])) == normal_form(-im*myδ(:i,:j)*Li(1,3)[:i])
    @test iszero(normal_form(comm(Li(1,2)[1], Li(2,3)[2])))

    # so(n): generators with disjoint indices commute, and the quadratic Casimir commutes with all generators
    for n in (4, 5)
        M = so(Symbol(:M,n), n)
        gens = [M(a,b) for a=1:n for b=a+1:n]
        @test iszero(normal_form(comm(M(1,2), M(3,4))))
        C = sum(G^2 for G in gens)
        @test all(G -> iszero(normal_form(comm(C, G))), gens)
        # Jacobi identity
        for _ = 1:20
            A, B, D = rand(gens, 3)
            @test iszero(normal_form(comm(A, comm(B, D)) + comm(B, comm(D, A)) + comm(D, comm(A, B))))
        end
    end

    # the same name cannot be used with different n (detected when the operators are reordered)
    @test_throws ArgumentError normal_form(so(:L, 4)(1,3)*L(1,2))
    # there is no vacuum state
    @test_throws ArgumentError vacExpVal(Lz)
end
