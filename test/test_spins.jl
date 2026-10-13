@testset "spins" begin
    @spin J
    @test string(Jp) == "J⁺()"
    @test Jp' == Jm && Jm' == Jp
    @test string(Jp*Jz) == "J⁺() Jᶻ()"
    # commutation relations of su(2) in the Cartan-Weyl basis
    @test normal_form(comm(Jz, Jp)) == Jp
    @test normal_form(comm(Jz, Jm)) == -Jm
    @test normal_form(comm(Jp, Jm)) == 2Jz
    Kp, Km, Kz = spin(:K; ndims=1)
    @test string(Kp) == "K⁺[:]"
    @test Kp'[:i] == Km[:i]
    @test latex(Km[:i]) == "{K}_{i}^{-}"
    @test normal_form(comm(Kp[:i], Km[:j])) == normal_form(2myδ(:i,:j)*Kz[:i])
    # without fixed j, products are not simplified
    @test normal_form(Jp*Jm) == Jp*Jm
    @test normal_form(Jm*Jp) == Jp*Jm - 2Jz
    @test_throws ArgumentError vacExpVal(Jz)

    # Casimir: J₊J₋ = j(j+1) - J_z² + J_z for fixed j
    @spin 3//2 L[:]
    @test normal_form(Lp[:i]*Lm[:i]) == normal_form(15//4 - Lz[:i]^2 + Lz[:i])
    @test normal_form(Lp[:i]*Lm[:j]) == Lp[:i]*Lm[:j]
    @test normal_form(Lp[:i]*Lp[:j]*Lm[:i]) == normal_form(Lp[:j]*(15//4 - Lz[:i]^2 + Lz[:i]))
    # J² = J₊J₋ + J_z² - J_z = j(j+1)
    @test normal_form(Lp[1]*Lm[1] + Lz[1]^2 - Lz[1]) == QuExpr(15//4)
    @test normal_form((Lp[1]*Lm[1] + Lm[1]*Lp[1])/2 + Lz[1]^2) == QuExpr(15//4)
    # vacuum: lowest-weight state
    @test vacExpVal(Lz[1]) == QuExpr(-3//2)
    @test vacExpVal(Lz[1], Lp[1]) == QuExpr(-1//2 * 3)
    @test iszero(vacExpVal(Lp[1]*Lm[1]))
    @test vacExpVal(Lm[1]*Lp[1]) == QuExpr(3)
    # different values of j for the same name are not allowed
    @test_throws ArgumentError normal_form(Lp[1]*spin(:L, 1; ndims=1)[2][1])

    # relations of the (2j+1)-dimensional representation: J₊ᵃ ∏_{m=-j}^{j-a} (J_z-m) = 0 and J₋ᵇ ∏_{m=b-j}^{j} (J_z-m) = 0
    @spin 1//2 S1[:]
    @test normal_form(S1z[1]^2) == QuExpr(1//4)
    @test iszero(normal_form(S1p[1]^2)) && iszero(normal_form(S1m[1]^2))
    @test normal_form(S1p[1]*S1z[1]) == -S1p[1]/2
    @test normal_form(S1z[1]*S1p[1]) == S1p[1]/2
    @test normal_form(S1m[1]*S1z[1]) == S1m[1]/2
    @test normal_form(S1p[1]*S1m[1]) == normal_form(1//2 + S1z[1])
    @test normal_form(S1z[1]^7) == S1z[1]/64
    @spin 1 K1
    @test normal_form(K1z^3) == K1z
    @test iszero(normal_form(K1p^3))
    @test normal_form(K1p*K1z^2) == normal_form(-K1p*K1z)
    @test normal_form(K1p^2*K1z) == normal_form(-K1p^2)
    @test normal_form(K1m^2*K1z) == normal_form(K1m^2)
    @test normal_form(Lz[1]^4) == normal_form(5//2*Lz[1]^2 - 9//16)
    # with operators in between that might act on the same site
    @test normal_form(S1p[:i]*S1m[:k]*S1z[:i]) == normal_form(-1//2*S1p[:i]*S1m[:k] + myδ(:i,:k)*(1//2 + S1z[:i]))
    # the normal form only depends on the operator: adding multiples of K_z³ - K_z does not change it
    Random.seed!(11)
    pool = QuExpr[K1p, K1m, K1z, Pr"g"*K1z, K1p*K1z]
    for _ = 1:100
        X, Y, Z = (prod(rand(pool, rand(0:3)); init=QuExpr(1)) for _ = 1:3)
        @test normal_form(X + Y*(K1z^3 - K1z)*Z) == normal_form(X)
        # and the products on one site are J₊ᵃ J_zⁿ and J₋ᵇ J_zⁿ with n ≤ 2j-a (n ≤ 2j-b)
        for t in keys(normal_form(X*Y*Z).terms)
            ops = collect(QuantumAlgebra.bares(t))
            a, n = count(O -> QuantumAlgebra.optype(O) != QuantumAlgebra.SpinZ_, ops), count(O -> QuantumAlgebra.optype(O) == QuantumAlgebra.SpinZ_, ops)
            @test n ≤ 2 - a
        end
    end

    # symbolic j (with Symbolics.jl)
    Symbolics.@variables jsym
    Sp, Sm, Sz = spin(:S, jsym)
    @test isequal(normal_form(Sp*Sm + Sz^2 - Sz), QuExpr(jsym*(jsym+1)))
    @spin j=jsym T
    @test isequal(normal_form(Tp*Tm + Tz^2 - Tz), QuExpr(jsym*(jsym+1)))

    # collective spin: one-axis twisting, d⟨J₋⟩/dt for H = χ J_z²
    H = Pr"χ"*Jz^2
    @test normal_form(heisenberg_eom(Jm,H)) == normal_form(1im*Pr"χ"*(Jm - 2Jm*Jz))
end
