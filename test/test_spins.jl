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
