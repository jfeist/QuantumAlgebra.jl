# the interface with operators as variables (@boson a[:,:] etc.), compared to the constructor functions
@testset "operator variables" begin
    bo, = boson_ops(:b)
    @boson b
    @test b == bo() && b' == bo'()
    @boson c[:,:]
    co, = boson_ops(:c)
    @test string(c) == "c[:,:]"
    @test c[:i,2] == co(:i,2)
    @test c'[:i,2] == co'(:i,2) == c[:i,2]'
    @test_throws ArgumentError c[:i]
    @test_throws ArgumentError c[:i,:j,:k]
    @test boson(:c; ndims=2)[:i,2] == c[:i,2]
    @test_throws ArgumentError boson(:c; ndims=-1)

    # indices
    @qindex i j[:]
    @test i == QuIndex(:i)
    @test j[2] == QuIndex(:j_2)
    @test string(j) == "j[:]"
    @test c[i,j[2]] == co(:i,:j_2)
    @test qindex(:k) == QuIndex(:k)
    @test qindex(:k; ndims=1)[3] == QuIndex(:k_3)
    @test_throws ArgumentError qindex(:k; ndims=2)

    @fermion f[:]
    fo, = fermion_ops(:f)
    @test f[1] == fo(1) && f'[:i] == fo'(:i)
    @test iszero(normal_form(f[1]*f[1]))
    @anticommuting_fermions c1[:] c2[:]
    @test iszero(normal_form(c1[1]*c2[2] + c2[2]*c1[1]))
    @test normal_form(c1[1]*c1'[1] + c1'[1]*c1[1]) == QuExpr(1)
    d1, d2 = anticommuting_fermions(:d1, :d2)
    @test anticommuting_fermions(:e1, :e2) isa Tuple
    @test iszero(normal_form(d1*d2 + d2*d1))

    # two-level systems
    @tls s[:]
    sold, = tlspm_ops(:s)
    @test s[1] == sold(1) && s'[1] == sold'(1)
    @tlsxyz t
    @test (tx, ty, tz) == map(C -> C(), tlsxyz_ops(:t))
    @test tlsxyz(:t; ndims=1)[3][:i] == tlsxyz_ops(:t)[3](:i)

    # parameters
    @qparam g[:]
    @test g[:i] == Pr"g_i"
    @qparam h::Complex
    @test h == Pc"h" && h' == (Pc"h")'
    @qparam k[:,:]::Complex
    @test k'[1,:j] == (Pc"k_1,j")'
    # several parameters at once
    @qparam ω χ[:] γ::Complex
    @test (ω, χ[1], γ) == (Pr"ω", Pr"χ_1", Pc"γ")
    @test qparam(:χ; ndims=1)[:i] == Pr"χ_i"
    @test qparam(:γ, Complex; ndims=2)'[1,2] == (Pc"γ_1,2")'
    @test_throws ArgumentError qparam(:x, Float64)

    # several operators at once
    @boson b1 b2[:]
    @test b1 == boson(:b1) && b2[1] == boson(:b2; ndims=1)[1]
    @test (@boson b3 b4) == (boson(:b3), boson(:b4))
    @fermion f1 f2[:,:]
    @test f2[1,:i] == fermion(:f2; ndims=2)[1,:i]
    @tls s1 s2[:]
    @test (s1', s2'[1]) == (tls(:s1)', tls(:s2; ndims=1)'[1])
    @tlsxyz t1 t2[:]
    @test (t1x, t2z[1]) == (tlsxyz(:t1)[1], tlsxyz(:t2; ndims=1)[3][1])
    @qudit 3 q1 q2[:]
    @test string(q1) == "q1 (3-level qudit)" && string(q2) == "q2[:] (3-level qudit)"
    @test q2(2,3)[:i] == qudit(:q2, 3; ndims=1)(2,3)[:i]
    @so 4 L1 L2[:]
    @test string(L1) == "L1 (so(4) generators)" && L2(1,2)[:i] == so(:L2, 4; ndims=1)(1,2)[:i]
    @spin J1 J2[:]
    @test (J1p, J2z[1]) == (spin(:J1)[1], spin(:J2; ndims=1)[3][1])
    @spin 1//2 S1 S2[:]
    @test normal_form(S1p*S1m) == normal_form(3//4 - S1z^2 + S1z)
    @test normal_form(S2p[1]*S2m[1]) == normal_form(3//4 - S2z[1]^2 + S2z[1])
    jval = 3//2
    @spin j=jval S3
    @test normal_form(S3p*S3m) == normal_form(15//4 - S3z^2 + S3z)

    # modes in vacuum are given by the operator variables (also for operators with indices)
    @boson va vb[:,:]
    @qudit 3 vq[:]
    @tls vσ
    X = va'*va*vb[1,2]*vb'[1,2]*vq(1,2)[3]
    @test Avac(X, vb) == Avac(X, (vb,)) == normal_form(va'*va*vq(1,2)[3])
    @test iszero(Avac(X, (vb, vq)))
    @test iszero(Avac(X, va))
    @test vacExpVal(X*vσ, 1, (va, vb, vq, vσ)) == QuExpr()
    @test vacA(vσ*vb'[1,1], vσ) == normal_form(vacA(vσ, vσ)*vb'[1,1])
    @test_throws ArgumentError Avac(X, vb[1,2])
    @test_throws ArgumentError Avac(X, qparam(:g; ndims=1))

    # printing in the syntax of the new interface (opt-in)
    @test string(c[:i,2]) == "c(i2)"
    QuantumAlgebra.use_bracket_display(true)
    try
        @test QuantumAlgebra.using_bracket_display()
        @test string(c[:i,2]) == "c[i,2]"
        @test string(c'[i,j[2]]) == "c†[i,j₂]"
        @test string(b'*b) == "b† b"
        @test string(normal_form(b*b')) == "1 + b† b"
        @test string(b^2) == "b²"
        @test string(Pr"g_i,j_2") == "g[i,j₂]"
        @test string(h') == "h*"
        @test string(expval(c'[:i,:j]*c[:j,1])) == "⟨c†[i,j] c[j,1]⟩"
        @test string(qudit(:q, 3; ndims=1)(2,3)[:i]) == "q²³[i]"
        @test string(qudit(:q, 3)(2,3)) == "q²³"
        @test string(∑(:i, c'[:i,1]*c[:i,1])) == "∑₁ c†[#₁,1] c[#₁,1]"
    finally
        QuantumAlgebra.use_bracket_display(false)
    end
end

# declarations with a name that is different from the variable name (e.g., not a valid variable name)
@testset "printed names" begin
    @boson au="a↑" ad[:]="a↓"
    @test string(au'*ad[:i]) == "a↑†() a↓(i)"
    @test au == boson("a↑") == boson(Symbol("a↑"))
    @fermion fu = "f↑"
    @anticommuting_fermions cu="c↑" cd="c↓"
    @test normal_form(cu*cd + cd*cu) == QuExpr()
    @test string(cu) == "c↑()"
    @tls t[:] = "τ̃"
    @tlsxyz s = "s̃"
    @test (t[1], sz) == (tls("τ̃"; ndims=1)[1], last(tlsxyz("s̃")))
    # Pauli matrices without name are printed as X, Y, Z
    X, Y, Z = tlsxyz(""; ndims=1)
    @test string(normal_form(comm(X[:i], Y[:j]))) == "2i δ(ij) Z(i)"
    @test latex(X[:i]*Y[1]) == "{X}_{i} {Y}_{1}"
    @tlsxyz P = ""
    @test (Px, Py, Pz) == tlsxyz("")
    @qudit 3 q = "q̃"
    @test string(q(2,3)) == "q̃²³()"
    @so 3 L[:] = "𝓛"
    @test L(1,2)[:i] == so("𝓛", 3; ndims=1)(1,2)[:i]
    # for @spin, `j = …` is the spin j, except with a string, which is the name of the spin `j`
    @spin j = "𝒥"
    @test string(jz) == "𝒥ᶻ()"
    jval = 1//2
    @spin j = jval K = "κ"
    @test normal_form(Kp*Km) == normal_form(3//4 + Kz - Kz^2)
    @test string(Kz) == "κᶻ()"
    # parameters, also with names that are given by a variable
    name = "x(t)"
    @qparam ft="f(t)" β[:]::Complex=raw"\tilde{\beta}" xt=name
    @test (ft, xt) == (Pr"f(t)", qparam("x(t)"))
    @test julia_expression(ft) == :(f(t))
    @test β'[:i] == param(Symbol("\\tilde{\\beta}"), 'c', :i)
    @test qparam("α", Complex; ndims=1)[1] == Pc"α_1"
end
