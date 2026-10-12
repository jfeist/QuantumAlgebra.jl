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
    so, = tlspm_ops(:s)
    @test s[1] == so(1) && s'[1] == so'(1)
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
        @test rstrip(string(Pr"g_i,j_2")) == "g[i,j₂]"
        @test rstrip(string(h')) == "h*"
        @test rstrip(string(expval(c'[:i,:j]*c[:j,1]))) == "⟨c†[i,j] c[j,1]⟩"
        @test string(qudit(:q, 3; ndims=1)(2,3)[:i]) == "q²³[i]"
        @test string(qudit(:q, 3)(2,3)) == "q²³"
        @test string(∑(:i, c'[:i,1]*c[:i,1])) == "∑₁ c†[#₁,1] c[#₁,1]"
    finally
        QuantumAlgebra.use_bracket_display(false)
    end
end
