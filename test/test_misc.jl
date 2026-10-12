# tests of functionality that is not covered elsewhere
@testset "approximate equality" begin
    @boson b
    @qparam g
    @test 1.0*b + 2g ≈ (1 + 1e-12)*b + 2g
    @test !(b ≈ 2b)
    # different terms (with coefficients that cannot be compared with numbers)
    Symbolics.@variables x
    @test !(x*b ≈ x*b')
    @test !(b ≈ b + g)
end
