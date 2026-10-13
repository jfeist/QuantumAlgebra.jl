# Check the prototype against QuantumAlgebra: products, sums and commutators of random expressions
using Test, Random, QuantumAlgebra
const QA = QuantumAlgebra
include(joinpath(@__DIR__, "ConcreteProto.jl"))
using .ConcreteProto

@boson a[:] b
@tls s[:]
@tlsxyz σ[:]
@qparam g[:] ω α::Complex

@testset "concrete representation" begin
    for pm in (false, true)
        QA.use_σpm(pm)
        pool = QuExpr[a[1], a'[1], a[2], a'[2], b, b', s[1], s'[1], s[2], s'[2], σx[1], σy[2], σz[3], σx[2], σy[1],
                      g[1], g[2]*a[1], α, α', ω*σz[1]]
        Random.seed!(1)
        randexpr() = sum(prod(rand(pool, rand(0:4)); init=QuExpr(1)) * rand([1, 2, -1, 1//2, 3im]) for _ in 1:rand(1:4))
        @testset "use_σpm($pm)" begin
            for _ in 1:500
                A, B = randexpr(), randexpr()
                CA, CB = concrete(A, B)
                @test QuExpr(CA) == normal_form(A)
                @test QuExpr(CA*CB) == normal_form(A*B)
                @test QuExpr(CA + CB) == normal_form(A + B)
                @test QuExpr(comm(CA, CB)) == comm(A, B)
            end
        end
    end
    QA.use_σpm(false)
end
