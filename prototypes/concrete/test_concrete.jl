# Check the prototype against QuantumAlgebra: products, sums and commutators of random expressions
using Test, Random, QuantumAlgebra
const QA = QuantumAlgebra
include(joinpath(@__DIR__, "ConcreteProto.jl"))
using .ConcreteProto

@boson a[:] b
@tls s[:]
@tlsxyz σ[:]
@fermion f[:]
@anticommuting_fermions c d
@spin 1//2 S[:]
@spin 1 K
@qparam g[:] ω α::Complex
q = qudit(:q, 3; ndims=1)
Jp, Jm, Jz = spin(:J; ndims=1)
L = so(:L, 4; ndims=1)

function check(pool, n)
    randexpr() = sum(prod(rand(pool, rand(0:3)); init=QuExpr(1)) * rand([1, 2, -1, 1//2, 3im]) for _ in 1:rand(1:3))
    for _ in 1:n
        A, B = randexpr(), randexpr()
        CA, CB = concrete(A, B)
        @test QuExpr(CA) == normal_form(A)
        @test QuExpr(CA*CB) == normal_form(A*B)
        @test QuExpr(CA + CB) == normal_form(A + B)
        @test QuExpr(comm(CA, CB)) == comm(A, B)
    end
end

@testset "concrete representation" begin
    for pm in (false, true)
        QA.use_σpm(pm)
        Random.seed!(1)
        pools = [
            "bosons, two-level systems, parameters" => QuExpr[a[1], a'[1], a[2], a'[2], b, b', s[1], s'[1], s[2], s'[2],
                σx[1], σy[2], σz[3], σx[2], σy[1], g[1], g[2]*a[1], α, α', ω*σz[1]],
            "fermions" => QuExpr[f[1], f'[1], f[2], f'[2], f[3], c, c', d, d', g[1]*f[2]],
            "qudits" => QuExpr[q(1,2)[1], q(2,1)[1], q(3,1)[1], q(2,3)[1], q(1,3)[2], q(3,3)[2], α*q(2,1)[2]],
            "spins" => QuExpr[Sp[1], Sm[1], Sp[2], Sz[1], Sz[2], Kp, Km, Kz, Jp[1], Jm[1], Jz[1], Jp[2]],
            "so(4)" => QuExpr[L(1,2)[1], L(2,3)[1], L(1,4)[1], L(3,4)[2], L(1,2)[2]],
            "all" => QuExpr[a[1], a'[1], s[1], s'[1], σx[2], f[1], f'[2], c', d, q(2,1)[1], q(1,3)[1], Sm[1], Sz[1], Kp, Jz[2], L(1,2)[1], g[1]],
        ]
        for (name, pool) in pools
            @testset "$name (use_σpm($pm))" begin
                check(pool, 300)
            end
        end
    end
    QA.use_σpm(false)
end
