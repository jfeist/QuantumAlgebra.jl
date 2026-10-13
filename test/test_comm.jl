# comm computes the commutator directly (with the Leibniz rule), check it against the definition
@testset "comm" begin
    Random.seed!(4321)
    QA.use_σpm(true)
    @boson_ops a
    @fermion_ops f
    @fermion_ops g
    @anticommuting_fermion_group c d
    @tlspm_ops σ
    τx, τy, τz = ((inds...) -> QuExpr(T(:τ, inds...)) for T in (QA.TLSx, QA.TLSy, QA.TLSz))
    @qudit 3 q[:]
    @spin J[:]
    @spin 1//2 K[:]
    @so 4 L[:]
    @un 3 S[:]
    # (Pauli matrices are only used with integer indices, since normal forms with symbolic indices are not unique)
    pool(s) = QuExpr[a(s), a'(s), f(s), f'(s), g(s), c(s), d'(s), σp(s), σm(s), q(2,1)[s], q(1,3)[s], q(2,3)[s],
                     Jp[s], Jm[s], Jz[s], Kp[s], Km[s], Kz[s], L(1,2)[s], L(2,3)[s], L(1,4)[s], S(2,1)[s], S(1,3)[s], S(2,2)[s],
                     Pr"g"*a(s), expval(a'(s))*f(s), (s isa Integer ? QuExpr[τx(s), τy(s), τz(s)] : QuExpr[])...]
    randprod(sites, n) = prod((rand(pool(rand(sites))) for _ in 1:n); init=QuExpr(1))
    randsum(idx) = rand() < 0.3 ? ∑(idx, prod((rand(pool(idx)) for _ in 1:rand(1:2)); init=QuExpr(1))) : QuExpr(0)
    for _ = 1:500
        sites = rand([(1,2), (:i,:j), (:i,1)])
        X = sum(randprod(sites, rand(0:3)) for _ in 1:rand(1:3)) + randsum(:k)
        Y = sum(randprod(sites, rand(0:3)) for _ in 1:rand(1:3)) + randsum(:l)
        @test comm(X, Y) == normal_form(X*Y - Y*X)
    end
    @test comm(a(1), 2) == comm(2, a(1)) == QuExpr()
    # the result is in normal form
    @test comm(a(:i), a'(:j)*a(:j)) == myδ(:i,:j)*a(:i)
    QA.use_σpm(false)
end
