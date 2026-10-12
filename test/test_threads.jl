# concurrent use from several threads (this only tests something if Julia runs with several threads,
# e.g., with JULIA_NUM_THREADS=4)
@testset "threads (nthreads = $(Threads.nthreads()))" begin
    # new names (and values of spins) can be created concurrently
    syms = [Symbol(:threadtest_, k) for k = 1:300]
    names = Vector{QA.QuOpName}(undef, length(syms))
    Threads.@threads for k in eachindex(syms)
        names[k] = QA.QuOpName(syms[k])
    end
    @test allunique([n.i for n in names])
    @test all(QA.sym(n) == s for (n,s) in zip(names, syms))
    @test sortperm(names) == sortperm(syms)
    js = [k//7 for k = 1:100]
    Jzs = Vector{QuExpr}(undef, length(js))
    Threads.@threads for k in eachindex(js)
        Jzs[k] = spin(:S, js[k]; ndims=1)[3][k]
    end
    @test all(QA.spin_j(only(QA.bares(only(keys(Jzs[k].terms))))) == js[k] for k in eachindex(js))

    # computations give the same results in parallel as in serial (the parallel ones are done first,
    # so that the caches are filled concurrently)
    QA.use_σpm(true)
    @boson_ops a
    @fermion_ops f
    @tlspm_ops σ
    @qudit 3 q[:]
    @spin J[:]
    Random.seed!(99)
    pool(s) = QuExpr[a(s), a'(s), f(s)*f'(s), σp(s), σm(s), q(2,3)[s], Jp[s], Jz[s], Pr"g"*a(s)]
    randexpr(sites) = sum(prod((rand(pool(rand(sites))) for _ in 1:rand(1:3)); init=QuExpr(1)) for _ in 1:rand(1:3))
    exprs = [(randexpr((:i,:j,1)), randexpr((:i,:j,1))) for _ = 1:200]
    H = ∑(:k, Pr"ω_k"*a'(:k)*a(:k)) + ∑(:k, Pr"g_k"*(a'(:k)*σm(:k) + σp(:k)*a(:k)))
    function calc(k)
        X, Y = exprs[k]
        # (including new parameter names, created concurrently)
        Z = normal_form(param(Symbol(:threadpar_, k), 'r') * X * Y)
        (comm(X, Y), Z, expval_as_corrs(Z), corr_as_expvals(corr(normal_form(X))), normal_form(heisenberg_eom(X, H, ((Pr"κ", a(1)),))))
    end
    par = Vector{Any}(undef, length(exprs))
    Threads.@threads for k in eachindex(exprs)
        par[k] = calc(k)
    end
    @test all(par[k] == calc(k) for k in eachindex(exprs))
    # equations of motion of correlators (which are cached)
    Hs = [Pr"ω"*a()'*a() + Pr"χ"*a()'*a()'*a()*a(), Pr"ω"*a()'*a() + Pr"g"*(a()'*σm() + σp()*a()) + Pr"Δ"*σp()*σm()]
    eqs = Vector{Any}(undef, 4)
    Threads.@threads for k = 1:4
        eqs[k] = heisenberg_eom_system(Corr, Hs[mod1(k,2)], 3, ((Pr"κ", a()),), (a(),))
    end
    for k = 1:2
        @test eqs[k].eqs == eqs[k+2].eqs == heisenberg_eom_system(Corr, Hs[k], 3, ((Pr"κ", a()),), (a(),)).eqs
    end
    QA.use_σpm(false)
end
