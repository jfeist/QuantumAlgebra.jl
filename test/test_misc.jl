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

# canonical numbering of sum indices (src/sum_indices.jl) for terms that need the full search
@testset "canonical sum indices" begin
    @qparam g[:,:]
    @boson b[:]
    # a directed 3-cycle and a directed 6-cycle of parameters g(i,j): colour refinement cannot distinguish the sum
    # indices (each appears once as the first and once as the second index of g), but no symmetry of the term maps
    # the cycles onto each other, so that the search has to compare different branches
    cyc(v) = prod(g[v[k], v[mod1(k+1,length(v))]] for k in eachindex(v))
    ks = [Symbol("k_$n") for n = 1:9]
    T(p, extra=(k -> 1)) = ∑(Tuple(ks), cyc(ks[p[1:3]]) * cyc(ks[p[4:9]]) * extra(ks[p[4]]))
    Tn = normal_form(T(1:9))
    @test length(Tn.terms) == 1
    # (with an operator on one of the sum indices, colour refinement distinguishes all of them)
    Tb(p) = T(p, k -> b[k])
    @test normal_form(Tb(randperm(9))) == normal_form(Tb(1:9))
    Random.seed!(7)
    for _ = 1:10
        # renaming the sum indices gives the same normal form
        p = randperm(9)
        @test normal_form(T(p)) == Tn
        @test iszero(normal_form(T(1:9) - T(p)))
    end
    # if the search exceeds its budget, the sum indices are numbered by their signatures only (so that equivalent
    # terms might not be combined, but the result is still correct)
    budget = QuantumAlgebra._MAX_CANONICAL_NUMBERINGS[]
    try
        QuantumAlgebra._MAX_CANONICAL_NUMBERINGS[] = 1
        X = T(randperm(9))
        Xn = normal_form(X)
        @test normal_form(Xn) == Xn
        @test normal_form(expand_sums(Xn, 1:2)) == normal_form(expand_sums(X, 1:2))
    finally
        QuantumAlgebra._MAX_CANONICAL_NUMBERINGS[] = budget
    end
end

@testset "equation systems" begin
    @boson a
    @qparam ω χ γ
    # dropcorr(n) neglects correlators of more than n operators
    @test dropcorr(1)(expval(a'*a)) == normal_form(expval(a')*expval(a))
    @test dropcorr(1, expval(a'*a*a)) == normal_form(expval(a')*expval(a)^2)
    @test_throws ArgumentError dropcorr(1, expval(a'*a)*a)
    H = ω*a'*a + χ*a'*a'*a*a
    Ls = ((γ, a),)
    eqs = heisenberg_eom_system(H, dropcorr(1), Ls, a)
    @test length(eqs) == 1
    @test eqs[a] == eqs[expval(a)] == eqs[1] == eqs[ExpVal(QuantumAlgebra.to_opprod(a))]
    @test eqs[a] == dropcorr(1)(expval(normal_form(heisenberg_eom(a, H, Ls))))
end

# LaTeX output (compared up to whitespace)
@testset "LaTeX output" begin
    nows(s) = strip(replace(string(s), r"\s+" => " "))
    @boson a
    @boson b[:]
    @qparam g[:] α::Complex
    @test nows(latex(g[:i]*α')) == "g_{i} {\\alpha}^{*}"
    @test nows(latex(normal_form(b[:i]*b'[:j]))) == "\\delta_{ij} + {b}_{j}^{\\dagger} {b}_{i}"
    @test nows(latex(expval(b'[:i]*b[1]))) == "\\langle {b}_{i}^{\\dagger} {b}_{1}\\rangle"
    @test nows(latex(corr(b'[:i]*b[1]))) == "\\langle {b}_{i}^{\\dagger} {b}_{1}\\rangle_{c}"
    @test occursin("{a}^{\\dagger} \\\\", repr(MIME("text/latex"), [a, a']))
    @test occursin("{a} & {a}^{\\dagger}", repr(MIME("text/latex"), [a a'; a' a]))
    @test occursin("g_{1}", repr(MIME("text/latex"), (a, g[1])))
    @qparam ω
    eqs = heisenberg_eom_system(ω*a'*a, 1, (), a)
    @test occursin("\\frac{\\mathrm{d}}{\\mathrm{d}t}\\langle {a}\\rangle", latex(eqs))
    # unicode names (with combining accents, subscripts, ...), which are enclosed in braces if they are not a
    # single character
    @qparam b̄[:] ω̄ₐ ω_c[:] κ²
    @boson Ω̃[:]
    @test nows(latex(b̄[:i]*ω̄ₐ*ω_c[:j]*κ²)) == "{\\bar{b}}_{i} {\\bar{\\omega}_a} {\\omega_c}_{j} {\\kappa^2}"
    @test nows(latex(Ω̃'[1]*Ω̃[:μ])) == "{{\\tilde{\\Omega}}}_{1}^{\\dagger} {{\\tilde{\\Omega}}}_{{\\mu}}"
end

@testset "vacuum" begin
    QuantumAlgebra.use_σxyz()
    @tls σ
    @tlsxyz σ
    # σᶻ is converted to jump operators for vacuum expectation values
    @test vacExpVal(σz) == QuExpr(-1)
    @test vacExpVal(σz, σ') == QuExpr(1)
    # (σᶻ between other operators is only converted at the end)
    @tls s[:]
    @tlsxyz s[:]
    @test vacExpVal(sx[1]*sz[2]*sx[3], 1//2*(1 + s'[1])*(1 + s'[3])) == QuExpr(-1)
    # spins with fixed j: J_z|0⟩ = -j|0⟩ (also from the left, and for only some modes in vacuum)
    @spin 3//2 J
    @boson b
    @test vacA(Jz) == QuExpr(-3//2)
    @test Avac(b'*Jz, Jp) == -3//2*b'
    @so 3 L
    @test_throws ArgumentError Avac(L(1,2))
    @test_throws ArgumentError Avac(b, :b)
    # a contraction of spin operators with commuting operators in between
    @spin 1//2 K
    @test normal_form(Jp*Km*Jm) == normal_form(Km*(15//4 - Jz^2 + Jz))
end

@testset "small functions" begin
    QA = QuantumAlgebra
    @boson b[:]
    @qparam g[:]
    @test one(QA.QuTerm) == QuExpr(1)
    @test *(b[1]) == b[1] && +(b[1]) == b[1]
    @test comm(2, 3) == 0
    @test julia_expression(2) == 2
    # "parameters" that are functions
    @test julia_expression(Pr"f(t)") == :(f(t))
    @test_logs (:warn, r"deprecated") QA.set_quindices_type("Vector")
    @test_throws ArgumentError QuIndex(:i_x)
    @test QA.QuOpName(QA.QuOpName(:b)) == QA.QuOpName(:b)
    @test repr("text/plain", first(boson_ops(:c))) == "c (QuExpr constructor)"
    @test string(QA.QuExprCstrctArrAlias(sin, 1)) == "sin[:]"
    @test_throws ArgumentError QA.parse_name_array(:(a+b))
    X = ∑(:i, g[:i]*b'[:i]) * ∑(:j, b[:j])
    @test QA.reorder_suminds()(X) == X
    # indices of the parts of a term, and the functions to split them by an index (used by other packages)
    t = only(keys(normal_form(g[:i]*expval(b'[:i]*b[1])*b[:k]).terms))
    @test QA.indices(t) == QuIndex.([:i, :i, 1, :k])
    @test QA.indices(QA.expvals(t)) == QA.indices(collect(QA.expvals(t))) == QuIndex.([:i, 1])
    @test QA.hasind(QuIndex(:i), only(QA.params(t))) && !QA.hasind(QuIndex(:k), only(QA.expvals(t)))
    tδ = only(t for t in keys(normal_form(b[:i]*b'[:j]).terms) if !isempty(QA.deltas(t)))
    @test QA.hasind(QuIndex(:j), only(QA.deltas(tδ)))
    @test QA.split_by_ind(QuIndex(:k), collect(QA.bares(t))) == (QA.BaseOperator[], collect(QA.bares(t)))
end

# the storage of terms (views into a vector of words)
@testset "term storage" begin
    QA = QuantumAlgebra
    @boson b[:]
    @qparam g[:]
    t = only(keys((g[1]*g[2]*expval(b[1])*corr(b'[2])*b'[1]*b[2]).terms))
    ps = QA.params(t)
    @test eltype(ps) == QA.Param && firstindex(ps) == 1 && lastindex(ps) == 2 && length(ps) == 2
    @test ps[2] == ps[end] == QA.Param(:g, 'r', 2)
    @test_throws BoundsError ps[3]
    t2 = only(keys((g[1]*g[2]*expval(b[1])*corr(b'[2])*b'[1]*b[2]).terms))
    @test ps == QA.params(t2) && hash(ps) == hash(QA.params(t2))
    @test QA.QuTerm(only(QA.expvals(t))) == only(keys(expval(b[1]).terms))
    @test QA.QuTerm(only(QA.corrs(t))) == only(keys(corr(b'[2]).terms))
    P = QA.BaseOpProduct(collect(QA.bares(t)))
    @test eltype(P) == QA.BaseOperator && firstindex(P) == 1 && collect(eachindex(P)) == [1, 2]
    @test collect(P) == collect(QA.bares(t)) && copy(P) == P
end
