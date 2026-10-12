# the tool for migrating code from the legacy interface (tools/migration/migrate.jl)
isdefined(Main, :QAMigration) || include(joinpath(pkgdir(QuantumAlgebra), "tools", "migration", "migrate.jl"))

@testset "migration tool" begin
    QuantumAlgebra.use_σxyz()
    old = raw"""
    using QuantumAlgebra
    @boson_ops b
    @tlspm_ops s
    @anticommuting_fermion_group c d
    H = Pr"ω"*a'()*a() + Pr"g_i"*(a'()*sm(:i) + sp(:i)*a()) + ∑(:k, Pc"η_k"*b'(:k)*b(:k)) + 2*bdag(1)*b(1)
    X = normal_form(comm(H, a()))
    Y = normal_form(c(1)*d'(2) + σx(:j)*σz(:j) + σy(:j)') + param(:α,'c',3)*f(1,2)
    Z = vacExpVal(a'()*a(), a'()^2)
    """
    new, msgs = QAMigration.migrate_source(old)
    @test new == raw"""
    using QuantumAlgebra
    # operators and parameters (migrated from the legacy QuantumAlgebra interface)
    @boson a
    @fermion f[:,:]
    @tlsxyz σ[:]
    @qparam ω
    @qparam g[:]
    @qparam η[:]::Complex
    @qparam α[:]::Complex
    @boson b[:]
    @tls s[:]
    @anticommuting_fermions c[:] d[:]
    H = ω*a'*a + g[:i]*(a'*s[:i] + s'[:i]*a) + ∑(:k, η[:k]*b'[:k]*b[:k]) + 2*b'[1]*b[1]
    X = normal_form(comm(H, a))
    Y = normal_form(c[1]*d'[2] + σx[:j]*σz[:j] + σy[:j]') + α'[3]*f[1,2]
    Z = vacExpVal(a'*a, a'^2)
    """
    @test all(m -> startswith(m, "migrated"), msgs)
    # both versions give the same results (each group of anticommuting fermions gets a unique internal
    # name, so results with them are compared as strings)
    mold, mnew = Module(), Module()
    Base.include_string(mold, old)
    Base.include_string(mnew, new)
    @test (mold.X, mold.Z) == (mnew.X, mnew.Z)
    @test string(mold.Y) == string(mnew.Y)

    # code that cannot be migrated is left unchanged and reported
    old = raw"""
    using QuantumAlgebra
    @boson_ops b
    x = b() + b(1)
    y = map(a, 1:3)
    g(x) = 2x
    z = Pr"g_i" + f'(2)
    """
    new, msgs = QAMigration.migrate_source(old)
    @test new == raw"""
    using QuantumAlgebra
    # operators and parameters (migrated from the legacy QuantumAlgebra interface)
    @fermion f[:]
    @boson_ops b
    x = b() + b(1)
    y = map(a, 1:3)
    g(x) = 2x
    z = Pr"g_i" + f'[2]
    """
    @test any(m -> occursin("`b` (declared at code:2): used with different numbers of indices (0, 1)", m), msgs)
    @test any(m -> occursin("`a` (default operator): used as a value", m), msgs)
    @test any(m -> occursin("parameter `g`: the new variable `g` would clash", m), msgs)

    # declarations in local scopes (e.g., in @testset, functions, let blocks) are only valid there
    old = raw"""
    using QuantumAlgebra, Test
    @testset "local" begin
        @boson_ops a
        @fermion_ops g
        global X = a(1)*a'(:i)*g(2) + Pr"h"*a(2)
    end
    function makeY()
        @tlspm_ops s
        sm(1)*sp(2)
    end
    Y = makeY()
    Z = a()*a'()
    """
    new, msgs = QAMigration.migrate_source(old)
    @test new == raw"""
    using QuantumAlgebra, Test
    # operators and parameters (migrated from the legacy QuantumAlgebra interface)
    @boson a
    @qparam h
    @testset "local" begin
        @boson a[:]
        @fermion g[:]
        global X = a[1]*a'[:i]*g[2] + h*a[2]
    end
    function makeY()
        @tls s[:]
        s[1]*s'[2]
    end
    Y = makeY()
    Z = a*a'
    """
    mold, mnew = Module(), Module()
    Base.include_string(mold, old)
    Base.include_string(mnew, new)
    @test (mold.X, mold.Y, mold.Z) == (mnew.X, mnew.Y, mnew.Z)

    # conditional declarations, clashes with local operators, invalid names, and quoted code
    old = raw"""
    using QuantumAlgebra
    if true
        @boson_ops b
    end
    let
        @fermion_ops d
        global W = d(1)*Pc"d_2"
    end
    for x in 1:2
        ex = Expr(:(=), :y, x)
    end
    V = b(1)*Pr"f(t)"*Pc"\bar{b}"*Pr"κ_k"
    """
    new, msgs = QAMigration.migrate_source(old)
    @test new == raw"""
    using QuantumAlgebra
    # operators and parameters (migrated from the legacy QuantumAlgebra interface)
    @qparam κ[:]
    if true
        @boson_ops b
    end
    let
        @fermion d[:]
        global W = d[1]*Pc"d_2"
    end
    for x in 1:2
        ex = Expr(:(=), :y, x)
    end
    V = b(1)*Pr"f(t)"*Pc"\bar{b}"*κ[:k]
    """
    @test any(m -> occursin("`b` (declared at code:3): declared conditionally", m), msgs)
    @test any(m -> occursin("parameter `d`: the new variable `d` would clash", m), msgs)
    @test any(m -> occursin("parameter `f(t)`: `f(t)` is not a valid variable name", m), msgs)
    @test any(m -> occursin("is not a valid variable name", m) && occursin("bar{b}", m), msgs)
    mold, mnew = Module(), Module()
    Base.include_string(mold, old)
    Base.include_string(mnew, new)
    @test (mold.V, mold.W) == (mnew.V, mnew.W)

    # `x = param(:x)` becomes the declaration of `x`, and indices can be variables or tuples
    old = raw"""
    using QuantumAlgebra
    α = param(:α,'n')
    gs(s) = param(:g,'r',s) + param(:κ,'r',(:i,:j))
    V = α*α'*gs(:i)*Pr"κ_k,l"
    """
    new, msgs = QAMigration.migrate_source(old)
    @test new == raw"""
    using QuantumAlgebra
    # operators and parameters (migrated from the legacy QuantumAlgebra interface)
    @qparam g[:]
    @qparam κ[:,:]
    @qparam α::Complex
    gs(s) = g[s] + κ[:i,:j]
    V = α*α'*gs(:i)*κ[:k,:l]
    """
    @test all(m -> startswith(m, "migrated"), msgs)
    mold, mnew = Module(), Module()
    Base.include_string(mold, old)
    Base.include_string(mnew, new)
    @test mold.V == mnew.V

    # a file is migrated together with the files it includes (with the includes pointing to the migrated files)
    mktempdir() do dir
        write(joinpath(dir, "main.jl"), "using QuantumAlgebra\ninclude(\"part.jl\")\nY = X + a(:j)\n")
        write(joinpath(dir, "part.jl"), "X = a(:i)*a'(:i)\n")
        msgs = QAMigration.migrate_file(joinpath(dir, "main.jl"))
        @test read(joinpath(dir, "main_migrated.jl"), String) ==
            "using QuantumAlgebra\n# operators and parameters (migrated from the legacy QuantumAlgebra interface)\n@boson a[:]\ninclude(\"part_migrated.jl\")\nY = X + a[:j]\n"
        @test read(joinpath(dir, "part_migrated.jl"), String) == "X = a[:i]*a'[:i]\n"
        run(file) = (m = Module(); Core.eval(m, :(include(p) = Base.include($m, p))); Base.include(m, file))
        @test run(joinpath(dir, "main.jl")) == run(joinpath(dir, "main_migrated.jl"))
        # operators used with different numbers of indices in different files are not migrated
        write(joinpath(dir, "part.jl"), "X = a()*a'(:i)\n")
        msgs = QAMigration.migrate_file(joinpath(dir, "main.jl"); inplace=true)
        @test any(m -> occursin("used with different numbers of indices (0, 1)", m), msgs)
        @test read(joinpath(dir, "main.jl"), String) == "using QuantumAlgebra\ninclude(\"part.jl\")\nY = X + a(:j)\n"
    end

    # new variables must not shadow names exported by QuantumAlgebra
    @test QAMigration.QA_EXPORTS == Set(filter(n -> !startswith(string(n), "@"), names(QuantumAlgebra)))

    # notebooks: the cells are treated as one program
    nb = """
    {"cells": [{"cell_type": "code", "metadata": {}, "outputs": [], "source": ["using QuantumAlgebra\\n", "QuantumAlgebra.use_σpm(true)"]},
               {"cell_type": "markdown", "metadata": {}, "source": ["text"]},
               {"cell_type": "code", "metadata": {}, "outputs": [], "source": ["H = Pr\\"ω\\"*a'(:i)*a(:i) + σp(:i)\\n", "normal_form(H)"]}],
     "metadata": {}, "nbformat": 4, "nbformat_minor": 5}
    """
    newnb, msgs = QAMigration.migrate_notebook(nb)
    cells = QAMigration.JSON.parse(newnb)["cells"]
    @test join(cells[1]["source"]) == "using QuantumAlgebra\n# operators and parameters (migrated from the legacy QuantumAlgebra interface)\n@boson a[:]\n@tls σ[:]\n@qparam ω\nQuantumAlgebra.use_σpm(true)"
    @test join(cells[3]["source"]) == "H = ω*a'[:i]*a[:i] + σ'[:i]\nnormal_form(H)"
    @test cells[2]["source"] == ["text"]
    QuantumAlgebra.use_σxyz()
end
