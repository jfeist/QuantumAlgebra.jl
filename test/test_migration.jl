# the tool for migrating code from the legacy interface (tools/migration/migrate.jl)
include(joinpath(@__DIR__, "..", "tools", "migration", "migrate.jl"))

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
    @test any(m -> occursin("`b` is used with different numbers of indices (0, 1)", m), msgs)
    @test any(m -> occursin("`a` is used as a value", m), msgs)
    @test any(m -> occursin("`g` would clash", m), msgs)

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
