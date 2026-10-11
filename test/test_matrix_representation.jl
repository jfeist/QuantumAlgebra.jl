# Check the algebraic rules of finite-dimensional operators (two-level systems, qudits, fermions)
# by representing expressions as explicit matrices on the tensor product space of all modes
# appearing in them, and comparing the matrices of random operator products before and after
# normal ordering.

using LinearAlgebra
using Random
const QA = QuantumAlgebra

_E(d,i,j) = (M = zeros(ComplexF64,d,d); M[i,j] = 1; M)

# a mode is identified by the operator name and its (integer) indices
_modekey(O::QA.BaseOperator) = (QA.sym(QA.opname(O)), Tuple(QA.inds(O)))
function _modedim(O::QA.BaseOperator)
    fam = QA.opfamily(O)
    fam == QA.QuditFamily ? QA.qudit_levels(O)[1] : fam in (QA.TLSFamily, QA.FermionFamily) ? 2 : error("no matrix representation for $O")
end
# local matrix of the operator (level 1 is the ground state / vacuum)
function _localmat(O::QA.BaseOperator)
    t = QA.optype(O)
    fam = QA.opfamily(O)
    if fam == QA.QuditFamily
        d, i, j = QA.qudit_levels(O)
        _E(d,i,j)
    elseif fam == QA.FermionFamily
        t == QA.FermionCreate_ ? _E(2,2,1) : _E(2,1,2)
    else
        t == QA.TLSCreate_ ? _E(2,2,1) : t == QA.TLSDestroy_ ? _E(2,1,2) :
        t == QA.TLSx_ ? _E(2,1,2) + _E(2,2,1) : t == QA.TLSy_ ? im*_E(2,1,2) - im*_E(2,2,1) : _E(2,2,2) - _E(2,1,1)
    end
end

_modes(A::QuExpr) = Dict(_modekey(O) => _modedim(O) for t in keys(A.terms) for O in QA.bares(t))

# matrix of operator O on the tensor product of modes (Jordan-Wigner strings for fermions of the same name)
function _fullmat(O::QA.BaseOperator, modes::Vector)
    key = _modekey(O)
    mats = map(modes) do (k,d)
        if k == key
            _localmat(O)
        elseif QA.opfamily(O) == QA.FermionFamily && k[1] == key[1] && k < key
            _E(2,1,1) - _E(2,2,2)
        else
            Matrix{ComplexF64}(I,d,d)
        end
    end
    foldl(kron, mats)
end

function matrep(A::QuExpr, modes::Vector)
    D = prod(last, modes; init=1)
    M = zeros(ComplexF64, D, D)
    for (t,s) in A.terms
        @assert isempty(QA.params(t)) && isempty(QA.expvals(t)) && isempty(QA.corrs(t)) && QA.nsuminds(t) == 0
        @assert all(d -> QA.isintindex(d.iA) && QA.isintindex(d.iB), QA.deltas(t))
        all(d -> d.iA == d.iB, QA.deltas(t)) || continue
        M .+= s .* foldl(*, (_fullmat(O,modes) for O in QA.bares(t)); init=Matrix{ComplexF64}(I,D,D))
    end
    M
end
function samematrix(A::QuExpr, B::QuExpr)
    modes = sort!(collect(merge(_modes(A), _modes(B))), by=first)
    matrep(A, modes) ≈ matrep(B, modes)
end

@testset "matrix representation" begin
    Random.seed!(1234)
    QA.use_σpm(true)
    q = qudit_ops(:q, 3)
    p = qudit_ops(:p, 2)
    @tlspm_ops σ
    @fermion_ops f
    # τ, ρ: Pauli matrices (note that tlsxyz_ops depends on use_σpm, so create them directly)
    τx, τy, τz = ((inds...) -> QuExpr(T(:τ, inds...)) for T in (QA.TLSx, QA.TLSy, QA.TLSz))
    ρx, ρy, ρz = ((inds...) -> QuExpr(T(:ρ, inds...)) for T in (QA.TLSx, QA.TLSy, QA.TLSz))

    # operators acting on sites given by `sites` (two-level systems either with
    # jump operators σ± or with Pauli matrices τˣʸᶻ, ρˣʸᶻ, since these cannot be mixed)
    function oppool(sites, pauli)
        ops = QuExpr[]
        for s in sites
            append!(ops, [q(i,j,s) for i=1:3, j=1:3])
            append!(ops, [p(i,j,s) for i=1:2, j=1:2])
            append!(ops, pauli ? [τx(s), τy(s), τz(s), ρx(s), ρz(s)] : [σp(s), σm(s)])
            append!(ops, [f(s), f'(s)])
        end
        ops
    end

    @testset "normal_form preserves operators (pauli=$pauli)" for pauli in (false,true)
        pool = oppool((1,2), pauli)
        for _ = 1:300
            X = prod(rand(pool, rand(1:6)))
            @test samematrix(normal_form(X), X)
        end
    end

    @testset "normal_form with symbolic indices (pauli=$pauli)" for pauli in (false,true)
        pool = oppool((:i,:j,1), pauli)
        substs = [(:i=>1, :j=>2), (:i=>1, :j=>1), (:i=>2, :j=>1), (:i=>2, :j=>2)]
        for _ = 1:500
            X = prod(rand(pool, rand(2:5)))
            Xn = normal_form(X)
            for sub in substs
                fsub = QA.replace_inds((QA.QuIndex(k) => QA.QuIndex(v) for (k,v) in sub)...)
                @test samematrix(normal_form(fsub(Xn)), fsub(X))
            end
        end
    end

    # the normal form should be unique, i.e., independent of the order in which products are
    # simplified (this does not hold for Pauli matrices with symbolic indices)
    @testset "unique normal form" begin
        pool = oppool((:i,:j,1), false)
        for _ = 1:300
            A, B, C = (prod(rand(pool, rand(1:2))) for _=1:3)
            @test normal_form(normal_form(A*B)*C) == normal_form(A*normal_form(B*C))
        end
    end
    QA.use_σpm(false)
end
