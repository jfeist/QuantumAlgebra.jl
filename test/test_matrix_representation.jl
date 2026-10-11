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
# spins without fixed j are represented with spin _ABSTRACT_SPIN_J[] (their algebra holds in any representation)
const _ABSTRACT_SPIN_J = Ref{Rational{Int}}(1)
_spinj(O::QA.BaseOperator) = something(QA.spin_j(O), _ABSTRACT_SPIN_J[])
function _modedim(O::QA.BaseOperator)
    fam = QA.opfamily(O)
    fam == QA.QuditFamily ? QA.qudit_levels(O)[1] : fam in (QA.TLSFamily, QA.FermionFamily) ? 2 :
    fam == QA.SpinFamily ? Int(2*_spinj(O)+1) : error("no matrix representation for $O")
end
# local matrix of the operator (level 1 is the ground state / vacuum)
function _localmat(O::QA.BaseOperator)
    t = QA.optype(O)
    fam = QA.opfamily(O)
    if fam == QA.QuditFamily
        d, i, j = QA.qudit_levels(O)
        _E(d,i,j)
    elseif fam == QA.SpinFamily
        # basis states m = -j,...,j (level 1 is the lowest-weight state)
        j = _spinj(O)
        ms = -j:j
        t == QA.SpinZ_ ? diagm(ComplexF64.(ms)) :
        t == QA.SpinPlus_ ? sum(sqrt(j*(j+1)-m*(m+1))*_E(length(ms),k+1,k) for (k,m) in enumerate(ms[1:end-1])) :
                            sum(sqrt(j*(j+1)-m*(m+1))*_E(length(ms),k,k+1) for (k,m) in enumerate(ms[1:end-1]))
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

# (arbitrary) numerical values for parameters with integer indices
function _paramval(p::QA.Param)
    x = 1.3 + 0.37*length(string(QA.opname(p))) + sum(i -> 0.731*i.num, QA.inds(p); init=0.0)
    v = QA.parstate(p) == 'r' ? complex(x) : complex(x, 0.41*x+0.2)
    QA.parstate(p) == 'c' ? conj(v) : v
end

function matrep(A::QuExpr, modes::Vector)
    D = prod(last, modes; init=1)
    M = zeros(ComplexF64, D, D)
    for (t,s) in A.terms
        @assert isempty(QA.expvals(t)) && isempty(QA.corrs(t)) && QA.nsuminds(t) == 0
        @assert all(d -> QA.isintindex(d.iA) && QA.isintindex(d.iB), QA.deltas(t))
        all(d -> d.iA == d.iB, QA.deltas(t)) || continue
        pv = prod(_paramval, QA.params(t); init=complex(1.0))
        M .+= (s*pv) .* foldl(*, (_fullmat(O,modes) for O in QA.bares(t)); init=Matrix{ComplexF64}(I,D,D))
    end
    M
end
# write out the sums in A explicitly, with sum indices taking the values in vals
function expand_sums(A::QuExpr, vals)
    B = QuExpr()
    for (t,s) in A.terms
        for assignment in Iterators.product(ntuple(_ -> vals, QA.nsuminds(t))...)
            f = QA.replace_inds((QA.sumindex(k) => QA.QuIndex(v) for (k,v) in enumerate(assignment))...)
            QA._add_sum_term!(B, f(t,0), s)
        end
    end
    B
end
function samematrix(A::QuExpr, B::QuExpr)
    modes = sort!(collect(merge(_modes(A), _modes(B))), by=first)
    # (with an absolute tolerance since the matrices can vanish up to rounding errors)
    isapprox(matrep(A, modes), matrep(B, modes); atol=1e-10)
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

    # operators acting on sites given by `sites` (two-level systems with jump operators σ±,
    # and if pauli == true, also other two-level systems with Pauli matrices τˣʸᶻ, ρˣʸᶻ)
    function oppool(sites, pauli)
        ops = QuExpr[]
        for s in sites
            append!(ops, [q(i,j,s) for i=1:3, j=1:3])
            append!(ops, [p(i,j,s) for i=1:2, j=1:2])
            append!(ops, [σp(s), σm(s)])
            pauli && append!(ops, [τx(s), τy(s), τz(s), ρx(s), ρz(s)])
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

    # sums (expanded explicitly with indices taking values 1 and 2), including the canonical numbering of sum indices
    @testset "sums" begin
        r = qudit_ops(:r, 3)
        sumpool(s) = QuExpr[f(s), f'(s), σp(s), σm(s), r(2,1,s), r(1,3,s), r(3,2,s), param(:g,'r',s)*f(s)]
        for _ = 1:150
            P1 = prod(rand(sumpool(:k), rand(1:2)))
            P2 = prod(rand([sumpool(:l); sumpool(1)], rand(1:3)))
            # P3 contains the free index l (that is not summed over)
            P3 = prod(rand([sumpool(:m); sumpool(:l)], rand(0:2)); init=QuExpr(1))
            X = ∑(:k, P1) * ∑(:l, P2) * ∑(:m, P3)
            Xn = normal_form(X)
            for lval in 1:2
                fl = QA.replace_inds(QA.QuIndex(:l) => QA.QuIndex(lval))
                @test samematrix(expand_sums(fl(Xn), 1:2), expand_sums(fl(X), 1:2))
            end
            # the order and nesting of sums does not matter
            @test Xn == normal_form(∑(:m, ∑(:k, P1 * ∑(:l, P2) * P3)))
        end
        # terms that only differ in the numbering of sum indices are combined
        @test iszero(normal_form(∑(:i,∑(:j,f(:i)*f(:j)))))
        @test iszero(normal_form(∑((:i,:j),Pr"g_i"*r(2,1,:i)*r(3,1,:j)) - ∑((:j,:i),Pr"g_i"*r(3,1,:j)*r(2,1,:i))))
    end

    @testset "spins (abstract: j=$jabs)" for jabs in (1, 3//2)
        _ABSTRACT_SPIN_J[] = jabs
        Jp, Jm, Jz = spin_ops(:J)
        Kp, Km, Kz = spin_ops(:K, 1//2)
        Lp, Lm, Lz = spin_ops(:L, 3//2)
        spinpool(s) = QuExpr[Jp(s), Jm(s), Jz(s), Kp(s), Km(s), Kz(s), Lp(s), Lm(s), Lz(s), f(s)]
        for _ = 1:200
            X = prod(rand(spinpool(1), rand(1:5)))
            @test samematrix(normal_form(X), X)
        end
        substs = [(:i=>1, :j=>2), (:i=>1, :j=>1)]
        for _ = 1:200
            X = prod(rand([spinpool(:i); spinpool(:j)], rand(2:4)))
            Xn = normal_form(X)
            for sub in substs
                fsub = QA.replace_inds((QA.QuIndex(k) => QA.QuIndex(v) for (k,v) in sub)...)
                @test samematrix(normal_form(fsub(Xn)), fsub(X))
            end
            A, B = (prod(rand([spinpool(:i); spinpool(:j)], rand(1:2))) for _=1:2)
            @test normal_form(normal_form(A*B)*X) == normal_form(A*normal_form(B*X))
        end
    end
    QA.use_σpm(false)
end
