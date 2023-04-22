using QuantumAlgebra
using QuantumAlgebra: δ, QuExpr, QuTerm, BaseOpProduct, BaseOperator, Param, QuIndex, _map_quexpr_ops, TLSx, TLSCreate, is_normal_form
using Test, Documenter

function myδ(i,j)
    iA,iB = QuIndex.((i,j))
    QuExpr(QuTerm([δ(min(iA,iB),max(iA,iB))],BaseOpProduct()))
end
scal(x) = x*one(QuExpr)

macro test_is_normal_form(x)
    x = esc(x)
    quote
        xn = normal_form($x)
        @test is_normal_form($x) == ($x == xn)
        @test is_normal_form(xn)
    end
end

@time @testset "QuantumAlgebra.jl" begin
    doctest(QuantumAlgebra)

    @testset "variables are consts" begin
        for name in names(QuantumAlgebra,all=true)
            startswith(string(name),"##meta#") && continue
            try
                QuantumAlgebra.eval(:($name = 5))
            catch err
                @test err isa ErrorException && startswith(err.msg,"invalid redefinition of constant")
            end
        end
    end

    @test isbitstype(QuantumAlgebra.QuIndex)
    @test isbitstype(QuantumAlgebra.QuOpName)
    if isbitstype(QuantumAlgebra.QuIndices)
        @test isbitstype(QuantumAlgebra.BaseOperator)
        @test isbitstype(QuantumAlgebra.Param)
    end
    @test isbitstype(QuantumAlgebra.δ)

    @testset "auto_normal_form($auto_norm)" for auto_norm in (false,true)
        QuantumAlgebra.auto_normal_form(auto_norm)

        @test QuantumAlgebra.using_auto_normal_form() == auto_norm

        @boson_ops b
        @test normal_form(b(:i)*b'(:j)) == b'(:j)*b(:i) + myδ(:i,:j)
        @test normal_form(b(:i)*a'(:j)) == a'(:j)*b(:i)

        @fermion_ops c
        @test normal_form(c(:i)*c'(:j) + c'(:j)*c(:i)) == myδ(:i,:j)
        # different fermionic species commute by default
        @test normal_form(c(:i)*f'(:j)) == f'(:j)*c(:i)

        # this creates two fermionic operators belonging to the same species (which anticommute, e.g., for itinerant and localized electrons)
        # NB: sort order is the order in which names are given (not alphabetic like normally)
        @anticommuting_fermion_group e d
        # fermionic operators of the same species anticommute
        @test iszero(normal_form(d()*e() + e()*d()))
        # ...but different fermionic species commute
        @test normal_form(f()*e() + e()*f()) == 2*e()*f()

        # all equal numbers should be equal scalars (ignore type)
        @test a() + 0 == a() + 0.0
        @test a() + 1 == a() + (1//1 + 0im)
        @test scal(1.5) == scal(3//2)

        # test params and their complex conjugation
        @test Pc"ω" == param(:ω,'n')
        @test Pc"ω_i" == param(:ω,'n',:i)
        @test Pc"ω_i,j,2" == param(:ω,'n',:i,:j,2)
        @test_throws ArgumentError param(:ω,'g')
        @test_throws ArgumentError param(:ω,'n',2,:i,"ag")
        @test_throws ArgumentError a("aa")

        @test adjoint(param(:g,'r')) == param(:g,'r')
        @test adjoint(param(:g,'n')) == param(:g,'c')
        @test adjoint(param(:g,'n',:i)) == param(:g,'c',:i)
        @test (Pc"g_i")' == param(:g,'c',:i)
        @test adjoint(param(:g,'c',:i,1,:m_3)) == param(:g,'n',:i,1,:m_3)

        tmp1 = param(:g,'n',2)*param(:g,'n',1)*param(:g,'c',3)*param(:b,'n')*param(:a,'n')*param(:d,'c',3)*param(:f,'r',1,:i)
        tmpc = param(:g,'c',2)*param(:g,'c',1)*param(:g,'n',3)*param(:b,'c')*param(:a,'c')*param(:d,'n',3)*param(:f,'r',1,:i)
        tmp2 = Pc"a"*Pc"b"*Pc"g_1"*Pc"g_2"*(Pc"g_3")' * (Pc"d_3")' * Pr"f_1,i"
        @test normal_form(tmp1) == normal_form(tmp2)
        @test normal_form(tmp1') == normal_form(tmpc)

        @test a()*5 == 5*a()
        @test a()+5 == 5+a()
        @test 5-a() == -a() + 5
        @test a()-5 == -5 + a()

        @test a()' == a'()
        @test normal_form(a(:m)*a'(:m)) == a'(:m)*a(:m) + 1
        @test normal_form((a(:m)*a(1))') == a'(1)*a'(:m)
        @test normal_form((a(1,2,:k)*a(:m))') == normal_form(a'(1,2,:k)*a'(:m))

        @test normal_form(f'(:a)*f(:b) + f(:b)*f'(:a)) == myδ(:a,:b)
        @test iszero(normal_form(f(:a)*f(:b) + f(:b)*f(:a)))
        @test iszero(normal_form(f'(:a)*f'(:b) + f'(:b)*f'(:a)))

        @test f()' == f'()
        @test f'()' == f()
        @test (f'()*f())' == f'()*f()
        @test comm(f'(),f()) == -comm(f(),f'())

        @test ∑(:j,a(:j))*∑(:i,a(:i)) == ∑(:i,∑(:j,a(:i)*a(:j)))
        @test ∑(:i,∑(:j,a(:i)*a(:j))) == ∑((:i,:j),a(:i)*a(:j))
        @test ∑(:i,∑(:j,∑(:k,a(:i)*a(:j)*a(:k)))) == ∑((:i,:j,:k),a(:i)*a(:j)*a(:k))

        tmp = ∑(:i,a'(:i)*a(:i))
        @test tmp' == tmp
        @test normal_form(a(:i)*∑(:i,a(:i))) == normal_form(∑(:i_1,a(:i_1)*a(:i)))
        @test normal_form(a'(:n)*tmp) == normal_form(∑(:i,a'(:n)*a'(:i)*a(:i)))
        @test normal_form(a(:n) *tmp) == normal_form(∑(:i,a'(:i)*a(:n)*a(:i)) + a(:n))
        @test normal_form(param(:g,'n',:i)*∑(:i,a(:i))) == normal_form(∑(:i_1,param(:g,'n',:i)*a(:i_1)))
        @test normal_form(param(:g,'n',:i_1)*a(:i)*∑(:i,a(:i))) == normal_form(∑(:i_2,param(:g,'n',:i_1)*a(:i)*a(:i_2)))
        @test normal_form(param(:g,'n',:n)*∑(:i,a(:i))) == normal_form(∑(:i,param(:g,'n',:n)*a(:i)))

        @test Pr"gz_i,μ" == param(:gz,'r',(:i,:μ))
        @test Pc"gz_i,μ" == param(:gz,'n',(:i,:μ))
        @test (Pc"gz_i,μ")' == param(:gz,'c',(:i,:μ))
        @test Pr"α_i_1,μ_2" == Pr"α_{i_1,μ_2}" == param(:α,'r',:i_1,:μ_2)
        @test Pr"α_◔,◔" == Pr"α_{◔,◔}" == param(:α,'r',:◔,:◔)

        @testset "index_functions" begin
            x = QuantumAlgebra.canon_inds_remember()
            if auto_norm
                @test QuantumAlgebra.canon_inds()(f(5)*f(2)*a(9,:k,:m) + a(10)) == -f(:i)*f(:j)*a(:k,:l,:m) + a(:i)
            else
                @test QuantumAlgebra.canon_inds()(f(5)*f(2)*a(9,:k,:m) + a(10)) == f(:i)*f(:j)*a(:k,:l,:m) + a(:i)
            end
            t1, _ = only((f(5)*f(2)*a(9,:k,:m)).terms)
            t2, _ = only((f(:i)*f(:j)*a(:k,:l,:m)).terms)
            indf = QuantumAlgebra.canon_inds_remember()
            @test indf(t1) == t2
            @test indf.replacements == Dict(QuIndex(:i) => auto_norm ? QuIndex(2) : QuIndex(5),
                                            QuIndex(:j) => auto_norm ? QuIndex(5) : QuIndex(2),
                                            QuIndex(:k) => QuIndex(9), QuIndex(:l) => QuIndex(:k),
                                            QuIndex(:m) => QuIndex(:m))
            @test_throws ErrorException indf(a())
    end

        @testset "normal_form with sums" begin
            # bug report from FJ Matute
            t1 = Pc"a_w"*Pc"B_x"*Pc"C_y"*Pc"D_z"*f(:w)*f(:x)'*f(:y)*f(:z)'
            X1 = ∑(:w,∑(:x,∑(:y,∑(:z,t1))))
            X2 = ∑(:w,∑(:x,∑(:y,∑(:z,normal_form(t1)))))
            @test normal_form(X1) == normal_form(X2)

            s = a(:i,:j,:k)*a'(:j,:k,:i)
            X1 = ∑((:i,:j),s)
            X2 = ∑((:i,:j),normal_form(s))
            @test normal_form(X1) == normal_form(X2)
        end

        @testset "use_σpm($with_σpm)" for with_σpm in (false,true)
            QuantumAlgebra.use_σpm(with_σpm)

            @test QuantumAlgebra.using_σpm() == with_σpm

            @test_throws ArgumentError σx(:i,"aa")
            @test σx(:i,:b) == σx((:i,:b))
            @test σx(:i,:j) != σx(:i,:b)

            α = param(:α)
            tmp = α' * a(1) + α * a(1)
            @test iszero(tmp - tmp)

            @test σx() == σp() + σm()
            @test σx(:i) == σp(:i) + σm(:i)
            @test σx(:i,:j) == σp((:i,:j)) + σm((:i,:j))

            @test normal_form(σy(:i)) == normal_form(-1im*(σp(:i) - σm(:i)))
            @test normal_form(σz(:i)) == normal_form(2*σp(:i)*σm(:i) - 1)

            for s=(σx,σy,σz)
                @test isone(normal_form(s(:i)*s(:i)))
                @test s(:i)' == s(:i)
            end
            @test normal_form(σx(:i)*σy(:i)) == 1im*σz(:i)
            @test normal_form(σy(:i)*σx(:i)) == -1im*σz(:i)
            @test normal_form(σy()*σz()) == 1im*σx()
            @test normal_form(σx(:i,:j)*σz(:i,:j)) == -1im*σy(:i,:j)
            @test normal_form(σx(:i,:j)*σz(:i,:j)) != -1im*σy(:i,:k)

            @test normal_form(a(1) * (σy(1) * a(1))') == normal_form(a(1) * (a'(1) * σy(1))) == normal_form(a'(1)*a(1)*σy(1) + σy(1))

            @test normal_form(σy(:i)*σy(:i)*σy(:i)*σy(:i)*σy(:i)*σy(:i)*σx(:i)) == σx(:i)
            @test normal_form(σy(:i)*σy(:i)*σy(:i)*σy(:i)*σy(:i)*σy(:i)*σx(:j)) == σx(:j)

            @test_throws ArgumentError normal_form(QuExpr(TLSx(:σ)) * QuExpr(TLSCreate(:σ)))

            tmp1 = normal_form(3*Pc"ω"*Pc"g"*expval(σp(:k))*σp(:k)*a'(5)*a(5))
            if QuantumAlgebra.using_σpm()
                @test length(tmp1.terms) == 1
                @test first(tmp1.terms)[1].bares == first((a'(5)*σp(:k)*a(5)).terms)[1].bares
            else
                @test length(tmp1.terms) == 4
            end

            x1 = ∑(:i,∑(:j,a(:i)*a(:j))) + f(:j,:k,:l)*f'(:j,:i,:m) + a(3,8,:i)*expval(a(:i)*a'(:k))*a'(:α,:μ,:ν)
            x2 = σp(:k) + σx(:m)
            @test normal_form(x1 + x2) == normal_form(x1) + normal_form(x2)
            @test normal_form(x1 * x2) == normal_form(normal_form(x1) * normal_form(x2))
            @test normal_form(x1 * x2 + x2) == normal_form(normal_form(x1) * normal_form(x2)) + normal_form(x2)

            @testset "is_normal_form" begin
                y = one(QuExpr)
                for x in (σx(:i)*σy(:i), a(1) * (σy(1) * a(1))', a(:d)*a'(:c), a(:i)*a'(:j)*σz(:α)*σy(:β)*σx(:α), expval(a'(:c)*a(:d)))
                    @test_is_normal_form x
                    y *= x
                    @test_is_normal_form y
                    @test_is_normal_form ∑(:c,∑(:i,∑(:α,y)))
                end
            end

            @testset "exponent" begin
                for x in (σx(:i)*σy(:i), a(1) * (σy(1) * a(1))', a(:d)*a'(:c), a(:i)*a'(:j)*σz(:α)*σy(:β)*σx(:α), expval(a'(:c)*a(:d)))
                    @test isone(x^0)
                    @test x^1 == x
                    @test x^2 == x*x
                    @test x^3 == x*x*x
                    @test x^4 == x*x*x*x
                    @test_throws ArgumentError x^-1
                    @test_throws ArgumentError x^-2
                end
            end

            @testset "Commutation inside ExpVal/Corr" begin
                for (x1,x2,s) in ((a(:d)*a'(:c), a(:i)*a'(:j)*σz(:α)*σy(:β)*σx(:α),  expval(a'(:c)*a(:d))),
                                (f(:d)*f'(:c), f(:i)*f'(:j)*σz(:α)*σy(:β)*σx(:α), -expval(f'(:c)*f(:d))))
                    x = ∑(:c,∑(:i,∑(:α,3*expval(x1)*x2)))
                    @test normal_form(expval(x1)) == myδ(:c,:d) + s
                    @test normal_form(3*expval(x1)) == expval(normal_form(3*x1))
                    @test normal_form(expval(x1)*expval(x1)) == normal_form(normal_form(expval(x1))*normal_form(expval(x1)))
                    @test normal_form(3*expval(x2)) == expval(normal_form(3*x2))
                    @test normal_form(3*expval(x1)*x2) == expval(3*normal_form(x1))*normal_form(x2)
                    @test normal_form(3*expval(x1)*corr(x2)) == normal_form(expval(3*normal_form(x1))*corr(normal_form(x2)))
                    @test normal_form(x) == normal_form(∑(:α,∑(:c,∑(:i,3*expval(normal_form(x1))*normal_form(x2)))))
                    xn = normal_form(x)
                    xsqn = normal_form(x*x)
                    @test xsqn == normal_form(xn*xn)

                    # do normal_form already on the first partial product, otherwise computation time explodes.
                    # NOTE: this test is broken not because of an error as far as I can tell, but because
                    # depending on the order of computation it produces equivalent terms that we cannot
                    # identify as equal yet, such as:
                    # ∑₁₂ σ⁺(#₁) σ⁺(#₂) σ⁺(β) σ⁻(#₁) and
                    # ∑₁₂ σ⁺(#₁) σ⁺(#₂) σ⁺(β) σ⁻(#₂)
                    #
                    # For use_σxyz, this is even worse, since even though non-adjacent operators are contracted,
                    # there are ambiguities in sums that even lead to different numbers of terms
                    if auto_norm
                        @test normal_form(xsqn*x) == normal_form(xn*xn*xn)
                    else
                        @test_broken normal_form(xsqn*x) == normal_form(xn*xn*xn)
                    end
                end
            end

            @test normal_form(myδ(:i,:k)*a(:k)) == normal_form(a(:i)*myδ(:k,:i))
            @test normal_form(myδ(:i,:k)*a(:k,:i)) == normal_form(a(:i,:i)*myδ(:k,:i))
            @test normal_form(myδ(:i,:k)*myδ(:i,:j)) == normal_form(myδ(:k,:i)*myδ(:j,:k))
            # k cannot be equal to 1 and 3 at the same time
            @test iszero(normal_form(myδ(1,:k)*myδ(:k,3)*σx(:k)))
            @test normal_form(myδ(1,:k)*myδ(:k,1)*σx(:k)) == myδ(1,:k)*σx(1)

            @test iszero(normal_form(comm(σx(5),σy(3))))
            @test iszero(normal_form(comm(σx(5),σx(5))))
            @test normal_form(comm(σx(1),σz(1))) == normal_form(-2im*σy(1))
            @test normal_form(comm(σx(:μ),σy(:ν))) == normal_form(2im*myδ(:μ,:ν)*σz(:ν))
            @test normal_form(1//2im * comm(σx(:m),σy(:m))) == normal_form(σz(:m))
            @test normal_form(σx(:a)*σy(:a)*σz(:a)) == 1im*one(QuExpr)

            @test iszero(normal_form(comm(Pc"g",a(5)+a(3))))
            @test iszero(normal_form(comm(Pc"g",a(5)*a(3))))
            @test iszero(normal_form(comm(a(5)+a(3),Pc"g")))
            @test iszero(normal_form(comm(a(5)*a(3),Pc"g")))

            @test normal_form(comm(a(5)+a(3),a'(5))) == normal_form(comm(a(5),a'(5))+ comm(a(3),a'(5)))
            @test normal_form(comm(a(5)+a(3),a'(5)*a(3))) == normal_form(comm(a(5),a'(5)*a(3))+ comm(a(3),a'(5)*a(3)))

            @test normal_form(a'(2)*σy(:i) - 14Pc"ω") == normal_form(-14Pc"ω" + σy(:i)*a'(2))
            @test normal_form(a'(2)*σy(:i)) == normal_form(σy(:i)*a'(2))

            @test normal_form(σp(1) * σm(1)) == normal_form(1//2*σz(1) + 1//2)
            @test normal_form(σm(1)*σp(1)) == normal_form(1 - σp(1)*σm(1))
            @test normal_form(σm(1)*σp(1)) == normal_form(σp(1)*σm(1) + comm(σm(1),σp(1)))
            @test normal_form(σz(1)) == normal_form(σp(1)*σm(1) - σm(1)*σp(1))

            @test iszero(comm(σp(1),σp(1)))
            @test normal_form(comm(σp(:n),σm(:n))) == normal_form(σz(:n))

            @test normal_form(comm(a(1),   a'(1)*a(1))) == a(1)
            @test normal_form(comm(a'(1),a'(1)*a(1))) == -a'(1)

            @test expval(scal(3)) == scal(3)
            @test corr(scal(3)) == scal(3)
            @test corr(3+a'(2)) == 3 + corr(a'(2))

            @test expval_as_corrs(scal(3)) == scal(3)
            @test expval_as_corrs(a(2)) == corr(a(2))
            @test expval_as_corrs(3Pc"g") == 3Pc"g"
            @test expval_as_corrs(3a(2)) == 3corr(a(2))
            @test expval_as_corrs(3+a'(2)) == 3 + corr(a'(2))

            @test expval_as_corrs(a(2)*a(2)) == corr(a(2)*a(2)) + corr(a(2))*corr(a(2))
            @test expval_as_corrs(a(2)*a(:m)) == corr(a(2)*a(:m)) + corr(a(2))*corr(a(:m))
            tmpas = a.(1:3)
            @test *(tmpas...) == a(1)*a(2)*a(3)
            tmpEVs = corr.(tmpas)
            @test expval_as_corrs(*(8,Pr"g",tmpas...)) == 8*Pr"g" * (corr(*(tmpas...)) + *(tmpEVs...) + tmpEVs[1]*corr(tmpas[2]*tmpas[3]) + tmpEVs[2]*corr(tmpas[1]*tmpas[3]) + tmpEVs[3]*corr(tmpas[1]*tmpas[2]))

            if QuantumAlgebra.using_σpm()
                @test expval_as_corrs(-param(:g,'r',1)*σp(1)) == -param(:g,'r',1)*corr(σp(1))
                @test expval_as_corrs(∑(:i,σp(:i)*σm(:n))) == ∑(:i,corr(σp(:i)*σm(:n))) + ∑(:i,corr(σp(:i))*corr(σm(:n)))
                @test corr_as_expvals(∑(:i,σp(:i)*σp(:j))) == ∑(:i,expval(σp(:i)*σp(:j))) - ∑(:i,expval(σp(:i))*expval(σp(:j)))
            else
                @test expval_as_corrs(-param(:g,'r',1)*σz(1)) == -param(:g,'r',1)*corr(σz(1))
                @test expval_as_corrs(∑(:i,σy(:i)*σy(:n))) == ∑(:i,corr(σy(:i)*σy(:n))) + ∑(:i,corr(σy(:i))*corr(σy(:n))) - corr(σy(:n))*corr(σy(:n))
                @test corr_as_expvals(∑(:i,σx(:i)*σy(:j))) == ∑(:i,expval(σx(:i)*σy(:j))) - ∑(:i,expval(σx(:i))*expval(σy(:j))) + expval(σx(:j))*expval(σy(:j))
            end

            @test corr_as_expvals(a'(:i)a(:j)) == expval(a'(:i)a(:j)) - expval(a'(:i))*expval(a(:j))
            tmpEVs = expval.(tmpas)
            @test corr_as_expvals(*(tmpas...)) == expval(*(tmpas...)) + 2 * *(tmpEVs...) - tmpEVs[1]*expval(tmpas[2]*tmpas[3]) - tmpEVs[2]*expval(tmpas[1]*tmpas[3]) - tmpEVs[3]*expval(tmpas[1]*tmpas[2])

            H = ∑(:i,param(:ω,'r',:i)*a'(:i)*a(:i))
            @test normal_form(comm(a(:i),H)) == param(:ω,'r',:i)*a(:i)
            @test normal_form(comm(a(:n),H)) == param(:ω,'r',:n)*a(:n)
            @test normal_form(comm(H,a(:n))) == -param(:ω,'r',:n)*a(:n)
            @test normal_form(comm(a'(:n),H)) == -param(:ω,'r',:n)*a'(:n)
            @test normal_form(comm(a'(:n)*a(:m),H)) == (param(:ω,'r',:m)-param(:ω,'r',:n))*a'(:n)*a(:m)

            @test normal_form(a()*H) == normal_form(∑(:i,param(:ω,'r',:i)*a'(:i)*a(:i)*a()))
            @test normal_form(a(:k)*H) == normal_form(param(:ω,'r',:k)*a(:k) + ∑(:i,param(:ω,'r',:i)*a'(:i)*a(:i)*a(:k)))
            HH = ∑(:i,param(:ω,'r',:i,:i)*a'(:i,:i)*a(:i,:i))
            @test normal_form(a(:k,:k)*HH) == normal_form(param(:ω,'r',:k,:k)*a(:k,:k) + ∑(:i,param(:ω,'r',:i,:i)*a'(:i,:i)*a(:i,:i)*a(:k,:k)))

            @test iszero(Avac(H))
            @test iszero(vacA(H))
            @test iszero(vacA(a'(3)*σp(1)*σm(1)))
            @test iszero(vacA(f'(:n)))
            @test vacA(f(:n)) == f(:n)
            @test vacA(a(:n)) == a(:n)
            @test Avac(f'(:n)) == f'(:n)
            @test iszero(Avac(f(:n)))
            @test iszero(Avac(a(3)*σp(1)*σm(1)))
            @test isone(Avac(σm(1)*σp(1)))
            @test iszero(Avac(σp(1)*σm(1)))
            @test Avac(σp(1)) == σp(1)
            @test vacA(σm(1)) == σm(1)
            if QuantumAlgebra.using_σpm()
                @test iszero(Avac(σm(1)))
                @test iszero(vacA(σp(1)))
            else
                @test Avac(σx(1)) == vacA(σx(1)) == σx(1)
            end
            @test iszero(vacExpVal(σx(1)))
            @test iszero(vacExpVal(σp(1)))
            @test iszero(vacExpVal(∑(:i,σp(:i))))
            @test iszero(vacExpVal(σp(:i)*σm(:k)))
            @test vacExpVal(σm(:j)*σp(:l)) == myδ(:j,:l)
            if QuantumAlgebra.using_σpm()
                @test normal_form(σp(:i)*σm(:k)) == σp(:i)*σm(:k)
            end
            @test normal_form(a(:n)*a'(:n)*a(:n)*a'(:n)) == scal(1) + 3a'(:n)*a(:n) + a'(:n)*a'(:n)*a(:n)*a(:n)

            S = (1/√(2*6))*a'(:n)*a'(:n)*a'(:n) + (1/√2)*a'(:m)
            for (A,val) in [(scal(1),scal(1)),
                            (a'(:n)*a(:n),scal(1.5) + 0.5 * myδ(:n,:m)),
                            (a'(:n)*a'(:n)*a(:n)*a(:n),scal(3))]
                @test vacExpVal(A,S) ≈ val
            end

            H = Pr"ωc"*a'()*a() + Pr"ωe"*σz() + Pr"g"*σx()*(a()+a'())
            Ls = ((Pr"κ",a()),(Pr"γ",σm()))
            @test normal_form(heisenberg_eom(a(),H,Ls)) == -1im*Pr"g"*σx() - 1im*Pr"ωc"*a() - 1//2*Pr"κ"*a()
            @test normal_form(heisenberg_eom(σz(),H,Ls)) == -Pr"γ"*(1+σz()) + 2Pr"g"*(σy()*a()+a'()*σy())

            H = ∑(:i,∑(:j,Pr"ω_i,j"*a'(:i)*a(:j)))
            Ls = ((:i,Pr"κ_i",a(:i)),)
            @test normal_form(heisenberg_eom(a(:i),H,Ls)) == -1im*∑(:j,Pr"ω_i,j"*a(:j)) - 1//2*Pr"κ_i"*a(:i)

            H = ∑((:i,:j,:k,:l),Pr"ω_i,j,k,l"*a'(:i,:j)*a(:k,:l))
            Ls = (((:k,:l),Pr"κ_k,l",a(:k,:l)),)
            @test normal_form(heisenberg_eom(a(:i,:j),H,Ls)) == -1im*∑((:k,:l),Pr"ω_i,j,k,l"*a(:k,:l)) - 1//2*Pr"κ_i,j"*a(:i,:j)

            tmp = ∑(:i,expval_as_corrs(a'(:n)*a(:i)))
            tmplatex = raw"\sum_{\#_{1}} \langle {a}_{n}^\dagger\rangle_{c} \langle {a}_{\#_{1}}\rangle_{c}  + \sum_{\#_{1}} \langle {a}_{n}^\dagger {a}_{\#_{1}}\rangle_{c} "
            @test latex(tmp) == tmplatex
            @test expval_as_corrs(tmp) == tmp
            @test sprint(show,"text/latex",tmp) == "\$$(tmplatex)\$"

            @test latex(a(:i_1,:i_2)) == raw"{a}_{i_{1}i_{2}}"

            @test latex(normal_form(e(:i)*d'(:j))) == " - {d}_{j}^\\dagger {e}_{i}"
            @test latex(normal_form(f(:i)*d'(:j))) == "{d}_{j}^\\dagger {f}_{i}"

            lσp = latex(σp())
            lσz = latex(σz())
            if QuantumAlgebra.using_σpm()
                @test lσp == "{{\\sigma}}^+"
                @test lσz == " - 1 + 2 {{\\sigma}}^+ {{\\sigma}}^-"
            else
                @test lσp == "\\frac{1}{2} {{\\sigma}}^x + \\frac{1}{2}\\mathit{i} {{\\sigma}}^y"
                @test lσz == "{{\\sigma}}^z"
            end

            @test QuantumAlgebra.symmetric_index_nums(a'(:i)*a'(:j)*a(:k)*a(:l)) == [2,2]

            @test string(normal_form(f(:i) * d'(:j))) == "d†(j) f(i)"
            @test string(normal_form(e(:i) * d'(:j))) == "-d†(j) e(i)"

            @test string(normal_form(∑(:i,a(:i)) * a'(:n))) == "1 + ∑₁ a†(n) a(#₁)"
            @test string(5a(:i) + 3a'(:j) - 3f(1)) == "3 a†(j) - 3 f(1) + 5 a(i)"
            @test string(5im*a(:i) + (3+2im)*a'(:j) - 3//2*f(1)) == "(3+2i) a†(j) - 3//2 f(1) + 5i a(i)"
            if QuantumAlgebra.using_σpm()
                @test string(normal_form(a(5)*a'(5)*σp(3)*expval_as_corrs(a'(5,:i)*a(5)))) == "⟨a†(5i)⟩c ⟨a(5)⟩c σ⁺(3) + ⟨a†(5i) a(5)⟩c σ⁺(3) + ⟨a†(5i)⟩c ⟨a(5)⟩c a†(5) σ⁺(3) a(5) + ⟨a†(5i) a(5)⟩c a†(5) σ⁺(3) a(5)"
            else
                @test string(normal_form(a(5)*a'(5)*σz(3)*expval_as_corrs(a'(5,:i)*a(5)))) == "⟨a†(5i)⟩c ⟨a(5)⟩c σᶻ(3) + ⟨a†(5i) a(5)⟩c σᶻ(3) + ⟨a†(5i)⟩c ⟨a(5)⟩c a†(5) σᶻ(3) a(5) + ⟨a†(5i) a(5)⟩c a†(5) σᶻ(3) a(5)"
            end

            @test julia_expression(QuExpr()) == 0
            # an empty QuTerm is the identity operator!
            @test julia_expression(QuTerm()) == 1

            @test julia_expression(expval(a'(:i)a(:j)),[:aᴴa]) == :(aᴴa[i, j])
            @test_throws ArgumentError julia_expression(expval(a'(:i)a(:j)),[:aa])

            x = ∑(:i,Pc"g_i,k"*a(:i,:j_2,:K)*a'(:i_1,:J,:k)*σp(:i))
            ex = julia_expression(expval(normal_form(x)))
            if QuantumAlgebra.using_σpm()
                @test ex == :(I[J, j₂] * I[K, k] * g[i₁, K] * σ⁺[i₁] + g[s̄₁, k] * aᴴσ⁺a[i₁, J, k, s̄₁, s̄₁, j₂, K])
            else
                # julia_expression interpolates complex numbers directly as complex, not as expressions. so make sure to do the same here
                halfim = 0.5im
                @test ex == :(0.5 * I[J, j₂] * I[K, k] * g[i₁, K] * σˣ[i₁] + $halfim * I[J, j₂] * I[K, k] * g[i₁, K] * σʸ[i₁] + 0.5 * g[s̄₁, k] * aᴴσˣa[i₁, J, k, s̄₁, s̄₁, j₂, K] + $halfim * g[s̄₁, k] * aᴴσʸa[i₁, J, k, s̄₁, s̄₁, j₂, K])
            end
        end
    end
end
