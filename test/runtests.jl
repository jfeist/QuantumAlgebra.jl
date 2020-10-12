using QuantumAlgebra
using QuantumAlgebra: δ, OpSum, OpTerm, BaseOpProduct, OpIndex, _map_opsum_ops #, prodtuples, prodtuple, sumtuple, distribute_indices!
using Test

function myδ(i,j)
    iA,iB = OpIndex.((i,j))
    OpSum(OpTerm([δ(min(iA,iB),max(iA,iB))],BaseOpProduct()))
end
scal(x) = OpSum(((OpTerm(),x),))
comm(A,B) = A*B - B*A
CorrOrExp(A::OpSum) = _map_opsum_ops(CorrOrExp,A)
CorrOrExp(A::OpTerm) = length(A)>1 ? corr(A) : expval(A)

@testset "QuantumAlgebra.jl" begin

    @boson_ops b
    @test normal_form(b(:i)*bdag(:j)) == bdag(:j)*b(:i) + myδ(:i,:j)
    @test normal_form(b(:i)*adag(:j)) == adag(:j)*b(:i)

    @fermion_ops c
    @test normal_form(c(:i)*cdag(:j) + cdag(:j)*c(:i)) == myδ(:i,:j)
    @test normal_form(c(:i)*fdag(:j)) == fdag(:j)*c(:i)

    # all equal numbers should be equal scalars (ignore type)
    @test a() + 0 == a() + 0.0
    @test a() + 1 == a() + (1//1 + 0im)
    #@test scal(1.5) == scal(3//2)

    @test QuantumAlgebra.SpatialIndex(QuantumAlgebra.x) == QuantumAlgebra.x

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

    @test a()' == adag()
    @test normal_form(a(:m)*adag(:m)) == adag(:m)*a(:m) + 1
    @test normal_form((a(:m)*a(1))') == adag(1)*adag(:m)
    @test normal_form((a(1,2,:k)*a(:m))') == normal_form(adag(1,2,:k)*adag(:m))

    @test normal_form(fdag(:a)*f(:b) + f(:b)*fdag(:a)) == myδ(:a,:b)
    @test normal_form(f(:a)*f(:b) + f(:b)*f(:a)) == scal(0)
    @test normal_form(fdag(:a)*fdag(:b) + fdag(:b)*fdag(:a)) == scal(0)

    @test f()' == fdag()
    @test fdag()' == f()
    @test (fdag()*f())' == fdag()*f()
    @test comm(fdag(),f()) == -comm(f(),fdag())

    @test ∑(:j,a(:j))*∑(:i,a(:i)) == ∑(:i,∑(:j,a(:i)*a(:j)))

    tmp = ∑(:i,adag(:i)*a(:i))
    @test tmp' == tmp
    @test normal_form(a(:i)*∑(:i,a(:i))) == normal_form(∑(:i_1,a(:i_1)*a(:i)))
    @test normal_form(adag(:n)*tmp) == normal_form(∑(:i,adag(:n)*adag(:i)*a(:i)))
    @test normal_form(a(:n)   *tmp) == normal_form(∑(:i,adag(:i)*a(:n)*a(:i)) + a(:n))
    @test normal_form(param(:g,'n',:i)*∑(:i,a(:i))) == normal_form(∑(:i_1,param(:g,'n',:i)*a(:i_1)))
    @test normal_form(param(:g,'n',:i_1)*a(:i)*∑(:i,a(:i))) == normal_form(∑(:i_2,param(:g,'n',:i_1)*a(:i)*a(:i_2)))
    @test normal_form(param(:g,'n',:n)*∑(:i,a(:i))) == normal_form(∑(:i,param(:g,'n',:n)*a(:i)))

    @test Pr"gz_i,μ" == param(:gz,'r',(:i,:μ))
    @test Pc"gz_i,μ" == param(:gz,'n',(:i,:μ))
    @test (Pc"gz_i,μ")' == param(:gz,'c',(:i,:μ))

    for with_σpm in (false,true)
        QuantumAlgebra.use_σpm(with_σpm)

        @test QuantumAlgebra.using_σpm() == with_σpm

        @test_throws ArgumentError σx(:i,"aa")
        @test σx(:i,:b) == σx((:i,:b))
        @test σx(:i,:j) != σx(:i,:b)

        α = param(:α)
        tmp = α' * a(1) + α * a(1)
        @test tmp - tmp == scal(0)

        @test σx() == σp() + σm()
        @test σx(:i) == σp(:i) + σm(:i)
        @test σx(:i,:j) == σp((:i,:j)) + σm((:i,:j))

        @test normal_form(σy(:i)) == normal_form(-1im*(σp(:i) - σm(:i)))
        @test normal_form(σz(:i)) == normal_form(2*σp(:i)*σm(:i) - 1)

        for s=(σx,σy,σz)
            @test normal_form(s(:i)*s(:i)) == scal(1)
            @test s(:i)' == s(:i)
        end
        @test normal_form(σx(:i)*σy(:i)) == 1im*σz(:i)
        @test normal_form(σy(:i)*σx(:i)) == -1im*σz(:i)
        @test normal_form(σy()*σz()) == 1im*σx()
        @test normal_form(σx(:i,:j)*σz(:i,:j)) == -1im*σy(:i,:j)
        @test normal_form(σx(:i,:j)*σz(:i,:j)) != -1im*σy(:i,:k)

        @test normal_form(a(1) * (σy(1) * a(1))') == normal_form(a(1) * (adag(1) * σy(1))) == normal_form(adag(1)*a(1)*σy(1) + σy(1))

        # @test prodtuple(a(5)) == (a(5),)
        # # tuples come out ordered!
        # @test prodtuple(a(5)*a(4)) == (a(4),a(5))
        # @test prodtuple(a(5)*a(4)) != (a(5),a(4))

        # @test sumtuple(a(5)) == (a(5),)
        # # tuples come out ordered!
        # @test sumtuple(a(5)+a(4)) == (a(4),a(5))
        # @test sumtuple(a(5)+a(4)) != (a(5),a(4))

        # @test_throws ArgumentError prodtuples(a(5)+a(4))
        tmp1 = normal_form(3*Pc"ω"*Pc"g"*expval(σp(:k))*σp(:k)*adag(5)*a(5))
        if QuantumAlgebra.using_σpm()
            @test length(tmp1.terms) == 1
            @test first(tmp1.terms)[1].bares == first((adag(5)*a(5)*σp(:k)).terms)[1].bares
        else
            @test length(tmp1.terms) == 4
        end

        @test normal_form(myδ(:i,:k)*a(:k)) == normal_form(a(:i)*myδ(:k,:i))
        @test normal_form(myδ(:i,:k)*a(:k,:i)) == normal_form(a(:i,:i)*myδ(:k,:i))
        @test normal_form(myδ(:i,:k)*myδ(:i,:j)) == normal_form(myδ(:k,:i)*myδ(:j,:k))
        # k cannot be equal to 1 and 3 at the same time
        @test normal_form(myδ(1,:k)*myδ(:k,3)*σx(:k)) == scal(0)
        @test normal_form(myδ(1,:k)*myδ(:k,1)*σx(:k)) == myδ(1,:k)*σx(1)

        @test normal_form(comm(σx(5),σy(3))) == scal(0)
        @test normal_form(comm(σx(5),σx(5))) == scal(0)
        @test normal_form(comm(σx(1),σz(1))) == normal_form(-2im*σy(1))
        @test normal_form(comm(σx(:μ),σy(:ν))) == normal_form(2im*myδ(:μ,:ν)*σz(:ν))
        @test normal_form(1//2im * comm(σx(:m),σy(:m))) == normal_form(σz(:m))
        @test normal_form(σx(:a)*σy(:a)*σz(:a)) == scal(1im)

        @test normal_form(comm(Pc"g",a(5)+a(3))) == scal(0)
        @test normal_form(comm(Pc"g",a(5)*a(3))) == scal(0)
        @test normal_form(comm(a(5)+a(3),Pc"g")) == scal(0)
        @test normal_form(comm(a(5)*a(3),Pc"g")) == scal(0)

        @test normal_form(comm(a(5)+a(3),adag(5))) == normal_form(comm(a(5),adag(5))+ comm(a(3),adag(5)))
        @test normal_form(comm(a(5)+a(3),adag(5)*a(3))) == normal_form(comm(a(5),adag(5)*a(3))+ comm(a(3),adag(5)*a(3)))

        @test normal_form(adag(2)*σy(:i) - 14Pc"ω") == normal_form(-14Pc"ω" + σy(:i)*adag(2))
        @test normal_form(adag(2)*σy(:i)) == normal_form(σy(:i)*adag(2))

        @test normal_form(σp(1) * σm(1)) == normal_form(1//2*σz(1) + 1//2)
        @test normal_form(σm(1)*σp(1)) == normal_form(1 - σp(1)*σm(1))
        @test normal_form(σm(1)*σp(1)) == normal_form(σp(1)*σm(1) + comm(σm(1),σp(1)))
        @test normal_form(σz(1)) == normal_form(σp(1)*σm(1) - σm(1)*σp(1))

        @test comm(σp(1),σp(1)) == scal(0)
        @test normal_form(comm(σp(:n),σm(:n))) == normal_form(σz(:n))

        @test normal_form(comm(a(1),   adag(1)*a(1))) == a(1)
        @test normal_form(comm(adag(1),adag(1)*a(1))) == -adag(1)

        @test expval(scal(3)) == scal(3)
        @test corr(scal(3)) == scal(3)
        @test corr(3+adag(2)) == 3 + corr(adag(2))

        # @test ascorr(scal(3)) == scal(3)
        # @test ascorr(a(2)) == expval(a(2))
        # @test ascorr(scal(3)*Pc"g") == scal(3)*Pc"g"
        # @test ascorr(scal(3)*a(2)) == scal(3)*ExpVal(a(2))
        # @test ascorr(scal(3)+adag(2)) == scal(3) + expval(adag(2))

        # @test ascorr(a(2)*a(2)) == corr(a(2)*a(2)) + expval(a(2))*ExpVal(a(2))
        # @test ascorr(a(2)*a(:m))' == corr(adag(2)*adag(:m)) + expval(adag(:m))*ExpVal(adag(2))
        # tmpas = a.(1:3)
        # @test *(tmpas...) == a(1)*a(2)*a(3)
        # tmpEVs = expval.(tmpas)
        # # multiply with scal(1) to trigger reordering (with prefactor, the ordering is different than without for now)
        # @test ascorr(*(8,Pr"g",tmpas...)) * scal(1) == 8*Pr"g" * (corr(*(tmpas...)) + *(tmpEVs...) + tmpEVs[1]*corr(tmpas[2]*tmpas[3]) + tmpEVs[2]*corr(tmpas[1]*tmpas[3]) + tmpEVs[3]*corr(tmpas[1]*tmpas[2]))

        # @test a(1) < ascorr(a(1)*a(2)*a(3)*a(4))
        # @test a(1) < ascorr(a(1)*a(2)*a(3)*a(4)*a(5))

        # if QuantumAlgebra.using_σpm()
        #     @test ascorr(scal(-1)*param(:g,'r',1)*σp(1)) == -param(:g,'r',1)*ExpVal(σp(1))
        #     @test ascorr(∑(:i,σp(:i)*σm(:n))) == ∑(:i,corr(σp(:i)*σm(:n))) + ∑(:i,ExpVal(σp(:i))*ExpVal(σm(:n)))
        # else
        #     @test ascorr(scal(-1)*param(:g,'r',1)*σz(1)) == -param(:g,'r',1)*ExpVal(σz(1))
        #     @test ascorr(∑(:i,σy(:i)*σy(:n))) == ∑(:i,corr(σy(:i)*σy(:n))) + ∑(:i,ExpVal(σy(:i))*ExpVal(σy(:n))) - expval(σy(:n))*ExpVal(σy(:n))
        # end

        @test CorrOrExp(a(5)) == expval(a(5))
        @test CorrOrExp(a(5)*a(:i)) == corr(a(5)*a(:i))

        H = ∑(:i,param(:ω,'r',:i)*adag(:i)*a(:i))
        @test normal_form(comm(a(:i),H)) == param(:ω,'r',:i)*a(:i)
        @test normal_form(comm(a(:n),H)) == param(:ω,'r',:n)*a(:n)
        @test normal_form(comm(H,a(:n))) == -param(:ω,'r',:n)*a(:n)
        @test normal_form(comm(adag(:n),H)) == -param(:ω,'r',:n)*adag(:n)
        @test normal_form(comm(adag(:n)*a(:m),H)) == (param(:ω,'r',:m)-param(:ω,'r',:n))*adag(:n)*a(:m)

        @test normal_form(a()*H) == normal_form(∑(:i,param(:ω,'r',:i)*adag(:i)*a(:i)*a()))
        @test normal_form(a(:k)*H) == normal_form(param(:ω,'r',:k)*a(:k) + ∑(:i,param(:ω,'r',:i)*adag(:i)*a(:i)*a(:k)))
        HH = ∑(:i,param(:ω,'r',:i,:i)*adag(:i,:i)*a(:i,:i))
        @test normal_form(a(:k,:k)*HH) == normal_form(param(:ω,'r',:k,:k)*a(:k,:k) + ∑(:i,param(:ω,'r',:i,:i)*adag(:i,:i)*a(:i,:i)*a(:k,:k)))

        # @test Avac(H) == scal(0)
        # @test vacA(H) == scal(0)
        # @test vacA(adag(3)*σp(1)*σm(1)) == scal(0)
        # @test vacA(fdag(:n)) == scal(0)
        # @test vacA(f(:n)) == f(:n)
        # @test vacA(a(:n)) == a(:n)
        # @test Avac(fdag(:n)) == fdag(:n)
        # @test Avac(f(:n)) == scal(0)
        # @test Avac(a(3)*σp(1)*σm(1)) == scal(0)
        # @test Avac(σm(1)*σp(1)) == scal(1)
        # @test Avac(σp(1)*σm(1)) == scal(0)
        # @test Avac(σp(1)) == σp(1)
        # @test vacA(σm(1)) == σm(1)
        # if QuantumAlgebra.using_σpm()
        #     @test Avac(σm(1)) == scal(0)
        #     @test vacA(σp(1)) == scal(0)
        # else
        #     @test Avac(σx(1)) == vacA(σx(1)) == σx(1)
        # end
        # @test vacExpVal(σx(1)) == scal(0)
        # @test vacExpVal(σp(1)) == scal(0)
        # @test vacExpVal(∑(:i,σp(:i))) == scal(0)
        # @test vacExpVal(σp(:i)*σm(:k)) == scal(0)

        @test normal_form(a(:n)*adag(:n)*a(:n)*adag(:n)) == scal(1) + 3adag(:n)*a(:n) + adag(:n)*adag(:n)*a(:n)*a(:n)

        # S = scal(1/√(2*6))*adag(:n)*adag(:n)*adag(:n) + scal(1/√2)*adag(:m)
        # for (A,val) in [(scal(1),scal(1)),
        #                 (adag(:n)*a(:n),scal(1.5) + 0.5 * myδ(:n,:m)),
        #                 (adag(:n)*adag(:n)*a(:n)*a(:n),scal(3))]
        #     @test vacExpVal(A,S) ≈ val
        # end

        # tmp = scal(1+2im)*∑(:i,a(:i)*adag(:i)*ascorr(adag(:n)*a(:m)))
        # @test latex(scal(1+2im)) == "(1+2i)"
        # @test latex(scal(1+2im//5)) == "\\left(1+\\frac{2}{5}i\\right)"

        # tmp = ∑(:i,ascorr(adag(:n)*a(:i)))
        # tmplatex = "\\sum_{i}\\langle {a}_{n}^\\dagger {a}_{i} \\rangle_{c} + \\sum_{i}\\langle {a}_{n}^\\dagger \\rangle \\langle {a}_{i} \\rangle"
        # @test latex(tmp) == tmplatex
        # @test ascorr(tmp) == tmp
        # @test sprint(show,"text/latex",tmp) == "\$$(tmplatex)\$"
        # if QuantumAlgebra.using_σpm()
        #     @test latex(σp()) == "\\sigma^+"
        # else
        #     @test latex(σz()) == "\\sigma_{z}"
        # end

        # inds = [:a,:b,:c,:d,:e,:f,:g,:h,:i,:j,:k,:l,:m,:n]
        # tmp1 = param(:ω,:y)*a(1)*adag(1)*a(3)*adag(4)*ExpVal(a(5))*corr(adag(5)*a(9))
        # tmp2 = param(:ω,:a)*ExpVal(a(:b))*corr(adag(:c)*a(:d))*adag(:e)*a(:f) + param(:ω,:g)*ExpVal(a(:h))*corr(adag(:i)*a(:j))*adag(:k)*adag(:l)*a(:m)*a(:n)
        # @test distribute_indices!(copy(inds),tmp1) == tmp2
        # if QuantumAlgebra.using_σpm()
        #     tmp1 = a(1,:n)*adag()*σm(1,:n)*σp()
        #     @test distribute_indices!(copy(inds),tmp1) == adag()*a(:a,:b)*σp()*σm(:c,:d)
        # else
        #     tmp1 = a(1,:n)*adag()*σz(1,:n)*σy()
        #     @test distribute_indices!(copy(inds),tmp1) == adag()*a(:a,:b)*σy()*σz(:c,:d)
        # end

        # @test_throws MethodError distribute_indices!(copy(inds),∑(:i,a(:i)))
        # @test_throws ArgumentError distribute_indices!([:a,:b],tmp1)

        # @test QuantumAlgebra.exchange_inds(adag(:j)*a(:k),:k,:j) == adag(:k)*a(:j)
        # @test QuantumAlgebra.extindices(∑(:i,adag(:i)*a(:k))) == [:k]
        # @test QuantumAlgebra.symmetric_index_nums(adag(:i)*adag(:j)*a(:k)*a(:l)) == [2,2]

        # @test string(∑(:i,a(:i)) * adag(:n)) == "1 + ∑_i a†(n) a(i)"
        # if QuantumAlgebra.using_σpm()
        #     @test string(a(5)*adag(5)*σp(3)*ascorr(adag(5,:i)*a(5))) == "⟨a†(5i)⟩ ⟨a(5)⟩ σ+(3) + ⟨a†(5i) a(5)⟩c σ+(3) + ⟨a†(5i)⟩ ⟨a(5)⟩ a†(5) a(5) σ+(3) + ⟨a†(5i) a(5)⟩c a†(5) a(5) σ+(3)"
        # else
        #     @test string(a(5)*adag(5)*σz(3)*ascorr(adag(5,:i)*a(5))) == "⟨a†(5i)⟩ ⟨a(5)⟩ σz(3) + ⟨a†(5i) a(5)⟩c σz(3) + ⟨a†(5i)⟩ ⟨a(5)⟩ a†(5) a(5) σz(3) + ⟨a†(5i) a(5)⟩c a†(5) a(5) σz(3)"
        # end
    end
end
