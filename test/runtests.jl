using QuantumAlgebra
using Test

@testset "QuantumAlgebra.jl" begin
    @test QuantumAlgebra.SpatialIndex(QuantumAlgebra.x) == QuantumAlgebra.x

    # all equal numbers should be equal scalars (ignore type)
    @test scal(0) == scal(0.0)
    @test scal(0) == scal(0im)
    @test scal(1.5) == scal(3//2)

    # test params and their complex conjugation
    @test param(:ω) == param(:ω,())
    @test param(:ω,(:i),'n') == param(:ω,:i)

    @test adjoint(param(:g,(),'r')) == param(:g,(),'r')
    @test adjoint(param(:g,(),'n')) == param(:g,(),'c')
    @test adjoint(param(:g,:i,'n')) == param(:g,:i,'c')
    @test adjoint(param(:g,(:i,1,:m),'c')) == param(:g,(:i,1,:m),'n')

    tmp1 = param(:g,2)*param(:g,1)*param(:g,3,'c')*param(:b)*param(:a)*param(:d,3,'c')*param(:f,(1,:i),'r')
    tmpc = param(:g,2,'c')*param(:g,1,'c')*param(:g,3)*param(:b,(),'c')*param(:a,(),'c')*param(:d,3,'n')*param(:f,(1,:i),'r')
    tmp2 = param(:a)*param(:b)*param(:g,1)*param(:g,2)*param(:g,3,'c')*param(:d,3,'c')*param(:f,(1,:i),'r')
    @test tmp1 == tmp2
    @test tmp1' == tmpc

    @test σx(:i) == σ(1,:i)
    @test σy(:i) == σ(2,:i)
    @test σz(:i) == σ(3,:i)
    @test σx(:i) == σ(QuantumAlgebra.x,:i)

    for ii=1:3
        @test σ(ii,:i)*σ(ii,:i) == scal(1)
        @test σ(ii,:i)' == σ(ii,:i)
    end
    @test σx(:i)*σy(:i) == scal(1im)*σz(:i)
    @test σy(:i)*σx(:i) == scal(-1im)*σz(:i)
    @test σy(:i)*σz(:i) == scal(1im)*σx(:i)
    @test σx(:i)*σz(:i) == scal(-1im)*σy(:i)

    @test a(:m)*adag(:m) == adag(:m)*a(:m) + scal(1)
    @test (a(:m)*a(1))' == adag(:m)*adag(1)
    @test (a(1)*a(:m))' == adag(1)*adag(:m)

    @test scal(1+1im) < scal(2-1im)
    @test scal(1+1im) > scal(1-1im)
    @test a(1) < σx(1)
    @test adag(1) < σx(1)
    @test !(a(1) < adag(1))
    @test adag(1) < a(1)
    @test a(5) < a(:i)
    @test !(adag(:m) < adag(2))

    @test a(1) * (σ(2,1) * a(1))' == a(1) * (adag(1) * σy(1)) == adag(1)*a(1)*σy(1) + σy(1)

    tmp = OpSumAnalytic(:i,adag(:i)*a(:i))
    @test tmp' == tmp
    @test_throws ErrorException a(:i)*OpSumAnalytic(:i,a(:i))
    @test adag(:n)*tmp == OpSumAnalytic(:i,adag(:n)*adag(:i)*a(:i))
    @test a(:n)   *tmp == OpSumAnalytic(:i,adag(:i)*a(:n)*a(:i)) + a(:n)
    @test_throws ErrorException param(:g,:i)*OpSumAnalytic(:i,a(:i))
    @test param(:g,:n)*OpSumAnalytic(:i,a(:i)) == OpSumAnalytic(:i,param(:g,:n)*a(:i))

    @test QuantumAlgebra.prodtuple(a(5)) == (a(5),)
    # tuples come out ordered!
    @test QuantumAlgebra.prodtuple(a(5)*a(4)) == (a(4),a(5))
    @test QuantumAlgebra.prodtuple(a(5)*a(4)) != (a(5),a(4))

    @test QuantumAlgebra.sumtuple(a(5)) == (a(5),)
    # tuples come out ordered!
    @test QuantumAlgebra.sumtuple(a(5)+a(4)) == (a(4),a(5))
    @test QuantumAlgebra.sumtuple(a(5)+a(4)) != (a(5),a(4))

    @test_throws ErrorException QuantumAlgebra.prodtuples(a(5)+a(4))
    # tuples come out ordered!
    tmp1 = scal(3)*param(:ω)*param(:g)*ExpVal(σz(:k))*σz(:k)*adag(5)*a(5)
    tmp2 = ( (scal(3),param(:g),param(:ω)), (ExpVal(σz(:k)),), (adag(5),a(5),σz(:k)) )
    @test QuantumAlgebra.prodtuples(tmp1) == tmp2

    @test comm(σ(1,5),σ(2,3)) == scal(0)
    @test comm(σ(1,5),σ(1,5)) == scal(0)
    @test comm(σx(1),σz(1)) == scal(-2im)*σy(1)
    @test comm(σx(:mu),σy(:muuu)) == scal(0)
    @test scal(1//2im)*comm(σx(:m),σy(:m)) == σz(:m)
    @test σx(:a)*σy(:a)*σz(:a) == scal(1im)

    @test comm(param(:g),a(5)+a(3)) == scal(0)
    @test comm(param(:g),a(5)*a(3)) == scal(0)
    @test comm(a(5)+a(3),param(:g)) == scal(0)
    @test comm(a(5)*a(3),param(:g)) == scal(0)

    @test comm(a(5)+a(3),adag(5)) == comm(a(5),adag(5))+ comm(a(3),adag(5))
    @test comm(a(5)+a(3),adag(5)*a(3)) == comm(a(5),adag(5)*a(3))+ comm(a(3),adag(5)*a(3))

    @test adag(2)*σy(:i) - scal(14)*param(:ω) == scal(-14)*param(:ω) + σy(:i)*adag(2)
    @test adag(2)*σy(:i) == σy(:i)*adag(2)

    @test σp(1) * σm(1) == scal(1//2)*σz(1) + scal(1//2)

    @test comm(σp(1),σp(1)) == scal(0)
    @test comm(σp(:n),σm(:n)) == σz(:n)

    @test comm(a(1),   adag(1)*a(1)) == a(1)
    @test comm(adag(1),adag(1)*a(1)) == -adag(1)

    @test ExpVal(scal(3)) == scal(3)
    @test Corr(scal(3)) == scal(3)
    @test Corr(scal(3)+adag(2)) == scal(3) + Corr(adag(2))

    @test ascorr(scal(3)) == scal(3)
    @test ascorr(a(2)) == ExpVal(a(2))
    @test ascorr(scal(3)*param(:g)) == scal(3)*param(:g)
    @test ascorr(scal(3)*a(2)) == scal(3)*ExpVal(a(2))
    @test ascorr(scal(3)+adag(2)) == scal(3) + ExpVal(adag(2))

    @test ascorr(a(2)*a(2)) == Corr(a(2)*a(2)) + ExpVal(a(2))*ExpVal(a(2))
    @test ascorr(a(2)*a(:m))' == Corr(adag(2)*adag(:m)) + ExpVal(adag(:m))*ExpVal(adag(2))
    tmpas = a.(1:3)
    @test *(tmpas...) == a(1)*a(2)*a(3)
    tmpEVs = ExpVal.(tmpas)
    @test ascorr(*(tmpas...)) == Corr(*(tmpas...)) + *(tmpEVs...) + tmpEVs[1]*Corr(tmpas[2]*tmpas[3]) + tmpEVs[2]*Corr(tmpas[1]*tmpas[3]) + tmpEVs[3]*Corr(tmpas[1]*tmpas[2])

    @test a(1) < ascorr(a(1)*a(2)*a(3)*a(4))
    @test_throws ErrorException ascorr(a(1)*a(2)*a(3)*a(4)*a(5))

    @test ascorr(scal(-1)*param(:g,1,'r')*σ(3,1)) == -param(:g,1,'r')*ExpVal(σz(1))
    @test ascorr(OpSumAnalytic(:i,σy(:i)*σy(:n))) == OpSumAnalytic(:i,Corr(σy(:i)*σy(:n))) + OpSumAnalytic(:i,ExpVal(σy(:i))*ExpVal(σy(:n))) - ExpVal(σy(:n))*ExpVal(σy(:n))

    @test CorrOrExp(a(5)) == ExpVal(a(5))
    @test CorrOrExp(a(5)*a(:i)) == Corr(a(5)*a(:i))

    H = OpSumAnalytic(:i,param(:ω,:i,'r')*adag(:i)*a(:i))
    # cannot commute with an operator with the same index as in the sum
    @test_throws ErrorException comm(a(:i),H)
    @test comm(a(:n),H) == param(:ω,:n,'r')*a(:n)
    @test comm(H,a(:n)) == -param(:ω,:n,'r')*a(:n)
    @test comm(adag(:n),H) == -param(:ω,:n,'r')*adag(:n)
    @test comm(adag(:n)*a(:m),H) == (param(:ω,:m,'r')-param(:ω,:n,'r'))*adag(:n)*a(:m)

    @test Avac(H) == scal(0)
    @test vacA(H) == scal(0)
    @test vacA(adag(3)*σp(1)*σm(1)) == scal(0)
    @test vacA(a(:n)) == a(:n)
    @test Avac(a(3)*σp(1)*σm(1)) == scal(0)
    @test Avac(σm(1)*σp(1)) == scal(1)
    @test Avac(σp(1)*σm(1)) == scal(0)
    @test Avac(σx(1)) == vacA(σx(1)) == σx(1)
    @test vacExpVal(σx(1)) == scal(0)
    @test vacExpVal(σp(1)) == scal(0)
    @test vacExpVal(OpSumAnalytic(:i,σp(:i))) == scal(0)
    @test vacExpVal(σp(:i)*σm(:k)) == scal(0)

    @test a(:n)*adag(:n)*a(:n)*adag(:n) == scal(1) + scal(3)*adag(:n)*a(:n) + adag(:n)*adag(:n)*a(:n)*a(:n)

    S = scal(1/√(2*6))*adag(:n)*adag(:n)*adag(:n) + scal(1/√2)*adag(:m)
    for (A,val) in [(scal(1),1),
                    (adag(:n)*a(:n),1.5),
                    (adag(:n)*adag(:n)*a(:n)*a(:n),3)]
        @test vacExpVal(A,S).v ≈ val
    end

    tmp = scal(1+2im)*OpSumAnalytic(:i,a(:i)*adag(:i)*ascorr(adag(:n)*a(:m)))
    @test latex(scal(1+2im)) == "(1+2i)"
    @test latex(scal(1+2im//5)) == "\\left(1+\\frac{2}{5}i\\right)"

    tmp = OpSumAnalytic(:i,ascorr(adag(:n)*a(:i)))
    tmplatex = "\\sum_{i}\\langle a_{n}^\\dagger a_{i} \\rangle_{c} + \\sum_{i}\\langle a_{n}^\\dagger \\rangle \\langle a_{i} \\rangle"
    @test latex(tmp) == tmplatex
    @test ascorr(tmp) == tmp
    @test sprint(show,"text/latex",tmp) == "\$$(tmplatex)\$"

    inds = [:a,:b,:c,:d,:e,:f,:g,:h,:i,:j,:k,:l,:m,:n]
    tmp1 = param(:ω,:y)*a(1)*adag(1)*a(3)*adag(:a)*ExpVal(a(:n))*Corr(adag(:n)*a(:n))
    tmp2 = param(:ω,:a)*ExpVal(a(:b))*Corr(adag(:c)*a(:d))*adag(:e)*a(:f) + param(:ω,:g)*ExpVal(a(:h))*Corr(adag(:i)*a(:j))*adag(:k)*adag(:l)*a(:m)*a(:n)
    @test tmp2 == QuantumAlgebra.distribute_indices!(copy(inds),tmp1)
    @test_throws MethodError QuantumAlgebra.distribute_indices!(copy(inds),OpSumAnalytic(:i,a(:i)))
    @test_throws ArgumentError QuantumAlgebra.distribute_indices!([:a,:b],tmp1)

    @test string(OpSumAnalytic(:i,a(:i)) * adag(:n)) == "1 + Σ_i a†(n) a(i)"
    @test string(a(5)*adag(5)*σz(3)*ascorr(adag(5)*a(5))) == "⟨a†(5)⟩ ⟨a(5)⟩ σz(3) + ⟨a†(5) a(5)⟩c σz(3) + ⟨a†(5)⟩ ⟨a(5)⟩ a†(5) a(5) σz(3) + ⟨a†(5) a(5)⟩c a†(5) a(5) σz(3)"
end
