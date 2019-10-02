for v in [0,1,2,3//2,3//1,1.5,1.5im,1im]
    scal(v)
end

param(:g)

σx(:i)*σy(:i)
σz(3)*σy(2)

a(:m)*adag(:n)*adag(3)*a(3)*adag(3)*a(:n)*adag(:m)

σy(1)'
σy(:m)'
(adag(1)*a(:n))' == adag(:n)*a(1)
a(1) * (σ(2,1) * a(1))'

let
    tmp = OpSumAnalytic(:i,adag(:i)*a(:i))
    adag(:n)*tmp
    a(:n)*tmp
    σz(:n)*tmp
    tmp*scal(2)
    tmp*param(:g,:n)
    tmp*a(:n)
    tmp*adag(:n)
    tmp*(adag(:n)*adag(:n))
    tmp*(adag(:n)*a(:n))
end
OpSumAnalytic(:i,σx(:i))*σy(:j)
OpSumAnalytic(:i,σx(:i))*σy(:a)
OpSumAnalytic(:i,σx(:i)*σy(:a))

comm(σ(1,5),σ(2,3))
comm(σ(1,5),σ(1,5))
comm(σ(1,1),σ(3,1))
scal(2)*comm(σ(1,5),σ(2,5))
comm(σ(2,5),comm(σ(1,5),σ(2,5)))

comm(scal(2//5)*param(:h)*σ(1,5),scal(3)*param(:g)*σ(2,5))

σ(1,1)*σ(2,1)
σ(1,3)*σ(2,2)*σ(2,2)

σ(1,5)*σ(2,5)*σ(3,5)
σ(1,5)*σ(2,5)*σ(1,5)

param(:g,2)*param(:g,1)*param(:g,'c',3)*param(:b)*param(:a)*param(:d,'c',3)*param(:f,'r',1)

a(2)*adag(2)*adag(3)*σ(2,1)*adag(2)*adag(4)*σ(3,1)*a(2)*adag(3)

adag(2)*a(2)*adag(2)*σ(2,1) + scal(2)*param(:g,'c',1)

scal(-1//2) + σp(1) * σm(1)

σp(2)
σm(2)
σp(1)*σm(1)

comm(σp(1),σp(1))
comm(σm(:n),σp(:n))
comm(a(2),adag(:m))

comm(a(1),adag(1)*a(1))
comm(a(1),adag(1)*a(2)*a(1))
comm(a(1),a(2)*adag(1)*a(1)*adag(2))

# ascorr(scal(3))
# ascorr(a(2))
# ascorr(scal(3)*a(2))

# ascorr(a(2)*a(2))
# ascorr(scal(2)*a(1)*a(2)*a(3))
# ascorr(a(1)*a(2)*a(3)*a(4))
# ascorr(adag(2)*a(1)*σ(3,1))

# ascorr(param(:g,'r')*adag(3)*a(2)*a(2))
# ascorr(scal(-1)*param(:g,'r',1)*σ(3,1))
# ascorr(adag(:n)*σy(:k)*a(:m)*adag(:m))

let
    H = OpSumAnalytic(:i,param(:ω,'r',:i)*adag(:i)*a(:i)) + OpSumAnalytic(:j,scal(1//2)*param(:ωe,'r',:j)*σ(3,:j)) +
        OpSumAnalytic(:i,OpSumAnalytic(:j,param(:g,'r',(:i,:j))*(adag(:i)+a(:i))*σ(1,:j)))
    comm(adag(:n)*a(:m),H)
    comm(adag(:n)*a(:n),H)
    for op in [a(:n),adag(:n),σ(1,:k),σ(2,:k),σ(3,:k)]
        latex(op),latex(comm(op,H))
    end
    for op in [adag(:n)*σz(:k),a(:n)*σz(:k)]
        latex(op),latex(comm(op,H))
    end

    ExpVal(H)
    #ascorr(H)
end

vacA(Avac(σp(1)*σm(1)))
vacA(Avac(σm(1)*σp(1)))

for i2 = (:n,:m)
    stateop = scal(1/√2)*adag(:n)*adag(:n)*adag(:n) + scal(1/√2)*adag(i2)
    for A in [scal(1),adag(:n)*a(:n),adag(:n)*adag(:n)*a(:n)*a(:n)]
        latex(A),latex(vacExpVal(A,stateop))
    end
end
