for v in [0,1,2,3//2,3//1,1.5,1.5im,1im]
    scal(v)
end

param(:g)'

Pr"ω_i,j,3"

σx(:i)*σy(:i)
σz(3)*σy(2)

a(:m)*adag(:n)*adag(3)*a(3)*adag(3)*a(:n)*adag(:m)

σy(1)'
σy(:m)'
(adag(1)*a(:n))' == adag(:n)*a(1)
a(1) * (σy(1) * a(1))'

let
    tmp = ∑(:i,adag(:i)*a(:i))
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
∑(:i,σx(:i))*σy(:j)
∑(:i,σx(:i))*σy(:a)
∑(:i,σx(:i)*σy(:a))

comm(σx(5),σy(3))
comm(σx(5),σx(5))
comm(σx(1),σz(1))
scal(2)*comm(σx(5),σy(5))
comm(σy(5),comm(σx(5),σy(5)))

comm(scal(2//5)*param(:h)*σx(5),scal(3)*param(:g)*σy(5))

σx(1)*σy(1)
σx(3)*σy(2)*σy(2)

σx(5)*σy(5)*σz(5)
σx(5)*σy(5)*σx(5)

param(:g,2)*param(:g,1)*param(:g,'c',3)*param(:b)*param(:a)*param(:d,'c',3)*param(:f,'r',1)

a(2)*adag(2)*adag(3)*σy(1)*adag(2)*adag(4)*σz(1)*a(2)*adag(3)

adag(2)*a(2)*adag(2)*σy(1) + scal(2)*param(:g,'c',1)

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

ascorr(scal(3))
ascorr(a(2))
ascorr(scal(3)*a(2))

ascorr(a(2)*a(2))
ascorr(scal(2)*a(1)*a(2)*a(3))
ascorr(a(1)*a(2)*a(3)*a(4))
ascorr(adag(2)*a(1)*σz(1))

ascorr(param(:g,'r')*adag(3)*a(2)*a(2))
ascorr(scal(-1)*param(:g,'r',1)*σz(1))
ascorr(adag(:n)*σy(:k)*a(:m)*adag(:m))

let
    H = ∑(:i,param(:ω,'r',:i)*adag(:i)*a(:i)) + ∑(:j,scal(1//2)*param(:ωe,'r',:j)*σz(:j)) +
        ∑(:i,∑(:j,param(:g,'r',(:i,:j))*(adag(:i)+a(:i))*σx(:j)))
    comm(adag(:n)*a(:m),H)
    comm(adag(:n)*a(:n),H)
    for op in [a(:n),adag(:n),σx(:k),σy(:k),σz(:k)]
        latex(op),latex(comm(op,H))
    end
    for op in [adag(:n)*σz(:k),a(:n)*σz(:k)]
        latex(op),latex(comm(op,H))
    end

    for op in [adag(:n)*σz(:k),a(:n)*σz(:k)]
        for x = 1:5
            op = comm(op,H)
        end
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
