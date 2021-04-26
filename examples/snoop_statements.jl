param(:g,'n',:k,3)'

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
    tmp*2
    2*tmp
    tmp*param(:g,'n',:n)
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
2*comm(σx(5),σy(5))
comm(σy(5),comm(σx(5),σy(5)))

comm((2//5)*Pc"h"*σx(5),3Pc"g"*σy(5))

σx(1)*σy(1)
σx(3)*σy(2)*σy(2)

σx(5)*σy(5)*σz(5)
σx(5)*σy(5)*σx(5)

Pc"g_2"*Pc"g_1"*(Pc"g_3")'*Pr"b"*Pc"a"*Pr"d_3"

normal_form(a(2)*adag(2)*adag(3)*σy(1)*adag(2)*adag(4)*σz(1)*a(2)*adag(3))

adag(2)*a(2)*adag(2)*σy(1) + 2*param(:g,'c',1)

-1//2 + σp(1) * σm(1)

σp(2)
σm(2)
σp(1)*σm(1)

comm(σp(1),σp(1))
comm(σm(:n),σp(:n))
comm(a(2),adag(:m))

comm(a(1),adag(1)*a(1))
comm(a(1),adag(1)*a(2)*a(1))
comm(a(1),a(2)*adag(1)*a(1)*adag(2))

expval_as_corrs(a(2))
expval_as_corrs(3*a(2))

expval_as_corrs(a(2)*a(2))
expval_as_corrs(2*a(1)*a(2)*a(3))
expval_as_corrs(a(1)*a(2)*a(3)*a(4))
expval_as_corrs(adag(2)*a(1)*σz(1))

expval_as_corrs(param(:g,'r')*adag(3)*a(2)*a(2))
expval_as_corrs(-1*param(:g,'r',1)*σz(1))
expval_as_corrs(adag(:n)*σy(:k)*a(:m)*adag(:m))

let
    H = ∑(:i,param(:ω,'r',:i)*adag(:i)*a(:i)) + ∑(:j,(1//2)*param(:ωe,'r',:j)*σz(:j)) +
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

    expval(H)
    #expval_as_corrs(H)
end

vacA(Avac(σp(1)*σm(1)))
vacA(Avac(σm(1)*σp(1)))

vacExpVal(adag(),a())

for i2 = (:n,:m)
    stateop = (1/√2)*adag(:n)*adag(:n)*adag(:n) + (1/√2)*adag(i2)
    for A in [adag(:n)*a(:n),adag(:n)*adag(:n)*a(:n)*a(:n)]
        latex(A),latex(vacExpVal(A,stateop))
    end
end

function dotest(H,A,n=5)
    for _ = 1:n
        A = normal_form(comm(H,A))
    end
    A
end

H = ∑(:n,∑(:m,∑(:k,Pr"ω_n,m"*adag(:n,:k)*a(:m,:k)))) + ∑(:i,1//2*Pr"ν_i"*σz(:i)) + ∑(:n,∑(:k,∑(:i,Pr"g_i,n,k"*σx(:i)*(adag(:n,:k)+a(:n,:k)))))
dotest(H,σz(:j))


allops = (a(1),f(:α),σz(:i),a(2),f(:β),σz(:j),a(3),f(:γ),σz(:k))
for N = 1:9
    ops = allops[1:N]
    A = normal_form(8*Pr"g_i"*prod(ops))
    expval_as_corrs(A)
end
