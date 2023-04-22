using PrecompileTools

@compile_workload let
    @static if !_DEFINE_DEFAULT_OPS
        @boson_ops a
        @fermion_ops f
        @tlspm_ops σ
        @tlsxyz_ops σ
    end

    param(:g,'n',:k,3)'
    Pr"ω_i,j,3"

    σx(:i)*σy(:i)
    σz(3)*σy(2)

    a(:m)*a'(:n)*a'(3)*a(3)*a'(3)*a(:n)*a'(:m)

    σy(1)'
    σy(:m)'
    (a'(1)*a(:n))' == a'(:n)*a(1)
    a(1) * (σy(1) * a(1))'

    tmp = ∑(:i,a'(:i)*a(:i))
    a'(:n)*tmp
    a(:n)*tmp
    σz(:n)*tmp
    tmp*2
    2*tmp
    tmp*param(:g,'n',:n)
    tmp*a(:n)
    tmp*a'(:n)
    tmp*(a'(:n)*a'(:n))
    tmp*(a'(:n)*a(:n))

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

    normal_form(a(2)*a'(2)*a'(3)*σy(1)*a'(2)*a'(4)*σz(1)*a(2)*a'(3))

    a'(2)*a(2)*a'(2)*σy(1) + 2*param(:g,'c',1)

    -1//2 + σp(1) * σm(1)

    σp(2)
    σm(2)
    σp(1)*σm(1)

    comm(σp(1),σp(1))
    comm(σm(:n),σp(:n))
    comm(a(2),a'(:m))

    comm(a(1),a'(1)*a(1))
    comm(a(1),a'(1)*a(2)*a(1))
    comm(a(1),a(2)*a'(1)*a(1)*a'(2))

    expval_as_corrs(a(2))
    expval_as_corrs(3*a(2))

    expval_as_corrs(a(2)*a(2))
    expval_as_corrs(2*a(1)*a(2)*a(3))
    expval_as_corrs(a(1)*a(2)*a(3)*a(4))
    expval_as_corrs(a'(2)*a(1)*σz(1))

    expval_as_corrs(param(:g,'r')*a'(3)*a(2)*a(2))
    expval_as_corrs(-1*param(:g,'r',1)*σz(1))
    expval_as_corrs(a'(:n)*σy(:k)*a(:m)*a'(:m))

    H = ∑(:i,param(:ω,'r',:i)*a'(:i)*a(:i)) + ∑(:j,(1//2)*param(:ωe,'r',:j)*σz(:j)) +
        ∑(:i,∑(:j,param(:g,'r',(:i,:j))*(a'(:i)+a(:i))*σx(:j)))
    comm(a'(:n)*a(:m),H)
    comm(a'(:n)*a(:n),H)
    for op in [a(:n),a'(:n),σx(:k),σy(:k),σz(:k)]
        comm(op,H)
    end
    for op in [a'(:n)*σz(:k),a(:n)*σz(:k)]
        comm(op,H)
    end

    for op in [a'(:n)*σz(:k),a(:n)*σz(:k)]
        for x = 1:3
            op = comm(op,H)
        end
    end

    expval(H)
    #expval_as_corrs(H)

    vacA(Avac(σp(1)*σm(1)))
    vacA(Avac(σm(1)*σp(1)))
    vacExpVal(a'(),a())

    for i2 = (:n,:m)
        stateop = (1/√2)*a'(:n)*a'(:n)*a'(:n) + (1/√2)*a'(i2)
        for A in [a'(:n)*a(:n),a'(:n)*a'(:n)*a(:n)*a(:n)]
            vacExpVal(A,stateop)
        end
    end

    let
        dotest = (H,A,n=3) -> begin
            for _ = 1:n
                A = normal_form(comm(H,A))
            end
            A
        end
        H = ∑(:n,∑(:m,∑(:k,Pr"ω_n,m"*a'(:n,:k)*a(:m,:k)))) + ∑(:i,1//2*Pr"ν_i"*σz(:i)) + ∑(:n,∑(:k,∑(:i,Pr"g_i,n,k"*σx(:i)*(a'(:n,:k)+a(:n,:k)))))
        dotest(H,σz(:j))
    end

    let allops = (a(1),f(:α),σz(:i),a(2),f(:β),σz(:j),a(3),f(:γ),σz(:k))
        for N = 1:5
            ops = allops[1:N]
            A = normal_form(8*Pr"g_i"*prod(ops))
            expval_as_corrs(A)
        end
    end
end