# Compare the prototype with QuantumAlgebra for nested commutators and powers of Hamiltonians
using QuantumAlgebra, Printf
const QA = QuantumAlgebra
include(joinpath(@__DIR__, "ConcreteProto.jl"))
using .ConcreteProto

function nested(H, A, n)
    for _ = 1:n
        A = comm(H, A)
    end
    A
end
powr(H, n) = (X = H; for _ = 2:n; X = X*H; end; X)
powr(H::QuExpr, n) = (X = H; for _ = 2:n; X = normal_form(X*H); end; X)

function compare(name, f, H, A, n)
    f(H, A, 2)
    CH, CA = concrete(H, A)
    f(CH, CA, 2)
    tqa = @elapsed (R = f(H, A, n))
    tto = @elapsed ((CH, CA) = concrete(H, A))
    tc = @elapsed (CR = f(CH, CA, n))
    tback = @elapsed (Rc = QuExpr(CR))
    T = valtype(CR.terms)
    @printf("%-44s %8d terms   QA %8.4f s   concrete %8.4f s (%5.1f×, %s)   to %.4f s, back %.4f s   %s\n",
            name, length(R.terms), tqa, tc, tqa/tc, T <: Complex{<:Rational} ? "rational" : T <: Complex{Int} ? "integer" : string(T),
            tto, tback, Rc == normal_form(R) ? "equal" : "DIFFERENT")
    flush(stdout)
end

@boson a b[:]
@tls σ[:]
@tlsxyz τ[:]
@fermion c[:,:]
@spin 1//2 S[:]
@spin 1 K[:]
@spin 5//2 C
@qparam ω Ω[:] g[:] J Δ h[:] χ G F t U ωc δ
q = qudit(:q, 3; ndims=1)

for pm in (false, true)
    QA.use_σpm(pm)
    println("use_σpm($pm)")
    N = 4
    H_TC = ω*a'*a + sum(Ω[i]*σ'[i]*σ[i] + g[i]*(a'*σ[i] + σ'[i]*a) for i in 1:N)
    for n in (6, 8)
        compare("  Tavis-Cummings (N=$N), [H,a], depth $n", nested, H_TC, a, n)
    end
    compare("  Tavis-Cummings (N=$N), H^4", (H, A, n) -> powr(H, n), H_TC, a, 4)
    L = 10
    H_XXZ = sum(J*(τx[i]*τx[i+1] + τy[i]*τy[i+1]) + Δ*τz[i]*τz[i+1] for i in 1:L-1) + sum(h[i]*τz[i] for i in 1:L)
    for n in (pm ? (5, 6) : (6, 8))
        compare("  XXZ chain (L=$L), [H,τˣ₅], depth $n", nested, H_XXZ, τx[5], n)
    end
end
QA.use_σpm(false)

println("other operators")
M = 3
H_Kerr = sum(χ*b'[i]^2*b[i]^2 for i in 1:M) + sum(G*(b'[i]*b[i+1] + b'[i+1]*b[i]) for i in 1:M-1) + F*(b[1] + b'[1])
compare("  driven Kerr modes (M=$M), [H,b₁], depth 7", nested, H_Kerr, b[1], 7)
Ls = 4
H_FH = -t*sum(c'[i,s]*c[i+1,s] + c'[i+1,s]*c[i,s] for i in 1:Ls-1, s in 1:2) + U*sum(c'[i,1]*c[i,1]*c'[i,2]*c[i,2] for i in 1:Ls)
for n in (4, 6)
    compare("  Fermi-Hubbard chain (L=$Ls), [H,c₁↑], depth $n", nested, H_FH, c[1,1], n)
end
H_S = sum(J*(Sp[i]*Sm[i+1] + Sm[i]*Sp[i+1]) + Δ*Sz[i]*Sz[i+1] for i in 1:7) + sum(h[i]*Sz[i] for i in 1:8)
compare("  spin-1/2 XXZ chain (L=8), [H,Sᶻ₄], depth 6", nested, H_S, Sz[4], 6)
H_K = sum(J*(Kp[i]*Km[i+1] + Km[i]*Kp[i+1]) + Δ*Kz[i]*Kz[i+1] for i in 1:5) + sum(h[i]*Kz[i] for i in 1:6)
compare("  spin-1 XXZ chain (L=6), [H,Kᶻ₃], depth 5", nested, H_K, Kz[3], 5)
compare("  collective spin j=5/2, [χJz²+Ω(J₊+J₋),J₋], depth 12", nested, χ*Cz^2 + Ω[1]*(Cp + Cm), Cm, 12)
Nq = 3
H_Λ = ωc*a'*a + sum(Δ*q(3,3)[i] + δ*q(2,2)[i] + g[i]*(a'*q(1,3)[i] + q(3,1)[i]*a) + Ω[i]*(q(3,2)[i] + q(2,3)[i]) for i in 1:Nq)
for n in (4, 6)
    compare("  Λ-type qudits in cavity (N=$Nq), [H,q¹²₁], depth $n", nested, H_Λ, q(1,2)[1], n)
end
