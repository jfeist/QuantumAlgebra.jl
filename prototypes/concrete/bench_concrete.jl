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

function compare(name, f, H, A, n)
    f(H, A, 2)
    CH, CA = concrete(H, A)
    f(CH, CA, 2)
    tqa = @elapsed (R = f(H, A, n))
    tto = @elapsed ((CH, CA) = concrete(H, A))
    tc = @elapsed (CR = f(CH, CA, n))
    tback = @elapsed (Rc = QuExpr(CR))
    @printf("%-40s %8d terms   QA %9.4f s   concrete %9.4f s (%6.1f×)   conversion to %.4f s, back %.4f s   %s\n",
            name, length(R.terms), tqa, tc, tqa/tc, tto, tback, Rc == normal_form(R) ? "equal" : "DIFFERENT")
    flush(stdout)
end

QA.use_σpm(false)
@boson a b[:]
@tls σ[:]
@tlsxyz τ[:]
@qparam ω Ω[:] g[:] J Δ h[:] χ G F

N = 4
H_TC = ω*a'*a + sum(Ω[i]*σ'[i]*σ[i] + g[i]*(a'*σ[i] + σ'[i]*a) for i in 1:N)
for n in (4, 6, 8)
    compare("Tavis-Cummings (N=$N), [H,a], depth $n", nested, H_TC, a, n)
end
for n in (3, 4)
    compare("Tavis-Cummings (N=$N), H^$n", (H, A, n) -> powr(H, n), H_TC, a, n)
end

L = 10
H_XXZ = sum(J*(τx[i]*τx[i+1] + τy[i]*τy[i+1]) + Δ*τz[i]*τz[i+1] for i in 1:L-1) + sum(h[i]*τz[i] for i in 1:L)
for n in (4, 6, 8)
    compare("XXZ chain (L=$L), [H,τˣ₅], depth $n", nested, H_XXZ, τx[5], n)
end

M = 3
H_Kerr = sum(χ*b'[i]^2*b[i]^2 for i in 1:M) + sum(G*(b'[i]*b[i+1] + b'[i+1]*b[i]) for i in 1:M-1) + F*(b[1] + b'[1])
for n in (4, 6, 7)
    compare("driven Kerr modes (M=$M), [H,b₁], depth $n", nested, H_Kerr, b[1], n)
end
