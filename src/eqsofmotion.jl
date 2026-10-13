export heisenberg_eom

"""
    AbstractLindbladTerm

A Lindblad term of a master equation (see [`heisenberg_eom`](@ref)), which is applied to an operator
`A` by calling it, `L(A)`. The subtypes are [`SingleLindbladTerm`](@ref),
[`MixedLindbladTerm`](@ref) and [`SummedLindbladTerm`](@ref); they are created from the tuples that
describe them with [`lindbladterm`](@ref).
"""
abstract type AbstractLindbladTerm end

"""
    SingleLindbladTerm(γ, L)

The Lindblad term with rate `γ` and jump operator `L`, which gives
``γ (L^† A L - \\frac{1}{2} \\{L^† L, A\\})`` in the equation of motion of an operator ``A``.
"""
struct SingleLindbladTerm{T} <: AbstractLindbladTerm where T
    γ::T
    L::QuExpr
end
"""
    MixedLindbladTerm(γ, (L1, L2))

The "mixed" Lindblad term with rate `γ` and two operators `L1` and `L2`, which gives the symmetrized
``\\frac{1}{2} [γ (L_2^† A L_1 - \\frac{1}{2} \\{L_2^† L_1, A\\}) + γ^* (L_1^† A L_2 - \\frac{1}{2} \\{L_1^† L_2, A\\})]``
in the equation of motion of an operator ``A``.
"""
struct MixedLindbladTerm{T} <: AbstractLindbladTerm where T
    γ::T
    Ls::NTuple{2,QuExpr}
end
"""
    SummedLindbladTerm(inds, L)

The sum of the Lindblad term `L` (a [`SingleLindbladTerm`](@ref) or [`MixedLindbladTerm`](@ref)) over
the indices `inds`, e.g., ``\\sum_i γ (σ^+_i A σ^-_i - \\frac{1}{2} \\{σ^+_i σ^-_i, A\\})``.
"""
struct SummedLindbladTerm{N,T} <: AbstractLindbladTerm where {N,T<:Union{SingleLindbladTerm,MixedLindbladTerm}}
    inds::NTuple{N,QuIndex}
    L::T
end
# L₂† A L₁ - ½{L₂†L₁, A} = L₂† [A,L₁] + ½[L₂†L₁, A], where the commutators are computed directly
# (and the large terms of L₂† A L₁, L₂†L₁ A and A L₂†L₁ that cancel are never generated)
_apply_lindblad(A,γ,L1,L2) = γ*(L2'*comm(A,L1) + 1//2*comm(L2'*L1,A))
(LT::SingleLindbladTerm{T})(A::QuExpr) where T = _apply_lindblad(A,LT.γ,LT.L,LT.L)
_apply_mixedlindblad(A,γ,L1,L2) = _apply_lindblad(A,γ,L1,L2)
function (LT::MixedLindbladTerm{T})(A::QuExpr) where T
    L1, L2 = LT.Ls
    # to ensure that dynamical equations are Hermitian, use a symmetrized form
    1//2*(_apply_mixedlindblad(A,LT.γ,L1,L2) + _apply_mixedlindblad(A,LT.γ',L2,L1))
end
function (LT::SummedLindbladTerm{N,T})(A::QuExpr) where {N,T}
    inds = LT.inds
    # make sure the sum indices do not appear in A by replacing any occurences with
    # temporary indices first, and then replacing the temporary indices back afterwards
    tmpinds = tmpindex.(1:length(inds))
    f1 = replace_inds(inds .=> tmpinds)
    f2 = replace_inds(tmpinds .=> inds)
    ∑(inds,LT.L(f1(A))) |> f2
end

"""
    lindbladterm(args...)

Return the Lindblad term described by `args` (see [`heisenberg_eom`](@ref)): `lindbladterm(L)` and
`lindbladterm(γ, L)` give a [`SingleLindbladTerm`](@ref) (with `γ = 1` if it is not given),
`lindbladterm((L1, L2))` and `lindbladterm(γ, (L1, L2))` a [`MixedLindbladTerm`](@ref), and with an
index or a tuple of indices as the first argument, e.g., `lindbladterm(:i, γ, L)`, the sum of the term
over these indices ([`SummedLindbladTerm`](@ref)).
"""
lindbladterm(L::QuExpr) = SingleLindbladTerm(1,L)
lindbladterm(γ::Union{Number,QuExpr},L::QuExpr) = SingleLindbladTerm(γ,L)
lindbladterm(Ls::NTuple{2,QuExpr}) = MixedLindbladTerm(1,Ls)
lindbladterm(γ::Union{Number,QuExpr},Ls::NTuple{2,QuExpr}) = MixedLindbladTerm(γ,Ls)
lindbladterm(ind::Union{Symbol,QuIndex},args...) = lindbladterm((ind,),args...)
lindbladterm(inds::T,args...) where T<:NTuple{N,Union{Symbol,QuIndex}} where N = SummedLindbladTerm(QuIndex.(inds),lindbladterm(args...))

_lindbladterm(args) = lindbladterm(args...)
_lindbladterm(L::AbstractLindbladTerm) = L

"""
    heisenberg_eom(A, H, Ls=())

Return the right-hand side of the Heisenberg equation of motion of the operator `A` for the
Hamiltonian `H` and the Lindblad terms `Ls` (with ``\\hbar = 1``),
```math
\\frac{dA}{dt} = i [H,A] + \\sum_k γ_k \\left(L_k^† A L_k - \\frac{1}{2} \\{L_k^† L_k, A\\}\\right).
```
The result is not normal ordered (use [`normal_form`](@ref)). Each element of the tuple `Ls`
describes a Lindblad term (see [`lindbladterm`](@ref)): `(L,)` or `(γ, L)` for the jump operator `L`
with rate `γ` (`1` if it is not given), `((L1, L2),)` or `(γ, (L1, L2))` for a "mixed" term with two
operators (see [`MixedLindbladTerm`](@ref)), and `(:i, γ, L)` (or with a tuple of indices) for the
sum of the term over the index `i`.

# Examples
```julia
julia> @boson a; @qparam ω γ;

julia> normal_form(heisenberg_eom(a, ω*a'*a, ((γ, a),)))
-1//2 γ a() - 1i ω a()
```
"""
heisenberg_eom(A::QuExpr,H::QuExpr,Ls::Tuple{Vararg{Union{<:Tuple,AbstractLindbladTerm}}}=()) = _heisenberg_eom(A,H,_lindbladterm.(Ls)...)
_heisenberg_eom(A::QuExpr,H::QuExpr,Ls::AbstractLindbladTerm...) = mapfoldl(L->L(A),+,Ls;init=1im*comm(H,A))
