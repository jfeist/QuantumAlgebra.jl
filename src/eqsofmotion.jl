export heisenberg_eom

abstract type AbstractLindbladTerm end

struct SingleLindbladTerm{T} <: AbstractLindbladTerm where T
    γ::T
    L::QuExpr
end
struct MixedLindbladTerm{T} <: AbstractLindbladTerm where T
    γ::T
    Ls::NTuple{2,QuExpr}
end
struct SummedLindbladTerm{N,T} <: AbstractLindbladTerm where {N,T<:Union{SingleLindbladTerm,MixedLindbladTerm}}
    inds::NTuple{N,QuIndex}
    L::T
end
function (LT::SingleLindbladTerm{T})(A::QuExpr) where T
    L = LT.L
    LdL = L'*L
    LT.γ*(L'*A*L - 1//2*(LdL*A + A*LdL))
end
_apply_mixedlindblad(A,γ,L1,L2) = (L2dL1 = L2'*L1; γ*(L2'*A*L1 - 1//2*(L2dL1*A + A*L2dL1)))
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

lindbladterm(L::QuExpr) = SingleLindbladTerm(1,L)
lindbladterm(γ::Union{Number,QuExpr},L::QuExpr) = SingleLindbladTerm(γ,L)
lindbladterm(Ls::NTuple{2,QuExpr}) = MixedLindbladTerm(1,Ls)
lindbladterm(γ::Union{Number,QuExpr},Ls::NTuple{2,QuExpr}) = MixedLindbladTerm(γ,Ls)
lindbladterm(ind::Union{Symbol,QuIndex},args...) = lindbladterm((ind,),args...)
lindbladterm(inds::T,args...) where T<:NTuple{N,Union{Symbol,QuIndex}} where N = SummedLindbladTerm(QuIndex.(inds),lindbladterm(args...))

_lindbladterm(args) = lindbladterm(args...)
_lindbladterm(L::AbstractLindbladTerm) = L

"""`heisenberg_eom(A,H,Ls=())` calculates ``\\dot{A} = i [H,A] + \\sum_i (L_i^† A L_i - ½ \\{L_i^† L_i, A\\})``, where `Ls = (L_1,L_2,...)` is an iterable of Lindblad operators."""
heisenberg_eom(A::QuExpr,H::QuExpr,Ls::Tuple{Vararg{Union{<:Tuple,AbstractLindbladTerm}}}=()) = _heisenberg_eom(A,H,_lindbladterm.(Ls)...)
_heisenberg_eom(A::QuExpr,H::QuExpr,Ls::AbstractLindbladTerm...) = mapfoldl(L->L(A),+,Ls;init=1im*comm(H,A))
