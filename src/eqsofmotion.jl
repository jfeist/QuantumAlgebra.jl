using OrderedCollections

export heisenberg_eom
export EqSys, droplen, dropcorr

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


# equation systems (EqSys)

to_bareterm(A::QuExpr) = ((t,s) = only(A.terms); @assert isone(s); t)
to_opprod(A::QuExpr)::BaseOpProduct = to_opprod(to_bareterm(A))
to_opprod(A::QuTerm)::BaseOpProduct = (@assert A.nsuminds == 0 && isempty(A.δs) && isempty(A.params) && isempty(A.expvals) && isempty(A.corrs); A.bares)

_getECs(A::QuTerm,::Type{Corr})   = (@assert isempty(A.expvals); A.corrs)
_getECs(A::QuTerm,::Type{ExpVal}) = (@assert isempty(A.corrs); A.expvals)

expheis(A::ExpVal,args...) = expheis(QuExpr(QuTerm(A.ops)),args...)
expheis(A::QuExpr,H::QuExpr,Ls=()) = expval(normal_form(heisenberg_eom(A,H,Ls)));

corrheis(A::Corr,args...) = corrheis(A.ops,args...)
corrheis(A::QuExpr,args...) = corrheis(to_opprod(A),args...)
function corrheis(ops::BaseOpProduct,H::QuExpr,Ls=())::QuExpr
    f = canon_inds_remember()
    rhs = _corrheis_cached(f(ops),H,Ls)
    replace_inds(f.replacements)(rhs)
end

const _CORRHEISDICT = Dict()
function _corrheis_cached(ops::BaseOpProduct,H::QuExpr,Ls)::QuExpr
    # if haskey(_CORRHEISDICT,(ops,H,Ls))
    #     println("found operator '$ops' in _CORRHEISDICT")
    # else
    #     println("did not find operator '$ops' in _CORRHEISDICT")
    # end
    get!(_CORRHEISDICT,(ops,H,Ls)) do
        # original right-hand side
        A = QuExpr(QuTerm(ops))
        dA = normal_form(heisenberg_eom(A,H,Ls))
        cdA = expval_as_corrs(dA)

        # terms from left-hand side, i.e., -d_t δ<A>
        δcA = corr(A) - expval_as_corrs(A)
        for (t,s) in δcA.terms
            # there should be no more "bare" operators in these expressions
            @assert isempty(t.bares) && isempty(t.expvals)
            for ii = 1:length(t.corrs)
                tn = deepcopy(t)
                O = tn.corrs[ii]
                deleteat!(tn.corrs,ii)
                tt = normal_form(QuExpr((tn=>s,))*corrheis(O,H,Ls))
                for (t2,s2) in tt.terms
                    _add_sum_term!(cdA,t2,s2)
                end
            end
        end
        cdA
    end
end

droplen(n) = A -> droplen(n,A)
droplen(n,A::QuTerm) = (any(@. length(A.corrs) > n) || any(@. length(A.expvals) > n)) ? 0 : 1
droplen(n,A::QuExpr) = QuExpr((t,droplen(n,t)*s) for (t,s) in A.terms)

dropcorr(n) = A -> dropcorr(n,A)
function dropcorr(n,A::QuExpr)
    Anew = QuExpr()
    terms_to_do = collect(A.terms)
    while !isempty(terms_to_do)
        (t,s) = pop!(terms_to_do)
        if all(@. length(t.expvals) <= n)
            _add_sum_term!(Anew,t,s)
        else
            isempty(t.corrs) && isempty(t.bares) || throw(ArgumentError("t.corrs and t.bares should be empty"))
            iev = findfirst(ev -> length(ev)>n, t.expvals)
            Otmp = QuExpr(QuTerm(t.nsuminds, t.δs, t.params, t.expvals[1:length(t.expvals) .!== iev], t.corrs, t.expvals[iev].ops))
            # express <Otmp> - <Otmp>c in terms of expvals
            Onew = normal_form(expval(Otmp) - corr_as_expvals(Otmp))
            for (tn,sn) in Onew.terms
                push!(terms_to_do, tn => sn*s)
            end
        end
    end
    Anew
end

const EqDict{T<:Union{ExpVal,Corr}} = OrderedDict{T,QuExpr}

struct EqSys{T<:Union{ExpVal,Corr}}
    eqs::EqDict{T}
end

function get_opstodo(ops::Nothing,H)
    opstodo = BaseOperator[]
    for t in keys(H.terms)
        for A in t.bares.v
            # we want to only get the annihilation operators
            if A.t in (BosonCreate_,FermionCreate_,TLSCreate_)
                A = A'
            end
            push!(opstodo,canon_inds()(A))
        end
    end
    sort!(opstodo)
    unique!(opstodo)
    # we actually want to return an array of BaseOpProducts
    [BaseOpProduct([A]) for A in opstodo]
end
get_opstodo(ops::QuExpr,H) = get_opstodo((ops,),H)
get_opstodo(ops,H) = [canon_inds()(to_opprod(normal_form(A))) for A in ops]

EqSys{LHSfunc}(H,maxord::Integer,Ls=(),ops=nothing) where LHSfunc = EqSys{LHSfunc}(H,droplen(maxord),Ls,ops)
function EqSys{LHSfunc}(H,rhsfilter,Ls=(),ops=nothing) where LHSfunc
    Lops = QuantumAlgebra._lindbladterm.(Ls)
    RHSfunc = LHSfunc === ExpVal ? expheis : (LHSfunc === Corr ? corrheis : throw(ArgumentError("LHSfunc must be Corr or ExpVal")))
    opstodo = LHSfunc.(get_opstodo(ops,H))
    eqs = EqDict{LHSfunc}()
    while length(opstodo)>0
        A = pop!(opstodo)
        RHS = rhsfilter(RHSfunc(A,H,Lops))
        # println("$A = $RHS")
        eqs[A] = RHS
        for t in keys(RHS.terms)
            @assert isempty(t.bares)
            for B in _getECs(t,LHSfunc)
                # we will need to calculate either B.A or B.A' (to use <A> = <A^†>^*)
                Bcan = canon_inds()(B)
                # "normal order" Bcand, without worrying about commutation
                # (since this is just about whether we want it in LHS list)
                Bcand = B'
                sort!(Bcand.ops.v)
                Bcand = canon_inds()(Bcand)
                if !(haskey(eqs,Bcan) || haskey(eqs,Bcand) || Bcan in opstodo || Bcand in opstodo)
                    push!(opstodo, nicety(Bcan)>=nicety(Bcand) ? Bcan : Bcand)
                end
            end
        end
    end
    EqSys{LHSfunc}(sort!(eqs))
end

nicety(A::Union{ExpVal,Corr}) = count(A->A.t ∈ (BosonDestroy_,FermionDestroy_,TLSDestroy_), A.ops.v)

function Base.show(io::IO, ::MIME"text/latex", eqsys::EqSys)
    print(io,"\\begin{align}")
    for (A,dAdt) in eqsys.eqs
        print(io,"&\\frac{d}{dt}",latex(A)," = ",latex(dAdt),"\\\\")
    end
    print(io,"\\end{align}")
end
