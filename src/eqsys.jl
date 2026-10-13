export heisenberg_eom_system, droplen, dropcorr, ExpVal, Corr

# equation systems (QuEqSys)

to_bareterm(A::QuExpr) = ((t,s) = only(A.terms); @assert isone(s); t)
to_opprod(A::QuExpr) = to_opprod(to_bareterm(A))
to_opprod(A::QuTerm) = (@assert nsuminds(A) == 0 && isempty(deltas(A)) && isempty(params(A)) && isempty(expvals(A)) && isempty(corrs(A)); bares(A))

_getECs(A::QuTerm,::Type{Corr})   = (@assert isempty(expvals(A)); corrs(A))
_getECs(A::QuTerm,::Type{ExpVal}) = (@assert isempty(corrs(A)); expvals(A))

"""
    expheis(A, H, Ls=())

Return the right-hand side of the equation of motion of the expectation value ``⟨A⟩`` (for an
expectation value `A` or an expression of operators), i.e., the expectation value of the normal form
of [`heisenberg_eom`](@ref)`(A, H, Ls)`.
"""
expheis(A::ExpVal,args...) = expheis(QuExpr(QuTerm(ops(A))),args...)
expheis(A::QuExpr,H::QuExpr,Ls=()) = expval(normal_form(heisenberg_eom(A,H,Ls)));

"""
    corrheis(A, H, Ls=())

Return the right-hand side of the equation of motion of the correlator (cumulant) ``⟨A⟩_c`` of a
product of operators `A` (or for a correlator `A`), expressed through correlators. The results are
cached, so that the equations of the correlators that appear in the equations of higher-order
correlators are only calculated once.
"""
corrheis(A::Corr,args...) = corrheis(ops(A),args...)
corrheis(A::QuExpr,args...) = corrheis(to_opprod(A),args...)
function corrheis(ops::OpProduct,H::QuExpr,Ls=())::QuExpr
    f = canon_inds_remember()
    rhs = _corrheis_cached(f(ops),H,Ls)
    replace_inds(f.replacements)(rhs)
end

# cache of the equations of motion of correlators (the lock is only held while accessing the cache,
# not during the calculation, so that the same equation might be calculated twice by different threads)
# the keys store their hash, since hashing H can be expensive and the cache is accessed twice
struct _HashedKey{T}
    key::T
    hash::UInt
end
_HashedKey(key) = _HashedKey(key, hash(key))
Base.hash(k::_HashedKey, h::UInt) = hash(k.hash, h)
Base.isequal(a::_HashedKey, b::_HashedKey) = a.hash == b.hash && isequal(a.key, b.key)
const _CORRHEISDICT = Dict{_HashedKey,QuExpr}()
const _CORRHEISLOCK = ReentrantLock()
function _corrheis_cached(ops::BaseOpProduct,H::QuExpr,Ls)::QuExpr
    key = _HashedKey((ops,H,Ls))
    cdA = @lock _CORRHEISLOCK get(_CORRHEISDICT,key,nothing)
    cdA === nothing || return cdA
    cdA = _corrheis_calc(ops,H,Ls)
    @lock _CORRHEISLOCK get!(_CORRHEISDICT,key,cdA)
end
function _corrheis_calc(ops::BaseOpProduct,H::QuExpr,Ls)
    # original right-hand side
    A = QuExpr(QuTerm(ops))
    dA = normal_form(heisenberg_eom(A,H,Ls))
    cdA = expval_as_corrs(dA)

    # terms from left-hand side, i.e., -d_t δ<A>
    δcA = corr(A) - expval_as_corrs(A)
    for (t,s) in δcA.terms
        # there should be no more "bare" operators in these expressions
        @assert isempty(bares(t)) && isempty(expvals(t))
        tcorrs = collect(corrs(t))
        for ii in eachindex(tcorrs)
            O = tcorrs[ii]
            tn = QuTerm(nsuminds(t),deltas(t),params(t),expvals(t),tcorrs[eachindex(tcorrs) .!= ii],bares(t))
            tt = normal_form(QuExpr((tn=>s,))*corrheis(O,H,Ls))
            for (t2,s2) in tt.terms
                _add_sum_term!(cdA,t2,s2)
            end
        end
    end
    cdA
end

"""
    droplen(n)

Create a filter that drops correlators/expectation values with length greater than `n`.
Useful for truncating hierarchies of equations.
"""
droplen(n) = Base.Fix1(droplen,n)
droplen(n,A::QuTerm) = (any(C -> length(C) > n, corrs(A)) || any(E -> length(E) > n, expvals(A))) ? 0 : 1
droplen(n,A::QuExpr) = QuExpr((t,droplen(n,t)*s) for (t,s) in A.terms)

"""
    dropcorr(n::Int)

Returns a filter function that rewrites terms with correlations `corrs` or
expectation values `expvals` of length greater than `n` in terms of lower-order
expressions (up to order `n`) and higher-order correlators, and then drops those
higher-order correlations.

For example, `dropcorr(1)` will replace
``⟨a^† a⟩ = ⟨a^† a⟩_c + ⟨a^†⟩ ⟨a⟩ ≈ ⟨a^†⟩ ⟨a⟩``.
"""
dropcorr(n) = Base.Fix1(dropcorr,n)
function dropcorr(n,A::QuExpr)
    Anew = QuExpr()
    terms_to_do = collect(A.terms)
    while !isempty(terms_to_do)
        (t,s) = pop!(terms_to_do)
        if all(E -> length(E) <= n, expvals(t))
            _add_sum_term!(Anew,t,s)
        else
            isempty(corrs(t)) && isempty(bares(t)) || throw(ArgumentError("corrs(t) and bares(t) should be empty"))
            tevs = collect(expvals(t))
            iev = findfirst(ev -> length(ev)>n, tevs)
            Otmp = QuExpr(QuTerm(nsuminds(t), deltas(t), params(t), tevs[eachindex(tevs) .!= iev], corrs(t), ops(tevs[iev])))
            # express <Otmp> - <Otmp>c in terms of expvals
            Onew = normal_form(expval(Otmp) - corr_as_expvals(Otmp))
            for (tn,sn) in Onew.terms
                push!(terms_to_do, tn => sn*s)
            end
        end
    end
    Anew
end

function get_opstodo(ops::Nothing,H)
    opstodo = BaseOperator[]
    for t in keys(H.terms)
        for A in bares(t)
            # we want to only get the annihilation operators
            if ladder(A) > 0
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

"""
    heisenberg_eom_system(H, rhsfilt, Ls=(), ops=nothing)
    heisenberg_eom_system(T, H, rhsfilt, Ls=(), ops=nothing)

Return the system of Heisenberg equations of motion for the expectation values (for `T = ExpVal`,
the default) or the correlators (cumulants, for `T = Corr`) of operator products, for the
Hamiltonian `H` and the Lindblad terms `Ls` (see [`heisenberg_eom`](@ref)), as a
[`QuEqSys`](@ref).

Starting from the operators `ops` (an operator, or a tuple of operators or products of operators),
the equations of all expectation values (or correlators) that appear on the right-hand sides are
added until the system is closed. If `ops` is `nothing` (the default), the system starts from the
annihilation operators that appear in `H`. Since the equations of products of ``n`` operators
usually contain products of more operators, the system is truncated by the filter function
`rhsfilt`, which is applied to the right-hand side of each equation and removes higher-order terms
or approximates them by lower-order ones (see [`droplen`](@ref) and [`dropcorr`](@ref)). An integer
`n` stands for `droplen(n)`.

# Examples
```julia
julia> @boson a; @qparam ω χ γ;

julia> H = ω*a'*a + χ*a'*(a'+a)*a;

julia> heisenberg_eom_system(H, 2, ((γ, a),), a)
dₜ⟨a()⟩ = -1//2 γ ⟨a()⟩ - 1i ω ⟨a()⟩ - 2i χ ⟨a†() a()⟩ - 1i χ ⟨a()²⟩
dₜ⟨a†() a()⟩ = -γ ⟨a†() a()⟩
dₜ⟨a()²⟩ = -2i χ ⟨a()⟩ - γ ⟨a()²⟩ - 2i ω ⟨a()²⟩
```
"""
heisenberg_eom_system(H::QuExpr,rhsfilt,Ls=(),ops=nothing) = heisenberg_eom_system(ExpVal,H,rhsfilt,Ls,ops)
heisenberg_eom_system(::Type{LHSfunc},H::QuExpr,maxord::Integer,Ls=(),ops=nothing) where LHSfunc = heisenberg_eom_system(LHSfunc,H,droplen(maxord),Ls,ops)
heisenberg_eom_system(::Type{LHSfunc},H::QuExpr,rhsfilter,Ls=(),ops=nothing) where LHSfunc = QuEqSys{LHSfunc}(H,rhsfilter,Ls,ops)

function QuEqSys{LHSfunc}(H,rhsfilter,Ls=(),lhsops=nothing) where LHSfunc
    Lops = QuantumAlgebra._lindbladterm.(Ls)
    RHSfunc = LHSfunc === ExpVal ? expheis : (LHSfunc === Corr ? corrheis : throw(ArgumentError("LHSfunc must be Corr or ExpVal")))
    opstodo = LHSfunc.(get_opstodo(lhsops,H))
    eqs = EqDict{LHSfunc}()
    while length(opstodo)>0
        A = pop!(opstodo)
        RHS = rhsfilter(RHSfunc(A,H,Lops))
        # println("$A = $RHS")
        eqs[A] = RHS
        # go through the terms in sorted order so that the result does not depend on the (hash-dependent) iteration order
        for t in sort!(collect(keys(RHS.terms)))
            @assert isempty(bares(t))
            for B in _getECs(t,LHSfunc)
                # we will need to calculate either B.A or B.A' (to use <A> = <A^†>^*)
                Bcan = canon_inds()(B)
                # "normal order" Bcand, without worrying about commutation
                # (since this is just about whether we want it in LHS list)
                Bcand = LHSfunc(BaseOpProduct(sort!(collect(ops(B')))))
                Bcand = canon_inds()(Bcand)
                if !(haskey(eqs,Bcan) || haskey(eqs,Bcand) || Bcan in opstodo || Bcand in opstodo)
                    push!(opstodo, nicety(Bcan)>=nicety(Bcand) ? Bcan : Bcand)
                end
            end
        end
    end
    QuEqSys{LHSfunc}(sort!(eqs))
end

nicety(A::Union{ExpVal,Corr}) = count(O->ladder(O) < 0, ops(A))

Base.getindex(EQ::QuEqSys{LHSfunc},A::QuExpr) where LHSfunc = getindex(EQ,to_bareterm(A))
Base.getindex(EQ::QuEqSys{LHSfunc},A::QuTerm) where LHSfunc = begin
    EC = isempty(bares(A)) ? only(_getECs(A,LHSfunc)) : LHSfunc(to_opprod(A))
    getindex(EQ,EC)
end
Base.getindex(EQ::QuEqSys{LHSfunc},A::LHSfunc) where LHSfunc = EQ.eqs[A]
Base.getindex(EQ::QuEqSys{LHSfunc},i::Int) where LHSfunc = EQ.eqs[collect(keys(EQ.eqs))[i]]
Base.length(EQ::QuEqSys) = length(EQ.eqs)
