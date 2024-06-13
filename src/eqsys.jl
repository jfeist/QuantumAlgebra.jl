export droplen, dropcorr

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
