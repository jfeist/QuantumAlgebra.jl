export julia_expression

# to make a variable name from operators
function varname(A::BaseOperator)
    if A.t in (BosonDestroy_,FermionDestroy_)
        string(A.name)
    elseif A.t in (BosonCreate_,FermionCreate_)
        string(A.name)*"ᴴ"
    elseif A.t == TLSx_
        string(A.name)*"ˣ"
    elseif A.t == TLSy_
        string(A.name)*"ʸ"
    elseif A.t == TLSz_
        string(A.name)*"ᶻ"
    elseif A.t == TLSDestroy_
        string(A.name)*"⁻"
    elseif A.t == TLSCreate_
        string(A.name)*"⁺"
    else
        error("should never be reached")
    end
end
varname(A::BaseOpProduct) = Symbol(varname.(A.v)...)
varname(A::Union{ExpVal,Corr}) = varname(A.ops)

# for sum indices use a valid and short symbol that cannot be expressed as a single char, i.e., no risk of collisions with "real" index names
# (the "standard" symbol '#' is not a variable name, so it becomes var"#")
indexpr(ind::OpIndex) = issumindex(ind) ? Symbol(:s̄,subscript(ind.num)) : Symbol(ind)
indexpr(A) = indexpr.(indices(A))

# the varnames optional (and ignored) argument allows to define a closure getexpr(A) = julia_expression(A,my_names)
function julia_expression(A::Param,varnames=nothing)
    # the parsing here allows "parameters" that are functions
    # e.g., with names such as "f(t)"
    x = Meta.parse(string(A.name))
    # is it a function call?
    if x isa Expr
        @assert isempty(indices(A))
    else
        x = :( $x[$(indexpr(A)...)] )
    end
    A.state == 'c' ? :( conj($x) ) : x
end
julia_expression(A::δ,varnames=nothing) = :( I[$(indexpr(A)...)] )
function julia_expression(A::Union{ExpVal,Corr},varnames=nothing)
    x = varname(A)
    if varnames === nothing || x ∈ varnames
        :( $x[$(indexpr(A)...)] )
    else
        # <A> = <A^†>^*
        Ap = A'
        x = varname(Ap)
        x ∈ varnames || throw(ValueError("passed allowed varnames = $(varnames), but neither $(varname(A)) nor $x are valid."))
        # since conjugation changes the order, we have to use the indices of A after conjugation
        :( conj($x[$(indexpr(Ap)...)]) )
    end
end
function julia_expression(A::OpTerm,varnames=nothing,s=1)
    isempty(A.bares) || throw(ArgumentError("Cannot convert term $A with bare operators to a julia expression."))
    sexpr = isone(s) ? [] : :( $(isreal(s) ? Float64(s) : ComplexF64(s)) )
    exprs = [sexpr; julia_expression.(A.δs); julia_expression.(A.params);
             julia_expression.(A.expvals,(varnames,)); julia_expression.(A.corrs,(varnames,))]
    length(exprs)==1 ? exprs[1] : :( *($(exprs...)) )
end
function julia_expression(A::OpSum,varnames=nothing)
    exprs = [julia_expression(term,varnames,s) for (term,s) in sort!(Tuple.(collect(A.terms)))]
    length(exprs)==1 ? exprs[1] : :( +($(exprs...)) )
end
