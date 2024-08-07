export julia_expression

# to make a variable name from operators
varname(A::BaseOperator) = string(A.name,BaseOpType_expr[Int(A.t)])
varname(A::BaseOpProduct) = Symbol(varname.(A.v)...)
varname(A::Union{ExpVal,Corr}) = varname(A.ops)

# for sum indices use a valid and short symbol that cannot be expressed as a single char, i.e., no risk of collisions with "real" index names
# (the "standard" symbol '#' is not a variable name, so it becomes var"#")
indexpr(ind::QuIndex) = issumindex(ind) ? Symbol(:s̄,subscript(ind.num)) : (isintindex(ind) ? ind.num : Symbol(ind))
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
    A = unalias(A)
    x = varname(A)
    if varnames === nothing || x ∈ varnames
        :( $x[$(indexpr(A)...)] )
    else
        # <A> = <A^†>^*
        Ap = A'
        x = varname(Ap)
        x ∈ varnames || throw(ArgumentError("Neither $(repr(varname(A))) nor $(repr(x)) are in the list of allowed variable names, varnames = $(varnames)."))
        # since conjugation changes the order, we have to use the indices of A after conjugation
        :( conj($x[$(indexpr(Ap)...)]) )
    end
end
num_expr(s::Number) = :( $(isreal(s) ? Float64(s) : ComplexF64(s)) )
function julia_expression(A::QuTerm,varnames=nothing,s=1)
    isempty(A) && return num_expr(s)
    isempty(A.bares) || throw(ArgumentError("Cannot convert term $A with bare operators to a julia expression."))
    sexpr = isone(s) ? [] : num_expr(s)
    exprs = [sexpr; julia_expression.(A.δs); julia_expression.(A.params);
             julia_expression.(A.expvals,(varnames,)); julia_expression.(A.corrs,(varnames,))]
    length(exprs)==1 ? exprs[1] : :( *($(exprs...)) )
end
function julia_expression(A::QuExpr,varnames=nothing)
    isempty(A) && return 0
    exprs = [julia_expression(term,varnames,s) for (term,s) in sort!(Tuple.(collect(A.terms)))]
    length(exprs)==1 ? exprs[1] : :( +($(exprs...)) )
end
julia_expression(A::Number,varnames=nothing) = julia_expression(QuExpr(A),varnames)
