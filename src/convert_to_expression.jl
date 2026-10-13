export julia_expression

# to make a variable name from operators
varname(A::BaseOperator) = string(opname(A),opsymbol(A,:expr))
varname(T::_QuditTransition) = string(opname(T.create),_qudit_symbol(_levels(T)...,:expr))
varname(A::OpProduct) = Symbol((varname(O) for O in _display_items(A))...)
varname(A::Union{ExpVal,Corr}) = varname(ops(A))

# for sum indices use a valid and short symbol that cannot be expressed as a single char, i.e., no risk of collisions with "real" index names
# (the "standard" symbol '#' is not a variable name, so it becomes var"#")
indexpr(ind::QuIndex) = issumindex(ind) ? Symbol(:s̄,subscript(ind.num)) : (isintindex(ind) ? ind.num : Symbol(ind))
indexpr(A) = indexpr.(indices(A))

"""
    julia_expression(A, varnames=nothing)

Return a Julia expression (an `Expr`) that computes the value of `A`, a `QuExpr` without bare
operators (e.g., with expectation values or correlators). Expectation values and correlators become
arrays with a name made from their operators, e.g., ``⟨b^†_i b_j⟩`` becomes `bᴴb[i,j]`, parameters
become arrays with their names (or, for names like `f(t)`, function calls), Kronecker deltas become
`I[i,j]` (from LinearAlgebra), and summation indices are named `s̄₁`, `s̄₂`, … (but the sums are not
written explicitly). If the collection of symbols `varnames` is given, expectation values (or
correlators) whose name is not in it are written as the complex conjugate of their adjoint (using
``⟨A⟩ = ⟨A^†⟩^*``), e.g., `julia_expression(expval(b[:i]), [:bᴴ])` gives `conj(bᴴ[i])`.

# Examples
```julia
julia> @boson b[:]; @qparam g[:];

julia> julia_expression(expval_as_corrs(g[:i]*b'[:i]*b[:j]))
:(g[i] * bᴴ[i] * b[j] + g[i] * bᴴb[i, j])
```
"""
function julia_expression(A::Param,varnames=nothing)
    # the parsing here allows "parameters" that are functions
    # e.g., with names such as "f(t)"
    x = Meta.parse(string(opname(A)))
    # is it a function call?
    if x isa Expr
        @assert isempty(indices(A))
    else
        x = :( $x[$(indexpr(A)...)] )
    end
    parstate(A) == 'c' ? :( conj($x) ) : x
end
julia_expression(A::δ,varnames=nothing) = :( I[$(indexpr(A)...)] )
function julia_expression(A::Union{ExpVal,Corr},varnames=nothing)
    A = unalias(A)
    x = varname(A)
    if varnames === nothing || x ∈ varnames
        :( $x[$(indexpr.(_display_indices(ops(A)))...)] )
    else
        # <A> = <A^†>^*
        Ap = A'
        x = varname(Ap)
        x ∈ varnames || throw(ArgumentError("Neither $(repr(varname(A))) nor $(repr(x)) are in the list of allowed variable names, varnames = $(varnames)."))
        # since conjugation changes the order, we have to use the indices of A after conjugation
        :( conj($x[$(indexpr.(_display_indices(ops(Ap)))...)]) )
    end
end
"""
    num_expr(x::Number)

The Julia expression for the coefficient `x` in [`julia_expression`](@ref), which is a `Float64` or a
`ComplexF64` by default. Packages that provide other number types (e.g., symbolic numbers) can add
methods for them.
"""
num_expr(s::Number) = :( $(isreal(s) ? Float64(s) : ComplexF64(s)) )
function julia_expression(A::QuTerm,varnames=nothing,s=1)
    isempty(A) && return num_expr(s)
    isempty(bares(A)) || throw(ArgumentError("Cannot convert term $A with bare operators to a julia expression."))
    sexpr = isone(s) ? [] : num_expr(s)
    exprs = [sexpr; julia_expression.(deltas(A)); julia_expression.(params(A));
             julia_expression.(expvals(A),(varnames,)); julia_expression.(corrs(A),(varnames,))]
    length(exprs)==1 ? exprs[1] : :( *($(exprs...)) )
end
function julia_expression(A::QuExpr,varnames=nothing)
    isempty(A) && return 0
    exprs = [julia_expression(term,varnames,s) for (term,s) in sort!(Tuple.(collect(A.terms)))]
    length(exprs)==1 ? exprs[1] : :( +($(exprs...)) )
end
julia_expression(A::Number,varnames=nothing) = julia_expression(QuExpr(A),varnames)
