# helper function used in some definitons, copied from ColorTypes.jl
#Base.@pure basetype(T::Type) = Base.typename(T).wrapper
#Base.@pure basetype(A) = basetype(typeof(A))

macro concrete(expr)
    @assert expr.head == :struct
    S = expr.args[2]
    return quote
        $(esc(expr))

        for n in fieldnames($S)
            if !isconcretetype(fieldtype($S, n))
                error("field $n is not concrete")
            end
        end
    end
end

simplify_number(x::Number) = x
simplify_number(x::AbstractFloat) = isinteger(x) ? Int(x) : x
simplify_number(x::Rational) = isone(denominator(x)) ? numerator(x) : x
simplify_number(x::Complex) = begin
    r, i = reim(x)
    iszero(i) ? simplify_number(r) : complex(simplify_number(r),simplify_number(i))
end
