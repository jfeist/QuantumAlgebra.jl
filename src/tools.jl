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
simplify_number(x::AbstractFloat) = isinteger(x) && typemin(Int) <= x <= typemax(Int) ? Int(x) : x
simplify_number(x::Rational) = isone(denominator(x)) ? numerator(x) : x
simplify_number(x::Complex) = begin
    r, i = reim(x)
    iszero(i) ? simplify_number(r) : complex(simplify_number(r),simplify_number(i))
end

## tables that are shared between threads

# A reference that can be read from several threads while it is replaced (by one thread at a time,
# under a lock). The values stored in it must not be modified after they have been stored. T should
# be a mutable type, since atomic fields of (larger) immutable types are implemented with a lock.
mutable struct _AtomicRef{T}
    @atomic x::T
end
Base.getindex(r::_AtomicRef) = @atomic :acquire r.x
Base.setindex!(r::_AtomicRef, x) = (@atomic :release r.x = x; x)

const _NOTFOUND = gensym(:notfound)

"""
    _COWDict{K,V}()

A dictionary that can be read from several threads without locking. It is never modified in place,
but replaced as a whole (under a lock) when entries are added, so it is only meant for small tables
(such as caches with few entries).
"""
struct _COWDict{K,V}
    d::_AtomicRef{Dict{K,V}}
    lock::ReentrantLock
end
_COWDict{K,V}(d::Dict{K,V}=Dict{K,V}()) where {K,V} = _COWDict{K,V}(_AtomicRef(d), ReentrantLock())
Base.get(c::_COWDict, k, default) = get(c.d[], k, default)
Base.setindex!(c::_COWDict, v, k) = @lock c.lock (c.d[] = setindex!(copy(c.d[]), v, k); v)
function Base.get!(f::Base.Callable, c::_COWDict{K,V}, k) where {K,V}
    v = get(c.d[], k, _NOTFOUND)
    v === _NOTFOUND || return v::V
    @lock c.lock begin
        v = get(c.d[], k, _NOTFOUND)
        v === _NOTFOUND || return v::V
        v = convert(V, f())::V
        # (f can add entries, so copy the current table only now)
        c.d[] = setindex!(copy(c.d[]), v, k)
        v
    end
end
