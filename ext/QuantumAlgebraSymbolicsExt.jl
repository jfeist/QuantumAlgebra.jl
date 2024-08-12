module QuantumAlgebraSymbolicsExt

using QuantumAlgebra
using Symbolics

@static if hasproperty(Symbolics.SymbolicUtils,:iscall)
    iscall = Symbolics.SymbolicUtils.iscall
else
    iscall = Symbolics.SymbolicUtils.istree
end

_issum(x::Num) = (v = Symbolics.value(x); iscall(v) && Symbolics.operation(v) === (+))

QuantumAlgebra.numstring(x::Num) = _issum(x) ? "($(string(x)))" : string(x)
QuantumAlgebra.numlatex(x::Num) = _issum(x) ? Expr(:latexifymerge,"\\left(",x, "\\right)") : x

QuantumAlgebra.num_expr(s::Union{Num,Complex{Num}}) = Symbolics.toexpr(s)

Symbolics.simplify(ex::QuExpr; kwargs...) = map_scalar_function(x -> Symbolics.simplify(x; kwargs...), ex)
# the definitions for Pair and Vector avoid method call ambiguities due to substitute(expr, s::{Pair,Vector})
for T in (Any,Pair,Vector)
    @eval Symbolics.substitute(ex::QuExpr, s::$T; kwargs...) = map_scalar_function(ex) do x
        Symbolics.substitute(x, s; kwargs...)
    end
end

end