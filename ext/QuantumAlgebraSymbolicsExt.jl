module QuantumAlgebraSymbolicsExt

using QuantumAlgebra
using Symbolics

_issum(x::Num) = (v = Symbolics.value(x); Symbolics.SymbolicUtils.iscall(v) && Symbolics.operation(v) === (+))

QuantumAlgebra.numstring(x::Num) = _issum(x) ? "($(string(x)))" : string(x)
QuantumAlgebra.numlatex(x::Num) = _issum(x) ? Expr(:latexifymerge,"\\left(",x, "\\right)") : x

QuantumAlgebra.num_expr(s::Union{Num,Complex{Num}}) = Symbolics.toexpr(s)

end