module QuantumAlgebraSymbolicsExt

using QuantumAlgebra
using Symbolics

QuantumAlgebra.numstring(x::Num) = iscall(Symbolics.value(x)) && operation(Symbolics.value(x))===(+) ? "($(string(x)))" : string(x)
QuantumAlgebra.numlatex(x::Num) = iscall(Symbolics.value(x)) && operation(Symbolics.value(x))===(+) ? Expr(:latexifymerge,"\\left(",x, "\\right)") : x

QuantumAlgebra.num_expr(s::Union{Num,Complex{Num}}) = Symbolics.toexpr(s)

end