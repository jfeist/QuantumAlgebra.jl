module QuantumAlgebraSymPyPythonCallExt

using QuantumAlgebra
using SymPyPythonCall

QuantumAlgebra.numstring(x::Sym) = Bool(x.o.is_Add) ? "($(string(x)))" : string(x)
QuantumAlgebra.numlatex(x::Sym) = Bool(x.o.is_Add) ? Expr(:latexifymerge, "\\left(", x, "\\right)") : x

QuantumAlgebra.num_expr(s::Sym) = convert(Expr, s)

end