module QuantumAlgebraSymPyExt

using QuantumAlgebra
using SymPy

QuantumAlgebra.numstring(x::Sym) = x.o.is_Add ? "($(string(x)))" : string(x)
QuantumAlgebra.numlatex(x::Sym) = x.o.is_Add ? Expr(:latexifymerge, "\\left(", x, "\\right)") : x

end