import Base: print
export latex

using Printf

print(io::IO,A::Union{BosonDestroy,FermionDestroy}) = print(io,"$(A.name)($(A.inds...))")
print(io::IO,A::Union{BosonCreate,FermionCreate}) = print(io,"$(A.name)†($(A.inds...))")
print(io::IO,A::σ) = print(io,"σ$(A.a)($(A.inds...))")
print(io::IO,A::σminus) = print(io,"σ-($(A.inds...))")
print(io::IO,A::σplus) = print(io,"σ+($(A.inds...))")
print(io::IO,A::scal) = print(io,imag(A.v)==0 ? @sprintf("%g",real(A.v)) : (real(A.v)==0 ? @sprintf("%gi",imag(A.v)) : @sprintf("(%g%+gi)",real(x),imag(x))))
print(io::IO,A::param) = print(io,string(A.name, A.state=='c' ? "'" : "",length(A.inds)==0 ? "" : "($(A.inds...))"))
print(io::IO,A::OpProd) = print(io,"$(A.A) $(A.B)")
print(io::IO,A::OpSum) = print(io,"$(A.A) + $(A.B)")
print(io::IO,A::ExpVal) = print(io,"⟨$(A.A)⟩")
print(io::IO,A::Corr) = print(io,"⟨$(A.A)⟩c")
print(io::IO,A::OpSumAnalytic) = print(io,"∑_$(A.ind) $(A.A)")

mystring(x::Number) = @sprintf "%g" x
mystring(x::Rational) = denominator(x)==1 ? "$(numerator(x))" : "\\frac{$(numerator(x))}{$(denominator(x))}"
mystring(x::Complex) = @sprintf "(%g%+gi)" real(x) imag(x)
mystring(x::Complex{Rational{T}}) where T = @sprintf "\\left(%s%s%si\\right)" mystring(real(x)) (imag(x)>=0 ? "+" : "-") mystring(abs(imag(x)))

Base.show(io::IO, ::MIME"text/latex", A::Operator) = print(io,"\$",latex(A),"\$")
latex(A::σ) = string("\\sigma_{$(A.a)",length(A.inds)>0 ? ",$(A.inds...)}" : "}")
latexindstr(inds) = length(inds)==0 ? "" : "_{$(inds...)}"
latex(A::δ) = "δ" * latexindstr(A.inds)
latex(A::Union{BosonDestroy,FermionDestroy}) = "{$(A.name)}$(latexindstr(A.inds))"
latex(A::Union{BosonCreate,FermionCreate}) = "{$(A.name)}$(latexindstr(A.inds))^\\dagger"
latex(A::σminus) = "\\sigma^-" * latexindstr(A.inds)
latex(A::σplus) = "\\sigma^+" * latexindstr(A.inds)
latex(A::scal) = imag(A.v)==0 ? mystring(real(A.v)) : (real(A.v)==0 ? mystring(imag(A.v))*"i" : mystring(A.v))
latex(A::param)  = string(A.name, latexindstr(A.inds), A.state=='c' ? "^*" : "")
latex(A::OpProd) = string(latex(A.A)," ",latex(A.B))
latex(A::OpSum) = string(latex(A.A)," + ",latex(A.B))
latex(A::ExpVal) = "\\langle $(latex(A.A)) \\rangle"
latex(A::Corr) = "\\langle $(latex(A.A)) \\rangle_{c}"
latex(A::OpSumAnalytic) = string("\\sum_{$(A.ind)}",latex(A.A))
