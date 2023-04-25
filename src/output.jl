export latex

using REPL
using Printf
using Latexify
using LaTeXStrings

const _SUPERSCRIPTNUMS = collect("⁰¹²³⁴⁵⁶⁷⁸⁹")
subscript(i::Integer) = (i<0 ? "₋" : "") * join('₀'+d for d in reverse(digits(abs(i))))
superscript(i::Integer) = (i<0 ? "⁻" : "") * join(_SUPERSCRIPTNUMS[d+1] for d in reverse(digits(abs(i))))

function printspaced(io::IO,v,trailing_space=false)
    isempty(v) && return
    print(io,v[1])
    iprinted = 1
    exponent = 1
    for i = 2:length(v)
        if v[i]==v[iprinted]
            exponent += 1
            # if the last element is repeated, print the exponent
            i == length(v) && print(io,superscript(exponent))
        else
            # print exponent for previous object
            exponent>1 && print(io,superscript(exponent))
            exponent = 1
            print(io," ")
            print(io,v[i])
            iprinted = i
        end
    end
    trailing_space && print(io," ")
end

Base.print(io::IO,A::BaseOpProduct) = printspaced(io,A.v)
function Base.print(io::IO,ii::QuIndex)
    if isnoindex(ii)
        return
    elseif isintindex(ii)
        print(io,ii.num)
    else
        print(io,ii.sym)
        ii.num != typemin(ii.num) && print(io,subscript(ii.num))
    end
end

Base.print(io::IO,A::BaseOperator) = (B = unalias(A); print(io,B.name,BaseOpType_sym[Int(B.t)],"(",B.inds...,")"))
Base.print(io::IO,A::δ) = print(io,"δ(",A.iA,A.iB,")")
Base.print(io::IO,A::Param) = print(io, A.name, A.state=='c' ? "*" : "", length(A.inds)==0 ? "" : "($(A.inds...))")
Base.print(io::IO,A::ExpVal) = print(io,"⟨", A.ops, "⟩")
Base.print(io::IO,A::Corr) = print(io,"⟨", A.ops, "⟩c")

function Base.print(io::IO,A::QuTerm)
    if A.nsuminds > 0
        print(io,"∑")
        for ii=1:A.nsuminds
            print(io,subscript(ii))
        end
        print(io," ")
    end
    isempty(A.δs) || printspaced(io,A.δs,true)
    isempty(A.params) || printspaced(io,A.params,true)
    isempty(A.expvals) || printspaced(io,A.expvals,true)
    isempty(A.corrs) || printspaced(io,A.corrs,true)
    isempty(A.bares) || print(io,A.bares)
end

function print_as_connector(io,str::String,print_connector::Bool)
    if print_connector
        if !isempty(str) && str[1] == '-'
            print(io," - ",str[2:end])
        else
            print(io," + ",str)
        end
    else
        print(io,str)
    end
end

function print_term_scalar(io::IO,t::QuTerm,s::Number, print_connector::Bool)
    if isempty(t)
        print_as_connector(io,string(s),print_connector)
        return
    elseif s == -one(s)
        str = "-"
    elseif s == 1
        str = ""
    else
        str = numstring(s)*" "
    end
    print_as_connector(io,str,print_connector)
    print(io,t)
end

numstring(x::Number) = @sprintf "%g" x
numstring(x::Rational) = denominator(x)==1 ? string(numerator(x)) : string(x)
numstring(x::Complex) = ((r,i) = reim(x); iszero(i) ? numstring(r) : (iszero(r) ? numstring(i)*"i" : "($(numstring(r))$(i<0 ? '-' : '+')$(numstring(abs(i)))i)"))

function Base.print(io::IO,A::QuExpr)
    if isempty(A.terms)
        print(io,0)
    else
        for (ii,(t,s)) in enumerate(sort!(Tuple.(collect(A.terms))))
            print_connector = ii != 1
            print_term_scalar(io,t,s, print_connector)
        end
    end
end

Base.show(io::IO, A::QuantumObject) = print(io,A)
Base.show(io::IO, m::MIME"text/latex", A::QuantumObject) = Base.show(io::IO, m, latexify(A))
Base.show(io::IO, m::MIME"text/latex", A::Vector{<:QuantumObject}) = Base.show(io::IO, m, latexify(A))
Base.show(io::IO, m::MIME"text/latex", A::Matrix{<:QuantumObject}) = Base.show(io::IO, m, latexify(A))
Base.show(io::IO, m::MIME"text/latex", A::NTuple{N,<:QuantumObject}) where N = Base.show(io::IO, m, latexify(A))

const _unicode_to_latex = Dict(v[1]=>"{$k}" for (k,v) in REPL.REPLCompletions.latex_symbols if length(v)==1)
unicode_to_latex(s::AbstractString) = join(map(c -> get(_unicode_to_latex,c,string(c)), collect(s)))
unicode_to_latex(s::QuOpName) = unicode_to_latex(string(s))

function _push_exponents!(ex,itr)
    prevO = nothing
    exponent = 0
    for O in itr
        if O==prevO
            exponent += 1
        else
            exponent>0 && push!(ex.args,exponent==1 ? prevO : LaTeXString("{$(latexify(prevO))}^{$exponent}"))
            exponent = 1
            prevO = O
        end
    end
    # last object
    exponent>0 && push!(ex.args,exponent==1 ? prevO : LaTeXString("{$(latexify(prevO))}^{$exponent}"))
end

@latexrecipe function f(ii::QuIndex)
    isnoindex(ii) && return
    isintindex(ii) && return LaTeXString("$(ii.num)")
    issumindex(ii) && return LaTeXString("\\#_{$(ii.num)}")
    s = unicode_to_latex(string(ii.sym))
    return LaTeXString(ii.num == typemin(ii.num) ? s : "$(s)_{$(ii.num)}")
end
latexindstr(inds) = isempty(inds) ? "" : "_{$(latexify.(inds)...)}"
@latexrecipe function f(A::BaseOperator)
    B = unalias(A)
    return LaTeXString("{$(unicode_to_latex(B.name))}$(latexindstr(B.inds))$(BaseOpType_latex[Int(B.t)])")
end
@latexrecipe function f(A::Param)
    return LaTeXString(string("$(unicode_to_latex(A.name))$(latexindstr(A.inds))", A.state=='c' ? "^{*}" : ""))
end
@latexrecipe function f(A::δ)
    return LaTeXString("\\delta$(latexindstr((A.iA,A.iB)))")
end
@latexrecipe function f(A::BaseOpProduct)
    isempty(A.v) && return
    ex = Expr(:call,:*)
    _push_exponents!(ex,A.v)
    return ex
end
@latexrecipe function f(A::ExpVal)
    return Expr(:latexifymerge, "\\langle ", A.ops, "\\rangle")
end
@latexrecipe function f(A::Corr)
    return Expr(:latexifymerge, "\\langle ", A.ops, "\\rangle_{c}")
end
@latexrecipe function f(A::QuTerm)
    isempty(A) && return
    ex = Expr(:call,:*)
    A.nsuminds > 0 && push!(ex.args,LaTeXString("\\sum_{$(latexify.(sumindex.(1:A.nsuminds))...)}"))
    _push_exponents!(ex,Iterators.flatten((A.δs,A.params,A.expvals,A.corrs,(A.bares,))))
    return ex
end
@latexrecipe function f(A::QuExpr)
    fmt --> FancyNumberFormatter()
    cdot --> false

    ex = 0
    # convert to tuples so sorting ignores the numbers if terms are different (which they are)
    for (t,s) in sort!(Tuple.(collect(A.terms)))
        sign, s = isreal(s) && s<0 ? (-1,-s) : (1,s)
        x = isone(s) ? (isempty(t) ? s : t) : :($s*$t)
        if ex === 0
            ex = isone(sign) ? x : :(-$x)
        else
            ex = isone(sign) ? :($ex + $x) : :($ex - $x)
        end
    end
    return ex isa Expr ? ex : Expr(:latexifymerge,ex)
end

latex(A::QuantumObject) = latexify(A,env=:raw).s