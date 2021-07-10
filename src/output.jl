export latex

using REPL
using Printf

function printspaced(io::IO,v,trailing_space=false)
    print(io,v[1])
    for i = 2:length(v)
        print(io," ")
        print(io,v[i])
    end
    trailing_space && print(io," ")
end

Base.print(io::IO,A::BaseOpProduct) = printspaced(io,A.v)
const _SUPERSCRIPTNUMS = collect("⁰¹²³⁴⁵⁶⁷⁸⁹")
subscript(i::Integer) = (i<0 ? "₋" : "") * join('₀'+d for d in reverse(digits(abs(i))))
superscript(i::Integer) = (i<0 ? "⁻" : "") * join(_SUPERSCRIPTNUMS[d+1] for d in reverse(digits(abs(i))))
function Base.print(io::IO,ii::OpIndex)
    if isnoindex(ii)
        return
    elseif isintindex(ii)
        print(io,ii.num)
    else
        print(io,ii.sym)
        ii.num != typemin(ii.num) && print(io,subscript(ii.num))
    end
end

Base.print(io::IO,A::BaseOperator) = print(io,A.name,OpType_sym[Int(A.t)],"(",A.inds...,")")
Base.print(io::IO,A::δ) = print(io,"δ($(A.iA)$(A.iB))")
Base.print(io::IO,A::Param) = print(io,string(A.name, A.state=='c' ? "*" : "",length(A.inds)==0 ? "" : "($(A.inds...))"))
Base.print(io::IO,A::ExpVal) = print(io,"⟨$(A.ops)⟩")
Base.print(io::IO,A::Corr) = print(io,"⟨$(A.ops)⟩c")

function Base.print(io::IO,A::OpTerm)
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

function print_term_scalar(io::IO,t::OpTerm,s::Number, print_connector::Bool)
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

function Base.print(io::IO,A::OpSum)
    if isempty(A.terms)
        print(io,0)
    else
        for (ii,(t,s)) in enumerate(sort!(Tuple.(collect(A.terms))))
            print_connector = ii != 1
            print_term_scalar(io,t,s, print_connector)
        end
    end
end

const QuantumObject = Union{OpIndex,OpName,BaseOperator,Param,BaseOpProduct,ExpVal,Corr,OpTerm,OpSum}
Base.show(io::IO, A::QuantumObject) = print(io,A)
Base.show(io::IO, ::MIME"text/latex", A::QuantumObject) = print(io,"\$",latex(A),"\$")

function latex(ii::OpIndex)
    if isnoindex(ii)
        ""
    elseif isintindex(ii)
        "$(ii.num)"
    elseif issumindex(ii)
        "\\#_{$(ii.num)}"
    elseif ii.num == typemin(ii.num)
        "$(ii.sym)"
    else
        "$(ii.sym)_{$(ii.num)}"
    end
end

latexjoin(args...) = join(latex.(args...))
latexindstr(inds::OpIndex...) = latexindstr(inds)
latexindstr(inds) = isempty(inds) ? "" : "_{$(latexjoin(inds))}"

numlatex(x::Number) = @sprintf "%g" x
numlatex(x::Rational) = denominator(x)==1 ? "$(numerator(x))" : "\\frac{$(numerator(x))}{$(denominator(x))}"
numlatex(x::Complex) = @sprintf "(%g%+gi)" real(x) imag(x)
numlatex(x::Complex{Rational{T}}) where T = @sprintf "\\left(%s%s%si\\right)" numlatex(real(x)) (imag(x)>=0 ? "+" : "-") numlatex(abs(imag(x)))

const _unicode_to_latex = Dict(v[1]=>"{$k}" for (k,v) in REPL.REPLCompletions.latex_symbols if length(v)==1)
unicode_to_latex(s::AbstractString) = join(map(c -> get(_unicode_to_latex,c,string(c)), collect(s)))
unicode_to_latex(s::OpName) = unicode_to_latex(string(s))

latex(x::Number) = imag(x)==0 ? numlatex(real(x)) : (real(x)==0 ? numlatex(imag(x))*"i" : numlatex(x))
latex(A::BaseOperator) = string("{", unicode_to_latex(A.name), "}", latexindstr(A.inds), OpType_latex[Int(A.t)])
latex(A::δ) = string("\\delta", latexindstr(A.iA,A.iB))
latex(A::Param) = string(unicode_to_latex(A.name), latexindstr(A.inds), A.state=='c' ? "^{*}" : "")
latex(A::BaseOpProduct) = join(latex.(A.v))
latex(A::ExpVal) = string("\\langle ", latex(A.ops), " \\rangle")
latex(A::Corr) = string("\\langle ", latex(A.ops), " \\rangle_{c}")

function latex(A::OpTerm)
    strings = String[]
    if A.nsuminds > 0
        push!(strings,"\\sum_{")
        append!(strings,latex.(sumindex.(1:A.nsuminds)))
        push!(strings,"} ")
    end
    append!(strings,latex.(A.δs))
    append!(strings,latex.(A.params))
    append!(strings,latex.(A.expvals))
    append!(strings,latex.(A.corrs))
    push!(strings,latex(A.bares))
    join(strings," ")
end

function latex(A::OpSum)
    isempty(A) && return "0"

    strings = String[]
    # convert to tuples so sorting ignores the numbers if terms are different (which they are)
    for (t,s) in sort!(Tuple.(collect(A.terms)))
        if isempty(t)
            push!(strings,latex(s))
        else
            ls = s==one(s) ? "" : (s==-one(s) ? "-" : latex(s))
            if !isempty(strings) && (ls=="" || (ls[1] ∉ ('-','+')))
                push!(strings," + ")
            end
            push!(strings,ls,latex(t))
        end
    end
    join(strings)
end
