export latex

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

function Base.print(io::IO,A::BaseOperator)
    A.t in (BosonDestroy_,FermionDestroy_) && return print(io,"$(A.name)($(A.inds...))")
    A.t in (BosonCreate_,FermionCreate_) && return print(io,"$(A.name)†($(A.inds...))")
    A.t == TLSx_ && return print(io,"$(A.name)ˣ($(A.inds...))")
    A.t == TLSy_ && return print(io,"$(A.name)ʸ($(A.inds...))")
    A.t == TLSz_ && return print(io,"$(A.name)ᶻ($(A.inds...))")
    A.t == TLSDestroy_ && return print(io,"$(A.name)⁻($(A.inds...))")
    A.t == TLSCreate_ && return print(io,"$(A.name)⁺($(A.inds...))")
    error("print should never reach this!")
end

Base.print(io::IO,A::δ) = print(io,"δ($(A.iA)$(A.iB))")
Base.print(io::IO,A::ExpVal) = print(io,"⟨$(A.ops)⟩")
Base.print(io::IO,A::Corr) = print(io,"⟨$(A.ops)⟩c")
Base.print(io::IO,A::Param) = print(io,string(A.name, A.state=='c' ? "*" : "",length(A.inds)==0 ? "" : "($(A.inds...))"))

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

function latex(A::BaseOperator)
    A.t in (BosonDestroy_,FermionDestroy_) && return "{$(A.name)}$(latexindstr(A.inds))"
    A.t in (BosonCreate_,FermionCreate_) && return "{$(A.name)}$(latexindstr(A.inds))^\\dagger"
    A.t == TLSx_ && return "\\sigma^x" * latexindstr(A.inds)
    A.t == TLSy_ && return "\\sigma^y" * latexindstr(A.inds)
    A.t == TLSz_ && return "\\sigma^z" * latexindstr(A.inds)
    A.t == TLSDestroy_ && return "\\sigma^-" * latexindstr(A.inds)
    A.t == TLSCreate_ && return "\\sigma^+" * latexindstr(A.inds)
    error("latex should not reach this!)")
end

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

latex(x::Number) = imag(x)==0 ? numlatex(real(x)) : (real(x)==0 ? numlatex(imag(x))*"i" : numlatex(x))
latex(A::δ) = string("\\delta", latexindstr(A.iA,A.iB))
latex(A::ExpVal) = "\\langle $(latex(A.ops)) \\rangle"
latex(A::Corr) = "\\langle $(latex(A.ops)) \\rangle_{c}"
latex(A::Param) = string(A.name, latexindstr(A.inds), A.state=='c' ? "^{*}" : "")
latex(A::BaseOpProduct) = join(latex.(A.v))

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
    if isempty(A.terms)
        "0"
    else
        strings = String[]
        # convert to tuples so sorting ignores the numbers if terms are different (which they are)
        for (ii,(t,s)) in enumerate(sort!(Tuple.(collect(A.terms))))
            if isempty(t)
                push!(strings,latex(s))
            else
                ls = s==one(s) ? "" : (s==-one(s) ? "-" : latex(s))
                if ii>1 && (ls=="" || (ls[1] ∉ ('-','+')))
                    push!(strings," + ")
                end
                push!(strings,ls,latex(t))
            end
        end
        join(strings)
    end
end
