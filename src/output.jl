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

subscript(i::Integer) = (i<0 ? "₋" : "") * join('₀'+d for d in reverse(digits(abs(i))))
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
    A.t == σ_ && return print(io,"σ$(A.a)($(A.inds...))")
    A.t == σminus_ && return print(io,"σ⁻($(A.inds...))")
    A.t == σplus_ && return print(io,"σ⁺($(A.inds...))")
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
        str = string(s)
    end
    print_as_connector(io,str,print_connector)
    print(io,t)
end

function Base.print(io::IO,A::OpSum)
    if isempty(A.terms)
        print(io,0)
    else
        for (ii,(t,s)) in enumerate(A.terms)
            print_connector = ii != 1
            print_term_scalar(io,t,s, print_connector)
        end
    end
end

mystring(x::Number) = @sprintf "%g" x
mystring(x::Rational) = denominator(x)==1 ? "$(numerator(x))" : "\\frac{$(numerator(x))}{$(denominator(x))}"
mystring(x::Complex) = @sprintf "(%g%+gi)" real(x) imag(x)
mystring(x::Complex{Rational{T}}) where T = @sprintf "\\left(%s%s%si\\right)" mystring(real(x)) (imag(x)>=0 ? "+" : "-") mystring(abs(imag(x)))

Base.show(io::IO, ::MIME"text/latex", A::Union{BaseOperator,BaseOpProduct,ExpVal,Corr,OpSum,OpTerm}) = print(io,"\$",latex(A),"\$")

function latex(A::BaseOperator)
    A.t in (BosonDestroy_,FermionDestroy_) && return "{$(A.name)}$(latexindstr(A.inds))"
    A.t in (BosonCreate_,FermionCreate_) && return "{$(A.name)}$(latexindstr(A.inds))^\\dagger"
    A.t == σ_ && return string("\\sigma_{$(A.a)",length(A.inds)>0 ? ",$(A.inds...)}" : "}")
    A.t == σminus_ && return "\\sigma^-" * latexindstr(A.inds)
    A.t == σplus_ && return "\\sigma^+" * latexindstr(A.inds)
    error("latex should not reach this!)")
end

function latex(ii::OpIndex)
    if isintindex(ii)
        "{$(ii.num)}"
    elseif issumindex(ii)
        "{\\#_{$(ii.num)}}"
    elseif ii.num == typemin(ii.num)
        "{$(ii.sym)}"
    else
        "{$(ii.sym)_{$(ii.num)}}"
    end
end

latexjoin(args...) = join(latex.(args...))
latexindstr(inds::OpIndex...) = latexindstr(inds)
latexindstr(inds) = isempty(inds) ? "" : "_{$(latexjoin(inds))}"

latex(x::Number) = imag(x)==0 ? mystring(real(x)) : (real(x)==0 ? mystring(imag(x))*"i" : mystring(x))
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
