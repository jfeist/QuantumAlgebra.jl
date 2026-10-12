export latex

using REPL
using Printf
using Latexify
using LaTeXStrings

const _SUPERSCRIPTNUMS = collect("⁰¹²³⁴⁵⁶⁷⁸⁹")
subscript(i::Integer) = (i<0 ? "₋" : "") * join('₀'+d for d in reverse(digits(abs(i))))
superscript(i::Integer) = (i<0 ? "⁻" : "") * join(_SUPERSCRIPTNUMS[d+1] for d in reverse(digits(abs(i))))

# print the elements of v separated by spaces, with repeated elements printed as powers
function printspaced(io::IO,v,trailing_space=false)
    isempty(v) && return
    prev = first(v)
    exponent = 0
    for x in v
        if exponent > 0 && x == prev
            exponent += 1
        else
            if exponent > 0
                # print exponent for previous object
                exponent>1 && print(io,superscript(exponent))
                print(io," ")
            end
            print(io,x)
            prev = x
            exponent = 1
        end
    end
    exponent>1 && print(io,superscript(exponent))
    trailing_space && print(io," ")
end

Base.print(io::IO,A::OpProduct) = printspaced(io,_display_items(A))

# For display, products |i⟩⟨1| |1⟩⟨j| of operators on the same qudit (see `qudit`) are shown as |i⟩⟨j|
struct _QuditTransition
    create::BaseOperator
    destroy::BaseOperator
end
Base.:(==)(A::_QuditTransition,B::_QuditTransition) = A.create == B.create && A.destroy == B.destroy
_levels(T::_QuditTransition) = (qudit_levels(T.create)[2], qudit_levels(T.destroy)[3])
Base.print(io::IO,T::_QuditTransition) = (print(io,opname(T.create),_qudit_symbol(_levels(T)...,:text)); _print_inds(io,inds(T.create),true))

# the items to display for a product of operators: operators, and pairs of operators on the same qudit that
# are merged into a transition operator (if the second one commutes with all operators between them)
function _display_items(ops)
    any(O -> optype(O) == QuditCreate_, ops) || return ops
    v = collect(ops)
    items = Union{BaseOperator,_QuditTransition}[]
    merged = falses(length(v))
    for k in eachindex(v)
        merged[k] && continue
        O = v[k]
        kp = optype(O) == QuditCreate_ ? findnext(P -> optype(P) == QuditDestroy_ && opname(P) == opname(O) && inds(P) == inds(O), v, k+1) : nothing
        if kp !== nothing && !merged[kp] && all(m -> commutes(v[m],v[kp]), k+1:kp-1)
            merged[kp] = true
            push!(items, _QuditTransition(O,v[kp]))
        else
            push!(items, O)
        end
    end
    items
end
# the indices of the displayed items
_display_indices(ops) = QuIndex[ind for x in _display_items(ops) for ind in (x isa _QuditTransition ? inds(x.create) : inds(x))]
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

# print the indices of an operator or parameter as, e.g., (ij₂k), or as [i,j₂,k] with use_bracket_display (where
# operators without indices are printed without brackets). With parentheses, these are only printed
# for no indices if `always` is true.
function _print_inds(io::IO,inds,always::Bool)
    if using_bracket_display()
        isempty(inds) && return
        print(io,"[")
        join(io,inds,",")
        print(io,"]")
    elseif always || !isempty(inds)
        print(io,"(",inds...,")")
    end
end
Base.print(io::IO,A::BaseOperator) = (B = unalias(A); print(io,opname(B),opsymbol(B,:text)); _print_inds(io,inds(B),true))
Base.print(io::IO,A::δ) = print(io,"δ(",A.iA,A.iB,")")
Base.print(io::IO,A::Param) = (print(io, opname(A), parstate(A)=='c' ? "*" : ""); _print_inds(io,assignedinds(inds(A)),false))
Base.print(io::IO,A::ExpVal) = print(io,"⟨", ops(A), "⟩")
Base.print(io::IO,A::Corr) = print(io,"⟨", ops(A), "⟩c")

function Base.print(io::IO,A::QuTerm)
    if nsuminds(A) > 0
        print(io,"∑")
        for ii=1:nsuminds(A)
            print(io,subscript(ii))
        end
        print(io," ")
    end
    isempty(deltas(A)) || printspaced(io,deltas(A),true)
    isempty(params(A)) || printspaced(io,params(A),true)
    isempty(expvals(A)) || printspaced(io,expvals(A),true)
    isempty(corrs(A)) || printspaced(io,corrs(A),true)
    isempty(bares(A)) || print(io,bares(A))
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
        print_as_connector(io,numstring(s),print_connector)
        return
    elseif isequal(s, -one(s))
        str = "-"
    elseif isequal(s, 1)
        str = ""
    else
        str = numstring(s)*" "
    end
    print_as_connector(io,str,print_connector)
    print(io,t)
end

numstring(x::Number) = @sprintf "%g" x
numstring(x::Rational) = denominator(x)==1 ? string(numerator(x)) : string(x)
numstring(x::Complex) = ((r,i) = reim(x); iszero(i) ? numstring(r) : (iszero(r) ? numstring(i)*"i" : "($(numstring(r))$((i<0)===true ? '-' : '+')$(numstring(abs(i)))i)"))

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

function Base.print(io::IO, eqsys::QuEqSys)
    for (A,dAdt) in eqsys.eqs
        println(io,"dₜ", A," = ", dAdt)
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

# this is just to serve as an easily overloaded function to deal with "strange" Number subtypes from, e.g., Symbolics or SymPy
numlatex(x) = x

function _push_exponents!(ex,itr)
    prevO = nothing
    exponent = 0
    for O in itr
        if O==prevO
            exponent += 1
        else
            exponent>0 && push!(ex.args,exponent==1 ? prevO : LaTeXString("{$(latexify(prevO,env=:raw))}^{$exponent}"))
            exponent = 1
            prevO = O
        end
    end
    # last object
    exponent>0 && push!(ex.args,exponent==1 ? prevO : LaTeXString("{$(latexify(prevO,env=:raw))}^{$exponent}"))
end

@latexrecipe function f(ii::QuIndex)
    isnoindex(ii) && return LaTeXString("")
    isintindex(ii) && return LaTeXString("$(ii.num)")
    issumindex(ii) && return LaTeXString("\\#_{$(ii.num)}")
    s = unicode_to_latex(string(ii.sym))
    return LaTeXString(ii.num == typemin(ii.num) ? s : "$(s)_{$(ii.num)}")
end
latexindstr(inds) = isempty(assignedinds(inds)) ? "" : "_{$(latexify.(inds)...)}"
@latexrecipe function f(A::BaseOperator)
    B = unalias(A)
    return LaTeXString("{$(unicode_to_latex(opname(B)))}$(latexindstr(inds(B)))$(opsymbol(B,:latex))")
end
@latexrecipe function f(A::Param)
    return LaTeXString(string("$(unicode_to_latex(opname(A)))$(latexindstr(inds(A)))", parstate(A)=='c' ? "^{*}" : ""))
end
@latexrecipe function f(A::δ)
    return LaTeXString("\\delta$(latexindstr(make_indices(A.iA,A.iB)))")
end
@latexrecipe function f(T::_QuditTransition)
    return LaTeXString("{$(unicode_to_latex(opname(T.create)))}$(latexindstr(inds(T.create)))$(_qudit_symbol(_levels(T)...,:latex))")
end
@latexrecipe function f(A::OpProduct)
    mult_symbol --> ""

    isempty(A) && return
    ex = Expr(:call,:*)
    _push_exponents!(ex,_display_items(A))
    return ex
end
@latexrecipe function f(A::ExpVal)
    mult_symbol --> ""

    return Expr(:latexifymerge, "\\langle ", ops(A), "\\rangle")
end
@latexrecipe function f(A::Corr)
    mult_symbol --> ""

    return Expr(:latexifymerge, "\\langle ", ops(A), "\\rangle_{c}")
end
@latexrecipe function f(A::QuTerm)
    mult_symbol --> ""

    isempty(A) && return
    ex = Expr(:call,:*)
    nsuminds(A) > 0 && push!(ex.args,LaTeXString("\\sum_{$(latexify.(sumindex.(1:nsuminds(A)))...)}"))
    _push_exponents!(ex,Iterators.flatten((deltas(A),params(A),expvals(A),corrs(A),(bares(A),))))
    return ex
end
@latexrecipe function f(A::QuExpr)
    fmt --> FancyNumberFormatter()
    mult_symbol --> ""

    ex = 0
    # convert to tuples so sorting ignores the numbers if terms are different (which they are)
    for (t,s) in sort!(Tuple.(collect(A.terms)))
        # wrap in try/catch to handle "strange" Number subtypes from, e.g., Symbolics or SymPy
        sign, s = try
            isreal(s) && (s<0) ? (-1,-s) : (1,s)
        catch
            (1,s)
        end
        # use numlatex in case we have to treat s specially (e.g., add parenthesis for Symbolics.Num sums)
        x = isone(s) ? (isempty(t) ? numlatex(s) : t) : Expr(:call,:*,numlatex(s),t)
        if ex === 0
            ex = isone(sign) ? x : :(-$x)
        else
            ex = isone(sign) ? :($ex + $x) : :($ex - $x)
        end
    end
    return ex isa Expr ? ex : Expr(:latexifymerge,ex)
end
@latexrecipe function f(eqsys::QuEqSys)
    env --> :align

    # for align, we can return two vectors for the lhs and rhs of each equation
    dAs = Expr.(:latexifymerge, "\\frac{\\mathrm{d}}{\\mathrm{d}t}", keys(eqsys.eqs))
    return dAs, values(eqsys.eqs)
end

"""
    latex(A)

Render any `QuantumObject` to a raw LaTeX string using `Latexify.jl` recipes.
Useful for pretty-printing expressions or embedding results in documents.
"""
latex(A::QuantumObject) = latexify(A,env=:raw).s
