# Migrate code from the legacy QuantumAlgebra interface (constructor functions, e.g., `@boson_ops a`
# and `a(:i)`) to the interface with operator variables (e.g., `@boson a[:]` and `a[:i]`).
#
# Usage: julia --project=tools/migration tools/migration/migrate.jl [--inplace | --dry-run] files...
#
# Julia files (.jl) and Jupyter notebooks (.ipynb) are supported. A file is migrated together with the
# files that it includes (with `include("file.jl")`), since they form one program. By default, the
# migrated code is written to new files `name_migrated.jl` (or `.ipynb`), and the includes are changed to
# point to them. Code that cannot be migrated automatically (e.g., an operator that is used with
# different numbers of indices, or that is passed as a function) is left unchanged and reported. Since
# both interfaces can be used together in QuantumAlgebra 1.x, the migrated code still works in that case.

module QAMigration

import JuliaSyntax as JS
import JSON
import Latexify
import REPL
using OrderedCollections: OrderedDict

export migrate_source, migrate_notebook, migrate_file

## parse units (a file or a notebook cell) and edits

struct Edit
    range::UnitRange{Int}  # byte range to replace (empty range = insertion before first(range))
    text::String
end

mutable struct Unit
    src::String
    label::String
    tree::Union{Nothing,JS.SyntaxNode}
    edits::Vector{Edit}
end
function Unit(src::String, label::String)
    tree = try
        JS.parseall(JS.SyntaxNode, src)
    catch
        nothing
    end
    Unit(src, label, tree, Edit[])
end

lineof(u::Unit, byte::Integer) = count(==(UInt8('\n')), view(codeunits(u.src), 1:min(byte, ncodeunits(u.src)))) + 1
location(u::Unit, n) = "$(u.label):$(lineof(u, first(JS.byte_range(n))))"
text(u::Unit, r::UnitRange) = String(view(codeunits(u.src), r))
text(u::Unit, n::JS.SyntaxNode) = text(u, JS.byte_range(n))

function apply_edits(u::Unit)
    s = codeunits(u.src)
    out = IOBuffer()
    pos = 1
    for e in sort(u.edits, by = e -> (first(e.range), last(e.range)))
        first(e.range) ≥ pos || error("overlapping edits in $(u.label)")
        write(out, view(s, pos:first(e.range)-1))
        write(out, e.text)
        pos = isempty(e.range) ? first(e.range) : last(e.range) + 1
    end
    write(out, view(s, pos:length(s)))
    String(take!(out))
end

## syntax helpers

kind(n) = JS.kind(n)
children(n) = JS.is_leaf(n) ? JS.SyntaxNode[] : JS.children(n)
isident(n, name=nothing) = kind(n) == JS.K"Identifier" && (name === nothing || n.val == name)
# the macro name of a macrocall (also for QuantumAlgebra.@name), or nothing
function macroname(n)
    kind(n) == JS.K"macrocall" || return nothing
    m = first(children(n))
    kind(m) == JS.K"." && (m = last(children(m)))
    kind(m) in (JS.K"MacroName", JS.K"StringMacroName") ? m.val : nothing
end
# a function call f(args...) (and not an operator call like a*b, -a or a')
isfuncall(n) = kind(n) == JS.K"call" && !JS.is_infix_op_call(n) && !JS.is_prefix_op_call(n) && !JS.is_postfix_op_call(n)
# `X'` (parsed as a postfix operator call)
isadjoint(n) = kind(n) == JS.K"call" && JS.is_postfix_op_call(n) && isident(last(children(n)), Symbol("'"))
# `x = …` is an assignment (and not a keyword argument or a named tuple element)
isassignment(n, parent) = kind(n) == JS.K"=" &&
    (parent === nothing || kind(parent) ∉ (JS.K"call", JS.K"parameters", JS.K"tuple", JS.K"vect", JS.K"braces", JS.K"ref", JS.K"curly", JS.K"macrocall", JS.K"dotcall"))
# the value of a character literal, or nothing
charval(n) = kind(n) == JS.K"Char" ? n.val :
    kind(n) == JS.K"char" && length(children(n)) == 1 && kind(only(children(n))) == JS.K"Char" ? only(children(n)).val : nothing
# the value of a plain string literal, or nothing
stringval(n) = kind(n) == JS.K"string" && length(children(n)) == 1 && kind(only(children(n))) == JS.K"String" ? only(children(n)).val : nothing
nindices(n) = n == 1 ? "1 index" : "$n indices"
islegacydecl(n) = (m = macroname(n); m in keys(LEGACYMACRO) || m == Symbol("@anticommuting_fermion_group"))

# walk the tree depth-first, calling f(n, stack) with the ancestors of n in stack (the parent is last);
# the children of n are skipped if f returns :skip
function walk(f, n, stack=JS.SyntaxNode[])
    f(n, stack) === :skip && return
    push!(stack, n)
    for c in children(n)
        walk(f, c, stack)
    end
    pop!(stack)
end
parentof(stack) = isempty(stack) ? nothing : last(stack)

## scopes

# constructs that introduce a local scope (short-form function definitions `f(x) = …` are also K"function")
const SCOPEKINDS = (JS.K"function", JS.K"->", JS.K"let", JS.K"for", JS.K"while", JS.K"try", JS.K"do",
                    JS.K"generator", JS.K"struct", JS.K"module", JS.K"macro")
# macros that wrap their arguments in a local scope
const SCOPEMACROS = (Symbol("@testset"), Symbol("@async"), Symbol("@spawn"))
isscope(n) = kind(n) in SCOPEKINDS || macroname(n) in SCOPEMACROS
# the innermost local scope that contains a node (nothing for global scope)
scopeof(stack) = (i = findlast(isscope, stack); i === nothing ? nothing : stack[i])
# whether a node is executed conditionally within its scope (e.g., inside `if`)
function isconditional(stack)
    i = something(findlast(isscope, stack), 0)
    any(n -> kind(n) in (JS.K"if", JS.K"elseif", JS.K"&&", JS.K"||", JS.K"?"), @view stack[i+1:end])
end
# the iteration specifications (`x in xs`, `x = xs`) of a for loop or a generator
function iterations(n)
    its = JS.SyntaxNode[]
    for c in children(n)
        kind(c) == JS.K"iteration" && append!(its, children(c))
        kind(c) == JS.K"filter" && foreach(i -> kind(i) == JS.K"iteration" && append!(its, children(i)), children(c))
    end
    filter(i -> kind(i) in (JS.K"in", JS.K"=") && !isempty(children(i)), its)
end

## the operators of the legacy interface

mutable struct Species
    family::Symbol          # :boson, :fermion, :tlspm, :tlsxyz
    name::Symbol            # species name, e.g., :b for `@boson_ops b`
    scope::Union{Nothing,JS.SyntaxNode}  # local scope of the declarations (nothing for global scope)
    decls::Vector{Any}      # (unit, macrocall node) of the legacy declarations (none for default operators)
    group::Vector{Symbol}   # names of a group of anticommuting fermions (declared together)
    ndims::Set{Int}         # numbers of indices in the uses
    uses::Vector{Any}       # (unit, call node, callee node, adjoint, index arguments) for each use
    problems::Vector{String}
end
Species(family, name, scope=nothing; group=Symbol[]) = Species(family, name, scope, Any[], group, Set{Int}(), Any[], String[])

# functions created by the legacy interface for a species: name => (variable name in the new interface, role),
# where role is :op (the operator), :dag (its adjoint) or :herm (a Hermitian operator, e.g., σx)
function legacy_functions(sp::Species)
    N = sp.name
    if sp.family in (:boson, :fermion)
        Dict(N => (N, :op), Symbol(N, :dag) => (N, :dag))
    elseif sp.family == :tlspm
        Dict(Symbol(N, :m) => (N, :op), Symbol(N, :p) => (N, :dag))
    else # :tlsxyz
        Dict(Symbol(N, c) => (Symbol(N, c), :herm) for c in (:x, :y, :z))
    end
end
newvars(sp::Species) = unique(var for (var, _) in values(legacy_functions(sp)))
function describe(sp::Species)
    N = sp.name
    names = sp.family == :tlspm ? "`$(N)p`/`$(N)m`" : sp.family == :tlsxyz ? "`$(N)x`/`$(N)y`/`$(N)z`" : "`$N`"
    where = isempty(sp.decls) ? "default operator" : "declared at $(location(first(sp.decls)...))"
    "$names ($where)"
end

const NEWMACRO = Dict(:boson => "@boson", :fermion => "@fermion", :tlspm => "@tls", :tlsxyz => "@tlsxyz")
const LEGACYMACRO = Dict(Symbol("@boson_ops") => :boson, Symbol("@fermion_ops") => :fermion,
                         Symbol("@tlspm_ops") => :tlspm, Symbol("@tlsxyz_ops") => :tlsxyz)
const DEFAULT_SPECIES = ((:boson, :a), (:fermion, :f), (:tlsxyz, :σ), (:tlspm, :σ))

# the names exported by QuantumAlgebra (except macros), Base and Core, which new variables must not shadow
const QA_EXPORTS = Set{Symbol}([:Avac, :Corr, :ExpVal, :QuExpr, :QuantumAlgebra, :a, :adag, :anticommuting_fermions,
    :boson, :boson_ops, :comm, :corr, :corr_as_expvals, :dropcorr, :droplen, :expval, :expval_as_corrs, :f, :fdag,
    :fermion, :fermion_ops, :gellmann, :heisenberg_eom, :heisenberg_eom_system, :julia_expression, :latex,
    :map_scalar_function, :normal_form, :param, :qindex, :qparam, :qudit, :so, :spin, :tls, :tlspm_ops, :tlsxyz,
    :tlsxyz_ops, :un, :vacA, :vacExpVal, :σm, :σp, :σx, :σy, :σz, :∑])
const RESERVED = union(QA_EXPORTS, names(Base), names(Core))

dimsstring(n) = n == 0 ? "" : "[" * join(fill(":", n), ",") * "]"

# the index arguments of a legacy call (a single tuple argument contains the indices, as in the legacy
# interface): (number of indices or nothing for a single splatted argument, nodes), or nothing if they
# cannot be converted
function indexargs(args)
    length(args) == 1 && kind(only(args)) == JS.K"tuple" && (args = children(only(args)))
    any(a -> kind(a) in (JS.K"parameters", JS.K"=", JS.K"..."), args) || return length(args), args
    length(args) == 1 && kind(only(args)) == JS.K"..." && return nothing, args
    nothing
end
# edits that replace the call `n` (with index arguments `args`) by `new[args]`, keeping the code of the arguments
# (so that edits inside them do not overlap)
function replace_call!(u::Unit, n, new, args)
    r = JS.byte_range(n)
    if isempty(args)
        push!(u.edits, Edit(r, new))
    else
        push!(u.edits, Edit(first(r):first(JS.byte_range(first(args)))-1, "$new["))
        push!(u.edits, Edit(last(JS.byte_range(last(args)))+1:last(r), "]"))
    end
end

## parameters

mutable struct Parameter
    name::Symbol
    var::Symbol             # the variable in the new interface (the name, or a unicode name for a name in LaTeX)
    complex::Set{Bool}
    ndims::Set{Int}
    uses::Vector{Any}       # (unit index, unit, node, conjugated, index arguments or index string)
    assigns::Vector{Any}    # (unit index, unit, assignment node) of statements `name = param(:name, …)` in global scope
    problems::Vector{String}
end
Parameter(name) = Parameter(name, something(unicode_name(string(name)), name), Set{Bool}(), Set{Int}(), Any[], Any[], String[])

# for a name in LaTeX (e.g., `\bar{b}`, `\tilde{\Omega}` or `\omega`), the unicode name with the same LaTeX output
# (`b̄`, `Ω̃`, `ω`) if that is a valid variable name, or nothing
function unicode_name(s::AbstractString)
    Base.isidentifier(s) && return nothing
    parts = String[]
    pos = 1
    while pos ≤ ncodeunits(s)
        # an accent on a character or a LaTeX symbol (e.g., `\bar{b}` or `\bar{\omega}`), a LaTeX symbol, or a character
        m = match(r"\G(?:\\([A-Za-z]+)\{(\\[A-Za-z]+|[^{}\\])\}|(\\[A-Za-z]+)|([^{}\\]))", s, pos)
        m === nothing && return nothing
        if m[1] !== nothing
            accent = get(REPL.REPLCompletions.latex_symbols, "\\" * m[1], nothing)
            base = startswith(m[2], "\\") ? get(REPL.REPLCompletions.latex_symbols, m[2], nothing) : m[2]
            (accent === nothing || base === nothing) && return nothing
            push!(parts, base * accent)
        elseif m[3] !== nothing
            c = get(REPL.REPLCompletions.latex_symbols, m[3], nothing)
            c === nothing && return nothing
            push!(parts, c)
        else
            push!(parts, m[4])
        end
        pos += ncodeunits(m.match)
    end
    u = Base.Unicode.normalize(join(parts), :NFC)
    # QuantumAlgebra converts names to LaTeX with Latexify
    Base.isidentifier(u) && Latexify.unicode2latex(u) == s ? Symbol(u) : nothing
end
# a string literal with the value s (raw"…" for strings with backslashes)
stringliteral(s) = occursin('\\', s) && !occursin('"', s) && !endswith(s, "\\") ? "raw\"$s\"" : repr(s)

# the name and indices of a parameter string `name_inds` (as in Pr"…"/Pc"…"), with the indices as code
function parse_paramstr(s::AbstractString)
    parts = split(s, "_"; limit=2)
    name = Symbol(parts[1])
    length(parts) == 1 && return name, String[]
    indstr = parts[2]
    startswith(indstr, "{") && endswith(indstr, "}") && (indstr = chop(indstr, head=1, tail=1))
    inds = map(split(indstr, ",")) do x
        x = strip(x)
        occursin(r"^\d+$", x) ? String(x) : (Base.isidentifier(x) ? ":" * x : nothing)
    end
    name, inds
end

## analysis

"""
    migrate_units!(units) -> messages

Analyze the units (the code of a program: a file, or the cells of a notebook, together with the files
they include) and add the edits that migrate them to the new interface. Returns a vector of messages
about what was done and what could not be migrated.
"""
function migrate_units!(units::Vector{Unit})
    msgs = String[]
    for u in units
        u.tree === nothing && push!(msgs, "$(u.label): could not be parsed, not migrated")
    end
    parsed = filter(u -> u.tree !== nothing, units)

    # names that are bound in the code (with the number of bindings), and the positions of all identifiers
    bindings = Dict{Symbol,Int}()
    idents = Dict{Symbol,Vector{Tuple{Int,Int}}}()
    function addbound!(n)
        if isident(n)
            bindings[n.val] = get(bindings, n.val, 0) + 1
        elseif kind(n) in (JS.K"=", JS.K"::", JS.K"where")
            # `x = …` (e.g., in `global x = …` or a keyword argument), `x::T`, `f(x) where T`
            isempty(children(n)) || addbound!(first(children(n)))
        elseif kind(n) in (JS.K"tuple", JS.K"parameters", JS.K"...", JS.K"call")
            # (a function signature for K"call")
            foreach(addbound!, children(n))
        end
    end
    for (k, u) in enumerate(parsed)
        walk(u.tree) do n, stack
            nk, parent = kind(n), parentof(stack)
            (nk == JS.K"quote" || islegacydecl(n)) && return :skip
            if isassignment(n, parent) || nk in (JS.K"function", JS.K"local", JS.K"global", JS.K"->", JS.K"do", JS.K"macro")
                addbound!(first(children(n)))
            elseif nk in (JS.K"for", JS.K"generator")
                foreach(it -> addbound!(first(children(it))), iterations(n))
            elseif nk == JS.K"let"
                foreach(c -> kind(c) == JS.K"=" || addbound!(c), children(first(children(n))))
            elseif nk == JS.K"catch" && isident(first(children(n)))
                addbound!(first(children(n)))
            end
            if isident(n) && !(parent !== nothing && kind(parent) == JS.K"." && last(children(parent)) === n)
                push!(get!(Vector{Tuple{Int,Int}}, idents, n.val), (k, first(JS.byte_range(n))))
            end
            nothing
        end
    end

    # legacy declarations, in their scopes
    species = Species[]
    function addspecies!(family, name, scope, decl; group=Symbol[])
        i = findfirst(s -> (s.family, s.name, s.group) == (family, name, group) && s.scope === scope, species)
        sp = i === nothing ? push!(species, Species(family, name, scope; group))[end] : species[i]
        push!(sp.decls, decl)
        sp
    end
    for u in parsed
        walk(u.tree) do n, stack
            kind(n) == JS.K"quote" && return :skip
            m = macroname(n)
            m === nothing && return nothing
            args = children(n)[2:end]
            if haskey(LEGACYMACRO, m) && length(args) == 1 && isident(only(args))
                new = [addspecies!(LEGACYMACRO[m], only(args).val, scopeof(stack), (u, n))]
            elseif m == Symbol("@anticommuting_fermion_group") && !isempty(args) && all(isident, args)
                names = Symbol[a.val for a in args]
                new = [addspecies!(:fermion, N, scopeof(stack), (u, n); group=names) for N in names]
            else
                return nothing  # look inside other macro calls (e.g., @testset)
            end
            isconditional(stack) && foreach(sp -> push!(sp.problems, "declared conditionally (e.g., inside `if`) at $(location(u, n))"), new)
            :skip
        end
    end
    declared = Dict{Symbol,Vector{Species}}()  # legacy function => species that declare it
    for sp in species, f in keys(legacy_functions(sp))
        push!(get!(Vector{Species}, declared, f), sp)
    end
    for (f, sps) in declared, sp in sps
        any(o -> o !== sp && o.scope === sp.scope, sps) && push!(sp.problems, "`$f` is declared as different operators in the same scope")
    end

    # default operators: used (called) in the code, but not bound by it
    defaults = Dict{Symbol,Species}()
    for (family, name) in DEFAULT_SPECIES
        sp = Species(family, name)
        any(f -> haskey(bindings, f), keys(legacy_functions(sp))) && continue
        push!(species, sp)
        foreach(f -> defaults[f] = sp, keys(legacy_functions(sp)))
    end

    # the species that a legacy function refers to at a position in the code: the declaration in the innermost
    # enclosing scope, or the default operator
    function resolve(f::Symbol, stack)
        sps = get(declared, f, Species[])
        for s in Iterators.flatten((Iterators.reverse(filter(isscope, stack)), (nothing,)))
            i = findfirst(sp -> sp.scope === s, sps)
            i === nothing || return sps[i]
        end
        get(defaults, f, nothing)
    end
    function bareuse!(u, n, stack)
        sp = resolve(n.val, stack)
        (sp === nothing || any(p -> startswith(p, "used as a value") || startswith(p, "imported"), sp.problems)) && return
        imported = any(s -> kind(s) in (JS.K"using", JS.K"import"), stack)
        push!(sp.problems, imported ? "imported explicitly at $(location(u, n))" :
                                      "used as a value (not called) at $(location(u, n)), e.g., passed as a function")
    end

    # uses of the legacy functions
    for u in parsed
        walk(u.tree) do n, stack
            (kind(n) == JS.K"quote" || islegacydecl(n)) && return :skip
            parent = parentof(stack)
            if isfuncall(n)
                callee = first(children(n))
                adj = isadjoint(callee)
                fnode = adj ? first(children(callee)) : callee
                isident(fnode) || return nothing
                sp = resolve(fnode.val, stack)
                sp === nothing && return nothing
                inds = indexargs(children(n)[2:end])
                if inds === nothing
                    push!(sp.problems, "the call `$(text(u, n))` at $(location(u, n)) cannot be converted automatically")
                else
                    first(inds) === nothing || push!(sp.ndims, first(inds))
                    push!(sp.uses, (u, n, fnode, adj, last(inds)))
                end
            elseif isadjoint(n) && isident(first(children(n)))
                # X' that is not called is used as a value
                called = parent !== nothing && isfuncall(parent) && first(children(parent)) === n
                called || bareuse!(u, first(children(n)), stack)
                return :skip
            elseif isident(n) && parent !== nothing
                called = isfuncall(parent) && first(children(parent)) === n
                field = kind(parent) == JS.K"." && last(children(parent)) === n
                called || field || bareuse!(u, n, stack)
            end
            nothing
        end
    end

    # the variables of the new interface must not clash with other names in the code: anywhere for variables in
    # global scope, and in the scope of the declaration for local variables
    unitindex = IdDict{Unit,Int}(u => k for (k, u) in enumerate(parsed))
    function clashes(v, sp)
        sp.scope === nothing && return haskey(idents, v) || v in RESERVED
        k, r = unitindex[first(first(sp.decls))], JS.byte_range(sp.scope)
        any(((kk, pos),) -> kk == k && pos in r, get(idents, v, ()))
    end
    for sp in species
        isempty(sp.decls) && isempty(sp.uses) && continue
        for v in newvars(sp)
            v ∉ keys(legacy_functions(sp)) && clashes(v, sp) &&
                push!(sp.problems, "the new variable `$v` would clash with a name that is used in the code")
        end
        length(sp.ndims) > 1 && push!(sp.problems, "used with different numbers of indices ($(join(sort!(collect(sp.ndims)), ", "))); use separate names for these operators")
        isempty(sp.ndims) && !isempty(sp.uses) && push!(sp.problems, "used with an unknown number of indices")
    end
    # a group of anticommuting fermions is migrated together
    for sp in species
        isempty(sp.group) && continue
        grp = filter(s -> s.group == sp.group && s.scope === sp.scope, species)
        dims = union((s.ndims for s in grp)...)
        if any(s -> !isempty(s.problems), grp) || length(dims) > 1
            isempty(sp.problems) && push!(sp.problems, "the group $(join(sp.group, ", ")) cannot be migrated together")
        end
    end

    # parameters
    params = OrderedDict{Symbol,Parameter}()
    for (k, u) in enumerate(parsed)
        walk(u.tree) do n, stack
            kind(n) == JS.K"quote" && return :skip
            m = macroname(n)
            if m in (Symbol("@Pr_str"), Symbol("@Pc_str"))
                s = stringval(last(children(n)))
                s === nothing && return :skip
                name, inds = parse_paramstr(s)
                p = get!(() -> Parameter(name), params, name)
                if any(isnothing, inds)
                    push!(p.problems, "the indices of `$(text(u, n))` at $(location(u, n)) cannot be converted automatically")
                else
                    addparamuse!(p, k, u, n, stack, m == Symbol("@Pc_str"), false, length(inds), join(inds, ","))
                end
                return :skip
            end
            if isfuncall(n) && isident(first(children(n)), :param)
                args = children(n)[2:end]
                (!isempty(args) && kind(args[1]) == JS.K"quote" && length(children(args[1])) == 1 && isident(only(children(args[1])))) || return nothing
                name = only(children(args[1])).val
                p = get!(() -> Parameter(name), params, name)
                rest = args[2:end]
                # the second argument is the state ('r', 'n' or 'c') if it is a character
                state = isempty(rest) ? 'n' : charval(rest[1])
                if state === nothing && kind(rest[1]) ∉ (JS.K"Integer", JS.K"quote", JS.K"tuple")
                    push!(p.problems, "`$(text(u, n))` at $(location(u, n)) cannot be converted automatically")
                    return nothing
                end
                state === nothing ? (state = 'n') : (rest = rest[2:end])
                inds = indexargs(rest)
                if inds === nothing || first(inds) === nothing
                    push!(p.problems, "`$(text(u, n))` at $(location(u, n)) cannot be converted automatically")
                else
                    addparamuse!(p, k, u, n, stack, state != 'r', state == 'c', first(inds), last(inds))
                end
            end
            nothing
        end
    end
    opnames = Set(v for sp in species if !(isempty(sp.decls) && isempty(sp.uses)) for v in newvars(sp))
    for p in values(params)
        length(p.complex) > 1 && push!(p.problems, "used both as a real and a complex parameter")
        length(p.ndims) > 1 && push!(p.problems, "used with different numbers of indices ($(join(sort!(collect(p.ndims)), ", ")))")
        clash = p.var in opnames || p.var in RESERVED || any(q -> q !== p && q.var == p.var, values(params))
        if isempty(p.assigns)
            clash |= haskey(idents, p.var) || haskey(bindings, p.var)
        else
            # `name = param(:name)` becomes the declaration `@qparam name`, so `name` must not be used before or
            # bound otherwise
            first_assign = minimum(((k, first(JS.byte_range(a))) for (k, _, a) in p.assigns))
            clash |= get(bindings, p.var, 0) > length(p.assigns) ||
                     any(<(first_assign), get(idents, p.var, ())) ||
                     any(((k, _, n),) -> (k, first(JS.byte_range(n))) < first_assign, p.uses)
        end
        if !Base.isidentifier(p.var)
            # suggest a declaration with a variable name (which has to be chosen by the user)
            var = filter(c -> Base.is_id_char(c), string(p.name))
            Base.isidentifier(var) || (var = "x")
            decl = "@qparam $var$(length(p.ndims) == 1 ? dimsstring(only(p.ndims)) : "")$(p.complex == Set([true]) ? "::Complex" : "")=$(stringliteral(string(p.name)))"
            push!(p.problems, "`$(p.name)` is not a valid variable name, declare it with a variable name of your choice, e.g., `$decl`")
        elseif clash
            push!(p.problems, "the new variable `$(p.var)` would clash with a name that is used in the code")
        end
    end

    # edits for operators
    newdecls = String[]
    tls_without_indices = false
    for sp in species
        if !isempty(sp.problems)
            push!(msgs, "not migrated: $(describe(sp)): " * join(unique(sp.problems), "; "))
            continue
        end
        isempty(sp.uses) && isempty(sp.decls) && continue
        nd = isempty(sp.ndims) ? 0 : only(sp.ndims)
        sp.family in (:tlspm, :tlsxyz) && nd == 0 && (tls_without_indices = true)
        if isempty(sp.decls)
            push!(newdecls, "$(NEWMACRO[sp.family]) $(sp.name)$(dimsstring(nd))")
        elseif isempty(sp.group)
            foreach(((u, n),) -> push!(u.edits, Edit(JS.byte_range(n), "$(NEWMACRO[sp.family]) $(sp.name)$(dimsstring(nd))")), sp.decls)
        elseif sp.name == first(sp.group)
            decl = "@anticommuting_fermions " * join(("$N$(dimsstring(nd))" for N in sp.group), " ")
            foreach(((u, n),) -> push!(u.edits, Edit(JS.byte_range(n), decl)), sp.decls)
        end
        for (u, n, fnode, adj, args) in sp.uses
            var, role = legacy_functions(sp)[fnode.val]
            dag = role == :dag
            adj && role != :herm && (dag = !dag)
            replace_call!(u, n, string(var) * (dag || (adj && role == :herm) ? "'" : ""), args)
        end
        push!(msgs, "migrated $(describe(sp)): $(nindices(nd)), $(length(sp.uses)) uses")
    end

    # edits for parameters
    for p in values(params)
        if !isempty(p.problems)
            push!(msgs, "not migrated: parameter `$(p.name)`: " * join(unique(p.problems), "; "))
            continue
        end
        cplx, nd = only(p.complex), only(p.ndims)
        decl = "@qparam $(p.var)$(dimsstring(nd))$(cplx ? "::Complex" : "")"
        isempty(p.assigns) && push!(newdecls, decl)
        foreach(((_, u, a),) -> push!(u.edits, Edit(JS.byte_range(a), decl)), p.assigns)
        for (_, u, n, conj, args) in p.uses
            any(((_, _, a),) -> last(children(a)) === n, p.assigns) && continue
            v = string(p.var) * (conj ? "'" : "")
            if args isa String
                push!(u.edits, Edit(JS.byte_range(n), isempty(args) ? v : "$v[$args]"))
            else
                replace_call!(u, n, v, args)
            end
        end
        renamed = p.var == p.name ? "" : " as `$(p.var)` (which gives the same LaTeX output)"
        push!(msgs, "migrated parameter `$(p.name)`$renamed: $(nindices(nd)), $(cplx ? "complex" : "real"), $(length(p.uses)) uses")
    end

    # declarations of default operators and parameters, inserted after `using QuantumAlgebra`
    if !isempty(newdecls)
        u, pos = insertion_point(parsed)
        decl = "# operators and parameters (migrated from the legacy QuantumAlgebra interface)\n" * join(newdecls, "\n")
        push!(u.edits, Edit(pos:pos-1, pos == 1 ? decl * "\n" : "\n" * decl))
    end
    # operators without indices are created when they are declared (instead of when they are used), so
    # the choice of basis for two-level systems (use_σpm) has to be made before that
    if tls_without_indices && any(u -> occursin(r"use_σ(pm|xyz)", u.src), parsed)
        push!(msgs, "warning: two-level system operators without indices are now created where they are declared, " *
                    "check that use_σpm/use_σxyz is called before that")
    end
    msgs
end

# record a use of a parameter (and whether it is the right-hand side of `name = …` in global scope)
function addparamuse!(p::Parameter, k, u, n, stack, cplx, conj, nd, args)
    push!(p.complex, cplx)
    push!(p.ndims, nd)
    push!(p.uses, (k, u, n, conj, args))
    length(stack) ≥ 2 || return
    a, ctx = stack[end], stack[end-1]
    if nd == 0 && !conj && isassignment(a, ctx) && last(children(a)) === n && isident(first(children(a)), p.var) &&
       kind(ctx) in (JS.K"toplevel", JS.K"block") && scopeof(stack) === nothing && !isconditional(stack)
        push!(p.assigns, (k, u, a))
    end
end

# where to insert new declarations: after the (first) statement `using QuantumAlgebra`, or at the beginning
function insertion_point(units)
    for u in units
        for st in children(u.tree)
            if kind(st) in (JS.K"using", JS.K"import") && occursin("QuantumAlgebra", text(u, st))
                return u, last(JS.byte_range(st)) + 1
            end
        end
    end
    first(units), 1
end

## files and programs

# a Julia file, or a notebook (whose code cells are the units)
mutable struct SourceFile
    path::String
    units::Vector{Unit}
    notebook::Any           # the parsed notebook, or nothing for Julia files
    cells::Vector{Any}      # the code cells of a notebook
    includes::Vector{Any}   # (unit, string node, included file) for each `include("file.jl")`
end
function SourceFile(src::String, path::String, label::String)
    if endswith(path, ".ipynb")
        nb = JSON.parse(src; dicttype=OrderedDict)
        cells = [c for c in nb["cells"] if c["cell_type"] == "code"]
        csrc(c) = c["source"] isa AbstractVector ? join(c["source"]) : c["source"]
        SourceFile(path, [Unit(String(csrc(c)), "$label, cell $k") for (k, c) in enumerate(cells)], nb, cells, Any[])
    else
        SourceFile(path, [Unit(src, label)], nothing, Any[], Any[])
    end
end
changed(f::SourceFile) = any(u -> !isempty(u.edits), f.units)
# the migrated source
function render(f::SourceFile)
    f.notebook === nothing && return apply_edits(only(f.units))
    for (c, u) in zip(f.cells, f.units)
        isempty(u.edits) && continue
        new = apply_edits(u)
        c["source"] = c["source"] isa AbstractVector ? [String(l) for l in eachline(IOBuffer(new); keep=true)] : new
    end
    io = IOBuffer()
    JSON.print(io, f.notebook, 1)
    String(take!(io)) * "\n"
end

# a file and the files that it includes (recursively), which form one program
function load_program(path::AbstractString)
    root = dirname(abspath(path))
    files = OrderedDict{String,SourceFile}()
    function load(p)
        haskey(files, p) && return files[p]
        f = files[p] = SourceFile(read(p, String), p, relpath(p, root))
        for u in f.units
            u.tree === nothing && continue
            walk(u.tree) do n, _
                kind(n) == JS.K"quote" && return :skip
                if isfuncall(n) && isident(first(children(n)), :include) && length(children(n)) == 2
                    s = children(n)[2]
                    inc = stringval(s)
                    if inc !== nothing && endswith(inc, ".jl") && isfile(joinpath(dirname(p), inc))
                        push!(f.includes, (u, only(children(s)), load(normpath(joinpath(dirname(p), inc)))))
                    end
                end
                nothing
            end
        end
        f
    end
    load(normpath(abspath(path)))
    collect(values(files))
end

## entry points

"""
    migrate_source(src::String; label="code") -> (newsrc, messages)

Migrate Julia code from the legacy QuantumAlgebra interface to the new interface.
"""
function migrate_source(src::String; label="code")
    f = SourceFile(src, "code.jl", label)
    msgs = migrate_units!(f.units)
    render(f), msgs
end

"""
    migrate_notebook(json::String; label="notebook") -> (newjson, messages)

Migrate the code cells of a Jupyter notebook (which are treated as one program). The outputs of the
cells are not changed (rerun the notebook to update them).
"""
function migrate_notebook(json::String; label="notebook")
    f = SourceFile(json, "notebook.ipynb", label)
    msgs = migrate_units!(f.units)
    render(f), msgs
end

"""
    migrate_program(path) -> (files, messages)

Migrate a Julia file or Jupyter notebook together with the files that it includes (with
`include("file.jl")`, recursively), which are treated as one program. Returns the files (with
`changed(f)` and the migrated source `render(f)`) and the messages, without writing anything.
"""
function migrate_program(path::AbstractString)
    files = load_program(path)
    files, migrate_units!(reduce(vcat, (f.units for f in files)))
end

migrated_path(p) = ((base, ext) = splitext(p); base * "_migrated" * ext)

"""
    migrate_file(path; inplace=false, dryrun=false)

Migrate a Julia file or Jupyter notebook together with the files that it includes. The results are
written to `name_migrated.ext` (with the includes pointing to the migrated files), or to the files
themselves with `inplace=true`, or not at all with `dryrun=true`. Returns the messages.
"""
function migrate_file(path::AbstractString; inplace=false, dryrun=false)
    files, msgs = migrate_program(path)
    length(files) > 1 && pushfirst!(msgs, "with the included files " * join((only(f.units).label for f in files[2:end] if f.notebook === nothing), ", "))
    dryrun && return msgs
    out = Set(filter(changed, files))
    if !inplace
        # files that include migrated files are also written, with the includes pointing to the migrated files
        while (new = filter(f -> f ∉ out && any(i -> i[3] in out, f.includes), files); !isempty(new))
            union!(out, new)
        end
        for f in out, (u, s, inc) in f.includes
            inc in out && push!(u.edits, Edit(JS.byte_range(s), migrated_path(s.val)))
        end
    end
    for f in files
        f in out || continue
        p = inplace ? f.path : migrated_path(f.path)
        write(p, render(f))
        push!(msgs, "wrote $p")
    end
    msgs
end

function main(args)
    inplace = "--inplace" in args
    dryrun = "--dry-run" in args
    paths = filter(a -> !startswith(a, "--"), args)
    isempty(paths) && (println("usage: julia --project=tools/migration tools/migration/migrate.jl [--inplace | --dry-run] files..."); return)
    # files that are included by another given file are migrated together with it
    included = Set{String}()
    for p in paths
        try
            union!(included, (f.path for f in load_program(p)[2:end]))
        catch
        end
    end
    failed = false
    for p in paths
        if normpath(abspath(p)) in included
            println(p, ": migrated together with the file that includes it")
            continue
        end
        println(p, ":")
        try
            foreach(m -> println("  ", m), migrate_file(p; inplace, dryrun))
        catch e
            println("  error: ", sprint(showerror, e))
            failed = true
        end
    end
    failed && exit(1)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    QAMigration.main(ARGS)
end
