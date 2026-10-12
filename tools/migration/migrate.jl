# Migrate code from the legacy QuantumAlgebra interface (constructor functions, e.g., `@boson_ops a`
# and `a(:i)`) to the interface with operator variables (e.g., `@boson a[:]` and `a[:i]`).
#
# Usage: julia --project=tools/migration tools/migration/migrate.jl [--inplace | --dry-run] files...
#
# Julia files (.jl) and Jupyter notebooks (.ipynb) are supported. By default, the migrated code is
# written to a new file `name_migrated.jl` (or `.ipynb`). Code that cannot be migrated automatically
# (e.g., an operator that is used with different numbers of indices, or that is passed as a function)
# is left unchanged and reported. Since both interfaces can be used together in QuantumAlgebra 1.x,
# the migrated code still works in that case.

module QAMigration

import JuliaSyntax as JS
import JSON
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
nindices(n) = n == 1 ? "1 index" : "$n indices"
islegacydecl(n) = (m = macroname(n); m in keys(LEGACYMACRO) || m == Symbol("@anticommuting_fermion_group"))

function walk(f, n, parent=nothing)
    f(n, parent) === :skip && return
    for c in children(n)
        walk(f, c, n)
    end
end

## the operators of the legacy interface

mutable struct Species
    family::Symbol          # :boson, :fermion, :tlspm, :tlsxyz
    name::Symbol            # species name, e.g., :b for `@boson_ops b`
    decl::Any               # (unit, macrocall node) of the legacy declaration, or nothing for default operators
    group::Vector{Symbol}   # names of a group of anticommuting fermions (declared together)
    ndims::Set{Int}         # numbers of indices in the uses
    uses::Vector{Any}       # (unit, call node, callee node, role) for each use
    problems::Vector{String}
end
Species(family, name, decl=nothing; group=Symbol[]) = Species(family, name, decl, group, Set{Int}(), Any[], String[])

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
const NEWMACRO = Dict(:boson => "@boson", :fermion => "@fermion", :tlspm => "@tls", :tlsxyz => "@tlsxyz")
const LEGACYMACRO = Dict(Symbol("@boson_ops") => :boson, Symbol("@fermion_ops") => :fermion,
                         Symbol("@tlspm_ops") => :tlspm, Symbol("@tlsxyz_ops") => :tlsxyz)
const DEFAULT_SPECIES = ((:boson, :a), (:fermion, :f), (:tlsxyz, :σ), (:tlspm, :σ))

dimsstring(n) = n == 0 ? "" : "[" * join(fill(":", n), ",") * "]"

## parameters

mutable struct Parameter
    name::Symbol
    complex::Set{Bool}
    ndims::Set{Int}
    uses::Vector{Any}       # (unit, node, conjugated, index strings)
    problems::Vector{String}
end
Parameter(name) = Parameter(name, Set{Bool}(), Set{Int}(), Any[], String[])

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

Analyze the units (a file, or the cells of a notebook, which are treated as one program) and add the
edits that migrate them to the new interface. Returns a vector of messages about what was done and
what could not be migrated.
"""
function migrate_units!(units::Vector{Unit})
    msgs = String[]
    for u in units
        u.tree === nothing && push!(msgs, "$(u.label): could not be parsed, not migrated")
    end
    parsed = filter(u -> u.tree !== nothing, units)

    # legacy declarations
    species = Species[]
    for u in parsed
        walk(u.tree) do n, _
            m = macroname(n)
            m === nothing && return nothing
            args = children(n)[2:end]
            if haskey(LEGACYMACRO, m) && length(args) == 1 && isident(only(args))
                push!(species, Species(LEGACYMACRO[m], only(args).val, (u, n)))
            elseif m == Symbol("@anticommuting_fermion_group") && all(isident, args)
                names = Symbol[a.val for a in args]
                append!(species, Species(:fermion, N, (u, n); group=names) for N in names)
            end
            :skip
        end
    end

    # names that are bound in the code (assigned, defined as functions, loop variables, arguments, ...)
    bound = Set{Symbol}()
    bare = Dict{Symbol,Any}()  # identifiers that are used as values (not called), first location
    function addbound!(n)
        if isident(n)
            push!(bound, n.val)
        elseif kind(n) in (JS.K"tuple", JS.K"parameters", JS.K"::", JS.K"=", JS.K"...")
            foreach(addbound!, children(n))
        elseif kind(n) == JS.K"call" && !isempty(children(n))
            # function signature (or a parenthesized operator definition): the name and the arguments
            foreach(addbound!, children(n))
        end
    end
    markbare!(u, n) = haskey(bare, n.val) || (bare[n.val] = (u, n))
    for u in parsed
        walk(u.tree) do n, parent
            k = kind(n)
            (k == JS.K"quote" || islegacydecl(n)) && return :skip
            if isassignment(n, parent) || k in (JS.K"function", JS.K"local", JS.K"global", JS.K"->")
                addbound!(first(children(n)))
            elseif k in (JS.K"for", JS.K"let", JS.K"generator")
                # iteration specifications / let bindings: `x = …` or `x in …`
                walk(n) do c, _
                    kind(c) in (JS.K"=", JS.K"in") && addbound!(first(children(c)))
                    nothing
                end
            end
            if isadjoint(n) && isident(first(children(n)))
                # X' that is not called is used as a value
                called = parent !== nothing && isfuncall(parent) && first(children(parent)) === n
                called || markbare!(u, first(children(n)))
                return :skip
            end
            if isident(n) && parent !== nothing
                called = isfuncall(parent) && first(children(parent)) === n
                field = kind(parent) == JS.K"." && last(children(parent)) === n
                called || field || markbare!(u, n)
            end
            nothing
        end
    end

    # default operators: used (called) in the code, but not declared or bound
    declared = Set(s.name for s in species)
    allfuncs = Dict{Symbol,Tuple{Species,Symbol,Symbol}}()
    for (family, name) in DEFAULT_SPECIES
        sp = Species(family, name)
        any(f -> f in declared || f in bound, keys(legacy_functions(sp))) && continue
        push!(species, sp)
    end
    for sp in species, (f, (var, role)) in legacy_functions(sp)
        allfuncs[f] = (sp, var, role)
    end

    # uses of the legacy functions
    for u in parsed
        walk(u.tree) do n, _
            (kind(n) == JS.K"quote" || islegacydecl(n)) && return :skip
            isfuncall(n) || return nothing
            callee = first(children(n))
            adj = isadjoint(callee)
            fnode = adj ? first(children(callee)) : callee
            isident(fnode) && haskey(allfuncs, fnode.val) || return nothing
            sp, _, _ = allfuncs[fnode.val]
            args = children(n)[2:end]
            if any(a -> kind(a) in (JS.K"...", JS.K"parameters", JS.K"="), args)
                kind(last(args)) == JS.K"..." && length(args) == 1 || push!(sp.problems, "$(location(u, n)): `$(text(u, n))` cannot be converted automatically")
            else
                push!(sp.ndims, length(args))
            end
            push!(sp.uses, (u, n, fnode, adj))
            nothing
        end
    end
    for sp in species, f in keys(legacy_functions(sp))
        if haskey(bare, f)
            u, n = bare[f]
            push!(sp.problems, "$(location(u, n)): `$f` is used as a value (not called), e.g., passed as a function")
        end
    end

    # the variables of the new interface must not clash with other names in the code
    newvars(sp) = unique(var for (var, _) in values(legacy_functions(sp)))
    for sp in species
        sp.decl === nothing && isempty(sp.uses) && continue
        for v in newvars(sp)
            v ∉ keys(legacy_functions(sp)) && (v in bound || haskey(bare, v)) &&
                push!(sp.problems, "the new variable `$v` would clash with a name that is used in the code")
        end
        length(sp.ndims) > 1 && push!(sp.problems, "`$(sp.name)` is used with different numbers of indices ($(join(sort!(collect(sp.ndims)), ", "))); use separate names for these operators")
        isempty(sp.ndims) && !isempty(sp.uses) && push!(sp.problems, "the number of indices of `$(sp.name)` cannot be determined")
    end
    # a group of anticommuting fermions is migrated together
    for sp in species
        isempty(sp.group) && continue
        grp = filter(s -> s.group == sp.group, species)
        dims = union((s.ndims for s in grp)...)
        if any(s -> !isempty(s.problems), grp) || length(dims) > 1
            isempty(sp.problems) && push!(sp.problems, "the group $(join(sp.group, ", ")) cannot be migrated together")
        end
    end

    # parameters
    params = OrderedDict{Symbol,Parameter}()
    for u in parsed
        walk(u.tree) do n, parent
            m = macroname(n)
            if m in (Symbol("@Pr_str"), Symbol("@Pc_str"))
                str = last(children(n))
                s = kind(str) == JS.K"string" && length(children(str)) == 1 && kind(only(children(str))) == JS.K"String" ? only(children(str)).val : nothing
                s === nothing && return :skip
                name, inds = parse_paramstr(s)
                p = get!(() -> Parameter(name), params, name)
                if any(isnothing, inds)
                    push!(p.problems, "$(location(u, n)): the indices of `$(text(u, n))` cannot be converted automatically")
                end
                push!(p.complex, m == Symbol("@Pc_str"))
                push!(p.ndims, length(inds))
                push!(p.uses, (u, n, false, inds))
                return :skip
            end
            if isfuncall(n) && isident(first(children(n)), :param)
                args = children(n)[2:end]
                (!isempty(args) && kind(args[1]) == JS.K"quote" && length(children(args[1])) == 1 && isident(only(children(args[1])))) || return nothing
                name = only(children(args[1])).val
                p = get!(() -> Parameter(name), params, name)
                state = length(args) ≥ 2 ? charval(args[2]) : nothing
                inds = args[(state === nothing ? 2 : 3):end]
                if any(a -> !(kind(a) == JS.K"Integer" || (kind(a) == JS.K"quote" && length(children(a)) == 1 && isident(only(children(a))))), inds)
                    push!(p.problems, "$(location(u, n)): `$(text(u, n))` cannot be converted automatically")
                end
                st = state === nothing ? 'n' : state
                push!(p.complex, st != 'r')
                push!(p.ndims, length(inds))
                push!(p.uses, (u, n, st == 'c', [text(u, a) for a in inds]))
            end
            nothing
        end
    end
    opnames = Set(v for sp in species for v in newvars(sp))
    for p in values(params)
        length(p.complex) > 1 && push!(p.problems, "parameter `$(p.name)` is used both as a real and a complex parameter")
        length(p.ndims) > 1 && push!(p.problems, "parameter `$(p.name)` is used with different numbers of indices ($(join(sort!(collect(p.ndims)), ", ")))")
        (p.name in bound || haskey(bare, p.name) || p.name in opnames) && push!(p.problems, "the new variable `$(p.name)` would clash with a name that is used in the code")
    end

    # edits for operators
    newdecls = String[]
    tls_without_indices = false
    for sp in species
        if !isempty(sp.problems)
            append!(msgs, "not migrated: " .* sp.problems)
            continue
        end
        isempty(sp.uses) && sp.decl === nothing && continue
        nd = isempty(sp.ndims) ? 0 : only(sp.ndims)
        sp.family in (:tlspm, :tlsxyz) && nd == 0 && (tls_without_indices = true)
        if sp.decl === nothing
            push!(newdecls, "$(NEWMACRO[sp.family]) $(sp.name)$(dimsstring(nd))")
        elseif isempty(sp.group)
            u, n = sp.decl
            push!(u.edits, Edit(JS.byte_range(n), "$(NEWMACRO[sp.family]) $(sp.name)$(dimsstring(nd))"))
        elseif sp.name == first(sp.group)
            u, n = sp.decl
            push!(u.edits, Edit(JS.byte_range(n), "@anticommuting_fermions " * join(("$N$(dimsstring(nd))" for N in sp.group), " ")))
        end
        for (u, n, fnode, adj) in sp.uses
            var, role = legacy_functions(sp)[fnode.val]
            dag = role == :dag
            adj && role != :herm && (dag = !dag)
            newname = string(var) * (dag || (adj && role == :herm) ? "'" : "")
            args = children(n)[2:end]
            if isempty(args)
                push!(u.edits, Edit(JS.byte_range(n), newname))
            else
                # keep the original text of the arguments
                inner = text(u, first(JS.byte_range(first(args))):last(JS.byte_range(last(args))))
                push!(u.edits, Edit(JS.byte_range(n), "$newname[$inner]"))
            end
        end
        push!(msgs, "migrated `$(sp.name)` ($(sp.decl === nothing ? "default operator" : "declared"), $(nindices(nd)), $(length(sp.uses)) uses)")
    end

    # edits for parameters
    for p in values(params)
        if !isempty(p.problems)
            append!(msgs, "not migrated: " .* p.problems)
            continue
        end
        cplx, nd = only(p.complex), only(p.ndims)
        push!(newdecls, "@param $(p.name)$(dimsstring(nd))$(cplx ? "::Complex" : "")")
        for (u, n, conj, inds) in p.uses
            v = string(p.name) * (conj ? "'" : "")
            push!(u.edits, Edit(JS.byte_range(n), isempty(inds) ? v : "$v[$(join(inds, ","))]"))
        end
        push!(msgs, "migrated parameter `$(p.name)` ($(nindices(nd)), $(cplx ? "complex" : "real"), $(length(p.uses)) uses)")
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

## entry points

"""
    migrate_source(src::String; label="code") -> (newsrc, messages)

Migrate Julia code from the legacy QuantumAlgebra interface to the new interface.
"""
function migrate_source(src::String; label="code")
    u = Unit(src, label)
    msgs = migrate_units!([u])
    apply_edits(u), msgs
end

"""
    migrate_notebook(json::String; label="notebook") -> (newjson, messages)

Migrate the code cells of a Jupyter notebook (which are treated as one program). The outputs of the
cells are not changed (rerun the notebook to update them).
"""
function migrate_notebook(json::String; label="notebook")
    nb = JSON.parse(json; dicttype=OrderedDict)
    cells = [c for c in nb["cells"] if c["cell_type"] == "code"]
    src(c) = c["source"] isa AbstractVector ? join(c["source"]) : c["source"]
    units = [Unit(String(src(c)), "$label, cell $k") for (k, c) in enumerate(cells)]
    msgs = migrate_units!(units)
    for (c, u) in zip(cells, units)
        isempty(u.edits) && continue
        new = apply_edits(u)
        c["source"] = c["source"] isa AbstractVector ? [String(l) for l in eachline(IOBuffer(new); keep=true)] : new
    end
    io = IOBuffer()
    JSON.print(io, nb, 1)
    String(take!(io)) * "\n", msgs
end

"""
    migrate_file(path; inplace=false, dryrun=false)

Migrate a Julia file or Jupyter notebook. The result is written to `name_migrated.ext` (or to the file
itself with `inplace=true`, or not at all with `dryrun=true`). Returns the messages.
"""
function migrate_file(path::AbstractString; inplace=false, dryrun=false)
    src = read(path, String)
    new, msgs = endswith(path, ".ipynb") ? migrate_notebook(src; label=basename(path)) : migrate_source(src; label=basename(path))
    if !dryrun && new != src
        base, ext = splitext(path)
        out = inplace ? path : base * "_migrated" * ext
        write(out, new)
        push!(msgs, "wrote $out")
    end
    msgs
end

function main(args)
    inplace = "--inplace" in args
    dryrun = "--dry-run" in args
    files = filter(a -> !startswith(a, "--"), args)
    isempty(files) && (println("usage: julia --project=tools/migration tools/migration/migrate.jl [--inplace | --dry-run] files..."); return)
    for f in files
        println(f, ":")
        foreach(m -> println("  ", m), migrate_file(f; inplace, dryrun))
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    QAMigration.main(ARGS)
end
