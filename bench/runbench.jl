module Benchmarks

using Chairmarks
using Printf
using Random
using QuantumAlgebra

@boson_ops a
@boson_ops b
@tlsxyz_ops σ

function benchmark_runner(f, args; seconds::Real)
    sample = @b $f($args...) seconds=seconds
    return Dict("time_ns" => sample.time * 1e9,
                "allocs" => sample.allocs,
                "bytes" => sample.bytes)
end

function nested_comm(Hamiltonian, A, n::Integer = 5)
    for _ = 1:n
        A = normal_form(comm(Hamiltonian, A))
    end
    return A
end

function corrheis_nocache(op, H, Ls)
    empty!(QuantumAlgebra._CORRHEISDICT)
    return QuantumAlgebra.corrheis(op, H, Ls)
end

function make_cases()
    suites = []

    H = ∑(:n, ∑(:m, ∑(:k, Pr"ω_n,m" * a(:n, :k)' * a(:m, :k)))) +
        ∑(:i, 1 // 2 * Pr"ν_i" * σz(:i)) +
        ∑(:n, ∑(:k, ∑(:i, Pr"g_i,n,k" * σx(:i) * (a(:n, :k)' + a(:n, :k)))))

    Acomm = σz(:j)

    push!(suites, (name = "commutator", cases = [
        (name = "commutator_n2", f = nested_comm, args = (H, Acomm, 2)),
        (name = "commutator_n4", f = nested_comm, args = (H, Acomm, 4)),
    ]))
        
    Random.seed!(4)
    # Shuffle so no operator appears twice (as σz() * σz() == 1).
    randops = shuffle!([a.(1:50)..., σz.(1:50)..., adjoint.(a.(1:50))...])

    push!(suites, (name = "prodcorr", cases = [
        begin
            ops = Tuple(sort(randops[1:n], by = op -> first(op.terms)[1]))
            expr = normal_form(8 * Pr"g_i" * prod(ops))
            (name = "prodcorr_n$(n)", f = expval_as_corrs, args = (expr,))
        end for n in 2:9
    ]))

    Avec1 = normal_form((a()' + a()) * (σx() + σz()))
    Svec1 = normal_form(1 + Pr"η" * a()' + Pr"β" * σx())
    Avec2 = normal_form((a()' + a() + b()' + b())^2 + σx() * (a()' * b() + b()' * a()))
    Svec2 = normal_form(1 + Pr"α" * a()' + Pr"β" * b()' + Pr"χ" * σx())

    vacxpval_cases = Any[(name = "vacexpval_single_mode", f = vacExpVal, args = (Avec1, Svec1)),
                         (name = "vacexpval_multi_mode", f = vacExpVal, args = (Avec2, Svec2))]
    if hasmethod(vacExpVal, (QuExpr, QuExpr, Tuple))
        push!(vacxpval_cases, (name = "vacexpval_restricted_modes", f = vacExpVal, args = (Avec2, Svec2, (a(),))))
    end
    push!(suites, (name = "vacexpval", cases = vacxpval_cases))

    xsum1 = a(:d) * a(:c)'
    xsum2 = a(:i) * a(:j)' * σz(:α) * σy(:β) * σx(:α)
    Anf_sum = ∑(:c, ∑(:i, ∑(:α, 3 * expval(xsum1) * xsum2)))
    Anf_tls = prod(fill(σx(:i) * σy(:i) * σz(:i), 4))
    Anf_poly = (a(:i) * a(:j)' + b(:j) * b(:i)' + σx(:k) * σy(:k) + expval(a(:i)' * b(:j)))^3

    push!(suites, (name = "normal_form_stress", cases = [
        (name = "normal_form_sum_nested", f = normal_form, args = (Anf_sum,)),
        (name = "normal_form_tls_chain", f = normal_form, args = (Anf_tls,)),
        (name = "normal_form_polynomial_power", f = normal_form, args = (Anf_poly,)),
    ]))

    if isdefined(QuantumAlgebra, :corrheis)
        Hcorr = Pr"ωa" * a()' * a() + Pr"ωb" * b()' * b() + Pr"g" * (a()' + a()) * σx()
        Lcorr = ((Pr"κ", a()), (Pr"γ", σx()))
        Ocorr = a()' * σz() * b()

        push!(suites, (name = "corrheis", cases = [
            (name = "corrheis_1", f = corrheis_nocache, args = (Ocorr, Hcorr, Lcorr)),
            (name = "corrheis_2", f = corrheis_nocache, args = (Ocorr, Hcorr * Hcorr, Lcorr)),
        ]))
    end


    if isdefined(QuantumAlgebra, :heisenberg_eom_system)
        Heq = Pr"ωa" * a()' * a() +
            Pr"ωb" * b()' * b() +
            Pr"gab" * (a()' * b() + b()' * a()) +
            1 // 2 * Pr"ωe" * σz()
        Leq = ((Pr"κa", a()), (Pr"κb", b()), (Pr"γe", σx()))

        push!(suites, (name = "eqsys", cases = [
            (name = "eqsys_expval_ord2_small",  f = heisenberg_eom_system, args = (ExpVal, Heq, 2, Leq, (a(),))),
            (name = "eqsys_expval_ord3_medium", f = heisenberg_eom_system, args = (ExpVal, Heq, 3, Leq, (a(), b(), σz()))),
            (name = "eqsys_corr_ord2_small",    f = heisenberg_eom_system, args = (Corr, Heq, 2, Leq, (a(),))),
            (name = "eqsys_corr_ord3_medium",   f = heisenberg_eom_system, args = (Corr, Heq, 3, Leq, (a() * b(), σz()))),
            ]))

        Hmixed = Pr"ωa" * a()' * a() + Pr"ωb" * b()' * b() + Pr"gm" * (a()' * b() + b()' * a())
        Lmixed = ((Pr"κa", a()), (Pr"κb", b()), (Pr"η", (a(), b())))

        push!(suites, (name = "mixed_lindblad", cases = [
            (name = "mixed_lindblad_eqsys_ord2", f = heisenberg_eom_system, args = (ExpVal, Hmixed, 2, Lmixed, (a(), b()))),
            (name = "mixed_lindblad_eqsys_ord3", f = heisenberg_eom_system, args = (Hmixed, 3, Lmixed, (a(), b(), a()' * b()))),
        ]))
    end

    return suites
end

function run_benchmarks(; seconds::Real = 1.0, io::IO = stdout)
    suites = make_cases()

    run_case = case -> begin
        metrics = benchmark_runner(case.f, case.args; seconds)
        @printf(io, "  %-40s time_us=%12.2f  allocs=%10d  bytes=%12d\n", case.name,
                metrics["time_ns"] / 1e3, metrics["allocs"], metrics["bytes"])
        Dict("name" => case.name, "metrics" => metrics)
    end
    run_suite = suite -> begin
        println(io, "[", suite.name, "]")
        Dict("name" => suite.name, "cases" => run_case.(suite.cases))
    end
    return Dict("suites" => run_suite.(suites))
end

end


using .Benchmarks
using ArgParse
using Dates
using TOML

function parse_cli_args(args)
    s = ArgParseSettings(autofix_names = true)

    @add_arg_table s begin
        "--output", "-o"
            help = "Output TOML file path for benchmark results"
            arg_type = Union{Nothing, String}
        "--seconds", "-s"
            help = "Target benchmark sampling duration per case"
            arg_type = Float64
            default = 1.0
    end

    return ArgParse.parse_args(args, s; as_symbols = true)
end

function main(args)
    opts = parse_cli_args(args)
    payload = Benchmarks.run_benchmarks(seconds = opts[:seconds])
    if !isnothing(opts[:output])
        output = abspath(opts[:output])
        mkpath(dirname(output))
        open(output, "w") do io
            TOML.print(io, payload)
        end
        println("Saved benchmark results to: ", output)
    end

    return nothing
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
