module BenchCases

using Chairmarks
using Printf
using Random
using QuantumAlgebra

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

function make_cases()
    @boson_ops a
    @tlsxyz_ops σ

    H = ∑(:n, ∑(:m, ∑(:k, Pr"ω_n,m" * a(:n, :k)' * a(:m, :k)))) +
        ∑(:i, 1 // 2 * Pr"ν_i" * σz(:i)) +
        ∑(:n, ∑(:k, ∑(:i, Pr"g_i,n,k" * σx(:i) * (a(:n, :k)' + a(:n, :k)))))

    Acomm = σz(:j)

    commutator_cases = [
        (name = "commutator_n2", f = nested_comm, args = (H, Acomm, 2)),
        (name = "commutator_n4", f = nested_comm, args = (H, Acomm, 4)),
    ]
        
    Random.seed!(4)
    # Shuffle so no operator appears twice (as σz() * σz() == 1).
    randops = shuffle!([a.(1:50)..., σz.(1:50)..., adjoint.(a.(1:50))...])
    
    prodcorr_cases = [
        begin
            ops = Tuple(sort(randops[1:n], by = op -> first(op.terms)[1]))
            expr = normal_form(8 * Pr"g_i" * prod(ops))
            (name = "prodcorr_n$(n)", f = expval_as_corrs, args = (expr,))
        end for n in 2:9
    ]

    return [
        (name = "commutator", cases = commutator_cases),
        (name = "prodcorr", cases = prodcorr_cases),
    ]
end

function run_benchmarks(; seconds::Real = 1.0, io::IO = stdout)
    suites = make_cases()

    run_case = case -> begin
        metrics = benchmark_runner(case.f, case.args; seconds)
        time_us = metrics["time_ns"] / 1e3
        allocs = round(Int, metrics["allocs"])
        bytes = round(Int, metrics["bytes"])
        @printf(io, "  %-20s time_us=%12.2f  allocs=%10d  bytes=%12d\n", case.name, time_us, allocs, bytes)
        Dict("name" => case.name, "metrics" => metrics)
    end
    run_suite = suite -> begin
        println(io, "[", suite.name, "]")
        Dict("name" => suite.name, "cases" => run_case.(suite.cases))
    end
    return run_suite.(suites)
end

end

if abspath(PROGRAM_FILE) == @__FILE__
    using .BenchCases
    BenchCases.run_benchmarks()
end