using ArgParse
using Dates
using TOML

include("bench_cases.jl")

function parse_cli_args(args)
    s = ArgParseSettings(autofix_names = true)

    @add_arg_table s begin
        "--output", "-o"
            help = "Output TOML file path for benchmark results"
            arg_type = String
            default = joinpath(@__DIR__, "results", "latest.toml")
        "--seconds", "-s"
            help = "Target benchmark sampling duration per case"
            arg_type = Float64
            default = 1.0
    end

    return ArgParse.parse_args(args, s; as_symbols = true)
end

function benchmark_payload(; seconds)
    suites = BenchCases.run_benchmarks(seconds = seconds)

    commit = get(ENV, "GITHUB_SHA", "")
    branch = get(ENV, "GITHUB_REF_NAME", "")

    return Dict(
        "meta" => Dict(
            "generated_at" => string(now(UTC)),
            "julia_version" => string(VERSION),
            "hostname" => gethostname(),
            "commit" => commit,
            "branch" => branch,
            "seconds" => seconds,
        ),
        "suites" => suites,
    )
end

function main(args)
    opts = parse_cli_args(args)
    payload = benchmark_payload(seconds = opts[:seconds])
    output = abspath(opts[:output])
    mkpath(dirname(output))
    open(output, "w") do io
        TOML.print(io, payload)
    end
    println("Saved benchmark results to: ", output)

    return nothing
end

main(ARGS)