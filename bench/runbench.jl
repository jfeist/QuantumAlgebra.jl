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
        "--baseline", "-b"
            help = "Baseline TOML file path to compare against"
            arg_type = String
        "--current", "-c"
            help = "Current TOML file path (skip benchmark execution and compare only)"
            arg_type = String
        "--seconds", "-s"
            help = "Target benchmark sampling duration per case"
            arg_type = Float64
            default = 1.0
        "--threshold", "-t"
            help = "Regression threshold ratio (e.g. 0.10 means 10% slower)"
            arg_type = Float64
            default = 0.10
        "--fail-on-regression"
            help = "Exit with error if regressions above threshold are detected"
            action = :store_true
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

function flatten_suites(suites)
    rows = Vector{Tuple{String, Float64}}()
    for suite in suites
        sname = String(suite["name"])
        for case in suite["cases"]
            cname = String(case["name"])
            metrics = case["metrics"]
            push!(rows, ("$(sname).$(cname)", Float64(metrics["time_ns"])))
        end
    end
    return rows
end

function compare_with_baseline(current_suites, baseline_suites; threshold = 0.10)
    regressions = Vector{NamedTuple{(:name, :baseline, :current, :ratio), Tuple{String, Float64, Float64, Float64}}}()

    now_rows = flatten_suites(current_suites)
    base_rows = Dict(flatten_suites(baseline_suites))

    for (name, new_val) in now_rows
        if !haskey(base_rows, name)
            continue
        end
        old_val = base_rows[name]
        if old_val <= 0
            continue
        end
        ratio = new_val / old_val
        if ratio > 1 + threshold
            push!(regressions, (name = name, baseline = old_val, current = new_val, ratio = ratio))
        end
    end

    return regressions
end

function print_regressions(regressions)
    if isempty(regressions)
        println("No regressions above threshold detected.")
        return
    end

    println("Regressions detected:")
    for r in regressions
        percent = (r.ratio - 1.0) * 100
        println(
            "  ",
            r.name,
            ": ",
            round(r.baseline; digits = 2),
            " ns -> ",
            round(r.current; digits = 2),
            " ns (",
            round(percent; digits = 1),
            "% slower)",
        )
    end
end

function main(args)
    opts = parse_cli_args(args)
    threshold = opts[:threshold]
    fail_on_regression = opts[:fail_on_regression]

    current_opt = get(opts, :current, nothing)
    baseline_opt = get(opts, :baseline, nothing)

    current_payload = if !isnothing(current_opt)
        current_path = abspath(current_opt)
        if !isfile(current_path)
            error("current results file not found: $(current_path)")
        end
        TOML.parsefile(current_path)
    else
        seconds = opts[:seconds]
        payload = benchmark_payload(seconds = seconds)
        output = abspath(opts[:output])
        mkpath(dirname(output))
        open(output, "w") do io
            TOML.print(io, payload)
        end
        println("Saved benchmark results to: ", output)
        payload
    end

    current_suites = current_payload["suites"]

    if !isnothing(baseline_opt)
        baseline_path = abspath(baseline_opt)
        if !isfile(baseline_path)
            error("baseline file not found: $(baseline_path)")
        end
        baseline = TOML.parsefile(baseline_path)
        regressions = compare_with_baseline(current_suites, baseline["suites"]; threshold = threshold)
        print_regressions(regressions)

        if fail_on_regression && !isempty(regressions)
            error("benchmark regressions detected")
        end
    end

    return nothing
end

main(ARGS)