using ArgParse
using Dates
using TOML

const BENCH_SCRIPT = abspath(joinpath(@__DIR__, "runbench.jl"))
const BENCH_PROJECT = abspath(joinpath(@__DIR__, "Project.toml"))
const DEFAULT_SCRATCH_DIR = joinpath(@__DIR__, ".bench_over_commits")

function parse_cli_args(args)
    s = ArgParseSettings(autofix_names = true)

    @add_arg_table s begin
        "revspec"
            help = "Git revision spec passed to `git rev-list --reverse` (e.g. main~20..main)"
            arg_type = String
            nargs = '*'
            default = ["HEAD"]
        "--seconds", "-s"
            help = "Benchmark sampling duration per case"
            arg_type = Float64
            default = 1.0
        "--database", "-d"
            help = "CSV output path (defaults to machine+Julia specific file)"
            arg_type = String
        "--rerun-existing"
            help = "Re-run commits already present in database"
            action = :store_true
        "--scratch-dir"
            help = "Temporary project-local directory containing the git worktree and benchmark environment"
            arg_type = String
            default = DEFAULT_SCRATCH_DIR
    end

    return ArgParse.parse_args(args, s; as_symbols = true)
end

function sanitize_slug(s::AbstractString)
    return replace(s, r"[^A-Za-z0-9._-]" => "_")
end

function default_database_path()
    host = sanitize_slug(gethostname())
    julia_id = sanitize_slug(string(VERSION))
    db_dir = joinpath(@__DIR__, "results")
    mkpath(db_dir)
    return joinpath(db_dir, "bench_history_$(host)_julia_$(julia_id).csv")
end

function git_output(repo::AbstractString, args::Vector{String})
    cmd = Cmd(vcat(["git", "-C", String(repo)], args))
    return chomp(read(cmd, String))
end

function collect_commits(repo::AbstractString, revspec::Vector{String})
    args = ["rev-list", "--reverse", revspec...]
    out = git_output(repo, args)
    if isempty(out)
        return String[]
    end
    return split(out, '\n')
end

function commit_metadata(repo::AbstractString, commit::AbstractString)
    out = git_output(repo, ["show", "-s", "--format=%H,%ct,%cI", String(commit)])
    parts = split(out, ',')
    if length(parts) != 3
        error("failed to parse commit metadata for $(commit)")
    end
    return (hash = parts[1], unix = parts[2], iso = parts[3])
end

function read_existing_commits(db_path::AbstractString)
    if !isfile(db_path)
        return Set{String}()
    end
    commits = Set{String}()
    open(db_path, "r") do io
        first_line = true
        for line in eachline(io)
            if first_line
                first_line = false
                continue
            end
            isempty(line) && continue
            cols = split(line, ',')
            isempty(cols) && continue
            push!(commits, cols[1])
        end
    end
    return commits
end

function ensure_database_header(db_path::AbstractString)
    mkpath(dirname(db_path))
    if !isfile(db_path)
        open(db_path, "w") do io
            println(io, "commit,commit_unix,commit_iso,julia_version,hostname,suite,case,time_ns,allocs,bytes")
        end
    end
end

function append_benchmark_rows(db_path::AbstractString, meta, payload)
    suites = payload["suites"]
    open(db_path, "a") do io
        for suite in suites
            suite_name = String(suite["name"])
            for case in suite["cases"]
                case_name = String(case["name"])
                metrics = case["metrics"]
                println(
                    io,
                    string(
                        meta.hash, ",",
                        meta.unix, ",",
                        meta.iso, ",",
                        VERSION, ",",
                        gethostname(), ",",
                        suite_name, ",",
                        case_name, ",",
                        Float64(metrics["time_ns"]), ",",
                        Int(round(metrics["allocs"])), ",",
                        Int(round(metrics["bytes"])),
                    ),
                )
            end
        end
    end
end

function ensure_worktree(repo::AbstractString, worktree::AbstractString, commit::AbstractString)
    run(`git -C $(repo) worktree prune`)
    if isdir(worktree)
        rm(worktree; recursive = true, force = true)
    end
    run(`git -C $(repo) worktree add --force --detach $(worktree) $(commit)`)
end

function checkout_in_worktree(worktree::AbstractString, commit::AbstractString)
    run(`git -C $(worktree) checkout --detach --force $(commit)`)
end

function setup_benchmark_env(env_dir::AbstractString, worktree::AbstractString)
    if isdir(env_dir)
        rm(env_dir; recursive = true, force = true)
    end
    mkpath(env_dir)
    cp(BENCH_PROJECT, joinpath(env_dir, "Project.toml"); force = true)
    if isfile(joinpath(env_dir, "Manifest.toml"))
        rm(joinpath(env_dir, "Manifest.toml"); force = true)
    end
    setup_expr = "using Pkg; Pkg.develop(PackageSpec(path=\"$(worktree)\")); Pkg.instantiate()"
    setup_cmd = `$(Base.julia_cmd()) --project=$(env_dir) -e $(setup_expr)`
    run(setup_cmd)
end

function run_benchmark_in_worktree(seconds::Real, env_dir::AbstractString)
    output_file, io = mktemp()
    close(io)
    cmd = `$(Base.julia_cmd()) --project=$(env_dir) $(BENCH_SCRIPT) --output $(output_file) --seconds $(seconds)`
    run(cmd)
    payload = TOML.parsefile(output_file)
    rm(output_file; force = true)
    return payload
end

function main(args)
    opts = parse_cli_args(args)

    repo = abspath(joinpath(@__DIR__, ".."))
    db_opt = get(opts, :database, nothing)
    db_path = isnothing(db_opt) ? default_database_path() : abspath(String(db_opt))
    rerun_existing = Bool(opts[:rerun_existing])
    seconds = opts[:seconds]
    scratch_dir = abspath(String(opts[:scratch_dir]))
    worktree = joinpath(scratch_dir, "worktree")
    bench_env = joinpath(scratch_dir, "env")

    revspec = Vector{String}(opts[:revspec])
    commits = collect_commits(repo, revspec)
    isempty(commits) && error("No commits matched revision spec: $(join(revspec, " "))")

    ensure_database_header(db_path)
    existing = read_existing_commits(db_path)

    println("Database: ", db_path)
    println("Commits matched: ", length(commits))

    mkpath(scratch_dir)
    ensure_worktree(repo, worktree, commits[1])
    try
        for commit in commits
            if !rerun_existing && (commit in existing)
                println("Skipping existing commit: ", commit)
                continue
            end

            checkout_in_worktree(worktree, commit)
            meta = commit_metadata(repo, commit)
            println("Running benchmarks for ", meta.hash, " (", meta.iso, ")")

            setup_benchmark_env(bench_env, worktree)
            payload = run_benchmark_in_worktree(seconds, bench_env)
            append_benchmark_rows(db_path, meta, payload)
            push!(existing, commit)
        end
    finally
        run(`git -C $(repo) worktree remove --force $(worktree)`)
        rm(scratch_dir; recursive = true, force = true)
    end

    println("Done.")
end

main(ARGS)
