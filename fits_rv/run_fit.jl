# Fit one system:  julia --project=. fits_rv/run_fit.jl <system_id> [options]
#
#   --config PATH        run configuration (default configs/pilot.toml)
#   --catalog PATH       stage-1 catalog (default sample/stage1_catalog.csv)
#   --rv-dir PATH        per-system RV tables (default rv_data/systems)
#   --out PATH           output root (default fits_rv/results)
#   --n-companions N     override the literature companion count
#   --test               short run (a few minutes): 2 chains x 200 draws, RVs only,
#                        no G23H download, escalation off
#
# <system_id> is the host's Gaia DR3 source ID (the key of rv_data/systems/<id>_*.csv).
# Outputs go to <out>/<display name>/ (e.g. HD111232): chain.csv, convergence.csv, attempts.csv, jitter.csv.

include(joinpath(@__DIR__, "pipeline.jl"))

"Readable name for a system (sample/system_names.csv), else the source ID itself."
function display_name(system_id)
    f = joinpath(@__DIR__, "..", "sample", "system_names.csv")
    isfile(f) || return system_id
    t = CSV.read(f, DataFrame; types = String)
    i = findfirst(==(system_id), t.system_id)
    return i === nothing ? system_id : t.display_name[i]
end

function parse_args(args)
    opts = Dict("config" => "configs/pilot.toml", "catalog" => "sample/stage1_catalog.csv",
                "rv-dir" => "rv_data/systems", "out" => "fits_rv/results",
                "n-companions" => nothing, "test" => false)
    positional = String[]
    i = 1
    while i <= length(args)
        a = args[i]
        if a == "--test"
            opts["test"] = true
        elseif startswith(a, "--")
            opts[a[3:end]] = args[i += 1]
        else
            push!(positional, a)
        end
        i += 1
    end
    length(positional) == 1 || error("usage: run_fit.jl <system_id> [options]")
    return positional[1], opts
end

function main(args)
    system_id, opts = parse_args(args)
    cfg = TOML.parsefile(opts["config"])
    if opts["test"]
        cfg["astrometry"]["mode"] = "none"
        cfg["sampler"]["adaptation"] = 200
        cfg["sampler"]["iterations"] = 200
        cfg["sampler"]["n_chains"] = 2
        cfg["sampler"]["max_escalations"] = 0
        ## Short chains cannot reach production ESS; the test checks that the machinery runs
        cfg["convergence"]["ess_bulk_min"] = 50
        cfg["convergence"]["ess_tail_min"] = 50
        cfg["convergence"]["indep_per_chain_min"] = 10
    end
    files = filter(f -> startswith(f, system_id * "_") && endswith(f, ".csv"),
                   readdir(opts["rv-dir"]))
    length(files) == 1 || error("expected one RV table for $system_id in $(opts["rv-dir"]), found $files")
    rvs = load_rvs(joinpath(opts["rv-dir"], files[1]))
    star = stellar_info(opts["catalog"], system_id)
    n = opts["n-companions"] === nothing ? star.n_companions : parse(Int, opts["n-companions"])
    name = replace(files[1], ".csv" => "")
    fit_system(rvs, star; cfg, name, n_companions = n,
               outdir = joinpath(opts["out"], display_name(system_id)))
    return nothing
end

main(ARGS)
