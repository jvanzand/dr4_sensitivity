# Octofitter pipeline for the "blind" ground-truth fits (PLAN.md, Step 3): one model per
# system, built from rv_data/systems/ and the stage-1 catalog, sampled with NUTS, with
# convergence checks after every fit. Loaded by run_fit.jl.
#
# Octofitter v9. The priors come from configs/*.toml; no literature orbital parameter ever
# enters a prior or a starting value (PLAN.md, Section 2.6). The number of companions is
# the count the literature reports (the census), or whatever the blind RV search
# (blind_search/rv_search.py) found; their orbits are blind.

using Octofitter
using OctofitterRadialVelocity
using Distributions
using CSV
using DataFrames
using Statistics
using Random
using TOML

const MCMCChains = Octofitter.MCMCChains
const MJD_OFFSET = 2400000.5          # BJD - MJD

# ---------------------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------------------

"""
    load_rvs(path) -> DataFrame

RVs of one system (rv_data/systems/*.csv): junk-flagged rows are dropped, times are BJD_TCB
converted to MJD (TCB), and the columns Octofitter expects are added (`epoch`, `rv`, `σ_rv`).
"""
function load_rvs(path)
    t = CSV.read(path, DataFrame; missingstring=["", "NA"], types=Dict(:isjunk => Bool))
    t = t[(.!t.isjunk) .& .!ismissing.(t.rv_ms) .& .!ismissing.(t.rv_err_ms), :]
    return DataFrame(
        dataset = String.(t.dataset),
        epoch = Float64.(t.bjd_tcb) .- MJD_OFFSET,
        rv = Float64.(t.rv_ms),
        σ_rv = Float64.(t.rv_err_ms),
    )
end

"Valid Julia identifier-like name for a dataset (the chain's column names embed it)."
function safe_name(s, taken)
    base = replace(String(s), r"[^A-Za-z0-9]+" => "_")
    base = strip(base, '_')
    name, k = base, 1
    while name in taken
        k += 1
        name = "$(base)_$k"
    end
    return name
end

"""
    stellar_info(catalog_csv, system_id) -> (; mass, mass_err, n_companions, gaia_id)

Host mass and its error (a published error, or the stand-in of census.py; the census records
which), and the number of companions the census holds for the system (the literature count,
used only to set how many Keplerians to fit).
"""
function stellar_info(catalog_csv, system_id)
    c = CSV.read(catalog_csv, DataFrame; types=Dict(:system_id => String), missingstring=["", "NA"])
    s = c[c.system_id .== system_id, :]
    nrow(s) == 0 && error("system $system_id is not in $catalog_csv")
    gaia_id = tryparse(Int64, system_id)
    return (; mass = Float64(s.st_mass[1]), mass_err = Float64(s.st_mass_err[1]),
            n_companions = nrow(s), gaia_id)
end

# ---------------------------------------------------------------------------------------
# Model
# ---------------------------------------------------------------------------------------

"""
    build_model(rvs, star; n_companions, cfg, name) -> (model, info)

Host body A and `n_companions` Keplerian bodies b, c, ..., with blind priors from `cfg`.
Each RV dataset gets its own jitter (log-uniform) and an analytically marginalized offset
(`MarginalizedRVObs`). With `cfg["astrometry"]["mode"] == "g23h"` the G23H likelihood is
added, the companions' inclinations are isotropic (`Sine()`), and Ω is free; for RVs alone
the orbits are edge-on (i = π/2), so `mass` is M sin i.
"""
function build_model(rvs, star; n_companions, cfg, name)
    pri = cfg["priors"]
    mode = cfg["astrometry"]["mode"]
    use_g23h = mode == "g23h"
    pmin = pri["period_min_known_d"]
    pmax = pri["period_max_yr"] * 365.25
    qmin, qmax, emax = pri["q_min"], pri["q_max"], pri["e_max"]
    jmin, jmax = pri["jitter_min"], pri["jitter_max"]
    tref = mean(rvs.epoch)

    A = Body(
        name = "A",
        variables = (@variables begin
            mass = system.M_pri
            flux_G = 1.0
            flux_Hp = 1.0
        end),
    )

    letters = collect('b':'z')
    n_companions <= length(letters) || error("too many companions")
    companions = map(1:n_companions) do k
        nm = string(letters[k])
        if use_g23h
            Body(
                name = nm, about = A,
                variables = (@variables begin
                    P ~ LogUniform(pmin, pmax)
                    q ~ LogUniform(qmin, qmax)
                    mass = q * system.M_pri
                    e ~ Uniform(0, emax)
                    ω ~ Uniform(0, 2pi)
                    i ~ Sine()
                    Ω ~ Uniform(0, 2pi)
                    M0 ~ Uniform(0, 2pi)
                    epoch = $tref
                    flux_G = 0.0
                    flux_Hp = 0.0
                end),
            )
        else
            Body(
                name = nm, about = A,
                variables = (@variables begin
                    P ~ LogUniform(pmin, pmax)
                    q ~ LogUniform(qmin, qmax)
                    mass = q * system.M_pri
                    e ~ Uniform(0, emax)
                    ω ~ Uniform(0, 2pi)
                    i = pi/2
                    Ω = 0.0
                    M0 ~ Uniform(0, 2pi)
                    epoch = $tref
                    flux_G = 0.0
                    flux_Hp = 0.0
                end),
            )
        end
    end

    ## One likelihood per RV dataset
    taken = String[]
    rvlikes = []
    for (ds, df) in pairs(groupby(rvs, :dataset))
        nm = safe_name(ds.dataset, taken)
        push!(taken, nm)
        tab = Table(epoch = df.epoch, rv = df.rv, σ_rv = df.σ_rv)
        push!(rvlikes, MarginalizedRVObs(
            tab; target = A, ref = Barycentre, name = nm,
            variables = (@variables begin
                jitter ~ LogUniform(jmin, jmax)
            end),
        ))
    end

    M, σM = star.mass, star.mass_err
    observations = Any[rvlikes...]
    if use_g23h
        star.gaia_id === nothing && error("G23H needs a Gaia DR3 source ID as the system ID")
        ## Local copies of the G23H tables (data/g23h/, git-ignored) are used when the config
        ## names them; otherwise Octofitter's DataDeps would download them to ~/.julia.
        local_path(key) = haskey(cfg["astrometry"], key) ?
            normpath(joinpath(@__DIR__, "..", cfg["astrometry"][key])) : nothing
        g23h = G23HObs(; gaia_id = star.gaia_id, target = A, blends = Tuple(companions),
                       ref = Barycentre, freeze_epochs = cfg["astrometry"]["freeze_epochs"],
                       catalog = local_path("g23h_catalog"),
                       dr2_transits_catalog = local_path("g23h_dr2_transits"))
        cat = g23h.catalog
        gid = star.gaia_id
        sys = System(
            name = name, bodies = [A, companions...],
            observations = [g23h, observations...],
            variables = (@variables begin
                M_pri ~ truncated(Normal(M, σM), lower = 0.05)
                plx ~ gaia_plx(; gaia_id = gid)
                ra = $(cat.ra)
                dec = $(cat.dec)
                pmra ~ Normal(cat.pmra_dr3, 10)
                pmdec ~ Normal(cat.pmdec_dr3, 10)
                rv = 0.0
                ref_epoch = 57388.5      # Gaia DR3 reference epoch, J2016.0 (MJD)
            end),
        )
    else
        sys = System(
            name = name, bodies = [A, companions...], observations = observations,
            variables = (@variables begin
                M_pri ~ truncated(Normal(M, σM), lower = 0.05)
            end),
        )
    end
    model = Octofitter.LogDensityModel(sys; verbosity = 0)
    return model, (; n_params = model.D, datasets = taken, tref)
end

# ---------------------------------------------------------------------------------------
# Sampling and convergence
# ---------------------------------------------------------------------------------------

"Run `n_chains` independent NUTS chains (different seeds) and join them."
function run_chains(model, scfg; adaptation, iterations, verbosity = 0)
    chains = map(1:scfg["n_chains"]) do k
        rng = Random.Xoshiro(scfg["seed"] + 1000k)
        octofit(rng, model, scfg["target_accept"]; adaptation, iterations,
                max_depth = scfg["max_depth"], verbosity)
    end
    return length(chains) == 1 ? chains[1] : MCMCChains.chainscat(chains...)
end

"Energy fraction of missing information for one chain (low = poor energy exploration)."
function ebfmi(E)
    return sum(abs2, diff(E)) / sum(abs2, E .- mean(E))
end

"""
    check_convergence(chain, ccfg) -> (; ok, report)

R-hat, bulk/tail ESS, effective independent draws per chain, divergences, tree-depth
saturation, and E-BFMI, against the thresholds in the config. Only the sampled parameters
are tested (not derived or internal columns).
"""
function check_convergence(chain, ccfg, max_depth)
    internals = Set([:n_steps, :is_accept, :acceptance_rate, :hamiltonian_energy,
                     :hamiltonian_energy_error, :max_hamiltonian_energy_error, :tree_depth,
                     :numerical_error, :step_size, :nom_step_size, :is_adapt, :loglike,
                     :logprior, :logpost])
    pars = [p for p in MCMCChains.names(chain, :parameters) if p ∉ internals]
    sub = chain[:, pars, :]
    rh = MCMCChains.rhat(sub)
    eb = MCMCChains.ess(sub; kind = :bulk)
    et = MCMCChains.ess(sub; kind = :tail)
    rhat_vals = filter(isfinite, rh[:, :rhat])
    ess_b = filter(isfinite, eb[:, :ess])
    ess_t = filter(isfinite, et[:, :ess])
    m = size(chain, 3)
    ndraw = size(chain, 1)
    worst_rhat = isempty(rhat_vals) ? NaN : maximum(rhat_vals)
    min_ess_b = isempty(ess_b) ? NaN : minimum(ess_b)
    min_ess_t = isempty(ess_t) ? NaN : minimum(ess_t)
    indep = min_ess_b / m
    nerr = vec(Array(chain[:numerical_error]))
    td = vec(Array(chain[:tree_depth]))
    divergent = mean(nerr)
    saturated = mean(td .== max_depth)
    E = Array(chain[:hamiltonian_energy])
    min_ebfmi = minimum(ebfmi(E[:, c]) for c in 1:size(E, 2))
    checks = [
        ("R-hat <= $(ccfg["rhat_max"])", worst_rhat <= ccfg["rhat_max"], worst_rhat),
        ("bulk ESS >= $(ccfg["ess_bulk_min"])", min_ess_b >= ccfg["ess_bulk_min"], min_ess_b),
        ("tail ESS >= $(ccfg["ess_tail_min"])", min_ess_t >= ccfg["ess_tail_min"], min_ess_t),
        ("independent draws per chain >= $(ccfg["indep_per_chain_min"])",
         indep >= ccfg["indep_per_chain_min"], indep),
        ("divergent fraction <= $(ccfg["max_divergent_frac"])",
         divergent <= ccfg["max_divergent_frac"], divergent),
        ("max-tree-depth fraction <= $(ccfg["max_treedepth_frac"])",
         saturated <= ccfg["max_treedepth_frac"], saturated),
        ("E-BFMI >= $(ccfg["ebfmi_min"])", min_ebfmi >= ccfg["ebfmi_min"], min_ebfmi),
    ]
    return (; ok = all(c[2] for c in checks), checks, n_chains = m, n_draws = ndraw,
            worst_rhat, min_ess_bulk = min_ess_b, min_ess_tail = min_ess_t)
end

"""
    jitter_report(chain, rvs, datasets, pcfg) -> Vector{NamedTuple}

Posterior-median jitter per dataset, against the median formal error and the prior bound.
A jitter far above the formal errors, or at the upper bound, often means the model fits the
data poorly (a missing companion, a wrong period), so such datasets are flagged.
"""
function jitter_report(chain, rvs, datasets, pcfg)
    out = NamedTuple[]
    groups = collect(groupby(rvs, :dataset))
    for (nm, g) in zip(datasets, groups)
        col = Symbol("$(nm)_jitter")
        if col ∉ MCMCChains.names(chain)
            continue
        end
        j = vec(Array(chain[col]))
        med = median(j)
        ferr = median(g.σ_rv)
        flag = med > pcfg["jitter_flag_factor"] * ferr || quantile(j, 0.84) > 0.95 * pcfg["jitter_max"]
        push!(out, (; dataset = nm, n = nrow(g), jitter_median = med,
                    jitter_p16 = quantile(j, 0.16), jitter_p84 = quantile(j, 0.84),
                    median_formal_error = ferr, flagged = flag))
    end
    return out
end

"""
    fit_system(rvs, star; cfg, name, outdir, test=false)

Build the model, sample with NUTS, check convergence, and on failure rerun with doubled
iterations (up to `max_escalations`). Writes the chain and a diagnostics table to `outdir`.
"""
function fit_system(rvs, star; cfg, name, outdir, n_companions = star.n_companions,
                    verbosity = 1)
    mkpath(outdir)
    model, info = build_model(rvs, star; n_companions, cfg, name)
    scfg, ccfg = cfg["sampler"], cfg["convergence"]
    verbosity > 0 && println("$name: $(nrow(rvs)) RVs in $(length(info.datasets)) datasets, ",
                             "$n_companions companion(s), $(info.n_params) parameters")
    attempts = NamedTuple[]
    chain, conv = nothing, nothing
    for k in 0:scfg["max_escalations"]
        iters = scfg["iterations"] * 2^k
        t0 = time()
        chain = run_chains(model, scfg; adaptation = scfg["adaptation"], iterations = iters)
        conv = check_convergence(chain, ccfg, scfg["max_depth"])
        push!(attempts, (; attempt = k + 1, iterations = iters, seconds = time() - t0,
                         ok = conv.ok, worst_rhat = conv.worst_rhat,
                         min_ess_bulk = conv.min_ess_bulk))
        if verbosity > 0
            println("  attempt $(k+1): $(iters) iterations x $(conv.n_chains) chains, ",
                    "$(round(attempts[end].seconds; digits = 1)) s, ",
                    conv.ok ? "converged" : "NOT converged")
            for (label, ok, val) in conv.checks
                println("    ", ok ? "ok  " : "FAIL", " ", label, "   (value ", round(val; sigdigits = 3), ")")
            end
        end
        conv.ok && break
    end
    jit = jitter_report(chain, rvs, info.datasets, cfg["priors"])
    if verbosity > 0
        for r in jit
            println("  jitter ", r.dataset, ": median ", round(r.jitter_median; digits = 2),
                    " m/s (formal error ", round(r.median_formal_error; digits = 2), ")",
                    r.flagged ? "  <-- FLAGGED" : "")
        end
    end
    ## Octofitter v9's `savehdf5` is an orbitize export and needs a `plx` column, which
    ## RV-only models lack; a plain CSV of the chain keeps one layout for every mode
    CSV.write(joinpath(outdir, "chain.csv"), DataFrame(chain))
    CSV.write(joinpath(outdir, "attempts.csv"), DataFrame(attempts))
    CSV.write(joinpath(outdir, "jitter.csv"), DataFrame(jit))
    CSV.write(joinpath(outdir, "convergence.csv"),
              DataFrame(check = [c[1] for c in conv.checks], ok = [c[2] for c in conv.checks],
                        value = [c[3] for c in conv.checks]))
    return (; chain, conv, attempts, jitter = jit, model)
end
