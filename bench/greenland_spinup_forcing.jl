# Validate `forcing_climatology`'s two methods against REAL ERA5-Land forcing at glacierized
# Greenland sites.
#
#   julia --project=bench --threads=auto bench/greenland_spinup_forcing.jl
#
# The `method=:representative` conclusions in `src/forcing_climatology.jl` were measured on 21
# synthetic sites derived from a single parent site, with only two sites in the critical
# melt/accumulation band and two in ablation. That is thin evidence for regime boundaries, and
# the band it is thinnest in — melt approaching accumulation — is exactly Greenland's percolation
# zone. This script re-runs the comparison on real forcing at sites spanning the dry-snow,
# percolation and ablation zones, and reports whether the synthetic conclusions hold.
#
# What it measures per site, on both cycles (`:average` and `:representative`):
#   * annual melt on the cycle vs on the full record — does averaging destroy it?
#   * equilibrated firn air content and column-mean density — do the two agree?
#   * cycles to convergence.
#
# Requires a CDS API key (ERA5-Land access). Downloads are cached under `CACHE_DIR` so a re-run
# is offline; the cache is keyed by site and time range.

using GEMB
using GEMB_ClimateForcing
using Statistics
using Printf
using Serialization
using Dates

const DD = GEMB.DimensionalData

const YEARS = (2000, 2019)          # 20 complete years: enough for a climatology and a block
const CACHE_DIR = joinpath(@__DIR__, "greenland_forcing_cache")
const MAX_CYCLES = 800

# Convergence criteria for the equilibrium comparison. All four, because this script compares
# *equilibria* between two forcing cycles, and a step-only test cannot support that: a column
# creeping steadily at just under the tolerance passes it while still drifting, so two columns
# declared converged may simply be two columns stopped at different points on their way down.
#
# That is not hypothetical here. An earlier run of this script used `delta_density = 1e-2` alone
# and reported EastGRIP converging in 18 cycles under `:average` against 214 under
# `:representative`, with a 20.8% difference in equilibrated firn air content at a site with
# 0.2 kg m-2 yr-1 of melt. An 18-cycle exit on a step-only criterion is far more likely a false
# positive than an equilibrium, which would make that 20.8% a measurement of non-convergence
# rather than of the method.
#
# The FAC criteria matter for the same reason they exist (see `gemb_spinup`): a density tolerance
# admits a FAC residual proportional to column depth, so on a 200 m Greenland column a density
# test alone is loose in exactly the quantity this script tabulates.
const CONV_DELTA_DENSITY = 1e-3     # [kg m-3]
const CONV_DRIFT_DENSITY = 1e-3     # [kg m-3 per cycle]
const CONV_DELTA_FAC = 1e-4         # [m of air] — 0.1 mm
const CONV_DRIFT_FAC = 1e-4         # [m of air per cycle]

"""
    sites() -> Vector{NamedTuple}

Glacierized Greenland locations spanning the three regimes the synthetic fleet split into.

Coordinates are long-running PROMICE/GC-Net station locations, which are on the ice sheet by
construction and have published mass-balance regimes to compare against — the point of using
real sites rather than more synthetic ones. ERA5-Land is on a ~9 km grid, so these are the
reanalysis cell containing the station, not the station microclimate.
"""
sites() = [
    # --- Dry-snow zone: high, cold, little or no melt ---
    (id="Summit",        lat=72.58, lon=-38.46, zone="dry-snow"),
    (id="NEEM",          lat=77.45, lon=-51.06, zone="dry-snow"),
    (id="EastGRIP",      lat=75.63, lon=-35.99, zone="dry-snow"),
    # --- Percolation zone: melt that refreezes in firn. The critical band. ---
    (id="DYE-2",         lat=66.48, lon=-46.28, zone="percolation"),
    (id="Saddle",        lat=66.00, lon=-44.50, zone="percolation"),
    (id="CP1",           lat=69.88, lon=-46.99, zone="percolation"),
    (id="KAN-U",         lat=67.00, lon=-47.03, zone="percolation"),
    (id="SwissCamp",     lat=69.57, lon=-49.30, zone="percolation"),
    # --- Ablation zone: melt exceeds accumulation, bare ice in summer ---
    (id="KAN-M",         lat=67.07, lon=-48.84, zone="ablation"),
    (id="KAN-L",         lat=67.10, lon=-49.95, zone="ablation"),
    (id="QAS-L",         lat=61.03, lon=-46.85, zone="ablation"),
    (id="THU-L",         lat=76.40, lon=-68.27, zone="ablation"),
]

"""
    cds_token() -> String

The CDS API key, from `CDS_API_KEY`/`CDSAPI_KEY` or from `~/.cdsapirc`.

Reading `~/.cdsapirc` is the point: it is where the `cdsapi` client already keeps the key, so
this needs no environment setup. Never logged — only its length is reported.
"""
function cds_token()
    for var in ("CDS_API_KEY", "CDSAPI_KEY")
        v = get(ENV, var, "")
        isempty(v) || return String(strip(v))
    end
    path = expanduser("~/.cdsapirc")
    if isfile(path)
        for line in eachline(path)
            m = match(r"^\s*key\s*:\s*(\S+)\s*$", line)
            m === nothing || return String(m.captures[1])
        end
    end
    error("no CDS API key: set CDS_API_KEY or put `key: <key>` in ~/.cdsapirc")
end

# Downloaded forcing, cached on disk. ERA5-Land point extraction is the slow step and it is
# perfectly reproducible, so a re-run of the analysis should not repeat it.
function forcing_for(site, token)
    mkpath(CACHE_DIR)
    path = joinpath(CACHE_DIR, "$(site.id)_$(YEARS[1])_$(YEARS[2]).jls")
    if isfile(path)
        return deserialize(path)
    end
    ds = climate_forcing(:era5land, site.lat, site.lon;
        time_range=(DateTime(YEARS[1], 1, 1), DateTime(YEARS[2], 12, 31, 23)),
        token=token)
    serialize(path, ds)
    return ds
end

mean_density(p) = sum(parent(p[:density]) .* parent(p[:dz])) / sum(parent(p[:dz]))

# Annual melt [kg m-2 yr-1] integrated over `forcing`, on the column `profile`. The rate rather
# than the total, so a one-year cycle, an n-year block and a 20-year record are comparable.
function melt_rate(profile, forcing, mp_last)
    out = gemb(profile, forcing, mp_last; thermal_workspace=GEMB.ThermalWorkspace())
    years = length(DD.dims(forcing, DD.Ti)) * forcing.time_step / GEMB.SECONDS_PER_YEAR
    return years > 0 ? sum(parent(out[:melt])) / years : NaN
end

# Spin up on `cycle` and report the equilibrated column.
function equilibrate(profile, cycle, mp)
    g = gemb_spinup(profile, cycle, mp;
        max_iterations=MAX_CYCLES,
        convergence_delta_density=CONV_DELTA_DENSITY,
        convergence_drift_density=CONV_DRIFT_DENSITY,
        convergence_delta_fac=CONV_DELTA_FAC,
        convergence_drift_fac=CONV_DRIFT_FAC,
        thermal_workspace=GEMB.ThermalWorkspace())
    md = DD.metadata(g)
    # `years` is what the two methods must be compared on: cycle counts are not comparable
    # because `:average` integrates one year per cycle and `:representative` integrates
    # `n_years`, so a longer cycle mechanically lowers the count without lowering the work.
    n_per_cycle = length(DD.dims(cycle, DD.Ti)) * cycle.time_step / GEMB.SECONDS_PER_YEAR
    return (
        fac = firn_air_content(parent(g[:dz]), parent(g[:density]), mp.density_ice),
        density = mean_density(g),
        temperature_deep = parent(g[:temperature])[end],
        cycles = md[:spinup_cycles],
        years = md[:spinup_cycles] * n_per_cycle,
        converged = md[:spinup_converged],
    )
end

function analyze(site, token, mp)
    ds = forcing_for(site, token)
    cf = initialize_forcing(ds)

    mp_last = GEMB.ModelParameters(;
        (f => getfield(mp, f) for f in fieldnames(GEMB.ModelParameters)
         if f != :output_frequency)..., output_frequency=:last)

    profile = initialize_profile(mp, cf)
    cs = GEMB.initialize_climate_summary(cf, mp)

    avg = forcing_climatology(cf)
    rep = forcing_climatology(cf; method=:representative, model_parameters=mp,
                              n_years=3, verbose=false)

    # Melt on the real record is the reference every cycle is judged against.
    melt_record = melt_rate(profile, cf, mp_last)

    # How each cycle's own climate differs from the record's. This separates the two ways a
    # cycle can move the equilibrium: through melt, or through the temperature and accumulation
    # that drive densification directly. Densification is Arrhenius, so a fraction of a kelvin
    # shifts the equilibrium density profile with no melt involved at all — which is why the
    # dry-snow sites are not expected to be identical under the two methods even though neither
    # melts.
    cs_avg = GEMB.initialize_climate_summary(avg, mp)
    cs_rep = GEMB.initialize_climate_summary(rep, mp)

    return (
        id = site.id, zone = site.zone,
        elevation = get(DD.metadata(ds), "elevation", NaN),
        accumulation = cs.accumulation,
        melt_record = melt_record,
        ratio = cs.accumulation > 0 ? melt_record / cs.accumulation : NaN,
        # Cycle-minus-record climate offsets, per method.
        dT_avg = cs_avg.temperature_air_mean - cs.temperature_air_mean,
        dT_rep = cs_rep.temperature_air_mean - cs.temperature_air_mean,
        dacc_avg = cs.accumulation > 0 ? cs_avg.accumulation / cs.accumulation - 1 : NaN,
        dacc_rep = cs.accumulation > 0 ? cs_rep.accumulation / cs.accumulation - 1 : NaN,
        melt_avg = melt_rate(profile, avg, mp_last),
        melt_rep = melt_rate(profile, rep, mp_last),
        rep_years = DD.metadata(rep)[:climatology_representative_year],
        eq_avg = equilibrate(profile, avg, mp),
        eq_rep = equilibrate(profile, rep, mp),
    )
end

#=============================================================================
# Run
=============================================================================#

const MP = GEMB.ModelParameters(output_frequency=:last)
token = cds_token()
@printf("CDS token loaded (%d chars). Threads: %d\n", length(token), Threads.nthreads())
@printf("Sites: %d, years %d-%d, cache: %s\n\n", length(sites()), YEARS[1], YEARS[2], CACHE_DIR)

# Downloads are serial (one shared remote store, and the cache makes a re-run free); the
# analysis per site is independent but each spinup is single-threaded, so sites run concurrently.
all_sites = sites()
println("Fetching forcing...")
for s in all_sites
    try
        t = @elapsed forcing_for(s, token)
        @printf("  %-12s ok (%.1f s)\n", s.id, t)
    catch err
        @printf("  %-12s FAILED: %s\n", s.id, first(split(sprint(showerror, err), '\n')))
    end
end

results = Any[]
lk = ReentrantLock()
Threads.@threads for s in all_sites
    try
        r = analyze(s, token, MP)
        lock(lk) do
            push!(results, r)
            @printf("  analyzed %-12s ratio=%.2f\n", r.id, r.ratio)
        end
    catch err
        lock(lk) do
            @printf("  analyzed %-12s FAILED: %s\n", s.id,
                    first(split(sprint(showerror, err), '\n')))
        end
    end
end

sort!(results, by = r -> -r.ratio)

println("\n", "="^118)
println("MELT: does averaging destroy it, and does a representative block recover it?")
println("="^118)
@printf("%-12s %-12s %7s %8s %8s %6s | %9s %9s %9s\n",
        "site", "zone", "elev", "accum", "melt_rec", "m/a", "melt_avg", "melt_rep", "rep_yrs")
for r in results
    @printf("%-12s %-12s %7.0f %8.1f %8.1f %6.2f | %9.1f %9.1f %9s\n",
            r.id, r.zone, r.elevation, r.accumulation, r.melt_record, r.ratio,
            r.melt_avg, r.melt_rep, string(r.rep_years))
end

println("\n", "="^118)
println("EQUILIBRATED COLUMN: do the two cycles agree?")
println("(`conv` flags whether BOTH spinups met every criterion; a `no` invalidates that row's")
println(" comparison, since two unconverged columns can differ for reasons unrelated to method.)")
println("="^118)
@printf("%-12s %6s %9s | %8s %8s %7s | %8s %8s | %6s %6s | %4s\n",
        "site", "m/a", "melt_rec", "FAC_avg", "FAC_rep", "dFAC%", "rho_avg", "rho_rep",
        "yr_avg", "yr_rep", "conv")
for r in results
    dfac = r.eq_avg.fac > 0 ? 100 * (r.eq_rep.fac - r.eq_avg.fac) / r.eq_avg.fac : NaN
    both = r.eq_avg.converged && r.eq_rep.converged
    @printf("%-12s %6.2f %9.1f | %8.3f %8.3f %7.1f | %8.2f %8.2f | %6.0f %6.0f | %4s\n",
            r.id, r.ratio, r.melt_record, r.eq_avg.fac, r.eq_rep.fac, dfac,
            r.eq_avg.density, r.eq_rep.density, r.eq_avg.years, r.eq_rep.years,
            both ? "yes" : "NO")
end

println("\n", "="^118)
println("WHY THEY DIFFER: each cycle's own climate against the record's")
println("(A no-melt site can still shift, because densification is Arrhenius in temperature.)")
println("="^118)
@printf("%-12s %9s | %8s %8s | %9s %9s | %8s %8s\n",
        "site", "melt_rec", "dT_avg", "dT_rep", "dacc_avg", "dacc_rep", "Tdeep_avg", "Tdeep_rep")
for r in results
    @printf("%-12s %9.1f | %+8.2f %+8.2f | %+8.1f%% %+8.1f%% | %8.2f %8.2f\n",
            r.id, r.melt_record, r.dT_avg, r.dT_rep,
            100 * r.dacc_avg, 100 * r.dacc_rep,
            r.eq_avg.temperature_deep, r.eq_rep.temperature_deep)
end

println("\n", "="^118)
println("VERDICT vs the synthetic conclusions")
println("="^118)

melters = [r for r in results if r.melt_record > 1.0]
if !isempty(melters)
    tot = sum(r.melt_record for r in melters)
    @printf("Melt recovery, aggregate over %d melting sites: :average = %.0f%%, :representative = %.0f%%\n",
            length(melters), 100 * sum(r.melt_avg for r in melters) / tot,
            100 * sum(r.melt_rep for r in melters) / tot)
    @printf("  (synthetic fleet reported 0%% and 116%%)\n")
    # The aggregate is dominated by the highest-melt sites and hides the sites where averaging
    # fails completely, so report the per-site distribution too. This is the same
    # cancelling-errors trap the synthetic `:estimate` ranking fell into.
    ret = sort([100 * r.melt_avg / r.melt_record for r in melters])
    @printf("  :average per-site retention: median %.0f%%, and %d of %d sites below 10%%\n",
            median(ret), count(<(10), ret), length(ret))
    @printf("  :representative per-site: median %.0f%%, worst %.0f%%\n",
            median(sort([100 * r.melt_rep / r.melt_record for r in melters])),
            minimum(100 * r.melt_rep / r.melt_record for r in melters))
end

# The claim under test: below ~0.3 the two methods agree, above it they diverge, and a cycle
# whose ratio approaches 1 can lose the firn a `:average` spinup keeps.
for (label, sel) in (("ratio <= 0.3", r -> r.ratio <= 0.3),
                     ("0.3 < ratio <= 1.0", r -> 0.3 < r.ratio <= 1.0),
                     ("ratio > 1.0", r -> r.ratio > 1.0))
    grp = [r for r in results if isfinite(r.ratio) && sel(r)]
    isempty(grp) && continue
    dfacs = [abs(r.eq_avg.fac > 0 ? 100 * (r.eq_rep.fac - r.eq_avg.fac) / r.eq_avg.fac : NaN)
             for r in grp]
    dfacs = filter(isfinite, dfacs)
    collapses = count(r -> r.eq_rep.fac < 1.0 && r.eq_avg.fac >= 1.0, grp)
    @printf("%-20s n=%2d  median |dFAC| = %5.1f%%  max = %5.1f%%  firn lost by :representative: %d\n",
            label, length(grp), isempty(dfacs) ? NaN : median(dfacs),
            isempty(dfacs) ? NaN : maximum(dfacs), collapses)
end

nonconv = [r.id for r in results if !r.eq_avg.converged || !r.eq_rep.converged]
isempty(nonconv) || @printf("\nDid not converge in %d cycles: %s\n",
                            MAX_CYCLES, join(nonconv, ", "))
