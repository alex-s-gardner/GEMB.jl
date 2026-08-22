# Shared fixtures for the `bench/greenland_*.jl` harnesses.
#
# Each of those scripts asks a different question of the same 12 ERA5-Land sites, so the site
# list, the download cache contract, and the "how much did this cycle melt" measurement have to
# be identical between them or the comparisons are not comparable. They were copy-pasted at
# first, which put five copies of the cache path expression and four of the site list in the
# tree — a site added to `sites()` would have been downloaded and then never analyzed.
#
# `include("greenland_common.jl")` from each harness.

using GEMB
using GEMB_ClimateForcing
using Statistics
using Printf
using Serialization
using Dates

const DD = GEMB.DimensionalData

const YEARS = (2000, 2019)          # 20 complete years: enough for a climatology and a block
const CACHE_DIR = joinpath(@__DIR__, "greenland_forcing_cache")

"""
    sites() -> Vector{NamedTuple}

Glacierized Greenland locations spanning the three mass-balance regimes.

Coordinates are long-running PROMICE/GC-Net station locations, which are on the ice sheet by
construction and have published regimes to compare against — the point of using real sites
rather than more synthetic ones. ERA5-Land is on a ~9 km grid, so these are the reanalysis cell
containing the station, not the station microclimate.
"""
sites() = [
    # --- Dry-snow zone: high, cold, little or no melt ---
    (id="Summit",        lat=72.58, lon=-38.46, zone="dry-snow"),
    (id="NEEM",          lat=77.45, lon=-51.06, zone="dry-snow"),
    (id="EastGRIP",      lat=75.63, lon=-35.99, zone="dry-snow"),
    # --- Percolation zone: melt that refreezes in firn ---
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

site_ids() = [s.id for s in sites()]

cache_path(id) = joinpath(CACHE_DIR, "$(id)_$(YEARS[1])_$(YEARS[2]).jls")

"""
    forcing_for(site, token) -> DimStack

Downloaded ERA5-Land forcing for one site, cached on disk. Point extraction is the slow step and
is perfectly reproducible, so re-running an analysis must not repeat it.
"""
function forcing_for(site, token)
    mkpath(CACHE_DIR)
    path = cache_path(site.id)
    isfile(path) && return deserialize(path)
    ds = climate_forcing(:era5land, site.lat, site.lon;
        time_range=(DateTime(YEARS[1], 1, 1), DateTime(YEARS[2], 12, 31, 23)),
        token=token)
    serialize(path, ds)
    return ds
end

# Read-only accessor for the harnesses that only ever consume what `forcing_for` wrote, so they
# neither need a CDS key nor can silently miss the cache through a diverging path expression.
function cached_forcing(id)
    path = cache_path(id)
    isfile(path) ||
        error("no cached forcing for $id; run bench/greenland_spinup_forcing.jl first")
    return deserialize(path)
end

# `mp` with output forced to a single final step, which is what every per-year and per-cycle
# measurement here wants. Mirrors what `gemb_spinup` does internally.
last_step_parameters(mp) = GEMB.ModelParameters(;
    (f => getfield(mp, f) for f in fieldnames(GEMB.ModelParameters)
     if f != :output_frequency)..., output_frequency=:last)

"""
    melt_rate(profile, forcing, mp_last) -> kg m-2 yr-1
    smb_rate(profile, forcing, mp_last) -> m of ice per year

Annual melt and surface mass balance over `forcing`, integrated on `profile`.

Both delegate to `GEMB._year_scores`, which is what `forcing_climatology(method=:representative)`
ranks candidate years on. That matters: these harnesses exist to check whether that ranking
reproduces the record, so measuring it with a re-implementation would let the two drift and
quietly invalidate the comparison.
"""
melt_rate(profile, forcing, mp_last) =
    GEMB._year_scores(profile, forcing, mp_last, GEMB.ThermalWorkspace())[1]

smb_rate(profile, forcing, mp_last) =
    GEMB._year_scores(profile, forcing, mp_last, GEMB.ThermalWorkspace())[2]

# Both reductions from one integration, for callers that want the pair.
melt_and_smb(profile, forcing, mp_last) =
    GEMB._year_scores(profile, forcing, mp_last, GEMB.ThermalWorkspace())

# Mass-weighted column mean density, via the same helper `gemb_spinup` converges on, so a
# harness cannot report a different mean density than the spinup was judged by.
mean_density(profile) = GEMB._column_mean_density(profile)

"""
    equilibrate(profile, cycle, mp; kwargs...) -> NamedTuple

Spin up on `cycle` and report the equilibrated column: `fac`, `density`, `temperature_deep`,
`cycles`, `years` and `converged`.

`years` is the field to compare methods on, not `cycles`: `:average` integrates one year per
cycle while `:representative` integrates `n_years`, so a longer cycle mechanically lowers the
cycle count without lowering the work.
"""
function equilibrate(profile, cycle, mp;
                     max_iterations=800,
                     convergence_delta_density=1e-3,
                     convergence_drift_density=1e-3,
                     convergence_delta_fac=1e-4,
                     convergence_drift_fac=1e-4)
    g = gemb_spinup(profile, cycle, mp;
        max_iterations=max_iterations,
        convergence_delta_density=convergence_delta_density,
        convergence_drift_density=convergence_drift_density,
        convergence_delta_fac=convergence_delta_fac,
        convergence_drift_fac=convergence_drift_fac,
        thermal_workspace=GEMB.ThermalWorkspace())
    md = DD.metadata(g)
    yr_per_cycle = length(DD.dims(cycle, DD.Ti)) * cycle.time_step / GEMB.SECONDS_PER_YEAR
    return (fac = GEMB._column_fac(g, mp),
            density = mean_density(g),
            temperature_deep = parent(g[:temperature])[end],
            cycles = md[:spinup_cycles],
            years = md[:spinup_cycles] * yr_per_cycle,
            converged = md[:spinup_converged])
end

# Equal-weight mean of whichever of the melt and accumulation errors are defined. Equal weight is
# a choice, not a derivation — an altimetry application would weight accumulation higher — so
# callers print the per-variable columns alongside it.
joint_error(melt_err, acc_err) =
    (parts = filter(isfinite, [melt_err, acc_err]); isempty(parts) ? NaN : mean(parts))

# Same, for a score NamedTuple carrying `melt_err` and `acc_err` fields.
joint_error_of(s) = joint_error(s.melt_err, s.acc_err)

# Run `f(id)` over every site concurrently, collecting results and reporting failures rather than
# aborting the sweep. Each site is an independent column, but a spinup is single-threaded.
function map_sites(f, ids=site_ids())
    results = Any[]
    lk = ReentrantLock()
    Threads.@threads for id in ids
        try
            r = f(id)
            lock(lk) do; push!(results, r); end
        catch err
            lock(lk) do
                @printf("  %-12s FAILED: %s\n", id, first(split(sprint(showerror, err), '\n')))
            end
        end
    end
    return results
end
