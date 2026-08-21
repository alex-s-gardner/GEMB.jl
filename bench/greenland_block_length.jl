# Does a single average-melt year beat a multi-year block, on real Greenland forcing?
#
# !! ABANDONED — DOES NOT COMPLETE. !!
# This runs a full strict-criteria spinup per site per variant (12 sites x 5 variants), and was
# killed after more than an hour and 4000+ model-years without finishing: a 5-year block at the
# 800-cycle cap is 4000 model-years for one site-variant, and 6 of the 12 sites cannot satisfy
# these criteria at all. `greenland_year_selection.jl` supersedes it by answering the same
# question spinup-free, in seconds per site. Kept only as the record of what the design cost;
# do not run it expecting output.
#
#   julia --project=bench --threads=auto bench/greenland_block_length.jl
#
# `greenland_spinup_forcing.jl` compared `:average` against `:representative` with `n_years=3`
# and found the block's weak point: it selects on melt alone, so its *accumulation* can be far
# off the record (-22.1% at EastGRIP, +15.3% at NEEM), and at no-melt sites that mismatch is the
# only thing distinguishing the two methods. Accumulation drives the firn column directly, so a
# cycle that gets melt right and accumulation wrong is not obviously better than one that gets
# melt wrong and accumulation right.
#
# This script isolates block length. `n_years=1` selects the single year whose melt is closest
# to the record mean; `n_years=3` and `5` select consecutive blocks whose mean melt is. Averaging
# `n` years of accumulation should pull the block's accumulation toward the record mean as
# ~1/sqrt(n), so the naive expectation is that longer blocks match accumulation better — but the
# 3-year block's -22.1% at EastGRIP suggests the sample is small enough that this does not
# reliably hold, which is what the accumulation column here tests.
#
# Reuses the forcing cache written by `greenland_spinup_forcing.jl`, so it needs no CDS key.

include("greenland_common.jl")
const BLOCK_LENGTHS = (1, 3, 5)

# Melt below which the hybrid `auto` variant uses `:average`. 5 kg m-2 yr-1 is well under the
# lowest site where averaging measurably lost melt (Saddle, 23.5) and well above the dry-snow
# sites (0.0-1.2), so it separates them without sitting on top of any observed value.
const AUTO_MELT_THRESHOLD = 5.0

# Same criteria as `greenland_spinup_forcing.jl`, so the numbers are directly comparable to that
# run's `:average` and 3-year columns. Six of twelve sites failed to converge within 800 cycles
# there, so `conv` is reported per row and a `NO` marks a comparison that is not an equilibrium.
const MAX_CYCLES = 800
const CONV_DELTA_DENSITY = 1e-3
const CONV_DRIFT_DENSITY = 1e-3
const CONV_DELTA_FAC = 1e-4
const CONV_DRIFT_FAC = 1e-4

function analyze(id, mp)
    cf = initialize_forcing(cached_forcing(id))
    mp_last = last_step_parameters(mp)
    profile = initialize_profile(mp, cf)
    cs = GEMB.initialize_climate_summary(cf, mp)
    melt_record = melt_rate(profile, cf, mp_last)

    variants = Dict{String,Any}()
    avg = forcing_climatology(cf)
    variants["avg"] = (cycle=avg, summary=GEMB.initialize_climate_summary(avg, mp))
    for n in BLOCK_LENGTHS
        c = forcing_climatology(cf; method=:representative, model_parameters=mp,
                                n_years=n, verbose=false)
        variants["b$n"] = (cycle=c, summary=GEMB.initialize_climate_summary(c, mp))
    end

    # The hybrid: `:representative` only where melt is significant, `:average` otherwise. This
    # is the rule the earlier runs point to — averaging preserves accumulation to +0.1% and is
    # cheaper, so it should be used wherever melt is not the thing being got right, and
    # `:representative` should be reserved for where averaging destroys the melt.
    #
    # The switch is on absolute melt, not on melt/accumulation: averaging's failure tracked
    # melt magnitude and elevation (0% retention at 23-152 kg m-2 yr-1, 86-115% at 777-3222),
    # not the ratio, so the ratio is the wrong discriminator for this choice.
    auto_is_rep = melt_record >= AUTO_MELT_THRESHOLD
    variants["auto"] = auto_is_rep ? variants["b1"] : variants["avg"]

    out = Dict{String,Any}()
    for (k, v) in variants
        out[k] = (
            melt = melt_rate(profile, v.cycle, mp_last),
            dacc = cs.accumulation > 0 ? v.summary.accumulation / cs.accumulation - 1 : NaN,
            eq = equilibrate(profile, v.cycle, mp),
            years_sel = get(DD.metadata(v.cycle), :climatology_representative_year, nothing),
        )
    end
    return (id=id, accumulation=cs.accumulation, melt_record=melt_record,
            ratio=cs.accumulation > 0 ? melt_record / cs.accumulation : NaN, v=out)
end

const MP = GEMB.ModelParameters(output_frequency=:last)
@printf("Block-length comparison, %d Greenland sites, threads=%d\n\n",
        length(site_ids()), Threads.nthreads())

results = Any[]
lk = ReentrantLock()
Threads.@threads for id in site_ids()
    try
        r = analyze(id, MP)
        lock(lk) do
            push!(results, r)
            @printf("  done %-12s (ratio %.2f)\n", id, r.ratio)
        end
    catch err
        lock(lk) do
            @printf("  done %-12s FAILED: %s\n", id,
                    first(split(sprint(showerror, err), '\n')))
        end
    end
end
sort!(results, by = r -> r.melt_record)

const KEYS = ["avg", ("b$n" for n in BLOCK_LENGTHS)..., "auto"]

println("\n", "="^112)
println("MELT recovered, as % of the real record (>100% overshoots)")
println("="^112)
@printf("%-12s %9s |%8s %8s %8s %8s %8s | selected years (n=1 / n=3 / n=5)\n",
        "site", "melt_rec", "avg", "n=1", "n=3", "n=5", "auto")
for r in results
    pct(k) = r.melt_record > 1.0 ? 100 * r.v[k].melt / r.melt_record : NaN
    @printf("%-12s %9.1f |%7.0f%% %7.0f%% %7.0f%% %7.0f%% %7.0f%% | %s / %s / %s\n",
            r.id, r.melt_record, pct("avg"), pct("b1"), pct("b3"), pct("b5"), pct("auto"),
            r.v["b1"].years_sel, r.v["b3"].years_sel, r.v["b5"].years_sel)
end

println("\n", "="^112)
println("ACCUMULATION offset of the cycle from the record — the block's weak point")
println("="^112)
@printf("%-12s %9s |%8s %8s %8s %8s %8s\n", "site", "accum",
        "avg", "n=1", "n=3", "n=5", "auto")
for r in results
    @printf("%-12s %9.1f |%+7.1f%% %+7.1f%% %+7.1f%% %+7.1f%% %+7.1f%%\n", r.id, r.accumulation,
            (100 * r.v[k].dacc for k in KEYS)...)
end

println("\n", "="^112)
println("EQUILIBRATED firn air content [m], and whether the spinup converged")
println("="^112)
@printf("%-12s %6s |%8s %8s %8s %8s %8s | %s\n", "site", "m/a",
        "avg", "n=1", "n=3", "n=5", "auto", "converged? avg/1/3/5/auto")
for r in results
    @printf("%-12s %6.2f |%8.3f %8.3f %8.3f %8.3f %8.3f | %s\n", r.id, r.ratio,
            (r.v[k].eq.fac for k in KEYS)...,
            join((r.v[k].eq.converged ? "y" : "N" for k in KEYS), "/"))
end

println("\n", "="^112)
println("SUMMARY")
println("="^112)
melters = [r for r in results if r.melt_record > 1.0]
for k in KEYS
    rets = sort([100 * r.v[k].melt / r.melt_record for r in melters])
    daccs = [abs(100 * r.v[k].dacc) for r in results if isfinite(r.v[k].dacc)]
    yrs = sum(r.v[k].eq.years for r in results)
    nconv = count(r -> r.v[k].eq.converged, results)
    @printf("%-4s melt: median %6.0f%%  <10%%: %d/%d  |  |dacc|: median %5.1f%%  max %5.1f%%  |  " *
            "model-years %6.0f  converged %2d/%2d\n",
            k, median(rets), count(<(10), rets), length(rets),
            median(daccs), maximum(daccs), yrs, nconv, length(results))
end
