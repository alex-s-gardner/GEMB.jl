# Test a regime-aware cycle-selection rule against the fixed strategies, on real Greenland
# forcing.
#
#   julia --project=bench --threads=auto bench/greenland_regime_rule.jl
#
# The rule under test, from the measurements in `greenland_year_selection.jl`:
#
#   no melt              -> `:average`                     (mean climatology)
#   moderate melt        -> `:representative`, SMB-ranked, n_years=3
#   ablation (m/a > 1)   -> `:representative`, melt-ranked, n_years=1
#
# The rationale for each leg is that no single fixed strategy won everywhere. `:average`
# preserves accumulation to ~0.1% at dry sites and costs less than half the model-years, but
# destroys melt where melt lives in the tail of the distribution. SMB-ranking matches
# accumulation to ~2.5% but under-recovers melt where accumulation dominates SMB. Melt-ranking
# matches melt to ~102% but misses accumulation by up to 38%. If the regimes really do want
# different things, a composite should beat every fixed choice; if it does not, that is evidence
# the regime split is not the right lever.
#
# Spinup-free by design: melt and accumulation fidelity are properties of the forcing, and the
# equilibrated-column comparison was abandoned as unaffordable (see
# `greenland_block_length.jl`). So this measures how well each cycle reproduces the record's melt
# and accumulation, NOT how well the spun-up column matches reality.
#
# Reuses the forcing cache written by `greenland_spinup_forcing.jl`.

using GEMB
using GEMB_ClimateForcing
using Statistics
using Printf
using Serialization
using Dates

const DD = GEMB.DimensionalData

const YEARS = (2000, 2019)
const CACHE_DIR = joinpath(@__DIR__, "greenland_forcing_cache")

# Regime thresholds. Both are read off the measured distribution rather than fitted:
#
# * `MELT_FLOOR` separates "no melt" from "moderate": the dry-snow sites melt 0.0-1.2 kg m-2 yr-1
#   and the lowest site where averaging measurably lost melt is Saddle at 23.5, so 5 sits in the
#   gap without landing on any observed value.
# * `ABLATION_RATIO` is melt/accumulation = 1, which is the definition of the ablation zone
#   rather than a tuned number: above it the column cannot retain a firn layer year on year.
const MELT_FLOOR = 5.0
const ABLATION_RATIO = 1.0

const SITE_IDS = ["Summit", "NEEM", "EastGRIP", "DYE-2", "Saddle", "CP1",
                  "KAN-U", "SwissCamp", "KAN-M", "KAN-L", "QAS-L", "THU-L"]

function cached_forcing(id)
    path = joinpath(CACHE_DIR, "$(id)_$(YEARS[1])_$(YEARS[2]).jls")
    isfile(path) || error("no cached forcing for $id; run bench/greenland_spinup_forcing.jl first")
    return deserialize(path)
end

function melt_rate(profile, forcing, mp_last)
    out = gemb(profile, forcing, mp_last; thermal_workspace=GEMB.ThermalWorkspace())
    years = length(DD.dims(forcing, DD.Ti)) * forcing.time_step / GEMB.SECONDS_PER_YEAR
    return years > 0 ? sum(parent(out[:melt])) / years : NaN
end

"""
    regime(melt, ratio) -> Symbol

Which leg of the composite rule a site falls in: `:no_melt`, `:moderate` or `:ablation`.

Two variables rather than one because they are not redundant: absolute melt decides whether melt
matters at all (a site melting 1 kg m-2 yr-1 needs no melt fidelity regardless of its ratio),
while the ratio decides whether a firn column can persist.
"""
function regime(melt, ratio)
    melt < MELT_FLOOR && return :no_melt
    isfinite(ratio) && ratio > ABLATION_RATIO && return :ablation
    return :moderate
end

# The candidate strategies, each a function of (cf, mp) returning a cycle. Named to match the
# columns of `greenland_year_selection.jl` so the two runs can be read together.
const STRATEGIES = [
    "avg" => (cf, mp) -> forcing_climatology(cf),
    "m1"  => (cf, mp) -> forcing_climatology(cf; method=:representative, model_parameters=mp,
                                             n_years=1, rank_by=:model, verbose=false),
    "m3"  => (cf, mp) -> forcing_climatology(cf; method=:representative, model_parameters=mp,
                                             n_years=3, rank_by=:model, verbose=false),
    "s3"  => (cf, mp) -> forcing_climatology(cf; method=:representative, model_parameters=mp,
                                             n_years=3, rank_by=:smb, verbose=false),
    # "s3 median" read literally: SMB-ranked, 3-year, targeting the median rather than the mean of
    # the candidate years' SMB. Tested alongside the mean form because the request was ambiguous
    # and the distinction is cheap to measure.
    "s3med" => (cf, mp) -> forcing_climatology(cf; method=:representative, model_parameters=mp,
                                               n_years=3, rank_by=:smb, statistic=:median,
                                               verbose=false),
]

function analyze(id, mp)
    cf = initialize_forcing(cached_forcing(id))
    mp_last = GEMB.ModelParameters(;
        (f => getfield(mp, f) for f in fieldnames(GEMB.ModelParameters)
         if f != :output_frequency)..., output_frequency=:last)
    profile = initialize_profile(mp, cf)
    cs = GEMB.initialize_climate_summary(cf, mp)
    melt_record = melt_rate(profile, cf, mp_last)
    ratio = cs.accumulation > 0 ? melt_record / cs.accumulation : NaN

    scores = Dict{String,Any}()
    for (name, build) in STRATEGIES
        cycle = build(cf, mp)
        cs_c = GEMB.initialize_climate_summary(cycle, mp)
        m = melt_rate(profile, cycle, mp_last)
        scores[name] = (
            # Melt error as a percentage of the record. Undefined where the record barely melts,
            # since a ratio against ~0 is meaningless — those sites are scored on accumulation
            # alone, which is the only thing that can be wrong there.
            melt_err = melt_record > MELT_FLOOR ? 100 * abs(m / melt_record - 1) : NaN,
            melt_pct = melt_record > 1.0 ? 100 * m / melt_record : NaN,
            acc_err = cs.accumulation > 0 ? 100 * abs(cs_c.accumulation / cs.accumulation - 1) : NaN,
        )
    end

    return (id=id, melt_record=melt_record, accumulation=cs.accumulation, ratio=ratio,
            regime=regime(melt_record, ratio), scores=scores)
end

# Joint error of one strategy at one site: the mean of its melt and accumulation errors, or
# whichever is defined. Equal weight is a choice, not a derivation — an altimetry application
# would weight accumulation higher — so the per-variable columns are printed too.
function joint(s)
    parts = filter(isfinite, [s.melt_err, s.acc_err])
    return isempty(parts) ? NaN : mean(parts)
end

const MP = GEMB.ModelParameters(output_frequency=:last)
@printf("Regime-rule test, %d sites, threads=%d\n", length(SITE_IDS), Threads.nthreads())

results = Any[]
lk = ReentrantLock()
Threads.@threads for id in SITE_IDS
    try
        r = analyze(id, MP)
        lock(lk) do; push!(results, r); end
    catch err
        lock(lk) do
            @printf("  %-12s FAILED: %s\n", id, first(split(sprint(showerror, err), '\n')))
        end
    end
end
sort!(results, by = r -> r.melt_record)

const NAMES = first.(STRATEGIES)

# The composite: which strategy each regime uses. `:moderate` is tested with both SMB forms.
const RULE = Dict(:no_melt => "avg", :moderate => "s3", :ablation => "m1")
const RULE_MED = Dict(:no_melt => "avg", :moderate => "s3med", :ablation => "m1")

println("\n", "="^108)
println("REGIME ASSIGNMENT and what the rule picks")
println("="^108)
@printf("%-12s %9s %9s %7s | %-10s %-8s\n",
        "site", "melt_rec", "accum", "m/a", "regime", "rule->")
for r in results
    @printf("%-12s %9.1f %9.1f %7.2f | %-10s %-8s\n",
            r.id, r.melt_record, r.accumulation, r.ratio, r.regime, RULE[r.regime])
end

println("\n", "="^108)
println("JOINT ERROR per site: composite rule vs every fixed strategy (lower is better)")
println("="^108)
@printf("%-12s %-10s | %s | %8s %8s\n", "site", "regime",
        join((lpad(n, 7) for n in NAMES), " "), "RULE", "RULEmed")
for r in results
    js = [joint(r.scores[n]) for n in NAMES]
    jr = joint(r.scores[RULE[r.regime]])
    jm = joint(r.scores[RULE_MED[r.regime]])
    @printf("%-12s %-10s | %s | %7.1f%% %7.1f%%\n", r.id, r.regime,
            join((@sprintf("%6.1f%%", j) for j in js), " "), jr, jm)
end

println("\n", "="^108)
println("SUMMARY: is the composite better than any single fixed choice?")
println("="^108)
function summarize(label, per_site)
    vals = filter(isfinite, per_site)
    @printf("%-9s median %6.1f%%   mean %6.1f%%   worst %6.1f%%\n",
            label, median(vals), mean(vals), maximum(vals))
end
for n in NAMES
    summarize(n, [joint(r.scores[n]) for r in results])
end
summarize("RULE", [joint(r.scores[RULE[r.regime]]) for r in results])
summarize("RULEmed", [joint(r.scores[RULE_MED[r.regime]]) for r in results])

# An oracle that picks the best strategy per site is the floor no rule can beat. If the composite
# is close to it, the regime split captures most of what is available; if a fixed strategy is
# already close, the split is not buying much.
summarize("ORACLE", [minimum(filter(isfinite, [joint(r.scores[n]) for n in NAMES]))
                     for r in results])

println("\nPer-regime breakdown (median joint error):")
@printf("%-10s %3s | %s | %8s\n", "regime", "n", join((lpad(n, 7) for n in NAMES), " "), "RULE")
for reg in (:no_melt, :moderate, :ablation)
    grp = [r for r in results if r.regime === reg]
    isempty(grp) && continue
    meds = [median(filter(isfinite, [joint(r.scores[n]) for r in grp])) for n in NAMES]
    ruled = median(filter(isfinite, [joint(r.scores[RULE[reg]]) for r in grp]))
    @printf("%-10s %3d | %s | %7.1f%%\n", reg, length(grp),
            join((@sprintf("%6.1f%%", m) for m in meds), " "), ruled)
end

# Does the rule pick the best available strategy at each site? This is the direct test of the
# regime hypothesis: the rule is only justified if its choice is at or near the per-site optimum.
println("\nAgreement with the per-site optimum:")
for (label, rule) in (("RULE", RULE), ("RULEmed", RULE_MED))
    agree = 0
    for r in results
        js = Dict(n => joint(r.scores[n]) for n in NAMES)
        best = argmin(n -> isfinite(js[n]) ? js[n] : Inf, NAMES)
        chosen = rule[r.regime]
        # "Agrees" if the rule's pick is the optimum, or within 2 percentage points of it — a
        # tighter test would be measuring noise in a 12-site sample.
        isfinite(js[chosen]) && js[chosen] <= js[best] + 2.0 && (agree += 1)
    end
    @printf("  %-8s picks the best (or within 2 points) at %d of %d sites\n",
            label, agree, length(results))
end
