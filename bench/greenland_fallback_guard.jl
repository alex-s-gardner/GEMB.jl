# Does an SMB-bias fallback rescue the two-regime rule?
#
#   julia --project=bench --threads=auto bench/greenland_fallback_guard.jl
#
# The two-regime rule that the measurements support is:
#
#     melt < MELT_FLOOR  ->  `:average`
#     otherwise          ->  `:representative`, melt-ranked, n_years=3   ("m3")
#
# The worry this script tests: m3 ranks on melt alone, so nothing stops the selected 3-year block
# from being unrepresentative in *accumulation* — and accumulation drives the firn column
# directly. The proposed guard is a cheap post-selection sanity check: if the chosen block's mass
# balance is badly biased against the record, re-select with SMB ranking ("s3"), which matched
# accumulation to ~2.5% median across these sites.
#
# The guard is cheap because both quantities are already computed during selection: the per-year
# scores are in hand, so the check is arithmetic on numbers the selection loop already has. No
# extra `gemb` run, no spinup.
#
# What is measured, per candidate trigger variable and threshold: how often the guard fires, and
# whether firing improves or degrades accumulation error, melt error, and their joint mean.
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
const MELT_FLOOR = 5.0          # kg m-2 yr-1; below this the rule uses `:average`

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

# SMB [m ice per year] over `forcing`, through the same `_smb_rate` the spinup and the `:smb`
# ranking use, so all three agree on the definition.
function smb_rate(profile, forcing, mp_last)
    out = gemb(profile, forcing, mp_last; thermal_workspace=GEMB.ThermalWorkspace())
    years = length(DD.dims(forcing, DD.Ti)) * forcing.time_step / GEMB.SECONDS_PER_YEAR
    return GEMB._smb_rate(out, years, mp_last.density_ice)
end

function analyze(id, mp)
    cf = initialize_forcing(cached_forcing(id))
    mp_last = GEMB.ModelParameters(;
        (f => getfield(mp, f) for f in fieldnames(GEMB.ModelParameters)
         if f != :output_frequency)..., output_frequency=:last)
    profile = initialize_profile(mp, cf)
    cs = GEMB.initialize_climate_summary(cf, mp)

    melt_rec = melt_rate(profile, cf, mp_last)
    smb_rec = smb_rate(profile, cf, mp_last)

    build(rank_by) = forcing_climatology(cf; method=:representative, model_parameters=mp,
                                        n_years=3, rank_by=rank_by, verbose=false)
    function score(cycle)
        cs_c = GEMB.initialize_climate_summary(cycle, mp)
        m = melt_rate(profile, cycle, mp_last)
        s = smb_rate(profile, cycle, mp_last)
        return (
            melt_err = melt_rec > MELT_FLOOR ? 100 * abs(m / melt_rec - 1) : NaN,
            acc_err = cs.accumulation > 0 ? 100 * abs(cs_c.accumulation / cs.accumulation - 1) : NaN,
            # Signed SMB bias, as a fraction of the record's SMB magnitude. The magnitude in the
            # denominator, not the value: SMB changes sign between accumulating and ablating
            # sites, and a relative error against a near-zero or negative SMB is not
            # interpretable otherwise.
            smb_bias = abs(smb_rec) > 1e-6 ? 100 * (s - smb_rec) / abs(smb_rec) : NaN,
        )
    end

    m3, s3 = build(:model), build(:smb)
    return (id=id, melt_rec=melt_rec, accumulation=cs.accumulation, smb_rec=smb_rec,
            ratio=cs.accumulation > 0 ? melt_rec / cs.accumulation : NaN,
            uses_average = melt_rec < MELT_FLOOR,
            m3=score(m3), s3=score(s3))
end

joint(s) = (p = filter(isfinite, [s.melt_err, s.acc_err]); isempty(p) ? NaN : mean(p))

const MP = GEMB.ModelParameters(output_frequency=:last)
@printf("SMB-bias fallback test, %d sites, threads=%d\n", length(SITE_IDS), Threads.nthreads())

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
sort!(results, by = r -> r.melt_rec)

# Only sites the rule sends to m3 can be rescued; the `:average` leg is already optimal there.
candidates = [r for r in results if !r.uses_average]

println("\n", "="^112)
println("The m3 block's bias per site, and what s3 would give instead")
println("="^112)
@printf("%-12s %8s %8s | %9s %9s | %8s %8s | %8s %8s\n",
        "site", "melt_rec", "smb_rec", "m3_smbbias", "s3_smbbias",
        "m3_acc", "s3_acc", "m3_joint", "s3_joint")
for r in candidates
    @printf("%-12s %8.1f %8.3f | %+8.1f%% %+8.1f%% | %7.1f%% %7.1f%% | %7.1f%% %7.1f%%\n",
            r.id, r.melt_rec, r.smb_rec, r.m3.smb_bias, r.s3.smb_bias,
            r.m3.acc_err, r.s3.acc_err, joint(r.m3), joint(r.s3))
end

println("\n", "="^112)
println("Would a bias-triggered fallback help? Sweeping the trigger variable and threshold.")
println("(`fires` = sites where the guard would switch m3 -> s3.)")
println("="^112)

# Two candidate trigger variables. `smb` is the user's proposal; `acc` is the direct measure of
# the thing being worried about, included because SMB and accumulation are not interchangeable at
# a melting site — SMB folds melt back in, so a large SMB bias need not mean a large accumulation
# bias, and vice versa.
const TRIGGERS = ["smb_bias" => (r -> abs(r.m3.smb_bias)),
                  "acc_err"  => (r -> r.m3.acc_err)]

@printf("%-9s %6s %6s | %-28s | %-28s\n", "trigger", "thresh", "fires",
        "median joint (m3 -> guarded)", "median acc_err (m3 -> guarded)")
for (tname, tfun) in TRIGGERS, thresh in (5.0, 10.0, 15.0, 20.0, 30.0)
    fires = [r for r in candidates if isfinite(tfun(r)) && tfun(r) > thresh]
    guarded_joint = [(isfinite(tfun(r)) && tfun(r) > thresh) ? joint(r.s3) : joint(r.m3)
                     for r in candidates]
    guarded_acc = [(isfinite(tfun(r)) && tfun(r) > thresh) ? r.s3.acc_err : r.m3.acc_err
                   for r in candidates]
    base_joint = [joint(r.m3) for r in candidates]
    base_acc = [r.m3.acc_err for r in candidates]
    f(v) = median(filter(isfinite, v))
    @printf("%-9s %5.0f%% %6d | %10.1f%% -> %10.1f%%     | %10.1f%% -> %10.1f%%\n",
            tname, thresh, length(fires), f(base_joint), f(guarded_joint),
            f(base_acc), f(guarded_acc))
end

println("\nPer-site effect of firing (acc_err and joint, m3 -> s3):")
for r in candidates
    Δacc = r.s3.acc_err - r.m3.acc_err
    Δj = joint(r.s3) - joint(r.m3)
    @printf("  %-12s smb_bias %+7.1f%%  acc %5.1f%% -> %5.1f%% (%+5.1f)  joint %5.1f%% -> %5.1f%% (%+5.1f) %s\n",
            r.id, r.m3.smb_bias, r.m3.acc_err, r.s3.acc_err, Δacc,
            joint(r.m3), joint(r.s3), Δj, Δj < 0 ? "BETTER" : "worse")
end
