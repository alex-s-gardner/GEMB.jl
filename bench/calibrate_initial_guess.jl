# Calibration harness for the steady-state initial guess (`initialize_profile`).
#
#   julia --project=bench --threads=auto bench/calibrate_initial_guess.jl truth
#   julia --project=bench --threads=auto bench/calibrate_initial_guess.jl report
#   julia --project=bench --threads=auto bench/calibrate_initial_guess.jl fit
#
# `initialize_profile` exists to cut the number of `gemb_spinup` cycles needed to reach
# density convergence. This measures whether it does, over a fleet of synthetic sites, and
# fits the parts of the guess that are free.
#
# ## Why the ground truth is not simply "the converged profile"
#
# The bottom cell is a Dirichlet reservoir — `calculate_temperature.jl` returns
# `temperature[m]` bit-unchanged by construction — so the *initialized* deep temperature is a
# permanent boundary condition that no amount of spinup forgets. Spinup therefore has no
# unique temperature attractor to converge to: it inherits whatever the initializer chose,
# and the converged column depends on it (measured at `test_1`: initializing the deep column
# at 247 K vs 256.4 K changes converged mass-weighted mean density by 11 kg m-3).
#
# So temperature ground truth is defined by *self-consistency* instead: the deep temperature
# `T_deep` whose converged column has no jump across the frozen bottom cell,
# `T[end] - T[end-1] = 0`. That residual is smooth and monotone in `T_deep`, so a secant
# root-find on it is well posed. This is the expensive part of the harness, and the reason
# the truth stage is cached.
#
# ## Why candidates are scored against the cached truth, not by re-spinning
#
# The deep temperature a candidate formula predicts is a pure function of the climate
# summary, so once the per-site self-consistent `T_deep` is known, any number of candidate
# formulas can be scored instantly against it. Full spinups are spent only twice: building
# the truth, and validating the winner. Holdout is by construction axis rather than at
# random, so it tests extrapolation rather than interpolation between neighbours.

using GEMB
using GEMB_ClimateForcing
using Statistics
using Printf
using Serialization
using Dates

const DD = GEMB.DimensionalData

# Convergence criteria for every spinup here: both the step and the trend test, so a column
# creeping steadily under the step tolerance does not pass. `max_iterations` is a backstop —
# any site that reaches it is dropped rather than scored, since a capped run reports the cap
# and not a cycle count.
const CONV_DELTA = 1e-3
const CONV_DRIFT = 1e-3
const MAX_CYCLES = 800

const TIME_STEP_HOURS = 3
const CACHE_FILE = joinpath(@__DIR__, "calibrate_initial_guess_truth.jls")

#=============================================================================
# The fleet
=============================================================================#

"""
    fleet() -> Vector{NamedTuple}

Site definitions: an elevation offset applied to the `test_1` synthetic site, times a
precipitation scale, times a wind scale.

Elevation is the physically coherent knob — `climate_adjust_for_elevation` moves temperature,
pressure, vapor pressure and longwave together under a lapse rate — so the ladder spans
heavy-melt through dry-cold as one family rather than as independent perturbations of
single variables. The precipitation and wind ladders are then applied on top to separate the
accumulation and turbulent-transfer dependence from the temperature dependence.

`split` is assigned **by construction axis, not at random**: whole elevation rungs are held
out, so a holdout site differs from every training site in its temperature climate rather than
sitting between two of them. A random split over a smooth grid would only ever test
interpolation between neighbours, which a badly overfitted formula passes. The two extreme
rungs (`-1500`, warmest/heavy-melt, and `+1200`, coldest/driest) are deliberately *among* the
held-out ones, so the holdout tests extrapolation past both ends of the training range.

Fleet size is set by measurement, not by taste. `probe` mode reports one spinup per corner: the
cost per site spans 8 s to 46 s (119 to 730 cycles), because a cold low-accumulation site builds
a 250 m column *and* starts far from its attractor. With the root-find spending up to
`T_DEEP_MAX_EVALS` spinups on top of the baseline one, a site costs up to ~5 minutes, so the
grid is kept to 36 sites — ~15 minutes on 8 threads. A denser grid would only add sites
between rungs that are already smooth, which is not what the holdout is testing.
"""
function fleet()
    # Relative to test_1's own 2200 m. Negative is warmer/lower.
    elevations = [-1500.0, -900.0, -300.0, 0.0, 600.0, 1200.0]
    precips    = [0.3, 1.0, 3.0]
    winds      = [0.6, 1.4]

    # Held out as whole temperature climates, including both ends of the ladder.
    holdout_elevations = Set([-1500.0, -300.0, 1200.0])

    sites = NamedTuple[]
    for dz in elevations, ps in precips, ws in winds
        holdout = dz in holdout_elevations
        push!(sites, (
            id = @sprintf("z%+05d_p%03d_w%03d", round(Int, dz), round(Int, 100ps), round(Int, 100ws)),
            delta_elevation = dz,
            precip_scale = ps,
            wind_scale = ws,
            split = holdout ? :holdout : :train,
        ))
    end
    return sites
end

"""
    floor_longwave(stack; floor=LONGWAVE_FLOOR) -> DimStack

`stack` with its `longwave_downward` layer floored.

The `test_1` synthetic series reaches 0 W m-2 on 0.97% of its steps, which is an artifact of
its additive noise model rather than a sky: 0 W m-2 is a 0 K atmosphere. It is harmless in the
transient run, but `climate_adjust_for_elevation` re-derives longwave from the adjusted
temperature and vapor pressure while *preserving the cloud increment* — so a 0 W m-2 step
carries a large negative increment through the adjustment and the adjuster's own unit
validation (lower bound 50 W m-2) rejects the whole series.

The floor is 80 W m-2 rather than the validator's own 50, because the *cooling* rungs of the
elevation ladder push longwave down further: at +1200 m a 50 W m-2 step becomes 35, so
flooring at the validator's bound only makes the warming half of the fleet legal. 80 W m-2 is
a ~194 K effective sky, and the ladder's worst cooling factor (~0.70) leaves it at 56 —
above the bound with margin. Only the noise tail is touched; the 1st percentile of the base
series is already above it.
"""
const LONGWAVE_FLOOR = 80.0

function floor_longwave(stack; floor::Real=LONGWAVE_FLOOR)
    # Clamp in place on a copy of the layer's own backing array. Rebuilding the stack with a
    # broadcast result instead would nest a `DimArray` inside a layer slot, which silently
    # yields mixed `OffsetVector`/`DimVector` layers after the adjuster round-trips it.
    out = deepcopy(stack)
    lw = parent(out[:longwave_downward])
    @. lw = max(lw, floor)
    return out
end

"""
    build_forcing(base_stack, site) -> ClimateForcing

Realize one fleet site as a `ClimateForcing`.

The elevation adjustment is delegated to `GEMB_ClimateForcing.climate_adjust_for_elevation`,
but the resulting stack is **not** handed to `initialize_forcing(::DimStack)`. That path would
read the adjuster's `precipitation_mean` metadata, which is written as `mean(precip)*8760` —
correct only for an hourly series, and 3x low for the 3-hourly stack used here. The layers are
extracted and passed to the vector method of `initialize_forcing` with the means computed from
the series itself, which is exact at any time step.
"""
function build_forcing(base_stack, site)
    # `base_stack` is expected already floored (see `floor_longwave`) — callers floor once and
    # reuse, rather than paying a deepcopy of the whole series per site.
    adjusted = GEMB_ClimateForcing.climate_adjust_for_elevation(
        base_stack, site.delta_elevation)

    time = collect(DD.lookup(DD.dims(adjusted, DD.Ti)))
    # `vec(Array(...))` rather than `parent(...)`: the adjuster's `_rebuild_forcing` can return
    # layers backed by an `OffsetVector`, whose `parent` is 0-indexed. Every consumer here
    # (and `initialize_forcing`'s own `@assert length(...) == n`) assumes 1-based, so
    # normalize once at the boundary instead of carrying the offset inward.
    layer(k) = vec(Array(adjusted[k]))
    T = layer(:temperature_air)
    P = layer(:pressure_air)
    e = layer(:vapor_pressure)
    LW = layer(:longwave_downward)
    SW = layer(:shortwave_downward)
    precip = layer(:precipitation) .* site.precip_scale
    wind = layer(:wind_speed) .* site.wind_scale

    # Annual rate from the series and its own step length, rather than an assumed one.
    dt_seconds = Dates.value(time[2] - time[1]) / 1000.0
    steps_per_year = GEMB.SECONDS_PER_YEAR / dt_seconds

    return initialize_forcing(time, T, P, precip, wind, SW, LW, e;
        temperature_air_mean = mean(T),
        wind_speed_mean = mean(wind),
        precipitation_mean = mean(precip) * steps_per_year,
        temperature_observation_height = 2.0,
        wind_observation_height = 10.0)
end

#=============================================================================
# Ground truth
=============================================================================#

# Mass-weighted mean density of a column, the quantity `gemb_spinup` converges on.
mean_density(p) = sum(parent(p[:density]) .* parent(p[:dz])) / sum(parent(p[:dz]))

spinup(p, cfc, mp; ws) = gemb_spinup(p, cfc, mp;
    max_iterations = MAX_CYCLES,
    convergence_delta_density = CONV_DELTA,
    convergence_drift_density = CONV_DRIFT,
    thermal_workspace = ws)

cycles(g) = DD.metadata(g)[:spinup_cycles]
converged(g) = DD.metadata(g)[:spinup_converged]

"""
    shift_temperature(profile, T_deep) -> profile

Copy of `profile` with its temperature profile rigidly shifted so the bottom cell sits at
`T_deep`, clamped to the melt point as every initialization path is.

A rigid shift rather than a rebuild: the shape of the guess (the damped annual wave) is not
what is being probed here, only the level the deep column is pinned at.
"""
function shift_temperature(profile, T_deep)
    p = deepcopy(profile)
    T = parent(p[:temperature])
    T .= min.(T .+ (T_deep - T[end]), GEMB.CtoK)
    return p
end

"""
    self_consistent_deep_temperature(profile, cfc, mp; ws) -> (T_deep, jump, evals)

The deep temperature at which the frozen bottom cell stops being a discontinuity: the root of
`r(T_deep) = T[end] - T[end-1]` after spinup to convergence.

Secant iteration, bracketed to the physically admissible range. `r` is monotone increasing in
`T_deep` (a warmer pinned base leaves a larger positive jump against a column set by the
surface climate), so the iteration is well behaved; returns as soon as the jump is inside
`T_JUMP_TOLERANCE`, which is the resolution the deep column can be said to have.
"""
const T_JUMP_TOLERANCE = 0.02   # [K] — below the ~0.13 K discretization floor of the column
const T_DEEP_MAX_EVALS = 6

function self_consistent_deep_temperature(profile, cfc, mp; ws)
    evals = 0
    function residual(T_deep)
        evals += 1
        g = spinup(shift_temperature(profile, T_deep), cfc, mp; ws=ws)
        converged(g) || return (NaN, g)
        T = parent(g[:temperature])
        return (T[end] - T[end - 1], g)
    end

    # Start from the guess's own deep temperature and one point 6 K below it — the sign of the
    # bias measured at `test_1`, so the bracket usually contains the root immediately.
    T1 = Float64(parent(profile[:temperature])[end])
    r1, _ = residual(T1)
    isnan(r1) && return (NaN, NaN, evals)
    abs(r1) < T_JUMP_TOLERANCE && return (T1, r1, evals)

    T0 = T1 - 6.0
    r0, _ = residual(T0)
    isnan(r0) && return (NaN, NaN, evals)

    while evals < T_DEEP_MAX_EVALS
        abs(r1) < T_JUMP_TOLERANCE && break
        denom = r1 - r0
        abs(denom) < 1e-12 && break
        T2 = clamp(T1 - r1 * (T1 - T0) / denom, 180.0, GEMB.CtoK)
        T0, r0 = T1, r1
        T1 = T2
        r1, _ = residual(T1)
        isnan(r1) && return (NaN, NaN, evals)
    end
    return (T1, r1, evals)
end

"""
    build_truth(mp; sites=fleet()) -> Dict{String,NamedTuple}

Per-site ground truth. For each site: the baseline cycles-to-convergence from the shipped
guess, the converged profile it reaches, and the self-consistent deep temperature.

Sites are independent columns, so they run concurrently. `mp` is a pure description and is
explicitly safe to share across threads (see `ModelParameters` in `src/types.jl`); the mutable
scratch is the `ThermalWorkspace`, so each thread gets its own.

A site whose spinup does not converge within `MAX_CYCLES` is recorded as such and excluded
from every score, rather than being credited with the cap.
"""
function build_truth(mp; sites=fleet())
    base = floor_longwave(simulate_climate_forcing("test_1", TIME_STEP_HOURS))
    results = Vector{Any}(undef, length(sites))
    # Progress is printed as sites land, not at the end: a site can cost minutes (measured up
    # to 730 spinup cycles), so a run with no output is indistinguishable from a hung one.
    done = Threads.Atomic{Int}(0)
    log_lock = ReentrantLock()

    Threads.@threads for i in eachindex(sites)
        site = sites[i]
        # One workspace per iteration rather than one per thread indexed by `threadid()`:
        # tasks may migrate between threads at any yield point, so a `threadid()`-indexed
        # pool can hand the same buffers to two concurrently-running columns. An allocation
        # per site is free against the hundreds of spinup cycles it serves.
        ws = GEMB.ThermalWorkspace()
        try
            cf = build_forcing(base, site)
            cfc = forcing_climatology(cf)
            cs = GEMB.initialize_climate_summary(cf, mp)
            guess = initialize_profile(mp, cf)

            g = spinup(guess, cfc, mp; ws=ws)
            ok = converged(g)
            T_deep, jump, evals = ok ?
                self_consistent_deep_temperature(guess, cfc, mp; ws=ws) : (NaN, NaN, 0)

            results[i] = (
                id = site.id, split = site.split,
                delta_elevation = site.delta_elevation,
                precip_scale = site.precip_scale, wind_scale = site.wind_scale,
                usable = ok && isfinite(T_deep),
                # Baseline: what the shipped guess costs.
                baseline_cycles = ok ? cycles(g) : MAX_CYCLES,
                baseline_deep_temperature = Float64(parent(guess[:temperature])[end]),
                baseline_jump = ok ? parent(g[:temperature])[end] - parent(g[:temperature])[end-1] : NaN,
                # Truth.
                deep_temperature = T_deep,
                deep_temperature_jump = jump,
                converged_mean_density = ok ? mean_density(g) : NaN,
                converged_density = ok ? copy(parent(g[:density])) : Float64[],
                converged_temperature = ok ? copy(parent(g[:temperature])) : Float64[],
                converged_dz = ok ? copy(parent(g[:dz])) : Float64[],
                # The climate scalars any candidate formula may read.
                summary = summary_scalars(cs),
                root_find_evals = evals,
            )
        catch err
            results[i] = (id = site.id, split = site.split,
                delta_elevation = site.delta_elevation,
                precip_scale = site.precip_scale, wind_scale = site.wind_scale,
                usable = false, error = sprint(showerror, err))
        end

        n = Threads.atomic_add!(done, 1) + 1
        r = results[i]
        lock(log_lock) do
            @printf("  [%2d/%2d] %-22s %s\n", n, length(sites), site.id,
                r.usable ? @sprintf("cycles=%d  T_deep=%.2f (guess %.2f)  %d spinups",
                                    r.baseline_cycles, r.deep_temperature,
                                    r.baseline_deep_temperature, r.root_find_evals + 1) :
                           "UNUSABLE")
            flush(stdout)
        end
    end

    return Dict(r.id => r for r in results)
end

load_truth() = deserialize(CACHE_FILE)
save_truth(t) = serialize(CACHE_FILE, t)

"""
    summary_scalars(cs) -> NamedTuple

The `ClimateSummary` reduced to a plain `NamedTuple`, which is what the cache stores.

Deliberately not the struct itself: `Serialization` records a struct's field layout, so a cache
written before a field is added to `ClimateSummary` cannot be read after — and adding
`temperature_surface_mean` is exactly what this harness exists to evaluate. A NamedTuple of
scalars survives that, and a candidate formula reading a field the cache predates gets a clean
`missing`-key error rather than a corrupt value.
"""
function summary_scalars(cs)
    return NamedTuple{fieldnames(typeof(cs))}(
        ntuple(i -> getfield(cs, i), fieldcount(typeof(cs))))
end

#=============================================================================
# Candidate deep-temperature formulas
=============================================================================#

"""
    candidates() -> Vector{Pair{String,Function}}

Candidate deep-temperature formulas, each a function of a `ClimateSummary`.

The search is a coarse sweep over *forms and scales*, not a descent from the shipped value,
so it is not anchored to today's defaults. Every candidate is a quantity the summary already
carries, so none of them requires a new pass over the forcing.
"""
function candidates()
    c = Pair{String,Function}[]
    push!(c, "shipped: T_air + latent" => cs -> cs.temperature_air_mean + cs.latent_warming)
    push!(c, "T_air" => cs -> cs.temperature_air_mean)
    # `temperature_surface_mean` does not exist in a cache built before it was added to
    # `ClimateSummary`, so these score as NaN rather than erroring — which is the honest
    # report for "this candidate cannot be evaluated against this cache".
    push!(c, "T_surface" => cs -> get(cs, :temperature_surface_mean, NaN))
    for f in (0.0, 0.25, 0.5, 0.75, 1.0)
        push!(c, @sprintf("T_surface + %.2f*latent", f) =>
            cs -> get(cs, :temperature_surface_mean, NaN) + f * cs.latent_warming)
    end
    return c
end

"""
    score(truth, predict; split) -> NamedTuple

Error statistics of a candidate deep-temperature formula against the self-consistent truth,
over the sites in `split`.

Reported as mean absolute error **and** p90, because a fit that improves the mean by getting
the easy majority right while leaving the hard tail worse is exactly the overfit this is meant
to catch. `bias` is signed, since a systematic offset is a different defect from scatter.
"""
function score(truth, predict; split)
    errs = Float64[]
    for r in values(truth)
        r.usable || continue
        (split === :all || r.split === split) || continue
        push!(errs, predict(r.summary) - r.deep_temperature)
    end
    isempty(errs) && return (n = 0, mae = NaN, p90 = NaN, bias = NaN, worst = NaN)
    a = abs.(errs)
    return (n = length(errs), mae = mean(a), p90 = quantile(a, 0.9),
            bias = mean(errs), worst = maximum(a))
end

#=============================================================================
# Reporting
=============================================================================#

function usable(truth)
    u = [r for r in values(truth) if r.usable]
    sort!(u, by = r -> (r.delta_elevation, r.precip_scale, r.wind_scale))
    return u
end

function report_truth(truth)
    all_sites = collect(values(truth))
    u = usable(truth)
    dropped = [r for r in all_sites if !r.usable]

    @printf("Fleet: %d sites, %d usable, %d dropped (train %d / holdout %d)\n",
        length(all_sites), length(u), length(dropped),
        count(r -> r.split === :train, u), count(r -> r.split === :holdout, u))
    for r in dropped
        # Truncated: a `MethodError` on a `DimArray` prints several hundred characters of
        # type parameters, which buries the rest of the report.
        msg = haskey(r, :error) ? first(split(r.error, '\n')) :
              "did not converge in $MAX_CYCLES cycles"
        @printf("  DROPPED %-22s %s\n", r.id, length(msg) > 110 ? msg[1:110] * "..." : msg)
    end

    @printf("\nBaseline cycles-to-convergence from the shipped guess: mean %.1f  p90 %.1f  max %d\n",
        mean(r.baseline_cycles for r in u), quantile([Float64(r.baseline_cycles) for r in u], 0.9),
        maximum(r.baseline_cycles for r in u))
    de = [r.baseline_deep_temperature - r.deep_temperature for r in u]
    @printf("Shipped deep-temperature error [K]: bias %+.2f  mae %.2f  worst %.2f\n",
        mean(de), mean(abs.(de)), maximum(abs.(de)))
    @printf("Root-find cost: %d spinups total across the fleet\n\n",
        sum(r.root_find_evals for r in u) + length(u))

    println("  site                    Δz  pscale  wscale  split    cycles   T_deep_true  T_deep_guess    err   jump")
    for r in u
        @printf("  %-22s %+6.0f %6.2f  %6.2f  %-8s %6d  %12.3f %12.3f %+7.2f %+6.3f\n",
            r.id, r.delta_elevation, r.precip_scale, r.wind_scale, r.split,
            r.baseline_cycles, r.deep_temperature, r.baseline_deep_temperature,
            r.baseline_deep_temperature - r.deep_temperature, r.baseline_jump)
    end
end

"""
    report_melt(truth, mp)

Score `initialize_climate_summary`'s surface-energy-balance melt estimate against the melt the
model itself produces on the equilibrated column.

This is the diagnostic that matters for the *density* half of the guess. The estimate feeds
`cs.balance`, whose **sign is the only regime discriminator** in `steady_state_profile`: where it
is non-positive the march returns solid ice at every depth. So an over-estimate of melt at a
marginal site does not merely shift the guess — it replaces a graded firn column with a block of
917 kg m-3 ice, which is the largest single error in the initial guess (measured: +296 kg m-3 of
mass-weighted mean density at `z-1500_p030_w060`, against ±12 K of deep-temperature error
elsewhere).

The reference melt is run on the site's own **converged** profile over one climatological cycle,
not on the guess: the estimate is a claim about the equilibrated site's climate, and running it
from the guess would fold the guess's own error back into the number being compared against.
Ratio rather than difference, because the failure mode is order-of-magnitude
(estimate/actual ~ 30x at the worst site), not a few percent.
"""
function report_melt(truth, mp)
    base = floor_longwave(simulate_climate_forcing("test_1", TIME_STEP_HOURS))
    u = usable(truth)
    rows = Vector{Any}(undef, length(u))

    Threads.@threads for i in eachindex(u)
        r = u[i]
        site = (id=r.id, delta_elevation=r.delta_elevation,
                precip_scale=r.precip_scale, wind_scale=r.wind_scale, split=r.split)
        cf = build_forcing(base, site)
        cfc = forcing_climatology(cf)
        cs = GEMB.initialize_climate_summary(cf, mp)

        # Melt on the converged column, measured on **both** forcings, because they are not
        # interchangeable and the difference is itself a defect:
        #
        #   * `cfc` (the climatology) is what `gemb_spinup` integrates, so its melt is what
        #     actually sets the attractor the guess is trying to hit.
        #   * `cf` (the full record) is what `initialize_climate_summary` reads.
        #
        # Melt is convex in temperature, so averaging years into one climatological cycle
        # cancels the warm excursions that carry it — the `test_1` parameter set documents this
        # for its own site. Comparing the estimate against only one of the two would attribute
        # that mismatch to the wrong place.
        profile = rebuild_profile(r)
        mp_out = GEMB.ModelParameters(;
            (f => getfield(mp, f) for f in fieldnames(GEMB.ModelParameters)
             if f != :output_frequency)..., output_frequency=:last)

        melt_on(forcing) = begin
            out = gemb(profile, forcing, mp_out; thermal_workspace=GEMB.ThermalWorkspace())
            years = length(DD.dims(forcing, DD.Ti)) * forcing.time_step / GEMB.SECONDS_PER_YEAR
            (sum(parent(out[:melt])) / years, sum(parent(out[:runoff])) / years)
        end
        melt_clim, runoff_clim = melt_on(cfc)
        melt_full, runoff_full = melt_on(cf)

        rows[i] = (r.id, r.split, cs.accumulation, cs.melt, melt_clim, melt_full,
                   runoff_full, cs.balance)
    end

    println("\nSEB melt estimate vs the model's own melt on the converged column [kg m-2 yr-1].")
    println("`melt_clim` is on the climatology the spinup integrates; `melt_full` on the full")
    println("record the estimate reads. `balance <= 0` sends the guess to the solid-ice branch.\n")
    println("  site                   split      acc   est_melt  melt_clim  melt_full   runoff  est/full   balance  ice?")
    for (id, sp, acc, est, mc, mf, ro, bal) in rows
        @printf("  %-22s %-8s %8.1f %9.1f %10.1f %10.1f %8.1f %9s %9.1f  %s\n",
            id, sp, acc, est, mc, mf, ro,
            mf > 1.0 ? @sprintf("%.1fx", est / mf) : "n/a", bal, bal <= 0 ? "ICE" : "")
    end

    for sp in (:train, :holdout, :all)
        sel = [r for r in rows if sp === :all || r[2] === sp]
        isempty(sel) && continue
        ratios = [r[4] / r[6] for r in sel if r[6] > 1.0]
        n_ice = count(r -> r[8] <= 0, sel)
        n_clim_zero = count(r -> r[5] <= 1.0 && r[6] > 1.0, sel)
        @printf("\n%-8s n=%2d  ice branch: %d  |  melt-free climatology despite real melt: %d",
            sp, length(sel), n_ice, n_clim_zero)
        isempty(ratios) || @printf("\n         est/full melt: median %.2fx  worst %.2fx",
            median(ratios), maximum(ratios))
        println()
    end
end

"""
    rebuild_profile(record) -> DimStack

The cached converged column as a profile `gemb` accepts. Only the fields the cache stores are
available, so grain state and water are reset to the bare-ice values — acceptable here because
this is used to measure *melt over one cycle* on an equilibrated column, which the surface
energy balance and the density/temperature profile determine, not the deep grain state.
"""
function rebuild_profile(r)
    m = length(r.converged_dz)
    zdim = DD.Z(1:m)
    layers = (
        dz = DD.DimArray(copy(r.converged_dz), (zdim,)),
        temperature = DD.DimArray(copy(r.converged_temperature), (zdim,)),
        density = DD.DimArray(copy(r.converged_density), (zdim,)),
        water = DD.DimArray(zeros(m), (zdim,)),
        grain_radius = DD.DimArray(fill(GEMB.RE_NEW_SNOW, m), (zdim,)),
        grain_dendricity = DD.DimArray(zeros(m), (zdim,)),
        grain_sphericity = DD.DimArray(zeros(m), (zdim,)),
        age = DD.DimArray(zeros(m), (zdim,)),
    )
    return DD.DimStack(layers)
end

"""
    refresh_summaries(truth, mp) -> Dict

Each site's cached record with its `ClimateSummary` scalars recomputed against the current
`src/` code.

The cache stores the truth (self-consistent deep temperature, converged profiles, baseline
cycles) *and* the climate summary that candidate formulas read. The truth costs ~14 minutes of
spinups and is unaffected by a change to how the summary is computed; the summary itself costs
a few forcing sweeps and **is** affected. Refreshing only the summaries means adding a field to
`ClimateSummary` — or changing one — does not force the expensive stage to be rebuilt.

Returns a **new** `Dict` rather than updating in place: the deserialized one has a concrete
`NamedTuple` value type that pins the summary's old field list, so assigning a record with an
added field into it fails to convert.
"""
function refresh_summaries(truth, mp)
    base = floor_longwave(simulate_climate_forcing("test_1", TIME_STEP_HOURS))
    ids = sort(collect(keys(truth)))
    updated = Vector{Any}(undef, length(ids))

    Threads.@threads for i in eachindex(ids)
        r = truth[ids[i]]
        if !r.usable
            updated[i] = r
        else
            site = (id=r.id, delta_elevation=r.delta_elevation,
                    precip_scale=r.precip_scale, wind_scale=r.wind_scale, split=r.split)
            cs = GEMB.initialize_climate_summary(build_forcing(base, site), mp)
            updated[i] = merge(r, (summary = summary_scalars(cs),))
        end
    end
    return Dict{String,Any}(id => updated[i] for (i, id) in enumerate(ids))
end

function report_candidates(truth)
    println("\nCandidate deep-temperature formulas, scored against the self-consistent truth.")
    println("A candidate is only acceptable if it improves BOTH train and holdout.\n")
    @printf("  %-30s %18s %18s\n", "", "train (mae/p90)", "holdout (mae/p90)")
    for (name, f) in candidates()
        tr = score(truth, f; split=:train)
        ho = score(truth, f; split=:holdout)
        @printf("  %-30s   %7.3f / %7.3f    %7.3f / %7.3f   bias %+7.3f\n",
            name, tr.mae, tr.p90, ho.mae, ho.p90, ho.bias)
    end
end

"""
    validate(mp; sites) -> nothing

Re-measure cycles-to-convergence with the current `src/` code and compare against the cached
baseline, per site. This is the only place a change to the guess is judged on the quantity it
exists to reduce; everything above judges the deep-temperature prediction that drives it.
"""
function validate(mp, truth)
    base = floor_longwave(simulate_climate_forcing("test_1", TIME_STEP_HOURS))
    u = usable(truth)
    now = Vector{Int}(undef, length(u))
    jumps = Vector{Float64}(undef, length(u))

    done = Threads.Atomic{Int}(0)
    log_lock = ReentrantLock()

    Threads.@threads for i in eachindex(u)
        r = u[i]
        ws = GEMB.ThermalWorkspace()   # per iteration; see the note in `build_truth`
        site = (id=r.id, delta_elevation=r.delta_elevation,
                precip_scale=r.precip_scale, wind_scale=r.wind_scale, split=r.split)
        cf = build_forcing(base, site)
        g = spinup(initialize_profile(mp, cf), forcing_climatology(cf), mp; ws=ws)
        now[i] = converged(g) ? cycles(g) : MAX_CYCLES
        T = parent(g[:temperature])
        jumps[i] = converged(g) ? T[end] - T[end-1] : NaN

        n = Threads.atomic_add!(done, 1) + 1
        lock(log_lock) do
            @printf("  [%2d/%2d] %-22s %d -> %d cycles\n",
                n, length(u), r.id, r.baseline_cycles, now[i])
            flush(stdout)
        end
    end

    println("\nValidation: cycles-to-convergence, current code vs cached baseline\n")
    println("  site                   split     baseline   now    delta   jump_base   jump_now")
    regressions = 0
    for (i, r) in enumerate(u)
        d = now[i] - r.baseline_cycles
        d > 0 && (regressions += 1)
        @printf("  %-22s %-8s %8d %6d %+7d %11.3f %10.3f\n",
            r.id, r.split, r.baseline_cycles, now[i], d, r.baseline_jump, jumps[i])
    end

    for sp in (:train, :holdout, :all)
        idx = [i for (i, r) in enumerate(u) if sp === :all || r.split === sp]
        isempty(idx) && continue
        b = [Float64(u[i].baseline_cycles) for i in idx]
        n = [Float64(now[i]) for i in idx]
        @printf("\n%-8s n=%3d  cycles mean %.1f -> %.1f   p90 %.1f -> %.1f   max %.0f -> %.0f\n",
            sp, length(idx), mean(b), mean(n), quantile(b, 0.9), quantile(n, 0.9),
            maximum(b), maximum(n))
        jb = filter(isfinite, [u[i].baseline_jump for i in idx])
        jn = filter(isfinite, [jumps[i] for i in idx])
        @printf("         self-consistency |jump| [K] mean %.4f -> %.4f   worst %.4f -> %.4f\n",
            mean(abs.(jb)), mean(abs.(jn)), maximum(abs.(jb)), maximum(abs.(jn)))
    end
    @printf("\n%d of %d sites regressed in cycles.\n", regressions, length(u))
end

#=============================================================================
# Entry point
=============================================================================#

const MODE = isempty(ARGS) ? "report" : ARGS[1]
const MP = ModelParameters(output_frequency=:last)

@printf("mode=%s  threads=%d\n\n", MODE, Threads.nthreads())

if MODE == "probe"
    # Cost of one spinup at the corners of the fleet, so the full run can be sized from
    # measurement. The cost per site is dominated by column depth and by how far the guess
    # starts from the attractor, both of which vary by an order of magnitude across the grid.
    base = floor_longwave(simulate_climate_forcing("test_1", TIME_STEP_HOURS))
    println("  dz        pscale  cells  depth   cycles  conv   seconds  s/cycle")
    for dz in (-1500.0, 0.0, 1200.0), ps in (0.3, 1.0, 3.0)
        site = (id="probe", delta_elevation=dz, precip_scale=ps, wind_scale=1.0, split=:train)
        cf = build_forcing(base, site)
        profile = initialize_profile(MP, cf)
        cfc = forcing_climatology(cf)
        t = @elapsed g = spinup(profile, cfc, MP; ws=GEMB.ThermalWorkspace())
        @printf("  %+7.0f %6.1f  %5d %6.1f  %6d  %-5s %8.2f %8.3f\n",
            dz, ps, length(profile[:dz]), sum(profile[:dz]),
            cycles(g), converged(g), t, t / cycles(g))
    end
elseif MODE == "truth"
    t = @elapsed truth = build_truth(MP)
    save_truth(truth)
    @printf("Built ground truth in %.1f s -> %s\n\n", t, basename(CACHE_FILE))
    report_truth(truth)
elseif MODE == "report"
    isfile(CACHE_FILE) || error("no cached truth; run with `truth` first")
    truth = refresh_summaries(load_truth(), MP)
    report_truth(truth)
    report_candidates(truth)
elseif MODE == "melt"
    isfile(CACHE_FILE) || error("no cached truth; run with `truth` first")
    report_melt(refresh_summaries(load_truth(), MP), MP)
elseif MODE == "fit"
    isfile(CACHE_FILE) || error("no cached truth; run with `truth` first")
    truth = refresh_summaries(load_truth(), MP)
    report_candidates(truth)
    validate(MP, truth)
else
    error("unknown mode $MODE; expected one of probe, truth, report, melt, fit")
end
